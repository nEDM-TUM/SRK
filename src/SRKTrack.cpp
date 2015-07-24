#include "SRKTrack.h"

#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TPolyLine3D.h"
#include "TView3D.h"
#include "TCanvas.h"
#include "TGeoTube.h"

#include <iostream>
#include <cfloat>
#include <vector>

using namespace std;

//#define SRKTRACKDEBUG 1

enum CYLSURF
{
	CYLSURF_RADIAL, CYLSURF_TOP, CYLSURF_BOTTOM
};

SRKTrack::SRKTrack()
{
	// TODO Auto-generated constructor stub
	trackID = 0;
	trackTree = NULL;
	trackFile = NULL;
	reflectionLimit = 1000000;
	use2D = true;
	useGravity = false;
	mass = 3.30e-025; //Hg mass kg
	meanVel = 193; //in m/s
	radius = 0.24;
	height = 0.1;
	cylinder.SetTubeDimensions(0, radius, height);
	timeLimit = 100;
	diffuseReflectionProb = 100;
	numTracks = 0;
	posTree = NULL;
	velTree = NULL;
	totalReflections=0;

}

SRKTrack::~SRKTrack()
{
	closeTrackFile();
	delete posTree;
	delete velTree;
}

void SRKTrack::closeTrackFile()
{

	if(trackFile != NULL && trackFile->IsOpen())
	{
		trackFile->Close();
	}
	else
	{
		delete trackTree;
	}
	delete trackFile;
	trackFile = NULL;
	trackTree = NULL;


	trackID = 0;
	numTracks = 0;
	totalReflections=0;

}

void SRKTrack::writeTrackToFile(TString inpTrackFilePath)
{
	if(trackFile != NULL && trackFile->IsOpen())
	{
		trackFile->Close();
	}
	trackFile = new TFile(inpTrackFilePath, "RECREATE");

	if(trackFile->IsZombie() || !trackFile->IsOpen())
	{
		cout << "Error opening track file: " << inpTrackFilePath << endl;
		return;

	}
	trackFile->cd();

	//trackFile->WriteTObject(trackTree);

	trackTree->Write("", TObject::kOverwrite);

	trackFile->Close();
	delete trackFile;
	trackFile = NULL;
	trackTree = NULL;
}

bool SRKTrack::loadTrackFromFile(TString filePath)
{
	closeTrackFile(); //Close if already open
	trackFile = new TFile(filePath, "READ");
	if(trackFile->IsZombie() || !trackFile->IsOpen())
	{
		cout << "Error opening track file: " << filePath << endl;
		return false;

	}
	trackTree = (TTree*) trackFile->Get("trackTree");
	trackTree->SetBranchAddress("trackID", &trackID);
	trackTree->SetBranchAddress("lastTrack", &lastTrack);
	trackTree->SetBranchAddress("time", &currentTime);
	trackTree->SetBranchAddress("pos", &posTree);
	trackTree->SetBranchAddress("vel", &velTree);
	return true;

}

void SRKTrack::getTrackTreeEntry(int entry, TVector3& posOut, TVector3& velOut, double& currentTimeOut, int& trackIDOut, bool& lastTrackOut)
{
	trackTree->GetEntry(entry);

	if(posTree == NULL)
		posOut = pos;
	else
		posOut = *posTree;

	if(velTree == NULL)
		velOut = vel;
	else
		velOut = *velTree;

	currentTimeOut = currentTime;
	trackIDOut = trackID;
	lastTrackOut = lastTrack;
}

void SRKTrack::makeTracksCylinder(int numTracksToAdd, TString inpTrackFilePath)
{

	makeTracksCylinder(numTracksToAdd);
	writeTrackToFile(inpTrackFilePath);

}

void SRKTrack::makeTracksCylinder(int numTracksToAdd)
{
	if(trackTree == NULL)
	{
		trackTree = new TTree("trackTree", "Initial points and points of reflection for particle in trap");
		trackTree->Branch("trackID", &trackID, "trackID/I");
		trackTree->Branch("lastTrack", &lastTrack, "lastTrack/O");
		trackTree->Branch("time", &currentTime, "time/D");
		trackTree->Branch("pos", "TVector3", &pos);
		trackTree->Branch("vel", "TVector3", &vel);
	}

	for (int i = 0; i < numTracksToAdd; i++)
	{

		makeTrackCylinder(i);
		numTracks++;
		if(i + 1 % 1000 == 0 || i +1 == numTracksToAdd)
			cout << i+1 << " tracks made." << endl;
	}
	cout << "Average Reflections: " << (double) totalReflections/ (double) numTracksToAdd << endl;
}

void SRKTrack::getRandomDirectionAndPointInCylinder(TVector3& posOut, TVector3& velOut)
{
	//Begin with random points
	posOut = getRandomPointInCylinder(); //Eventually need to make this account for gravity based density
	velOut = getRandomDirection();
//	posOut=TVector3();
//	velOut=TVector3(1,0,0);

	//pos.Print();
	//vel.Print();

	if(use2D)
	{
		velOut.SetZ(0.);
		velOut.SetMag(1);
		posOut.SetZ(0.);
	}
	velOut *= meanVel;
}

void SRKTrack::makeTrackCylinder(int inpTrackID)
{

	trackID = inpTrackID;
	lastTrack = false;
	currentTime = 0;

	getRandomDirectionAndPointInCylinder(pos, vel);

	TVector3 posOut, velOut;
	bool timeLimitReached = false;

	trackTree->Fill(); //Record initial point

#ifdef SRKTRACKDEBUG
	cout << "Pos: "; pos.Print();
	cout << "Vel: "; vel.Print();
#endif

	for (int i = 0; i < reflectionLimit && !timeLimitReached; i++)
	{
		double timeStep = getReflectionCylinder(pos, vel, posOut, velOut);

		//If time limit happens earlier than when the next reflection has, stop it
		if(currentTime + timeStep >= timeLimit)
		{
			velOut = vel;
			posOut = pos + vel * (timeLimit - currentTime);
			currentTime = timeLimit;
			timeLimitReached = true;
			lastTrack = true;

		}
		else
		{
			currentTime += timeStep;
		}

		pos = posOut;
		vel = velOut;

#ifdef SRKTRACKDEBUG
		cout << "Pos: "; pos.Print();
		cout << "Vel: "; vel.Print();
#endif

		if(i == reflectionLimit - 1)
		{
			lastTrack = true;
		}

		trackTree->Fill();
		totalReflections++;
	}

	return;
}

double SRKTrack::getReflectionCylinder(TVector3 pos0, TVector3 vel0, TVector3& posOut, TVector3& velOut)
{

	double point[3] =
	{ pos0.X(), pos0.Y(), pos0.Z() };
	TVector3 directionVec = vel0;
	directionVec.SetMag(1);
	double direction[3] =
	{ directionVec.X(), directionVec.Y(), directionVec.Z() };
	double distance = cylinder.DistFromInside(point, direction);

	TVector3 step = vel0;
	step.SetMag(distance);

	posOut = pos0 + step;

	double minTime = distance / vel0.Mag();

	double normArray[3];
	point[0] = posOut.X();
	point[1] = posOut.Y();
	point[2] = posOut.Z();
	cylinder.ComputeNormal(point, direction, normArray);
	TVector3 norm(-normArray[0], -normArray[1], -normArray[2]);

	velOut = getReflectedVector(diffuseReflectionProb, vel0, norm);
	return minTime;

}

//Returns time to point
double SRKTrack::getTimeIntersectVecInCircle(TVector2 pos0, TVector2 vel0, double radius)
{

	TVector2 posTrans = pos0.Rotate(-vel0.Phi());

	posTrans.Set(sqrt(radius * radius - posTrans.Y() * posTrans.Y()), posTrans.Y());

	TVector2 posOut = posTrans.Rotate(vel0.Phi());

	double time = (posOut.X() - pos0.X()) / vel0.X(); //Could have used X or Y here.  Should be the same
	return time;
}

TVector3 SRKTrack::getRandomDirection()
{

	double phi = gRandom->Rndm() * 2. * TMath::Pi();
	double cosTheta = gRandom->Rndm() * 2. - 1.;
	double sinTheta = sqrt(1. - cosTheta * cosTheta);
	return TVector3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
}

TVector3 SRKTrack::getRandomPointInCylinder()
{
	double r = sqrt(gRandom->Rndm()) * radius;
	double theta = gRandom->Rndm() * 2. * TMath::Pi();
	double z = (gRandom->Rndm() - 0.5) * height;
	return TVector3(r * cos(theta), r * sin(theta), z);
}

TVector3 SRKTrack::getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal)
{
	if(DiffCoefficient > 0 && gRandom->Rndm() < DiffCoefficient)
	{

		TVector3 outVec = getDiffuseReflectedVector(normal);
		if(use2D) //Approximate 2D diffuse scattering by just setting z direction to zero
		{
			outVec.SetZ(0.);
			outVec.SetMag(currentDirection.Mag());
		}
		return outVec;
	}

	//Specular Reflection
	double NormDirection = currentDirection.Dot(normal);
	// proportional component normal to the surface
	return currentDirection - 2 * normal * NormDirection;

}

TVector3 SRKTrack::getDiffuseReflectedVector(const TVector3 normal)
{
	double theta = asin(sqrt(gRandom->Rndm()));
	double phi = gRandom->Rndm() * 2 * TMath::Pi() - TMath::Pi();
	TVector3 newDirection;
	newDirection.SetPtThetaPhi(1, theta, phi);
	newDirection.RotateUz(normal);

	return newDirection;
}

void SRKTrack::drawTrack(int trackIDToDraw)
{
	//First need to convert track to array
	TVector3* pos;
	int trackID;

	trackTree->SetBranchAddress("trackID", &trackID);
	trackTree->SetBranchAddress("pos", &pos);

	const int numEntries = trackTree->GetEntries();
	vector<double> xVector;
	vector<double> yVector;
	vector<double> zVector;

	bool found=false;
	for (int i = 0; i < numEntries; i++)
	{
		trackTree->GetEntry(i);
		if(trackIDToDraw == trackID)
		{
			xVector.push_back(pos->X());
			yVector.push_back(pos->Y());
			zVector.push_back(pos->Z());
			found=true;
		}
		else
		{
			if(found)
			{
				break;  //Break loop if we found them all (they should all be in a row)
			}
		}

	}

	TPolyLine3D* pl1 = new TPolyLine3D(xVector.size(), &xVector[0], &yVector[0], &zVector[0]);
	TCanvas* c1 = new TCanvas("c1");
	c1->cd();
	TView3D* view = (TView3D*) TView::CreateView(1);

	view->SetRange(-radius, -radius, -0.5 * height, radius, radius, 0.5 * height);
	TGeoTube* cylinder = new TGeoTube(0., radius, 0.5 * height);

	pl1->SetLineColor(kBlue);
	pl1->Draw();
	cylinder->Draw("same");
}
