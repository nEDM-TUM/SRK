#include "SRKMotionTracker.h"

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

SRKMotionTracker::SRKMotionTracker()
{
	// TODO Auto-generated constructor stub
	trackID = 0;
	trackFile = NULL;
	trackTree = NULL;
	reflectionLimit = 1000000;
	use2D = false;
	useGravity = false;
	mass = 3.30e-025; //Hg mass kg
	meanVel = 193; //in m/s
	radius = 0.235;
	height = 0.12;
	theShape = new TGeoTube(0, radius, height*.5);
	timeLimit = 100;
	diffuseReflectionProb = 100;
	numTracks = 0;
	posTree = NULL;
	velTree = NULL;
	totalReflections = 0;
	currentTime = 0.;
	lastTrack = false;
	currentEntry = 0;
	manualTracking=false;
	vel.SetXYZ(meanVel,0,0);
	additionalRandomVelZ=0.;
	velProfHistPath="";
	velProfHist=NULL;
	temperature=0;

}

SRKMotionTracker::~SRKMotionTracker()
{
	gROOT->cd();
	closeTrackFile();
	delete theShape;
	delete posTree;
	delete velTree;
	delete velProfHist;
}

void SRKMotionTracker::closeTrackFile()
{

	if(trackFile != NULL && trackFile->IsOpen())
	{
		trackFile->Close();
	}
	else
	{
		delete trackTree;
		trackTree = NULL;
	}
	delete trackFile;
	trackFile = NULL;
	trackTree = NULL;
	trackID = 0;
	numTracks = 0;
	totalReflections = 0;
	gROOT->cd();
	currentEntry = 0;

}

void SRKMotionTracker::openTrackFile(TString inpTrackFilePath)
{
	if(trackFile != NULL && trackFile->IsOpen())
	{
		closeTrackFile();
	}
	trackFile = new TFile(inpTrackFilePath, "RECREATE");
	trackFile->cd();
	if(trackTree == NULL)
	{
		makeTrackTree();
	}

}

void SRKMotionTracker::makeTrackTree()
{
	trackTree = new TTree("trackTree", "Initial points and points of reflection for particle in trap");
	trackTree->Branch("trackID", &trackID, "trackID/I");
	trackTree->Branch("lastTrack", &lastTrack, "lastTrack/O");
	trackTree->Branch("time", &currentTime, "time/D");
	trackTree->Branch("pos", "TVector3", &pos);
	trackTree->Branch("vel", "TVector3", &vel);
}

void SRKMotionTracker::writeTrackToFile()
{
	if(trackFile->IsZombie() || !trackFile->IsOpen())
	{
		cout << "Error track file not open." << endl;
		return;

	}
	trackFile->cd();

	trackTree->Write("", TObject::kOverwrite);
}

bool SRKMotionTracker::loadTrackFile(TString filePath)
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
	currentEntry = 0;
	return true;

}

void SRKMotionTracker::getNextTrackTreeEntry(TVector3& posOut, TVector3& velOut, double& currentTimeOut, int& trackIDOut, bool& lastTrackOut)
{
	trackTree->GetEntry(currentEntry);

	posOut = *posTree;

	velOut = *velTree;

	currentTimeOut = currentTime;
	trackIDOut = trackID;
	lastTrackOut = lastTrack;
	currentEntry++;
}

void SRKMotionTracker::makeTracks(int numTracksToAdd, TString inpTrackFilePath)
{
	bool useDynamic = false; // For testing whether it's better to keep the tree in memory or let ROOT buffer as it goes
	if(useDynamic)
	{
		if(trackTree == NULL)
		{
			makeTrackTree();
		}
	}
	else
	{
		openTrackFile(inpTrackFilePath);
	}
	makeTracks(numTracksToAdd);

	if(useDynamic)
	{
		openTrackFile(inpTrackFilePath);
	}

	writeTrackToFile();
	closeTrackFile();
}

void SRKMotionTracker::makeTracks(int numTracksToAdd)
{

	for (int i = 0; i < numTracksToAdd; i++)
	{

		makeTrack(i);
		numTracks++;
		if((i + 1) % 1000 == 0 || i + 1 == numTracksToAdd) cout << i + 1 << " tracks made." << endl;
	}
	cout << "Average Reflections: " << (double) totalReflections / (double) numTracksToAdd << endl;
}

void SRKMotionTracker::getRandomVelocityVectorAndPosition(TVector3& posOut, TVector3& velOut)
{

	//Begin with random points
	posOut = getRandomPointInCylinder(); //Eventually need to make this account for gravity based density
	velOut = getRandomDirection();

	if(use2D)
	{
		velOut.SetZ(0.);
		//velOut.SetMag(1);
		posOut.SetZ(0.);
	}

	if(velProfHist != NULL)
	{
		velOut.SetMag(velProfHist->GetRandom());
	}
	else
	{
		//use maxwell distribution with the rest of the path being the temp
		if(temperature > 0)
		{
			double maxwellVariance = TMath::Sqrt(TMath::K()*temperature/mass);
			for(int i = 0; i < 3; i++)
			{
				velOut[i]=gRandom->Gaus(0,maxwellVariance);
			}
		}
		else
		{
			velOut.SetMag(meanVel);
		}
	}

	if(additionalRandomVelZ != 0.)
	{
		velOut.SetZ(velOut.Z()+(2.*gRandom->Rndm()-1.)*additionalRandomVelZ);
	}

}

bool SRKMotionTracker::loadVelProfHist()
{
	delete velProfHist;
	velProfHist=NULL;
	temperature=0;
	if(velProfHistPath == "")
		return true;
	if(velProfHistPath[0] =='!')
	{
		TString temp = velProfHistPath;
		temp.Remove(0,1);
		temperature = temp.Atof();
		cout << "Using Maxwell distribution at temperature " << temperature << endl;
		return true;
	}
	TFile velProfFile(velProfHistPath, "READ");
	if(velProfFile.IsZombie() || !velProfFile.IsOpen())
	{
		cout << "Error opening track file: " << velProfHistPath << endl;
		return false;
	}

	velProfHist=(TH1D*) velProfFile.Get("VelProfHist");
	velProfHist->SetDirectory(0);
	velProfFile.Close();

	cout << "Loaded velocity profile histogram, " << velProfHist->GetTitle() << ", from " << velProfHistPath << endl;

	return true;
}


void SRKMotionTracker::makeTrack(int inpTrackID)
{

	trackID = inpTrackID;

	currentTime = 0;

	if(!manualTracking)
	{
		getRandomVelocityVectorAndPosition(pos, vel);
	}


	trackTree->Fill(); //Record initial point


#ifdef SRKTRACKDEBUG
	//Initial pos/vel
	cout << "Pos: "; pos.Print();
	cout << "Vel: "; vel.Print();
#endif

	lastTrack = false;
	for (currentEntry = 0; currentEntry < reflectionLimit && !lastTrack; currentEntry++)
	{

		getNextTrackingPoint(pos, vel, currentTime);
		trackTree->Fill();

	}


	return;
}

bool SRKMotionTracker::getNextTrackingPoint(TVector3& posIn, TVector3& velIn, double& timeIn)
{
	TVector3 posRef, velRef;
	double timeStep = getNextReflection(posIn, velIn, posRef, velRef);
	lastTrack = false;

	//If time limit happens earlier than when the next reflection has, stop it
	if(timeIn + timeStep >= timeLimit)
	{

		posIn += velIn * (timeLimit - timeIn);
		timeIn = timeLimit;
		lastTrack = true;
		totalReflections = 0;
	}
	else
	{
		timeIn += timeStep;
		posIn = posRef;
		velIn = velRef;
		lastTrack = false;
		totalReflections++;

		if(totalReflections >= reflectionLimit)
		{
			lastTrack = true;
			totalReflections = 0;
		}
	}

#ifdef SRKTRACKDEBUG
	cout << "Pos: "; posIn.Print();
	cout << "Vel: "; velIn.Print();
	cout << "Time: " << timeIn << endl;
#endif

	return lastTrack;
}

double SRKMotionTracker::getNextReflection(TVector3 pos0, TVector3 vel0, TVector3& posOut, TVector3& velOut)
{

	double point[3] = { pos0.X(), pos0.Y(), pos0.Z() };
	TVector3 directionVec = vel0;
	directionVec.SetMag(1);
	double direction[3] = { directionVec.X(), directionVec.Y(), directionVec.Z() };
	double distance = theShape->DistFromInside(point, direction);

	TVector3 step = vel0;
	step.SetMag(distance);

	posOut = pos0 + step;

	double minTime = distance / vel0.Mag();

	double normArray[3];
	point[0] = posOut.X();
	point[1] = posOut.Y();
	point[2] = posOut.Z();
	theShape->ComputeNormal(point, direction, normArray);
	TVector3 norm(-normArray[0], -normArray[1], -normArray[2]);

	velOut = getReflectedVector(diffuseReflectionProb, vel0, norm);
	velOut.SetMag(vel0.Mag());

	return minTime;

}

//Returns time to point
double SRKMotionTracker::getTimeIntersectVecInCircle(TVector2 pos0, TVector2 vel0, double radius)
{

	TVector2 posTrans = pos0.Rotate(-vel0.Phi());

	posTrans.Set(sqrt(radius * radius - posTrans.Y() * posTrans.Y()), posTrans.Y());

	TVector2 posOut = posTrans.Rotate(vel0.Phi());

	double time = (posOut.X() - pos0.X()) / vel0.X(); //Could have used X or Y here.  Should be the same
	return time;
}

TVector3 SRKMotionTracker::getRandomDirection()
{

	double phi = gRandom->Rndm() * 2. * TMath::Pi();
	double cosTheta = gRandom->Rndm() * 2. - 1.;
	double sinTheta = sqrt(1. - cosTheta * cosTheta);
	return TVector3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
}

TVector3 SRKMotionTracker::getRandomPointInCylinder()
{
	double r = sqrt(gRandom->Rndm()) * radius;
	double theta = gRandom->Rndm() * 2. * TMath::Pi();
	double z = (gRandom->Rndm() - 0.5) * height;
	return TVector3(r * cos(theta), r * sin(theta), z);
}

TVector3 SRKMotionTracker::getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal)
{
	TVector3 outVec;
	if(DiffCoefficient > 0 && gRandom->Rndm() < DiffCoefficient)
	{
		outVec = getDiffuseReflectedVector(normal);
		if(use2D) //Approximate 2D diffuse scattering by just setting z direction to zero
		{

			if(additionalRandomVelZ != 0.) //If we're in this special circumstance then specularly reflect the z direction only)
			{
				//Ensure that x y component's share is the same as before
				outVec.SetMag(currentDirection.Mag2() / (currentDirection.x() * currentDirection.x() + currentDirection.y() * currentDirection.y()));
				outVec.SetZ(currentDirection.z() * (1. - 2. * normal.z() * normal.z())); //reflect specularly the z direction
			}
			else
			{
				outVec.SetZ(0.);
			}
		}
	}
	else
	{

		//Specular Reflection
		double NormDirection = currentDirection.Dot(normal);
		// proportional component normal to the surface
		outVec=currentDirection - 2 * normal * NormDirection;
	}
	outVec.SetMag(1);
	return outVec;

}

TVector3 SRKMotionTracker::getDiffuseReflectedVector(const TVector3 normal)
{
	double theta = asin(sqrt(gRandom->Rndm()));
	double phi = gRandom->Rndm() * 2 * TMath::Pi() - TMath::Pi();
	TVector3 newDirection;
	newDirection.SetPtThetaPhi(1, theta, phi);
	newDirection.RotateUz(normal);

	return newDirection;
}

void SRKMotionTracker::drawTrack(int trackIDToDraw)
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

	bool found = false;
	for (int i = 0; i < numEntries; i++)
	{
		trackTree->GetEntry(i);
		if(trackIDToDraw == trackID)
		{
			xVector.push_back(pos->X());
			yVector.push_back(pos->Y());
			zVector.push_back(pos->Z());
			found = true;
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

	pl1->SetLineColor(kBlue);
	pl1->Draw();
	theShape->Draw("same");
}
