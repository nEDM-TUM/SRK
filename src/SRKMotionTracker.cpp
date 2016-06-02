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
#include <limits>

using namespace std;

//#define SRKTRACKDEBUG 1


SRKMotionTracker::SRKMotionTracker()
: theGeoManager("theManager", "SRK Simulation Geometry")
{
	trackID = 0;
	trackTree = nullptr;
	reflectionLimit = numeric_limits<int>::max();
	use2D = false;
	useGravity = false;
	mass = 3.30e-025; //Hg mass kg
	meanVel = 193; //in m/s
	chamberRadius = 0.235;
	chamberHeight = 0.12;
	chamberPhi = 0;
	chamberTheta = 0;
	chamberPsi = 0;
	timeLimit = 100;
	diffuseReflectionProb = 1;
	depolAtWallProb = 0;
	numTracks = 0;
	posForBranch = new TVector3();
	velForBranch = new TVector3();
	totalReflections = 0;
	lastTrack = false;
	currentEntry = 0;
	manualTracking = false;
	defaultState.vel.SetXYZ(meanVel, 0, 0);
	additionalRandomVelZ = 0.;
	velProfHistPath = "";
	velProfHist = nullptr;
	temperature = 0;
	meanFreePath = -1;
	totalGasCollisions = 0;
	theRotation.SetAngles(0, 0, 0);
//	theGeoManager = new TGeoManager("theManager", "SRK Simulation Geometry");
	vacMat = new TGeoMaterial("Vacuum", 0, 0, 0); //Removed with GeoManager
	vacMed = new TGeoMedium("Vacuum", 1, vacMat); //Removed with GeoManager
	safety = 1e-9;
	maxTrackSize = -1000; //defined as negative number for proper normal vector determination
	maxTimePerStep=numeric_limits<double>::max();

}

SRKMotionTracker::~SRKMotionTracker()
{
	gROOT->cd();
	closeTrackFile();
	delete posForBranch;
	delete velForBranch;
	delete velProfHist;

}

void SRKMotionTracker::makeCylinderGeometry()
{
	theGeoManager.ClearPhysicalNodes(true);
	double maxDist = 2.01 * sqrt(chamberRadius * chamberRadius + chamberHeight * chamberHeight);
	TGeoVolume* top = theGeoManager.MakeBox("world", vacMed, maxDist, maxDist, maxDist);

	TGeoRotation *r1 = new TGeoRotation();
	r1->SetAngles(chamberPhi, chamberTheta, chamberPsi); // all angles in degrees
	TGeoVolume* chamber = theGeoManager.MakeTube("chamber", vacMed, 0, chamberRadius, chamberHeight * 0.5);

	theGeoManager.SetTopVolume(top);
	top->AddNode(chamber, 1, r1);
	cout << "Cylinderical Geometry made with r: " << chamberRadius << " h: " << chamberHeight << " rotated by Euler angles (deg): " << chamberPhi << " " << chamberTheta << " " << chamberPsi << endl;
	theGeoManager.CloseGeometry();

}

void SRKMotionTracker::closeTrackFile()
{

	if(trackFile.IsOpen())
	{
		trackFile.Close();
	}
	else
	{
		delete trackTree;
		trackTree = nullptr;
	}
	trackTree = nullptr;
	trackID = 0;
	numTracks = 0;
	totalReflections = 0;
	totalGasCollisions = 0;
	gROOT->cd();
	currentEntry = 0;
}

void SRKMotionTracker::openTrackFile(TString inpTrackFilePath)
{
	if(trackFile.IsOpen())
	{
		closeTrackFile();
	}
	trackFile.Open(inpTrackFilePath, "RECREATE");
	trackFile.cd();
	if(trackTree == nullptr)
	{
		makeTrackTree();
	}

}

void SRKMotionTracker::makeTrackTree()
{
	trackTree = new TTree("trackTree", "Initial points and points of reflection for particle in trap");
	trackTree->Branch("trackID", &trackID, "trackID/I");
	trackTree->Branch("lastTrack", &lastTrack, "lastTrack/O");
	trackTree->Branch("time", &currentState.time, "time/D");
	trackTree->Branch("pos", "TVector3", &posForBranch);
	trackTree->Branch("vel", "TVector3", &velForBranch);
}

void SRKMotionTracker::writeTrackToFile()
{
	if(trackFile.IsZombie() || !trackFile.IsOpen())
	{
		cout << "Error track file not open." << endl;
		return;

	}
	trackFile.cd();

	trackTree->Write("", TObject::kOverwrite);
}

bool SRKMotionTracker::loadTrackFile(TString filePath)
{
	closeTrackFile(); //Close if already open
	trackFile.Open(filePath, "READ");
	if(trackFile.IsZombie() || !trackFile.IsOpen())
	{
		cout << "Error opening track file: " << filePath << endl;
		return false;
	}
	trackTree = (TTree*) trackFile.Get("trackTree");
	trackTree->SetBranchAddress("trackID", &trackID);
	trackTree->SetBranchAddress("lastTrack", &lastTrack);
	trackTree->SetBranchAddress("time", &currentState.time);
	trackTree->SetBranchAddress("pos", &posForBranch);
	trackTree->SetBranchAddress("vel", &velForBranch);
	currentEntry = 0;
	return true;

}

void SRKMotionTracker::getNextTrackTreeEntry(SRKMotionState& stateOut, int& trackIDOut, bool& lastTrackOut)
{
	trackTree->GetEntry(currentEntry);

	stateOut.pos = *posForBranch;
	stateOut.vel = *velForBranch;

	trackIDOut = trackID;
	lastTrackOut = lastTrack;
	currentEntry++;
}

void SRKMotionTracker::fillTrackTree()
{
	*posForBranch=currentState.pos;
	*velForBranch=currentState.vel;
	trackTree->Fill();
}

void SRKMotionTracker::makeTracks(int numTracksToAdd, TString inpTrackFilePath)
{
	bool useDynamicLoading = false; // For testing whether it's better to keep the tree in memory or let ROOT buffer as it goes
	if(useDynamicLoading)
	{
		if(trackTree == nullptr)
		{
			makeTrackTree();
		}
	}
	else
	{
		openTrackFile(inpTrackFilePath);
	}
	makeTracks(numTracksToAdd);

	if(useDynamicLoading)
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

TVector3 SRKMotionTracker::getRandomVelocityVector()
{
	TVector3 velOut = getRandomDirection();

	if(use2D)
	{
		velOut.SetZ(0.);
	}

	if(velProfHist != nullptr)
	{
		velOut.SetMag(velProfHist->GetRandom());
	}
	else
	{
		//use maxwell distribution with the rest of the path being the temp
		if(temperature > 0)
		{
			double maxwellVariance = TMath::Sqrt(TMath::K() * temperature / mass);
			for (int i = 0; i < 3; i++)
			{
				velOut[i] = gRandom->Gaus(0, maxwellVariance);
			}
		}
		else
		{
			velOut.SetMag(meanVel);
		}
	}

	if(additionalRandomVelZ != 0.)
	{
		velOut.SetZ(velOut.Z() + (2. * gRandom->Rndm() - 1.) * additionalRandomVelZ);
	}
	return velOut;
}

void SRKMotionTracker::getRandomVelocityVectorAndPosition(SRKMotionState& stateOut)
{

	//Begin with random points
	stateOut.pos = getRandomPointInCylinder(); //Eventually need to make this account for gravity based density
	stateOut.vel = getRandomVelocityVector();  //Also does not account for gravity

}

bool SRKMotionTracker::loadVelProfHist()
{
	delete velProfHist;
	velProfHist = nullptr;
	temperature = 0;
	if(velProfHistPath == "") return true;
	if(velProfHistPath[0] == '!')
	{
		TString temp = velProfHistPath;
		temp.Remove(0, 1);
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

	velProfHist = (TH1D*) velProfFile.Get("VelProfHist");
	velProfHist->SetDirectory(0);
	velProfFile.Close();

	cout << "Loaded velocity profile histogram, " << velProfHist->GetTitle() << ", from " << velProfHistPath << endl;

	return true;
}

void SRKMotionTracker::getInitialState(SRKMotionState& stateOut)
{
	if(!manualTracking)
	{
		getRandomVelocityVectorAndPosition(stateOut);
	}
	else
	{
		stateOut=defaultState;
	}
	stateOut.type=SRKStepPointType::START;

}

void SRKMotionTracker::makeTrack(int inpTrackID)
{

	trackID = inpTrackID;
	getInitialState(currentState);

	fillTrackTree(); //Record initial point

#ifdef SRKTRACKDEBUG
	//Initial pos/vel
	cout << "___________________________________________________________________________________________" << endl;
	cout << "Initial Position:"<< endl;
	currentState.print();
#endif

	lastTrack = false;
	SRKMotionState stateOut;
	for (currentEntry = 0; currentEntry < reflectionLimit && !lastTrack; currentEntry++)
	{

		getNextTrackingPoint(currentState,stateOut);
		currentState=stateOut;
		fillTrackTree();

	}

	return;
}

bool SRKMotionTracker::getNextTrackingPoint(const SRKMotionState& stateIn, SRKMotionState& stateOut)
{

#ifdef SRKTRACKDEBUG
	cout << "Before reflection: " << endl;
	stateIn.print();
#endif

	lastTrack = false;
	TVector3 posRef, velRef;
	double timeToWall = getNextReflection(stateIn,stateOut);

#ifdef SRKTRACKDEBUG
	cout << "Time till reflection: " << timeToReflection << endl;
#endif

	//Use Mean Free Path to calculate a time
	double meanFreePathTime = numeric_limits<double>::max();
	if(meanFreePath > 0)
	{
		double mfpDist = -meanFreePath * log(1. - gRandom->Rndm());  //Randomly determine when the next collision will happen
		meanFreePathTime = mfpDist / stateIn.vel.Mag();

#ifdef SRKTRACKDEBUG
		cout << "Mean Free Path Time: " << meanFreePathTime << endl;
#endif
	}

	//If time limit happens earlier than when the next reflection/meanfreepath has, stop it
	if(stateIn.time + timeToWall >= timeLimit && stateIn.time + meanFreePathTime >= timeLimit)
	{

		stateOut.pos += stateIn.vel * (timeLimit - stateIn.time);
		stateOut.time = timeLimit;
		stateOut.type=SRKStepPointType::TIMELIMIT;
		lastTrack = true;
		totalReflections = 0;
		totalGasCollisions = 0;
	}
	else
	{
		//Check to see if mean free path collision happens earlier
		if(meanFreePathTime < timeToWall)
		{
			stateOut.time += meanFreePathTime;
			stateOut.vel = getRandomVelocityVector(); //Either samples the velocity profile or does random direction
			stateOut.pos += stateIn.vel * meanFreePathTime;
			stateOut.type=SRKStepPointType::GASSCATTER;
			lastTrack = false;
			totalGasCollisions++;
		}
		else
		{
			stateOut.time += timeToWall;
			stateOut.vel = posRef;
			stateOut.pos = velRef;
			stateOut.type=SRKStepPointType::REFLECTION;
			lastTrack = false;
			totalReflections++;

			if(gRandom->Rndm() < depolAtWallProb)
			{
				stateOut.type=SRKStepPointType::DEPOLARIZED;
				lastTrack=true;
			}

			if(totalReflections >= reflectionLimit)
			{
				lastTrack = true;
				totalReflections = 0;
			}
		}
	}

#ifdef SRKTRACKDEBUG
	cout << "After reflection: " << endl;
	stateOut.print();
	cout << "___________________________________________________________________________________________" << endl;
#endif

	return lastTrack;
}

double SRKMotionTracker::getNextReflection(const SRKMotionState& stateIn, SRKMotionState& stateOut)
{

	double point[3] = { stateIn.pos.X(), stateIn.pos.Y(), stateIn.pos.Z() };

	TVector3 directionVec = stateIn.vel;
	directionVec.SetMag(1);
	double direction[3] = { directionVec.X(), directionVec.Y(), directionVec.Z() };

	theGeoManager.InitTrack(point, direction);

#ifdef SRKTRACKDEBUG
	cout << "Starting in: "; theGeoManager->FindNode(point[0],point[1],point[2])->Print();
	cout << "Currently in: " << theGeoManager->GetCurrentVolume()->GetName() << endl;
#endif
	theGeoManager.FindNextBoundary(maxTrackSize);
	double distance = theGeoManager.GetStep();

//	TGeoShape* theShape=theGeoManager->GetCurrentVolume()->GetShape();
//	double distance = theShape->DistFromInside(point, direction);
	distance -= safety;
	directionVec.SetMag(distance); //Now at stepsize

	stateOut.pos = stateIn.pos + directionVec;
#ifdef SRKTRACKDEBUG
	cout << "Distance to boundary: " << distance << endl;
#endif

	double minTime = distance / stateIn.vel.Mag();
#ifdef SRKTRACKDEBUG
	cout << "Time to Boundary: " << minTime << endl;
#endif
	if(minTime < 0) cout << "Error" << endl;

	point[0] = stateOut.pos.X();
	point[1] = stateOut.pos.Y();
	point[2] = stateOut.pos.Z();

	theGeoManager.SetCurrentPoint(point);
	//	theGeoManager->SetCurrentPoint(0,0,0);
	//	theGeoManager->SetCurrentDirection(0,1,0);
//	theGeoManager->Step(true,false);

#ifdef SRKTRACKDEBUG
	cout << "On Boundary: " << theGeoManager->IsOnBoundary() << endl;
	cout << theGeoManager->GetCurrentPoint()[0] << " " << theGeoManager->GetCurrentPoint()[1] << " " << theGeoManager->GetCurrentPoint()[2] << endl;
	cout << theGeoManager->GetCurrentDirection()[0] << " " << theGeoManager->GetCurrentDirection()[1] << " " << theGeoManager->GetCurrentDirection()[2] << endl;
#endif
	const double* normArray = nullptr;
//	theGeoManager->FindNextBoundary();
	normArray = theGeoManager.FindNormal();

//

//	double normArray[3];
//	theShape->ComputeNormal(point, direction, normArray);

	TVector3 norm(-normArray[0], -normArray[1], -normArray[2]);
#ifdef SRKTRACKDEBUG
	cout << "Inner Normal: "; norm.Print();
#endif

	stateOut.vel = getReflectedVector(diffuseReflectionProb, stateIn.vel, norm);
	stateOut.vel.SetMag(stateIn.vel.Mag());
#ifdef SRKTRACKDEBUG
	cout << "Reflected Vel: "; stateOut.vel.Print();
#endif

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
	double r = sqrt(gRandom->Rndm()) * chamberRadius;
	double theta = gRandom->Rndm() * 2. * TMath::Pi();
	double z = 0;
	if(!use2D)
	{
		z = (gRandom->Rndm() - 0.5) * chamberHeight;
	}

	TGeoRotation *r1 = new TGeoRotation();
	r1->SetAngles(chamberPhi, chamberTheta, chamberPsi); // all angles in degrees

	double local[3] = { r * cos(theta), r * sin(theta), z };
	double master[3];

	r1->LocalToMaster(local, master);

	return TVector3(master[0], master[1], master[2]);
}

TVector3 SRKMotionTracker::getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal)
{
	TVector3 outVec;
	if(DiffCoefficient > 0 && gRandom->Rndm() < DiffCoefficient)
	{
#ifdef SRKTRACKDEBUG
		cout << "Diffuse Reflection!" << endl;
#endif
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
#ifdef SRKTRACKDEBUG
		cout << "Specular Reflection" << endl;
#endif
		//Specular Reflection
		double NormDirection = currentDirection.Dot(normal);
		// proportional component normal to the surface
		outVec = currentDirection - 2 * normal * NormDirection;
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
	view->SetRange(-chamberRadius, -chamberRadius, -0.5 * chamberHeight, chamberRadius, chamberRadius, 0.5 * chamberHeight);

	pl1->SetLineColor(kBlue);
	pl1->Draw();
	theGeoManager.Draw("same");
}
