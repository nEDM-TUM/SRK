#ifndef SRKTRACK_H_
#define SRKTRACK_H_

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TGeoShape.h"
#include "TGeoTube.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TH1.h"

////////////////////////////////////////////////////////////////
/// class SRKMotionTracker
///
/// Tracks particle motion inside a geometry (i.e. chamber)
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////
class SRKMotionTracker
{
public:
	SRKMotionTracker();
	virtual ~SRKMotionTracker();

	void makeTracks(int numTracksToAdd, TString inpTrackFilePath);  /// Make tracks to a ROOT file
	void makeTracks(int numTracksToAdd); /// Make tracks to a local tree
	void drawTrack(int trackID);  /// Draw a track to a TCanvas

	void makeCylinderGeometry();  /// Make a cylindrical geometry
	double getNextReflection(TVector3 pos0, TVector3 vel0, TVector3& posOut, TVector3& velOut); //Return time till next reflection
	bool getNextTrackingPoint(TVector3& posIn, TVector3& velIn, double& timeIn); //Gets the next tracking point (typically point of reflection) returns true if it's the last track

	TVector3 getRandomDirection();  //Get's a random direction
	TVector3 getRandomPointInCylinder(); //Get's a randomly sampled point in the a cylinder
	void getRandomVelocityVectorAndPosition(TVector3& posOut, TVector3& velOut);  //gets a random velocity vector and position
	TVector3 getRandomVelocityVector();

	bool loadTrackFile(TString filePath); //returns true if successful
	void openTrackFile(TString inpTrackFilePath); //Open track file for trackTree
	void closeTrackFile(); //Close track file
	void writeTrackToFile(); //Writing track to file

	void getNextTrackTreeEntry(TVector3& posOut, TVector3& velOut, double& currentTimeOut, int& trackIDOut, bool& lastTrackOut);  //Getting data from trackTree based on currentEntry

	inline double getTimeLimit(){return timeLimit;}
	inline double getDiffuseReflectionProb(){return diffuseReflectionProb;}
	inline double getChamberRadius(){return chamberRadius;}
	inline double getChamberHeight(){return chamberHeight;}
	inline double getMeanVel(){return meanVel;}
	inline int getTrackTreeEntries(){return trackTree->GetEntries();}
	inline int getTotalReflections(){return totalReflections;}
	inline TVector3 getPos(){return pos;}
	inline TVector3 getVel(){return vel;}
	inline bool getManualTracking(){return manualTracking;}
	inline int getReflectionLimit(){return reflectionLimit;}
	inline bool getUse2D(){return use2D;}
	inline double getAdditionalRandomVelZ(){return additionalRandomVelZ;}
	inline const TString getVelProfHistPath(){return velProfHistPath;}
	inline double getMass(){return mass;}
	inline double getMeanFreePath(){return meanFreePath;}

	inline void setTimeLimit(double inp){timeLimit=inp;}
	inline void setDiffuseReflectionProb(double inp){diffuseReflectionProb=inp;}
	inline void setMeanVel(double inp){meanVel=inp;}
	inline void setReflectionLimit(int inp){reflectionLimit=inp;}
	inline void setUse2D(bool inp){use2D=inp;}
	inline void setAdditionalRandomVelZ(double inp){additionalRandomVelZ=inp;}
	inline void setChamberRadius(double inp){chamberRadius=inp; }
	inline void setChamberHeight(double inp){chamberHeight=inp; }
	inline void setChamberRotation(double phi,double theta, double psi){chamberPhi=phi;chamberTheta=theta;chamberPsi=psi;}
	inline void setPos(const TVector3& inp){ pos=inp;}
	inline void setVel(const TVector3& inp){ vel=inp;}
	inline void setManualTracking(const bool inp){ manualTracking=inp;}
	inline void setVelProfHistPath(const TString inp){velProfHistPath=inp; loadVelProfHist();}
	inline void setMass(const double inp){mass=inp;}
	inline void setMeanFreePath(const double inp){meanFreePath=inp;}

protected:

	double getTimeIntersectVecInCircle(TVector2 pos0, TVector2 vel0, double radius);  //Determine based on velocity and position and radius when particle would reach boundary.
	TVector3 getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal);  //Get the reflected vector in the geometry
	bool loadVelProfHist();  //Load the velocity profile histogram to sample from
	double getRandomVelFromProfHist();  //Randomly sample the velProfHist

	TVector3 getDiffuseReflectedVector(const TVector3 normal);  //Get a Lambert diffuse reflection given a normal surface vector

	void makeTrack(int inpTrackID);  //make a single track

	void makeTrackTree();  //Set up the track tree

	TTree* trackTree; //Track tree containing position and reflection information
	int currentEntry; //Current entry being tracked
	int numTracks;   //Number of tracks to simulate
	double timeLimit; //Stop after a time limit
	int reflectionLimit; //Stop after a numer of reflections
	bool use2D;  //Whether to simulate in 2D or 3D
	double additionalRandomVelZ; //Allows for a z velocity to be added independent from the determination of the x and y velocities
	bool useGravity;  //Whether gravity is enabled (currently not)
	bool manualTracking; //Instead of random
	double diffuseReflectionProb;  //0 = Specular 1=Diffuse Lambert

	//Particle
	double mass;
	double meanVel;
	double temperature; //Temperature of the particles (for Maxwell distributions)
	TH1* velProfHist;
	TString velProfHistPath; //the file path for velProfHist.  Note that if it begins with "!"  the number following it will be used as temperature in Kelvin for a Maxwell distribution
	double meanFreePath; //in meters, negative means no mean free path

	//Geometry
	double chamberRadius;
	double chamberHeight;
	double chamberPhi;  //Euler angle for chamber rotation
	double chamberTheta; //Euler angle for chamber rotation
	double chamberPsi; //Euler angle for chamber rotation
	TGeoRotation theRotation; //The chamber rotation
	TGeoManager* theGeoManager;  //ROOT based geomanager
	TGeoMaterial* vacMat;  //Vaccuum material
	TGeoMedium* vacMed; //Vacuum medium
	double safety;  //A small distance unimportant for the physical problem to mitigate floating point errors
	double maxTrackSize;  //Current implementation with TGeoManger requires a negative number

	TFile trackFile;

	//Tracking Variables
	TVector3 pos; //Current position
	TVector3 vel; //Current velocity
	double currentTime;
	int trackID; //current track id number
	bool lastTrack; //Whether this is the last track
	int totalReflections; //How many reflections with chamber walls
	int totalGasCollisions; //How many gas collisions

	//For branch addresses for trees...This is dumb...should probably replace with a state class
	TVector3* posForBranch;
	TVector3* velForBranch;

};

#endif /* SRKTRACK_H_ */
