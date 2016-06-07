#ifndef SRKTRACK_H_
#define SRKTRACK_H_

#include <array>
#include <vector>

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

#include "SRKMotionState.h"

enum class SRKPointType { TIMELIMIT, PERIODICSTOP, WALL, GASSCATTER, COUNT};

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
	bool getNextTrackingPoint(const SRKMotionState& stateIn, SRKMotionState& stateOut); // Gets the next tracking point (typically point of reflection) returns true if it's the last track
	void getNextReflection(const SRKMotionState& stateIn, SRKMotionState& stateOut); /// Return time till next reflection

	TVector3 getRandomDirection();  //G/ et's a random direction
	TVector3 getRandomPointInCylinder(); // /Get's a randomly sampled point in the a cylinder
	void getRandomVelocityVectorAndPosition(SRKMotionState& stateOut);  /// gets a random velocity vector and position
	TVector3 getRandomVelocityVector();
	void getInitialState(SRKMotionState& stateOut); //If manualTracking, then sets state with defaults, otherwise gets it randomly

	bool loadTrackFile(TString filePath); /// returns true if successful
	void openTrackFile(TString inpTrackFilePath); /// Open track file for trackTree
	void closeTrackFile(); /// Close track file
	void writeTrackToFile(); /// Writing track to file

	void getNextTrackTreeEntry(SRKMotionState& stateOut, int& trackIDOut, bool& lastTrackOut);  /// Getting data from trackTree based on currentEntry

	inline const double getTimeLimit(){return timeLimit;}
	inline const double getDiffuseReflectionProb(){return diffuseReflectionProb;}
	inline const double getChamberRadius(){return chamberRadius;}
	inline const double getChamberHeight(){return chamberHeight;}
	inline const double getMeanVel(){return meanVel;}
	inline const int getTrackTreeEntries(){return trackTree->GetEntries();}
	inline const int getTotalReflections(){return totalReflections;}
	inline const TVector3 getDefaultPos(){return defaultState.pos;}
	inline const TVector3 getDefaultVel(){return defaultState.vel;}
	inline const bool isManualTracking(){return manualTracking;}
	inline const int getReflectionLimit(){return reflectionLimit;}
	inline const bool isUse2D(){return use2D;}
	inline const double getAdditionalRandomVelZ(){return additionalRandomVelZ;}
	inline const TString getVelProfHistPath(){return velProfHistPath;}
	inline const double getMass(){return mass;}
	inline const double getMeanFreePath(){return meanFreePath;}
	inline const double getDepolAtWallProb(){return depolAtWallProb;}
	inline const double getPeriodicStopTime(){return periodicStopTime;}

	inline void setTimeLimit(double inp){timeLimit=inp;}
	inline void setDiffuseReflectionProb(double inp){diffuseReflectionProb=inp;}
	inline void setMeanVel(double inp){meanVel=inp;}
	inline void setReflectionLimit(int inp){reflectionLimit=inp;}
	inline void setUse2D(bool inp){use2D=inp;}
	inline void setAdditionalRandomVelZ(double inp){additionalRandomVelZ=inp;}
	inline void setChamberRadius(double inp){chamberRadius=inp; }
	inline void setChamberHeight(double inp){chamberHeight=inp; }
	inline void setChamberRotation(double phi,double theta, double psi){chamberPhi=phi;chamberTheta=theta;chamberPsi=psi;}
	inline void setDefaultPos(const TVector3& inp){ defaultState.pos=inp;}
	inline void setDefaultVel(const TVector3& inp){ defaultState.vel=inp;}
	inline void setManualTracking(const bool inp){ manualTracking=inp;}
	inline void setVelProfHistPath(const TString inp){velProfHistPath=inp; loadVelProfHist();}
	inline void setMass(const double inp){mass=inp;}
	inline void setMeanFreePath(const double inp){meanFreePath=inp;}
	inline void setDepolAtWallProb(const double inp){depolAtWallProb=inp;}
	inline void setPeriodicStopTime(const double inp){periodicStopTime=inp;}

protected:

	double getTimeIntersectVecInCircle(TVector2 pos0, TVector2 vel0, double radius);  /// Determine based on velocity and position and radius when particle would reach boundary.
	TVector3 getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal);  /// Get the reflected vector in the geometry
	bool loadVelProfHist();  /// Load the velocity profile histogram to sample from
	double getRandomVelFromProfHist();  /// Randomly sample the velProfHist
	void fillTrackTree();

	TVector3 getDiffuseReflectedVector(const TVector3 normal);  /// Get a Lambert diffuse reflection given a normal surface vector


	void makeTrack(int inpTrackID);  /// make a single track

	void makeTrackTree();  /// Set up the track tree

	TTree* trackTree; /// Track tree containing position and reflection information
	int currentEntry; /// Current entry being tracked
	int numTracks;   /// Number of tracks to simulate
	double timeLimit; /// Stop after a time limit
	int reflectionLimit; /// Stop after a numer of reflections
	bool use2D;  /// Whether to simulate in 2D or 3D
	double additionalRandomVelZ; /// Allows for a z velocity to be added independent from the determination of the x and y velocities
	bool useGravity;  /// Whether gravity is enabled (currently not)
	bool manualTracking; /// Instead of random
	double diffuseReflectionProb;  /// 0 = Specular 1=Diffuse Lambert
	double depolAtWallProb;  /// Probability of depolarization at the wall upon reflection

	//Particle
	double mass;
	double meanVel;
	double temperature; /// Temperature of the particles (for Maxwell distributions)
	TH1* velProfHist;
	TString velProfHistPath; /// the file path for velProfHist.  Note that if it begins with "!"  the number following it will be used as temperature in Kelvin for a Maxwell distribution
	double meanFreePath; /// in meters, negative means no mean free path

	//Geometry
	double chamberRadius;
	double chamberHeight;
	double chamberPhi;  /// Euler angle for chamber rotation
	double chamberTheta; /// Euler angle for chamber rotation
	double chamberPsi; /// Euler angle for chamber rotation
	TGeoRotation theRotation; /// The chamber rotation
	TGeoManager theGeoManager;  /// ROOT based geomanager
	TGeoMaterial* vacMat;  /// Vaccuum material
	TGeoMedium* vacMed; /// Vacuum medium
	double safety;  /// A small distance unimportant for the physical problem to mitigate floating point errors
	double maxTrackSize;  /// Current implementation with TGeoManger requires a negative number

	TFile* trackFile;

	//Tracking Variables
	SRKMotionState currentState;
	SRKMotionState defaultState;
	int trackID; /// current track id number
	bool lastTrack; /// Whether this is the last track
	int totalReflections; /// How many reflections with chamber walls
	int totalGasCollisions; /// How many gas collisions
	double periodicStopTime; /// Max time per track
	double nextPeriodicStop; // Last time a periodic stop was made
//	std::array<double,(int) int(SRKPointType::COUNT> timeChecks;
	std::vector<double> timeChecks;

	//For branch addresses for trees...This is dumb, but necessary
	TVector3* posForBranch;
	TVector3* velForBranch;
	char* typeForBranch; //Uses SRKStepPointType
};

#endif /* SRKTRACK_H_ */
