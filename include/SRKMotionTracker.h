#ifndef SRKTRACK_H_
#define SRKTRACK_H_

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TGeoShape.h"
#include "TGeoTube.h"

class SRKMotionTracker
{
public:
	SRKMotionTracker();
	virtual ~SRKMotionTracker();

	void makeTracks(int numTracksToAdd, TString inpTrackFilePath);
	void makeTracks(int numTracksToAdd);
	void drawTrack(int trackID);

	double getNextReflection(TVector3 pos0, TVector3 vel0, TVector3& posOut, TVector3& velOut); //Return time till next reflection
	bool getNextTrackingPoint(TVector3& posIn, TVector3& velIn, double& timeIn);

	TVector3 getRandomDirection();
	TVector3 getRandomPointInCylinder();
	void getRandomDirectionAndPointInCylinder(TVector3& posOut, TVector3& velOut);

	bool loadTrackFile(TString filePath); //returns true if successful
	void openTrackFile(TString inpTrackFilePath);
	void closeTrackFile();
	void writeTrackToFile();

	void getNextTrackTreeEntry(TVector3& posOut, TVector3& velOut, double& currentTimeOut, int& trackIDOut, bool& lastTrackOut);

	inline double getTimeLimit(){return timeLimit;}
	inline double getDiffuseReflectionProb(){return diffuseReflectionProb;}
	inline double getChamberRadius(){return radius;}
	inline double getChamberHeight(){return height;}
	inline double getMeanVel(){return meanVel;}
	inline int getTrackTreeEntries(){return trackTree->GetEntries();}
	inline int getTotalReflections(){return totalReflections;}
	inline TVector3 getPos(){return pos;}
	inline TVector3 getVel(){return vel;}
	inline bool getManualTracking(){return manualTracking;}
	inline int getReflectionLimit(){return reflectionLimit;}
	inline bool getUse2D(){return use2D;}
	inline double getAdditionalRandomVelZ(){return additionalRandomVelZ;}

	inline void setTimeLimit(double inp){timeLimit=inp;}
	inline void setDiffuseReflectionProb(double inp){diffuseReflectionProb=inp;}
	inline void setMeanVel(double inp){meanVel=inp;}
	inline void setReflectionLimit(int inp){reflectionLimit=inp;}
	inline void setUse2D(bool inp){use2D=inp;}
	inline void setAdditionalRandomVelZ(bool inp){additionalRandomVelZ=inp;}
	inline void setChamberRadius(double inp){delete theShape; radius=inp;theShape= new TGeoTube(0,radius,height*.5);}
	inline void setChamberHeight(double inp){delete theShape; height=inp;theShape= new TGeoTube(0,radius,height*.5);}
	inline void setPos(const TVector3& inp){ pos=inp;}
	inline void setVel(const TVector3& inp){ vel=inp;}
	inline void setManualTracking(const bool inp){ manualTracking=inp;}

protected:

	double getTimeIntersectVecInCircle(TVector2 pos0, TVector2 vel0, double radius);
	TVector3 getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal);

	TVector3 getDiffuseReflectedVector(const TVector3 normal);

	void makeTrack(int inpTrackID);

	void makeTrackTree();

	TTree* trackTree; //Track tree containing position and reflection information
	int currentEntry;
	int numTracks;
	double timeLimit; //Stop after a time limit
	int reflectionLimit; //Stop after a numer of reflections
	bool use2D;
	double additionalRandomVelZ;
	bool useGravity;
	bool manualTracking; //Instead of random
	double diffuseReflectionProb;  //0 = Specular 1=Diffuse Lambert

	//Particle
	double mass;
	double meanVel;

	//Geometry
	double radius;
	double height;
	TGeoShape* theShape;

	TFile* trackFile;

	//Tracking Variables
	TVector3 pos;
	TVector3 vel;
	double currentTime;
	int trackID;
	bool lastTrack;
	int totalReflections;

	//For branch addresses for trees...This is dumb...should probably replace with a state class
	TVector3* posTree;
	TVector3* velTree;


};

#endif /* SRKTRACK_H_ */
