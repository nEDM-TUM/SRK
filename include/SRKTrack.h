#ifndef SRKTRACK_H_
#define SRKTRACK_H_

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TGeoTube.h"

class SRKTrack
{
public:
	SRKTrack();
	virtual ~SRKTrack();

	double getReflectionCylinder(TVector3 pos0,TVector3 vel0,TVector3& posOut,TVector3& velOut); //Return true if error/limit
	void makeTracksCylinder(int numTracksToAdd);
	void makeTracksCylinder(int numTracksToAdd,TString inpTrackFilePath);
	void drawTrack(int trackID);

	TVector3 getRandomDirection();
	TVector3 getRandomPointInCylinder();
	void getRandomDirectionAndPointInCylinder(TVector3& posOut, TVector3& velOut);

	bool loadTrackFromFile(TString filePath); //returns true if succesful
	void closeTrackFile();
	void writeTrackToFile(TString outTrackFilePath);

	void getTrackTreeEntry(int entry,TVector3& posOut, TVector3& velOut, double& currentTimeOut, int& trackIDOut, bool& lastTrackOut);

	inline double getTimeLimit(){return timeLimit;}
	inline double getDiffuseReflectionProb(){return diffuseReflectionProb;}
	inline double getChamberRadius(){return radius;}
	inline double getChamberHeight(){return height;}
	inline double getMeanVel(){return meanVel;}
	inline int getTrackTreeEntries(){return trackTree->GetEntries();}
	inline int getTotalReflections(){return totalReflections;}

	inline void setTimeLimit(double inp){timeLimit=inp;}
	inline void setDiffuseReflectionProb(double inp){diffuseReflectionProb=inp;}
	inline void setMeanVel(double inp){meanVel=inp;}
	inline void setReflectionLimit(int inp){reflectionLimit=inp;}
	inline void setUse2D(bool inp){use2D=inp;};
	inline void setChamberRadius(double inp){radius=inp;cylinder.SetTubeDimensions(0,radius,height);}
	inline void setChamberHeight(double inp){height=inp;cylinder.SetTubeDimensions(0,radius,height);}


protected:


	double getTimeIntersectVecInCircle(TVector2 pos0,TVector2 vel0,double radius);
	TVector3 getReflectedVector(const double DiffCoefficient, const TVector3 currentDirection, const TVector3 normal);
	TVector3 getDiffuseReflectedVector(const TVector3 normal);

	void makeTrackCylinder(int inpTrackID);


	TTree* trackTree; //Track tree containing position and reflection information
	int numTracks;
	double timeLimit; //Stop after a time limit
	int reflectionLimit; //Stop after a numer of reflections
	bool use2D;
	bool useGravity;
	double diffuseReflectionProb;  //0 = Specular 1=Diffuse Lambert

	//Particle
	double mass;
	double meanVel;

	//Geometry
	double radius;
	double height;
	TGeoTube cylinder;

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
