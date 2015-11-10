#ifndef SRKMANAGER_H_
#define SRKMANAGER_H_

#include "SRKMotionTracker.h"
#include "SRKSpinTracker.h"
#include "SRKGlobalField.h"
#include "SRKMotionTracker.h"

#include "TString.h"
#include "TGraphErrors.h"
#include "TList.h"

const TString SRKHISTSDIR = "/home/mjbales/work/nedm/hists/";
const TString SRKGRAPHSDIR = "/home/mjbales/work/nedm/graphs/";
const TString SRKTRACKSDIR = "/home/mjbales/work/nedm/tracks/";
const TString SRKMACROSDIR = "/home/mjbales/work/nedm/macros/";

struct SRKRunStats
{
	SRKRunStats():
		numEvents(0),sZDetProb(0.),
		phiMean(0.),phiError(0.), phiStDev(0.), phiKurtosis(0.),phiKurtosisError(0), phiSkewness(0.), phiSkewnessError(0.), phiTsallisPower(0.), phiTsallisPowerError(0.),
		thetaMean(0.),thetaError(0.), thetaStDev(0.), thetaKurtosis(0.),thetaKurtosisError(0), thetaSkewness(0.), thetaSkewnessError(0.), thetaTsallisPower(0.), thetaTsallisPowerError(0.) { }
	int numEvents;
	double sZDetProb;

	double phiMean;
	double phiError;
	double phiStDev;
	double phiKurtosis;
	double phiKurtosisError;
	double phiSkewness;
	double phiSkewnessError;
	double phiTsallisPower;
	double phiTsallisPowerError;

	double thetaMean;
	double thetaError;
	double thetaStDev;
	double thetaKurtosis;
	double thetaKurtosisError;
	double thetaSkewness;
	double thetaSkewnessError;
	double thetaTsallisPower;
	double thetaTsallisPowerError;

};

class SRKManager
{
public:
	SRKManager();
	virtual ~SRKManager();

	bool trackSpins(int numTracks);
	void trackSpinsDeltaOmega(int numTracks);
	void loadParametersFromResultsFile(TString filePath);
	void outputDataForRIDs(TString rangeString); //Format of int int
	void makeTracks(int numTracks);
	//TGraphErrors* trackSpinsDeltaOmegaSteyerlPlot(int numTracksPerPoint, TString runNameString, int numOmega, double omegaStart, double omegaEnd, bool useLog = true, int approximateReflectionsFixedTime = 0);

	//Getters
	inline bool getRecordAllSteps(){return recordAllSteps;}
	inline bool getUseAltStepping(){ return useAltStepping;}
	inline bool getParallelFields(){ return parallelFields;}
	inline double getB0FieldStrength(){ return b0FieldStrength;}
	inline double getE0FieldStrength(){ return e0FieldStrength;}
	inline double getBGradFieldStrength(){ return bGradFieldStrength;}
	inline double getDipoleFieldStrength(){ return dipoleFieldStrength;}
	inline TVector3 getDipolePosition(){ return dipolePosition;}
	inline TVector3 getDipoleDirection(){ return dipoleDirection;}
	inline TString getTrackFilePath() { return trackFilePath;}
	inline TString getResultsFilePath() { return resultsFilePath;}
	inline TString getRunID() { return runID;}
	inline double getPhaseMean(){return phaseMean;}
	inline double getPhaseError(){return phaseError;}

	inline double getGyromagneticRatio(){return theSpinTracker->getGyromagneticRatio();}
	inline int getStepsTaken(){return theSpinTracker->getStepsTaken();}

	inline bool getConstStepper(){return theSpinTracker->getConstStepper();}

	inline double getTimeLimit(){return theMotionTracker->getTimeLimit();}
	inline double getDiffuseReflectionProb(){return theMotionTracker->getDiffuseReflectionProb();}
	inline double getChamberRadius(){return theMotionTracker->getChamberRadius();}
	inline double getChamberHeight(){return theMotionTracker->getChamberHeight();}
	inline double getMeanVel(){return theMotionTracker->getMeanVel();}
	inline int getReflectionLimit(){return theMotionTracker->getReflectionLimit();}
	inline bool getUse2D(){return theMotionTracker->getUse2D();}
	inline double getAdditionalRandomVelZ(){return theMotionTracker->getAdditionalRandomVelZ();}
	inline bool getManualTracking(){return theMotionTracker->getManualTracking();}
	inline TVector3 getPos(){return theMotionTracker->getPos();}
	inline TVector3 getVel(){return theMotionTracker->getVel();}
	inline int getRandomSeed(){ return randomSeed;}
	inline TString getVelProfHistPath(){return theMotionTracker->getVelProfHistPath();}
	inline double getEPSAbs(){return theSpinTracker->getEPSAbs();}
	inline double getEPSRel(){return theSpinTracker->getEPSRel();}
	inline double getZeta(){return bGradFieldStrength*getChamberRadius()/(2*b0FieldStrength);}
	inline double getEta(){return getChamberRadius()*getGyromagneticRatio()*e0FieldStrength/(299792458.*299792458.);}
	inline double getOmega0(){return getGyromagneticRatio()*b0FieldStrength;}
	inline double getMass(){return theMotionTracker->getMass();}



	//Setters
	inline void setRecordAllSteps(bool inp){recordAllSteps=inp;}
	inline void setUseAltStepping(bool inp){useAltStepping=inp;}
	inline void setParallelFields(bool inp){parallelFields=inp;}
	inline void setB0FieldStrength(double inp){b0FieldStrength=inp;}
	inline void setE0FieldStrength(double inp){e0FieldStrength=inp;}
	inline void setBGradFieldStrength(double inp){bGradFieldStrength=inp;}
	inline void setDipoleFieldStrength(double inp){dipoleFieldStrength=inp;}
	inline void setDipolePosition(TVector3 inp){dipolePosition=inp;}
	inline void setDipoleDirection(TVector3 inp){dipoleDirection=inp;}
	inline void setTrackFilePath(TString inp) {trackFilePath=inp;}
	inline void setResultsFilePath(TString inp) {resultsFilePath=inp;}
	inline void setDefaultResultsDir(TString inp) {defaultResultsDir=inp;}
	inline void setRunID(TString inp) {runID=inp;}

	inline void setGyromagneticRatio(double inp){theSpinTracker->setGyromagneticRatio(inp);}
	inline void setEPSAbs(double inp){theSpinTracker->setEPSAbs(inp);}
	inline void setEPSRel(double inp){theSpinTracker->setEPSRel(inp);}
	inline void setInitialStepSize(double inp){theSpinTracker->setInitialStepSize(inp);}
	inline void setConstStepper(bool inp){theSpinTracker->setConstStepper(inp);}

	inline void setTimeLimit(double inp){theMotionTracker->setTimeLimit(inp);}
	inline void setDiffuseReflectionProb(double inp){theMotionTracker->setDiffuseReflectionProb(inp);}
	inline void setMeanVel(double inp){theMotionTracker->setMeanVel(inp);}
	inline void setReflectionLimit(int inp){theMotionTracker->setReflectionLimit(inp);}
	inline void setUse2D(bool inp){theMotionTracker->setUse2D(inp);}
	inline void setAdditionalRandomVelZ(double inp){theMotionTracker->setAdditionalRandomVelZ(inp);}
	inline void setChamberRadius(double inp){theMotionTracker->setChamberRadius(inp);}
	inline void setChamberHeight(double inp){theMotionTracker->setChamberHeight(inp);}
	inline void setVelByOmegaSteyerl(double OmegaSteyerl){theMotionTracker->setMeanVel(fabs(OmegaSteyerl*getGyromagneticRatio()*b0FieldStrength*getChamberRadius()));}
	inline void setManualTracking(bool inp){ theMotionTracker->setManualTracking(inp);}
	inline void setPos(const TVector3& inp){ theMotionTracker->setPos(inp);}
	inline void setVel(const TVector3& inp){ theMotionTracker->setVel(inp);}
	inline void setRandomSeed(const int inp){ randomSeed=inp;}
	inline void setVelProfHistPath(const TString inp){theMotionTracker->setVelProfHistPath(inp);}
	inline void setMass(const double inp){theMotionTracker->setMass(inp);}



protected:
	void createResultsFile(TString resultsFilePath);
	void closeResultsFile();
	void setInitialState(SRKMotionState& initialState);
	void setFinalState(SRKMotionState& finalState);
	void writeEvent();
	void writeAllSteps(std::vector<SRKMotionState>* stepRecord, std::vector<double>* stepTimes);
	void loadFields();
	void calcDeltaPhaseMean(TString inpRunID);
	SRKRunStats calcResultsFileStats(TString filePath, bool useWrapping);
	double calculateSzDetectionProbability(double phi, double theta);
	bool fileExists(TString filePath);
	bool fileExistsAndNotZombie(TString strFileName);


	TVector3 pos0, pos, vel0, vel; //For recording
	double phi0, phi, theta0, theta; //For recording
	double time0, time;
	int trackID;
	int randomSeed;

	bool recordAllSteps;

	bool useAltStepping;
	bool parallelFields;

	TString trackFilePath;
	TString resultsFilePath;
	TString defaultResultsDir;
	TString runID;

	TFile* resultsFile;
	TTree* hitTree;

	double phaseMean, phaseError;
	double deltaPhaseMean,deltaPhaseError;

	SRKGlobalField* theGlobalField;
	SRKSpinTracker* theSpinTracker;
	SRKMotionTracker* theMotionTracker;

	//Magnetic field settings that eventually should be moved to SRKGlobalField
	double b0FieldStrength;
	double e0FieldStrength;
	double bGradFieldStrength;
	double dipoleFieldStrength;
	TVector3 dipolePosition;
	TVector3 dipoleDirection;

};

#endif /* SRKMANAGER_H_ */
