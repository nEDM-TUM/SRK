#ifndef SRKMANAGER_H_
#define SRKMANAGER_H_

#include "SRKMotionTracker.h"
#include "SRKSpinTracker.h"
#include "SRKGlobalField.h"
#include "SRKMotionTracker.h"

#include "TString.h"
#include "TGraphErrors.h"

const TString SRKRESULTSDIR = "/home/mjbales/work/nedm/results/";
const TString SRKHISTSDIR = "/home/mjbales/work/nedm/hists/";
const TString SRKGRAPHSDIR = "/home/mjbales/work/nedm/graphs/";
const TString SRKTRACKSDIR = "/home/mjbales/work/nedm/tracks/";

class SRKManager
{
public:
	SRKManager();
	virtual ~SRKManager();

	double trackSpins(int numTracks, TString trackFilePath, TString resultsFilePath);
	double trackSpinsDeltaOmega(int numTracks, TString runNameString, double& deltaOmegaError);
	double trackSpinsDeltaOmegaSteyerl(int numTracks, TString runNameString, double& deltaOmegaSteyerlError);
	void trackSpinsDeltaOmegaSteyerlPlot(int numTracksPerPoint, TString runNameString, int numOmega, double omegaStart, double omegaEnd, bool useLog = true, int approximateReflectionsFixedTime = 0);
	double trackSpinsFalseEDM(int numTracks, TString runNameString, double& falseEDMError);  //False EDM in e*cm

	//Getters
	inline double getMeanOmega(){return meanOmega;}
	inline double getErrorOmega(){return errorOmega;}
	inline double getGyromagneticRatio(){return theSpinTracker->getGyromagneticRatio();}
	inline double getTimeLimit(){return theMotionTracker->getTimeLimit();}
	inline double getDiffuseReflectionProb(){return theMotionTracker->getDiffuseReflectionProb();}
	inline double getChamberRadius(){return theMotionTracker->getChamberRadius();}
	inline double getChamberHeight(){return theMotionTracker->getChamberHeight();}
	inline double getMeanVel(){return theMotionTracker->getMeanVel();}
	inline double getZeta(){return bGradFieldStrength*getChamberRadius()/(2*b0FieldStrength);}
	inline double getEta(){return getChamberRadius()*getGyromagneticRatio()*e0FieldStrength/(299792458.*299792458.);}
	inline double getOmega0(){return getGyromagneticRatio()*b0FieldStrength;}
	inline int getStepsTaken(){return theSpinTracker->getStepsTaken();}
	//Setters
	inline void setRecordAllSteps(bool inp){recordAllSteps=inp;}
	inline void setUseAltStepping(bool inp){useAltStepping=inp;}
	inline void setParallelFields(bool inp){parallelFields=inp;}
	inline void setGyromagneticRatio(double inp){theSpinTracker->setGyromagneticRatio(inp);}
	inline void setTimeLimit(double inp){theMotionTracker->setTimeLimit(inp);}
	inline void setDiffuseReflectionProb(double inp){theMotionTracker->setDiffuseReflectionProb(inp);}
	inline void setMeanVel(double inp){theMotionTracker->setMeanVel(inp);}
	inline void setReflectionLimit(int inp){theMotionTracker->setReflectionLimit(inp);}
	inline void setUse2D(bool inp){theMotionTracker->setUse2D(inp);};
	inline void setChamberRadius(double inp){theMotionTracker->setChamberRadius(inp);}
	inline void setChamberHeight(double inp){theMotionTracker->setChamberHeight(inp);}
	inline void setVelByOmegaSteyerl(double OmegaSteyerl){theMotionTracker->setMeanVel(fabs(OmegaSteyerl*getGyromagneticRatio()*b0FieldStrength*getChamberRadius()));}
	inline void setB0FieldStrength(double inp){b0FieldStrength=inp;}
	inline void setE0FieldStrength(double inp){e0FieldStrength=inp;}
	inline void setBGradFieldStrength(double inp){bGradFieldStrength=inp;}
	inline void setDipoleFieldStrength(double inp){dipoleFieldStrength=inp;}
	inline void setDipolePosition(TVector3 inp){dipolePosition=inp;}
	inline void setDipoleDirection(TVector3 inp){dipoleDirection=inp;}
	inline void setPerStepError(double inp_eps_abs,double inp_eps_rel){theSpinTracker->setPerStepError(inp_eps_abs,inp_eps_rel);}
	inline void setInitialStepSize(double inp){theSpinTracker->setInitialStepSize(inp);}
	inline void setConstStepper(bool inp){theSpinTracker->setConstStepper(inp);}

protected:
	void createResultsFile(TString resultsFilePath);
	void closeResultsFile();
	void setInitialState(SRKMotionState& initialState);
	void setFinalState(SRKMotionState& finalState);
	void writeEvent();
	void writeAllSteps(std::vector<SRKMotionState>* stepRecord, std::vector<double>* stepTimes);
	void loadFields();
	double reducePeriodicNumber(double inp, double start, double end);

	TVector3 pos0, pos, vel0, vel; //For recording
	double phi0, phi, theta0, theta; //For recording
	double time0, time;
	int trackID;

	bool recordAllSteps;

	bool useAltStepping;
	bool parallelFields;

	TFile* resultsFile;
	TTree* hitTree;

	double meanOmega, errorOmega;

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
