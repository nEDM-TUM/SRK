#ifndef SRKMANAGER_H_
#define SRKMANAGER_H_

#include "SRKMotionTracker.h"
#include "SRKSpinTracker.h"
#include "SRKGlobalField.h"
#include "SRKMotionTracker.h"
#include "SRKRunStats.h"

#include "TString.h"
#include "TGraphErrors.h"
#include "TList.h"

const TString SRKHISTSDIR = "/home/mjbales/work/nedm/hists/";
const TString SRKGRAPHSDIR = "/home/mjbales/work/nedm/graphs/";
const TString SRKTRACKSDIR = "/home/mjbales/work/nedm/tracks/";
const TString SRKMACROSDIR = "/home/mjbales/work/nedm/macros/";


////////////////////////////////////////////////////////////////
/// SRKManager
///
/// A bit monolithic class which runs and manages the simulation
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

class SRKManager
{
public:
	SRKManager();
	~SRKManager();

	bool precessSpinsAlongTracks(int numTracks); /// Track/simulate/precess a number of tracks through the geometry
	void precessSpinsAlongTracksParAndAnti(int numTracks); /// Track/simulate a number of tracks through the geometry and then repeat with flipped electric field (for EDM)
	void loadParametersFromResultsFile(TString filePath);  ///Rerun the macro commands from a previously run results file
	void makeTracks(int numTracks);  /// Make a number of motion tracks through the geometry

	//Getters
	inline bool isRecordAllSteps(){return recordAllSteps;}
	inline bool isUseAltStepping(){ return useAltStepping;}
	inline bool isParallelFields(){ return parallelFields;}
	inline bool isUseDynamicTracking(){ return useDynamicTracking;}
	inline int getRandomSeed(){ return randomSeed;}
	inline double getB0FieldStrength(){ return b0FieldStrength;}
	inline double getE0FieldStrength(){ return e0FieldStrength;}
	inline double getBGradFieldStrength(){ return bGradFieldStrength;}
	inline double getEGradFieldStrength(){ return eGradFieldStrength;}
	inline double getDipoleFieldStrength(){ return dipoleFieldStrength;}
	inline double getPhaseMean(){return phaseMean;}
	inline double getPhaseError(){return phaseError;}
	inline double getPhiStart(){return phiStart;}
	inline double getThetaStart(){return thetaStart;}
	inline double getZeta(){return bGradFieldStrength*theMotionTracker.getChamberRadius()/(2*b0FieldStrength);}
	inline double getEta(){return theMotionTracker.getChamberRadius()*theSpinTracker.getGyromagneticRatio()*e0FieldStrength/(299792458.*299792458.);}
	inline double getOmega0(){return theSpinTracker.getGyromagneticRatio()*b0FieldStrength;}
	inline TVector3 getDipolePosition(){ return dipolePosition;}
	inline TVector3 getDipoleDirection(){ return dipoleDirection;}
	inline TVector3 getE0FieldDirection(){ return e0FieldDirection;}
	inline TVector3 getB0FieldDirection(){ return b0FieldDirection;}
	inline TString getTrackFilePath() { return trackFilePath;}
	inline TString getResultsFilePath() { return resultsFilePath;}
	inline TString getRunID() { return runID;}
	inline TString getVelProfHistPath(){return theMotionTracker.getVelProfHistPath();}
	inline SRKMotionTracker* getSRKMotionTracker(){return &theMotionTracker;}
	inline SRKSpinTracker* getSRKSpinTracker(){return &theSpinTracker;}

	//Setters
	inline void setRecordAllSteps(bool inp){recordAllSteps=inp;}
	inline void setUseAltStepping(bool inp){useAltStepping=inp;}
	inline void setParallelFields(bool inp){parallelFields=inp;}
	inline void setUseDynamicTracking(bool inp){useDynamicTracking=inp;}
	inline void setRandomSeed(const int inp){ randomSeed=inp;} /// If set to non zero, same seed will be used every simulation run until set again
	inline void setB0FieldStrength(double inp){b0FieldStrength=inp;}
	inline void setE0FieldStrength(double inp){e0FieldStrength=inp;}
	inline void setPhiStart(const double inp){phiStart=inp;}
	inline void setThetaStart(const double inp){thetaStart=inp;}
	inline void setBGradFieldStrength(double inp){bGradFieldStrength=inp;}
	inline void setEGradFieldStrength(double inp){eGradFieldStrength=inp;}
	inline void setDipoleFieldStrength(double inp){dipoleFieldStrength=inp;}
	inline void setDipolePosition(TVector3 inp){dipolePosition=inp;}
	inline void setDipoleDirection(TVector3 inp){dipoleDirection=inp;}
	inline void setTrackFilePath(TString inp) {trackFilePath=inp; if(trackFilePath=="!dynamic") useDynamicTracking=true;}
	inline void setResultsFilePath(TString inp) {resultsFilePath=inp;}
	inline void setDefaultResultsDir(TString inp) {defaultResultsDir=inp; resultsFilePath = defaultResultsDir + runID+".root";}
	inline void setRunID(TString inp) {runID=inp; resultsFilePath = defaultResultsDir + runID+".root";}
	inline void setE0FieldDirection(TVector3 inp){e0FieldDirection=inp;}
	inline void setB0FieldDirection(TVector3 inp){b0FieldDirection=inp;}


protected:
	void createResultsFile(TString resultsFilePath); /// Create and leave open the root file which stores the tree of results and the macrofile used to create them
	void closeResultsFile();  /// Closes the result file
	void setInitialState(SRKODEState& initialState); /// Takes an initial motion state and loads it into variables used by resultsTree
	void setFinalState(SRKODEState& finalState); /// Takes a final motion state and loads it into variables used by resultsTree
	void writeEvent();  /// Fills tree
	void writeAllSteps(std::vector<SRKODEState>* stepRecord, std::vector<double>* stepTimes);  /// Takes a vector of the results of each step and writes them to the tree
	void loadFields();  /// Load the electric and magnetic fields (could load ROOT files for interpolation)
	void calcDeltaPhaseMean(TString inpRunID); /// Given a parallell
	SRKRunStats calcResultsFileStats(TString filePath, bool useWrapping); /// Calculates some statistics for run file
	double calculateSzDetectionProbability(double phi, double theta); /// Calculates probability of detecting spin in the z+ direction
	bool fileExists(TString filePath);
	bool fileExistsAndNotZombie(TString strFileName);

	void precessSpinsAlongTracksDynamic(int numTracks);
	void precessSpinsAlongTracksWithTrackFile(int numTracks);

	TVector3 pos0, pos, vel0, vel; //For recording
	double phi0, phi, theta0, theta; //For recording

	double phiStart,thetaStart;
	double time0, time; //Initial time and current/final time
	int trackID;  //Current track/event
	int randomSeed;

	bool recordAllSteps; //Should all steps of the simulation be recorded to the results file

	bool useAltStepping; //Whether to use an alternate stepping method
	bool parallelFields; //Whether the electric field should be set parallel or antiparallel

	bool useDynamicTracking;
	TString trackFilePath; /// Path for tracks through geometry file
	TString resultsFilePath; /// Path for results of simulation
	TString defaultResultsDir; /// The default directory results are stored in
	TString runID; /// An identifying number/string for the simulation run

	TFile* resultsFile; /// File where results are stored
	TTree* resultsTree; /// Tree where results are stored

	double phaseMean, phaseError;
	double deltaPhaseMean,deltaPhaseError;

	SRKGlobalField theGlobalField; /// Electric and magnetic field manager
	SRKSpinTracker theSpinTracker; /// Tracks spins along a path
	SRKMotionTracker theMotionTracker; /// Creates motion tracks in a geometry

	/// Magnetic field settings that eventually should be moved to SRKGlobalField
	double b0FieldStrength;  /// Constant magnetic field
	double e0FieldStrength;  /// Constant electric field
	TVector3 b0FieldDirection; /// Direction for b0Field
	TVector3 e0FieldDirection; /// direction for e0Field
	double bGradFieldStrength; /// Strength of magnetic gradient field (z direction only for now)
	double eGradFieldStrength; /// Strength of electric gradient field (z direction only for now)
	double dipoleFieldStrength; /// Strength of magnetic dipole
	TVector3 dipolePosition; /// position of magnetic dipole
	TVector3 dipoleDirection; /// Direction of magnetic dipole

	std::vector<SRKODEState>* stepRecord = nullptr; //Used for recording all steps
	std::vector<double>* stepTimes = nullptr; //Used for recording all steps

};

#endif /* SRKMANAGER_H_ */
