#ifndef SRKSPINTRACKER_H_
#define SRKSPINTRACKER_H_

#include <iostream>
#include <vector>

#include "SRKEquationOfMotion.h"

#include "TString.h"

class SRKSpinTracker
{
public:
	SRKSpinTracker();
	SRKSpinTracker(SRKGlobalField* inpGlobalField);
	virtual ~SRKSpinTracker();
	void trackSpin(SRKMotionState& theState, double timeToTrack, std::vector<SRKMotionState>* stepRecord = NULL, std::vector<double>* stepTimes = NULL);
	void trackSpinAltA(SRKMotionState& theState, double timeToTrack,std::vector<SRKMotionState>* stepRecord=NULL, std::vector<double>* stepTimes=NULL);

	inline void resetStepsTaken(){stepsTaken=0;}
	inline void setEPSAbs(double inp){eps_abs =inp;}
	inline void setEPSRel(double inp){eps_rel =inp;}
	inline void setInitialStepSize(double inp){initialStepSize=inp;}
	inline void setConstStepper(bool inp){constStepper=inp;}
	inline void setGyromagneticRatio(double inp){theEquationOfMotion.setGyromagneticRatio(inp);}

	inline double getEPSAbs(){return eps_abs;}
	inline double getEPSRel(){return eps_rel;}
	inline int getStepsTaken(){return stepsTaken;}
	inline double getGyromagneticRatio(){return theEquationOfMotion.getGyromagneticRatio();}
	inline bool getConstStepper(){return constStepper;}

protected:
	double eps_abs; //absolute error in radians per step (applies to phi only)
	double eps_rel; //relative error per step (applies to phi only)
	double initialStepSize; //in time
	int stepsTaken;
	SRKEquationOfMotion theEquationOfMotion;

	bool constStepper;

};

#endif /* SRKSPINTRACKER_H_ */
