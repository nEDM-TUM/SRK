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
	inline void setPerStepError(double inp_eps_abs,double inp_eps_rel){eps_abs =inp_eps_abs; eps_rel =inp_eps_rel;}
	inline void setInitialStepSize(double inp){initialStepSize=inp;}
	inline void setConstStepper(bool inp){constStepper=inp;}
	inline void setGyromagneticRatio(double inp){theEquationOfMotion.setGyromagneticRatio(inp);}
	inline int getStepsTaken(){return stepsTaken;}
	inline double getGyromagneticRatio(){return theEquationOfMotion.getGyromagneticRatio();}

protected:
	double eps_abs; //absolute error in radians per step (applies to phi only)
	double eps_rel; //relative error per step (applies to phi only)
	double initialStepSize; //in time
	int stepsTaken;
	SRKEquationOfMotion theEquationOfMotion;

	bool constStepper;

};

#endif /* SRKSPINTRACKER_H_ */
