#ifndef SRKSPINTRACKER_H_
#define SRKSPINTRACKER_H_

#include <SRKEqOfMNonRelLinSpherical.h>
#include <SRKEqOfMHarmOscillation.h>
#include <iostream>
#include <vector>

#include "TString.h"

////////////////////////////////////////////////////////////////
/// class SRKSpinTracker
///
/// Tracks spin precession between two points using ODEINT
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

class SRKSpinTracker
{
public:
	SRKSpinTracker();
	SRKSpinTracker(SRKGlobalField* theField);
	virtual ~SRKSpinTracker();

	void trackSpin(SRKODEState& theState, double timeToTrack, std::vector<SRKODEState>* stepRecord = nullptr, std::vector<double>* stepTimes = nullptr);  //Primary spin tracker
	void trackSpinAltA(SRKODEState& theState, double timeToTrack,std::vector<SRKODEState>* stepRecord=nullptr, std::vector<double>* stepTimes=nullptr); //Alternate spin tracker (more manual)

	inline void resetStepsTaken(){stepsTaken=0;}
	inline void setEPSAbs(double inp){eps_abs =inp;}
	inline void setEPSRel(double inp){eps_rel =inp;}
	inline void setInitialStepSize(double inp){initialStepSize=inp;}
	inline void setConstStepper(bool inp){constStepper=inp;}
	inline void setGyromagneticRatio(double inp){theEquationOfMotion.setGyromagneticRatio(inp);}
	inline void setGlobalField(SRKGlobalField* inp){theEquationOfMotion.setGlobalField(inp);}

	inline double getEPSAbs(){return eps_abs;}
	inline double getEPSRel(){return eps_rel;}
	inline int getStepsTaken(){return stepsTaken;}
	inline double getGyromagneticRatio(){return theEquationOfMotion.getGyromagneticRatio();}
	inline bool isConstStepper(){return constStepper;}

protected:
	double eps_abs; //absolute error in radians per step (applies to phi only)
	double eps_rel; //relative error per step (applies to phi only)
	double initialStepSize; //in time
	int stepsTaken;  //Number of steps taken
	SRKEqOfMHarmOscillation theEquationOfMotion; //The equation of motion

	bool constStepper; //Whether to use the constant stepper algorithm

};

#endif /* SRKSPINTRACKER_H_ */
