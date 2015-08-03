#include "SRKSpinTracker.h"
#include <iostream>
#include <algorithm>
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
using namespace std;

//typedef runge_kutta_cash_karp54< SRKMotionState > error_stepper_type;
typedef runge_kutta_dopri5<SRKMotionState, SRKSpinFloat> error_stepper_type;

typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

//#define SRKSPINTRACKERDEBUG 1

SRKSpinTracker::SRKSpinTracker(SRKGlobalField* inpGlobalField) :
	theEquationOfMotion(inpGlobalField)
{
	eps_abs = 1.e-7;
	eps_rel = 1.e-7;
	initialStepSize = 0.01;
	stepsTaken = 0;
	constStepper = false;
}

SRKSpinTracker::~SRKSpinTracker()
{

}

void SRKSpinTracker::trackSpin(SRKMotionState& theState, double timeToTrack, std::vector<SRKMotionState>* stepRecord, std::vector<double>* stepTimes)
{
	runge_kutta4<SRKMotionState, SRKSpinFloat> stepper;
	controlled_stepper_type controlled_stepper(default_error_checker<SRKSpinFloat, range_algebra, default_operations>(SRKSpinFloat(eps_abs), SRKSpinFloat(eps_rel), SRKSpinFloat(1.0), SRKSpinFloat(1.0)));

#ifdef SRKSPINTRACKERDEBUG
	cout <<"PRESPIN TRACKING";
	printMotionState(theState);
#endif

	if(stepRecord != NULL)
	{
		if(constStepper)
		{
			integrate_const(stepper, theEquationOfMotion, theState, SRKSpinFloat(0.0), SRKSpinFloat(timeToTrack), SRKSpinFloat(initialStepSize), push_back_state_and_time(stepRecord, stepTimes));
		}
		else
		{
			integrate_adaptive(controlled_stepper, theEquationOfMotion, theState, theState[8], theState[8] + SRKSpinFloat(timeToTrack), SRKSpinFloat(initialStepSize), push_back_state_and_time(stepRecord, stepTimes));
		}

	}
	else
	{
		if(constStepper)
		{
			integrate_const(stepper, theEquationOfMotion, theState, SRKSpinFloat(0.0), SRKSpinFloat(timeToTrack), SRKSpinFloat(initialStepSize));
		}
		else
		{
			integrate_adaptive(controlled_stepper, theEquationOfMotion, theState, theState[8], theState[8] + SRKSpinFloat(timeToTrack), SRKSpinFloat(initialStepSize));
		}
	}
	theState[8] += timeToTrack;

#ifdef SRKSPINTRACKERDEBUG
	cout <<"PostSpin TRACKING";
	printMotionState(theState);
#endif

}

void SRKSpinTracker::trackSpinAltA(SRKMotionState& theState, double timeToTrack, std::vector<SRKMotionState>* stepRecord, std::vector<double>* stepTimes)
{
	SRKSpinFloat timeToTrackConv(timeToTrack);
	error_stepper_type theStepper;
	SRKSpinFloat dt = SRKSpinFloat(initialStepSize);
	SRKSpinFloat t0 = theState[8];  //Initial time simulation started at
	SRKSpinFloat deltaPhi;
	SRKSpinFloat val;
	SRKMotionState previousState(9);
	SRKMotionState stateError(9);
	bool potentialLastStep = false;
	for (;;)
	{

		previousState = theState;
		theStepper.do_step(theEquationOfMotion, theState, theState[8], dt, stateError);
		deltaPhi = theState[6] - previousState[6];
		if(deltaPhi > SRKSpinFloat(.785)) //if it rotates greater than pi slow it down!
		{
			dt = dt * SRKSpinFloat(.5) / deltaPhi;
			theState = previousState;
			potentialLastStep = false;
		}
		else
		{
			val = abs(stateError[6]);
//			val/=(eps_abs+eps_rel*abs(deltaPhi));
			val /= (SRKSpinFloat(eps_rel) * abs(deltaPhi));
			if(val > SRKSpinFloat(1.)) //too much error, slow down and repeat
			{
//				SRKSpinFloat a=(SRKSpinFloat(0.9)*SRKSpinFloat(pow(val,SRKSpinFloat(-1.)/(SRKSpinFloat(theStepper.error_order_value) - SRKSpinFloat(1.)))));
				SRKSpinFloat a = SRKSpinFloat(0);
				SRKSpinFloat b = SRKSpinFloat(0.2);
				dt = dt * max(a, b);
				theState = previousState;
				potentialLastStep = false;
			}
			else if(val < SRKSpinFloat(0.5))
			{
				theState[8] += dt;
				stepsTaken++;
				if(stepRecord != NULL && stepTimes != NULL)
				{
					stepRecord->push_back(theState);
					stepTimes->push_back(static_cast<double>(theState[8]));
				}
				if(potentialLastStep) //If potential last step worked we break!
				{
					break;
				}

//				SRKSpinFloat a=SRKSpinFloat(0.9)*SRKSpinFloat(pow(val,SRKSpinFloat(-1.)/SRKSpinFloat(theStepper.order_value )));
				SRKSpinFloat a = SRKSpinFloat(0);
				SRKSpinFloat b = SRKSpinFloat(5);
				dt = dt * max(a, b);

			}
			else
			{
				theState[8] += dt;
				stepsTaken++;
				if(stepRecord != NULL && stepTimes != NULL)
				{
					stepRecord->push_back(theState);
					stepTimes->push_back(static_cast<double>(theState[8]));
				}
				if(potentialLastStep) //If potential last step worked we break!
				{
					break;
				}
			}

		}

		if(theState[8] + dt > t0 + timeToTrackConv)
		{
			dt = t0 + timeToTrackConv - theState[8];
			potentialLastStep = true;
		}

	}

}
