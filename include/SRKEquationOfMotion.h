#ifndef SRKEQUATIONOFMOTION_H_
#define SRKEQUATIONOFMOTION_H_

#include "SRKODEState.h"
#include "SRKGlobalField.h"

////////////////////////////////////////////////////////////////
/// class SRKFieldSettings
/// Class of equation of motion for BOOST ODE.
/// This can be seen in the commented out code.
/// Also included is definitions of the SRKMotionState which defines
/// what state of the particle at a point in time.
/// This includes related helper functions
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

class SRKEquationOfMotion
{
public:
	SRKEquationOfMotion();
	SRKEquationOfMotion(SRKGlobalField* inpGlobalField);
	virtual ~SRKEquationOfMotion(){}
	virtual void operator()(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat /* t */); //THE equation of motion for use with BOOST odeint


	inline void setGlobalField(SRKGlobalField* inp){theGlobalField=inp;}
	inline void setGyromagneticRatio(double inp){ gyromagneticRatio = SRKSpinFloat(inp); }
	inline double getGyromagneticRatio(){ return static_cast<double> (gyromagneticRatio); }

protected:
	void SRKEqOfMNonRelLinearSpherical(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat /* t */);
	SRKGlobalField* theGlobalField;

	//For retrieval of field
	double fieldDouble[9];
	double posDouble[3];

	SRKSpinFloat gyromagneticRatio;
	SRKSpinFloat theField[9];

};




#endif /* SRKEQUATIONOFMOTION_H_ */
