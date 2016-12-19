#ifndef SRKEQOFMHARMOSCILLATION_H_
#define SRKEQOFMHARMOSCILLATION_H_

#include <SRKEquationOfMotion.h>

////////////////////////////////////////////////////////////////
/// class SRKEqOfMHarmOscillation
/// Class of equation of motion of harmonic oscillation.
///
/// Author: Eva Krägeloh
///////////////////////////////////////////////////////////////

class SRKEqOfMHarmOscillation: public SRKEquationOfMotion
{
public:
	using SRKEquationOfMotion::SRKEquationOfMotion;
	virtual ~SRKEqOfMHarmOscillation(){}

private:
	void SRKEqOfMotion(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat /* t */);

};




#endif /* SRKEQOFMHARMOSCILLATION_H_ */
