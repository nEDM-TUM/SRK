#ifndef SRKEQOFMNONRELLINSPHERICAL_H_
#define SRKEQOFMNONRELLINSPHERICAL_H_

#include <SRKEquationOfMotion.h>

////////////////////////////////////////////////////////////////
/// class SRKEqOfMNonRelLinSpherical
/// Class of equation of motion for magnetic fields resulting
/// from motional electric fields (non-relativistic).
///
/// Author: Matthew Bales, Eva Krägeloh
///////////////////////////////////////////////////////////////

class SRKEqOfMNonRelLinSpherical: public SRKEquationOfMotion
{
public:
	using SRKEquationOfMotion::SRKEquationOfMotion;
	virtual ~SRKEqOfMNonRelLinSpherical(){}

private:
	void SRKEqOfMotion(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat /* t */);

};




#endif /* SRKEQOFMNONRELLINSPHERICAL_H_ */
