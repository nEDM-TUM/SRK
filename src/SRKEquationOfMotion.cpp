#include "SRKEquationOfMotion.h"

#include <iostream>
#include <iomanip>

using namespace std;

//#define SRKEQOFMNONRELLINEARSPHEREICALDEBUG 1

SRKEquationOfMotion::SRKEquationOfMotion()
{
	theGlobalField=nullptr;
	gyromagneticRatio = -4.84578839927e7;  //Hg radians/s/T
}

SRKEquationOfMotion::SRKEquationOfMotion(SRKGlobalField* inpGlobalField)
{
	theGlobalField = inpGlobalField;
	gyromagneticRatio=0;
}

void SRKEquationOfMotion::operator()(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat t)
{
	SRKEqOfMNonRelLinearSpherical(x,dxdt,t);
}

void SRKEquationOfMotion::SRKEqOfMNonRelLinearSpherical(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat /* t */)
{
	posDouble[0] = static_cast<double> (x[0]);
	posDouble[1] = static_cast<double> (x[1]);
	posDouble[2] = static_cast<double> (x[2]);
	theGlobalField->getFieldValue(posDouble, fieldDouble);
	for(int i=0;i<9;i++)
	{
		theField[i]=fieldDouble[i];
	}
	const SRKSpinFloat motionalCoeff = 1. / (2.99792458e8 * 2.99792458e8);

	dxdt[0] = x[3];
	dxdt[1] = x[4];
	dxdt[2] = x[5];
	//acceleration
	dxdt[3] = 0.;
	dxdt[4] = 0.;
	dxdt[5] = 0.;

	//Motional magnetic field from electric field here.
	theField[0] = theField[0] - motionalCoeff * (x[4] * theField[5] - x[5] * theField[4]);
	theField[1] = theField[1] - motionalCoeff * (x[5] * theField[3] - x[3] * theField[5]);
	theField[2] = theField[2] - motionalCoeff * (x[3] * theField[4] - x[4] * theField[3]);

	SRKSpinFloat phi = x[6];
	SRKSpinFloat thetaPrime = x[7];

	SRKSpinFloat sinphi = sin(phi);
	SRKSpinFloat cosphi = cos(phi);

//	SRKSpinFloat dphi = gyromagneticRatio * (thetaPrime * (theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Presumes tan(theta)~=theta
	SRKSpinFloat dphi = gyromagneticRatio * (tan(thetaPrime) * (theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Exact equation (though issues if |theta| >pi/2 wraps around)
//	SRKSpinFloat dphi = gyromagneticRatio * ((theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Fake equation with constant
	SRKSpinFloat dthetaPrime = gyromagneticRatio * (theField[1] * cosphi - theField[0] * sinphi);
//	SRKSpinFloat dthetaPrime = 0;

	dxdt[6] = dphi;
	dxdt[7] = dthetaPrime;


#ifdef SRKEQOFMNONRELLINEARSPHEREICALDEBUG
	cout << "InitialState: ";
	printMotionState(x);
	cout << "Delta: ";
	printMotionState(dxdt);
#endif

}

