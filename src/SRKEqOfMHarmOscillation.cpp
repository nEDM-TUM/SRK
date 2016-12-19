#include <SRKEqOfMHarmOscillation.h>
#include <iostream>
#include <iomanip>

using namespace std;

//#define SRKEQOFMHARMOSCILLATIONDEBUG 1

void SRKEqOfMHarmOscillation::SRKEqOfMotion(const SRKODEState& x, SRKODEState& dxdt, const SRKSpinFloat /* t */)
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
	const double chargeToMassRatio = 1.6021766208e-19 / (105.6583745e6 * 1.783e-36);
	const double quadrupoleFactor_v = -6.002882e-3;    //from periods in g-2 paper -6.051949e-3;
	const double quadrupoleFactor_h = -0.952949e-3;    //from periods in g-2 paper -0.958701e-3;
	const double velocity_y = 0.9994*2.99792458e8;     //from g-2 paper 2.99619080538093934e8;

	dxdt[0] = x[3];
	dxdt[1] = x[4];
	dxdt[2] = x[5];
	//acceleration
	dxdt[3] = 29.3 * 29.3 * chargeToMassRatio * quadrupoleFactor_v * velocity_y * x[0];
//	dxdt[4] = chargeToMassRatio * (quadrupoleFactor_v * x[0] * x[3] + quadrupoleFactor_h * x[2] * x[5]);
	dxdt[4] = 0;
	dxdt[5] = 29.3 * 29.3 * chargeToMassRatio * quadrupoleFactor_h * velocity_y * x[2]; //gamma squared bc of omega into frame

	//Motional magnetic field from electric field here.
	theField[0] = theField[0] - motionalCoeff * (x[4] * theField[5] - x[5] * theField[4]);
	theField[1] = theField[1] - motionalCoeff * (x[5] * theField[3] - x[3] * theField[5]);
	theField[2] = theField[2] - motionalCoeff * (x[3] * theField[4] - x[4] * theField[3]);

	SRKSpinFloat phi = x[6];
	SRKSpinFloat thetaPrime = x[7];

	SRKSpinFloat sinphi = sin(phi);
	SRKSpinFloat cosphi = cos(phi);

//	SRKSpinFloat dphi = chargeToMassRatio * gyromagneticRatio * (thetaPrime * (theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Presumes tan(theta)~=theta
	SRKSpinFloat dphi = chargeToMassRatio * gyromagneticRatio * (tan(thetaPrime) * (theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Exact equation (though issues if |theta| >pi/2 wraps around)
//	SRKSpinFloat dphi = chargeToMassRatio * gyromagneticRatio * ((theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Fake equation with constant
	SRKSpinFloat dthetaPrime = chargeToMassRatio * gyromagneticRatio * (theField[1] * cosphi - theField[0] * sinphi);
//	SRKSpinFloat dthetaPrime = 0;

	dxdt[6] = dphi;
	dxdt[7] = dthetaPrime;


#ifdef SRKEQOFMHARMOSCILLATIONDEBUG
	cout << "InitialState: ";
	printMotionState(x);
	cout << "Delta: ";
	printMotionState(dxdt);
#endif

}
