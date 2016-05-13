#include "SRKEquationOfMotion.h"
#include <iostream>
#include <iomanip>

using namespace std;

SRKEquationOfMotion::SRKEquationOfMotion()
{
	theGlobalField=NULL;
	gyromagneticRatio = -4.84578839927e7;  //Hg radians/s/T
}

SRKEquationOfMotion::SRKEquationOfMotion(SRKGlobalField* inpGlobalField)
{
	theGlobalField = inpGlobalField;
	gyromagneticRatio=0;
}

SRKEquationOfMotion::~SRKEquationOfMotion()
{

}

void SRKEquationOfMotion::operator()(const SRKMotionState& x, SRKMotionState& dxdt, const SRKSpinFloat t)
{
	SRKEqOfMNonRelLinearSpherical(x,dxdt,t);
}

void SRKEquationOfMotion::SRKEqOfMNonRelLinearSpherical(const SRKMotionState& x, SRKMotionState& dxdt, const SRKSpinFloat /* t */)
{
	posDouble[0] = static_cast<double> (x[0]);
	posDouble[1] = static_cast<double> (x[1]);
	posDouble[2] = static_cast<double> (x[2]);
	theGlobalField->GetFieldValue(posDouble, fieldDouble);
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
	SRKSpinFloat dphi = gyromagneticRatio * (tan(thetaPrime) * (theField[0] * cosphi + theField[1] * sinphi) - theField[2]); //Exact equation (though issues if theta wraps around)
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

void setMotionState(SRKMotionState& outState, TVector3* pos, TVector3* vel, double phi, double theta)
{
	outState[0] = SRKSpinFloat(pos->X());
	outState[1] = SRKSpinFloat(pos->Y());
	outState[2] = SRKSpinFloat(pos->Z());
	outState[3] = SRKSpinFloat(vel->X());
	outState[4] = SRKSpinFloat(vel->Y());
	outState[5] = SRKSpinFloat(vel->Z());
	outState[6] = SRKSpinFloat(phi);
	outState[7] = SRKSpinFloat(theta);
}

void updateMotionStatePosVel(SRKMotionState& outState, TVector3 pos, TVector3 vel, double currentTime)
{
	outState[0] = pos.X();
	outState[1] = pos.Y();
	outState[2] = pos.Z();
	outState[3] = vel.X();
	outState[4] = vel.Y();
	outState[5] = vel.Z();
	outState[8] = currentTime;
}

void printMotionState(const SRKMotionState& theState)
{
	cout << "----------------------------" << endl;
	cout << "Pos: " << static_cast<double>(theState[0]) << ", " << static_cast<double>(theState[1]) << ", " << static_cast<double>(theState[2]) << endl;
	cout << "Vel: " << static_cast<double>(theState[3]) << ", " << static_cast<double>(theState[4]) << ", " << static_cast<double>(theState[5]) << endl;
	cout << "Phi: " << static_cast<double>(theState[6]) << "   Theta: " << static_cast<double>(theState[7]) << "   Time: " << static_cast<double>(theState[8]) << endl;
	cout << "----------------------------" << endl;
}
