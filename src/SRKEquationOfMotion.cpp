#include <SRKEquationOfMotion.h>
#include <iostream>
#include <iomanip>

using namespace std;

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
	SRKEqOfMotion(x,dxdt,t);
}

