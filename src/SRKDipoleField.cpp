#include "SRKDipoleField.h"

using namespace std;
SRKDipoleField::SRKDipoleField(FieldSettings inpFS):SRKField(inpFS)
{

}

void SRKDipoleField::addFieldValue(const double globalPoint[4], double *Bfield) const
{
	const double Point[] = { (globalPoint[0] - fs.centerPos[0]), (globalPoint[1] - fs.centerPos[1]), (globalPoint[2] - fs.centerPos[2]) };
	const double r = sqrt(Point[0] * Point[0] + Point[1] * Point[1]  + Point[2] * Point[2] );
	const double r3inv = 1./(r*r*r);
	const double r5inv = r3inv/(r*r);
	const double mdotr = 3.*r5inv*(fs.moment.x()*Point[0] + fs.moment.y()*Point[1]  + fs.moment.z()*Point[2]);
	Bfield[g4FieldX] += scalingFactorWUnits*(Point[0] *mdotr - fs.moment.x() * r3inv);
	Bfield[g4FieldY] += scalingFactorWUnits*(Point[1] *mdotr - fs.moment.y() * r3inv);
	Bfield[g4FieldZ] += scalingFactorWUnits*(Point[2] *mdotr - fs.moment.z() * r3inv);

}



