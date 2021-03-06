#include "SRKDipoleField.h"

using namespace std;
SRKDipoleField::SRKDipoleField(SRKFieldSettings inpFS) :
	SRKField(inpFS)
{

}

void SRKDipoleField::addFieldValue(const double globalPoint[4], double fieldValue[9])
{
	const double Point[] = { (globalPoint[0] - fs.centerPos[0]), (globalPoint[1] - fs.centerPos[1]), (globalPoint[2] - fs.centerPos[2]) };
	const double r = sqrt(Point[0] * Point[0] + Point[1] * Point[1] + Point[2] * Point[2]);
	const double r3inv = 1. / (r * r * r); //precalc for speed
	const double r5inv = r3inv / (r * r); //precal for speed
	const double mdotr = 3. * r5inv * (fs.moment.x() * Point[0] + fs.moment.y() * Point[1] + fs.moment.z() * Point[2]);
	fieldValue[g4FieldX] += scalingFactorWUnits * (Point[0] * mdotr - fs.moment.x() * r3inv);
	fieldValue[g4FieldY] += scalingFactorWUnits * (Point[1] * mdotr - fs.moment.y() * r3inv);
	fieldValue[g4FieldZ] += scalingFactorWUnits * (Point[2] * mdotr - fs.moment.z() * r3inv);

//	cout << "Point: " << Point[0] << "\t" << Point[1] << "\t" << Point[2] << endl;
//	cout << "Field: " << Bfield[0] << "\t" << Bfield[1] << "\t" << Bfield[2] << endl;

}

