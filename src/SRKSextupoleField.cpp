#include "SRKSextupoleField.h"
#include <iostream>

//#define SRKSEXTDEBUG 1

SRKSextupoleField::SRKSextupoleField(SRKFieldSettings inpFS) :
	SRKField(inpFS)
{

}

void SRKSextupoleField::addFieldValue(const double globalPoint[4], double fieldValue[9])
{

//	const TVector3 direction2 = fs.axisDirection.Cross(fs.direction);
	const double Point[3] =
	{ globalPoint[0] - fs.centerPos[0], globalPoint[1] - fs.centerPos[1], globalPoint[2] - fs.centerPos[2] };
	//Assuming axis in y, direction in x at the moment. Change to projection on given directions later (Point)
	fieldValue[g4FieldX] += scalingFactorWUnits*(Point[0]*Point[0] - Point[2]*Point[2] + 2*Point[0]*Point[2]);
//	fieldValue[g4FieldY] += -0.5 * Point[1] * scalingFactorWUnits;
	fieldValue[g4FieldZ] += scalingFactorWUnits*(Point[0]*Point[0] - Point[2]*Point[2] - 2*Point[0]*Point[2]);

#ifdef SRKSEXTDEBUG
			//Magnetic field values at position
			std::cout << "B Sext x:"<< double(scalingFactorWUnits*(Point[0]*Point[0] - Point[2]*Point[2] + 2*Point[0]*Point[2])) << std::endl;
			std::cout << "B Sext z:"<< double(scalingFactorWUnits*(Point[0]*Point[0] - Point[2]*Point[2] - 2*Point[0]*Point[2])) << std::endl;
#endif

}
