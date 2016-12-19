#include "SRKQuadrupoleField.h"
#include <iostream>

//#define SRKQUADDEBUG 1

SRKQuadrupoleField::SRKQuadrupoleField(SRKFieldSettings inpFS) :
	SRKField(inpFS)
{

}

void SRKQuadrupoleField::addFieldValue(const double globalPoint[4], double fieldValue[9])
{

//	const TVector3 direction2 = fs.axisDirection.Cross(fs.direction);
	const double Point[3] =
	{ globalPoint[0] - fs.centerPos[0], globalPoint[1] - fs.centerPos[1], globalPoint[2] - fs.centerPos[2] };
	//Assuming axis in y, direction in x at the moment. Change to projection on given directions later (Point)
	fieldValue[g4FieldX] += scalingFactorWUnits*Point[2];
//	fieldValue[g4FieldY] += -0.5 * Point[1] * scalingFactorWUnits;
	fieldValue[g4FieldZ] += scalingFactorWUnits*Point[0];

#ifdef SRKQUADDEBUG
			//Magnetic field values at position
			std::cout << "B Quad x:"<< double(scalingFactorWUnits*Point[2]) << std::endl;
			std::cout << "B Quad z:"<< double(scalingFactorWUnits*Point[0]) << std::endl;
#endif

}
