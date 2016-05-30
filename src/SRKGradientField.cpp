#include "SRKGradientField.h"

SRKGradientField::SRKGradientField(SRKFieldSettings inpFS) :
	SRKField(inpFS)
{

}

void SRKGradientField::addFieldValue(const double globalPoint[4], double fieldValue[9])
{

	const double Point[3] =
	{ globalPoint[0] - fs.centerPos[0], globalPoint[1] - fs.centerPos[1], globalPoint[2] - fs.centerPos[2] };
	fieldValue[g4FieldX] += -0.5 * Point[0] * scalingFactorWUnits;
	fieldValue[g4FieldY] += -0.5 * Point[1] * scalingFactorWUnits;
	fieldValue[g4FieldZ] += scalingFactorWUnits * Point[2];

}
