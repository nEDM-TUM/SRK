#include "SRKGradientField.h"

SRKGradientField::SRKGradientField(FieldSettings inpFS):SRKField(inpFS)
{

}

void SRKGradientField::addFieldValue(const double localPoint[4], double *Bfield) const
{

	const double Point[3] = { localPoint[0] - fs.centerPos[0], localPoint[1]- fs.centerPos[1], localPoint[2] - fs.centerPos[2] };
	Bfield[g4FieldX] += -0.5*Point[0] *scalingFactorWUnits;
	Bfield[g4FieldY] += -0.5*Point[1] *scalingFactorWUnits;
	Bfield[g4FieldZ] += scalingFactorWUnits * Point[2];

}
