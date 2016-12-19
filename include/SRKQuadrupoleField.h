#ifndef SRKQuadrupoleField_H
#define SRKQuadrupoleField_H 1

#include "SRKField.h"

////////////////////////////////////////////////////////////////
/// class SRKQuadrupoleField
/// Quadrupole field defined in the x and z directions
/// The field is zero at centerPos and diverges from there.
/// Uses scaleFactor as quadrupole constant K in [FIELDUNITS]/m
///
/// Author: Eva Kraegeloh
///////////////////////////////////////////////////////////////

class SRKQuadrupoleField: public SRKField
{
public:
	SRKQuadrupoleField(SRKFieldSettings inpFS);
	~SRKQuadrupoleField(){}
	void addFieldValue(const double globalPoint[4], double fieldValue[9]);
};
#endif
