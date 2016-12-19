#ifndef SRKSextupoleField_H
#define SRKSextupoleField_H 1

#include "SRKField.h"

////////////////////////////////////////////////////////////////
/// class SRKSextupoleField
/// Sextupole field defined in the x and z directions
/// The field is zero at centerPos and diverges from there.
/// Uses scaleFactor as sextupole constant Q in [FIELDUNITS]/m²
///
/// Author: Eva Kraegeloh
///////////////////////////////////////////////////////////////

class SRKSextupoleField: public SRKField
{
public:
	SRKSextupoleField(SRKFieldSettings inpFS);
	~SRKSextupoleField(){}
	void addFieldValue(const double globalPoint[4], double fieldValue[9]);
};
#endif
