#ifndef SRKGradientField_H
#define SRKGradientField_H 1

#include "SRKField.h"

////////////////////////////////////////////////////////////////
/// class SRKGradientField
/// Gradient field defined in the z direction (Later will add rotation)
/// The field is of the set field strength at centerPos and diverges from there.
/// Uses scaleFactor as gradient strength in [FIELDUNITS]/m
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

class SRKGradientField: public SRKField
{
public:
	SRKGradientField(SRKFieldSettings inpFS);
	~SRKGradientField(){}
	void addFieldValue(const double globalPoint[4], double fieldValue[9]);
};
#endif
