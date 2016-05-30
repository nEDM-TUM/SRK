#ifndef SRKDipoleField_H
#define SRKDipoleField_H 1

#include "SRKField.h"
////////////////////////////////////////////////////////////////
/// SRKDipoleField
/// Implements an SRKField which creates a point-like dipole field
/// where the center is at fs.centerPos and it's magnitude is based on
/// scalingValue.  It's direction is defined by fs.direction
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////
class SRKDipoleField: public SRKField
{
public:

	SRKDipoleField(SRKFieldSettings inpFS);
	~SRKDipoleField(){};
	void addFieldValue(const double globalPoint[4], double fieldValue[9]);

private:

};

#endif
