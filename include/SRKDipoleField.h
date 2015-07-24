#ifndef SRKDipoleField_H
#define SRKDipoleField_H 1


#include "SRKField.h"

//SRKDipoleField creates a dipole field using the field settings moment. Rotations and symmetry from the field setting are not used.
class SRKDipoleField : public SRKField
{
public:

	SRKDipoleField(FieldSettings inpFS);
	~SRKDipoleField(){};
	void  addFieldValue(const double Point[4], double* Bfield) const;

private:

};

#endif
