#ifndef SRKGradientField_H
#define SRKGradientField_H 1

#include "SRKField.h"

//Gradient field defined in the z direction (Later will add rotation)
//The field is of the set field strength at centerPos and diverges from there.
//Uses scaleFactor as gradient strength in [FIELDUNITS]/m

class SRKGradientField : public SRKField
{
public:
	SRKGradientField(FieldSettings inpFS);
  ~SRKGradientField(){};
  void  addFieldValue( const double Point[4], double *Bfield ) const;
};
#endif
