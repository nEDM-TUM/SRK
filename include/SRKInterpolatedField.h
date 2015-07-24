#ifndef SRKINTERPOLATEDFIELD_HH_
#define SRKINTERPOLATEDFIELD_HH_

#include "SRKField.h"

#include "SRKROOTField.h"
#include "SRKField.h"
#include "TVector3.h"

//Interpolated field using a 2D or 3D field from a ROOT File using UCNROOTField

class SRKInterpolatedField: public SRKField
{
public:
	SRKInterpolatedField(FieldSettings fs);
	virtual ~SRKInterpolatedField();

	void addFieldValue(const double Point[4],
            double* Bfield ) const;
private:
	SRKROOTField* theROOTField;
	TVector3 offset;
};

#endif /* SRKINTERPOLATEDFIELD_HH_ */
