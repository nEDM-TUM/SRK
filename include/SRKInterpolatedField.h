#ifndef SRKINTERPOLATEDFIELD_HH_
#define SRKINTERPOLATEDFIELD_HH_

#include "SRKField.h"

#include "SRKROOTField.h"
#include "SRKField.h"
#include "TVector3.h"

////////////////////////////////////////////////////////////////
/// class SRKInterpolatedField
///
/// A field stored in a ROOT file that interpolation is performed on.
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

class SRKInterpolatedField: public SRKField
{
public:
	SRKInterpolatedField(SRKFieldSettings fs);
	~SRKInterpolatedField(){}

	void addFieldValue(const double globalPoint[4], double fieldValue[9]);

private:
	SRKROOTField theROOTField;
	TVector3 offset;
};

#endif /* SRKINTERPOLATEDFIELD_HH_ */
