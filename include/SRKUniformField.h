#ifndef SRKUniformField_HH
#define SRKUniformField_HH

#include <TVector3.h>

#include "SRKField.h"

class SRKUniformField: public SRKField
{
public:
	SRKUniformField(SRKFieldSettings inpFS);

	virtual ~SRKUniformField();

	void addFieldValue(const double globalPoint[4], double fieldValue[9]);

private:
	double fieldComponents[3];
};

#endif
