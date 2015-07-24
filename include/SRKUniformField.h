#ifndef SRKUniformField_HH
#define SRKUniformField_HH

#include <TVector3.h>

#include "SRKField.h"

class SRKUniformField : public SRKField
{
public:
	SRKUniformField(FieldSettings inpFS);

	virtual ~SRKUniformField();

	void addFieldValue(const double pos[4], double *field) const;

private:
	double fieldComponents[3];
};

#endif
