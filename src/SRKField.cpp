#include "SRKField.h"

#include "TMath.h"

#include <iostream>

using namespace std;

SRKField::SRKField()
{
	scalingFactorWUnits = 0;
	g4FieldX = 0;
	g4FieldY = 0;
	g4FieldZ = 0;
}

SRKField::SRKField(SRKFieldSettings inpFS)
{
	fs = inpFS;
//	fs.direction = TVector3(0, 0, 1);
	fs.direction.SetMag(1.0);

	scalingFactorWUnits = fs.scalingValue;

	//presumed fs.scalingValue is the magnetic moment multiplied by \mu_0 in units of T*m^3 or the electric dipole moment divided by \epsilon_0 in units of V*m^2
	//Similarly for gravity
	if(fs.fieldClass == FIELDCLASS_DIPOLE)
	{
		//scalingFactorWUnits *= 1. / (4 * TMath::Pi());
		scalingFactorWUnits *= 1.; //Pignol and Roccia definition
	}
	else if(fs.fieldClass == FIELDCLASS_GRADIENT)
	{
		scalingFactorWUnits *= 1.;
	}

	if(fs.fieldType == FIELD_MAGNETIC)
	{
		g4FieldX = 0;
		g4FieldY = 1;
		g4FieldZ = 2;
	}
	else if(fs.fieldType == FIELD_ELECTRIC)
	{
		g4FieldX = 3;
		g4FieldY = 4;
		g4FieldZ = 5;
	}
	else if(fs.fieldType == FIELD_GRAVITY)
	{
		g4FieldX = 6;
		g4FieldY = 7;
		g4FieldZ = 8;
	}
	else
	{
		cout << "Field type not recognized!" << endl;
	}
	cout << defaultfloat;
	cout << "Loading " << fs.getFieldTypeString() << " " << fs.getFieldClassString() << "Field with scaling value of: " << fs.scalingValue << endl;

}
