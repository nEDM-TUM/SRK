//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//

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

SRKField::SRKField(FieldSettings inpFS)
{
	fs = inpFS;
	fs.direction = TVector3(0, 0, 1);
	fs.direction.SetMag(1.0);

	scalingFactorWUnits = fs.scalingValue;

	//presumed fs.scalingValue is the magnetic moment multiplied by \mu_0 in units of T*m^3 or the electric dipole moment divided by \epsilon_0 in units of V*m^2
	//Similarly for gravity
	if(fs.fieldClass == FIELDCLASS_DIPOLE)
	{
		//scalingFactorWUnits *= 1. / (4 * TMath::Pi());
		scalingFactorWUnits*=1.; //Believe this is how Pignol and Roccia define it.
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
	cout.unsetf(ios_base::floatfield); //Default "adaptive" floating point (cleaner in C++11....)
	cout << "Loading " << fs.getFieldTypeString() << " " << fs.getFieldClassString() << "Field with scaling value of: "  << fs.scalingValue << endl;

}

SRKField::~SRKField()
{

}
