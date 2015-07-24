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

#include <time.h>
#include <iostream>

#include "SRKGlobalField.h"

using namespace std;

SRKGlobalField::SRKGlobalField() :
	SRKField()
{
	theFieldList = new FieldList();

	fs.fieldType = FIELD_MAGNETIC; //Highest field type used;

	currentFieldSettingsToModify = -1;
}

SRKGlobalField::~SRKGlobalField()
{
	clear();
}

void SRKGlobalField::updateField()
{
	first = true;

	numFields = 0;
	theFieldArray = 0;

	clear();
	currentFieldSettingsToModify = -1;

	constructFields();

}


void SRKGlobalField::addFieldValue(const double* point, double* outField) const
{

	// protect against Geant4 bug that calls us with point[] NaN.
	if(point[0] != point[0]) return;

	// (can't use nfp or outField, as they may change)
	if(first) ((SRKGlobalField*) this)->setupArray();   // (cast away const)

	for (int i = 0; i < numFields; ++i)
	{

		theFieldArray[i]->addFieldValue(point, outField);

	}

}

//clears the field list NOT the fieldSettingsToLoad
void SRKGlobalField::clear()
{
	if(theFieldList)
	{
		if(theFieldList->size() > 0)
		{
			FieldList::iterator i;
			for (i = theFieldList->begin(); i != theFieldList->end(); ++i)
				delete* i;
			theFieldList->clear();
		}
	}

	if(theFieldArray) delete[] theFieldArray;

	first = true;

	numFields = 0;
	theFieldArray = NULL;
}

void SRKGlobalField::setupArray()
{
	first = false;
	numFields = theFieldList->size();
	theFieldArray = new SRKField*[numFields + 1]; // add 1 so it's never 0
	for (int i = 0; i < numFields; ++i)
		theFieldArray[i] = (*theFieldList)[i];
}

void SRKGlobalField::setCurrentFieldSettingsToModify(int inp)
{
	if(inp >= 0 && inp < MAX_SRK_NUM_FIELDS) //Is the number within a reasonable range
	{
		currentFieldSettingsToModify = inp;
		if(fieldSettingsToLoad.size() <= (unsigned int) inp) //Do we need to grow the vector
		{
			fieldSettingsToLoad.resize(inp + 1); //New FieldSettings that are created are not loaded because their scaling value is set to zero by default
		}
	}
	else
	{
		cout << "Field to modify not within appropriate range of 0 to " << MAX_SRK_NUM_FIELDS << "." << endl;
	}
}

void SRKGlobalField::constructFields()
{
	for (unsigned int i = 0; i < fieldSettingsToLoad.size(); i++)
	{
		if(fieldSettingsToLoad[i].scalingValue != 0) //Only load fields if their scaling values are not zero
		{
			if(fieldSettingsToLoad[i].fieldClass == FIELDCLASS_INTERPOLATION)
			{
				addField(new SRKInterpolatedField(fieldSettingsToLoad[i]));
			}
			else if(fieldSettingsToLoad[i].fieldClass == FIELDCLASS_DIPOLE)
			{
				addField(new SRKDipoleField(fieldSettingsToLoad[i]));
			}
			else if(fieldSettingsToLoad[i].fieldClass == FIELDCLASS_UNIFORM)
			{
				addField(new SRKUniformField(fieldSettingsToLoad[i]));
			}
			else if(fieldSettingsToLoad[i].fieldClass == FIELDCLASS_GRADIENT)
			{
				addField(new SRKGradientField(fieldSettingsToLoad[i]));
			}
		}
		if(fieldSettingsToLoad[i].fieldType == FIELD_GRAVITY)
		{
			fs.fieldType = FIELD_GRAVITY; //Highest field type used;
		}
		else if(fieldSettingsToLoad[i].fieldType == FIELD_ELECTRIC)
		{
			fs.fieldType = FIELD_ELECTRIC; //Highest field type used;
		}
	}
	return;
}
