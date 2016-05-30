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
// $Id: G4UniformElectricField.cc,v 1.13 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 
//
// Class for creation of uniform Electric Field
//
// 30.1.97 V.Grichine
//
// -------------------------------------------------------------------

#include "SRKUniformField.h"

SRKUniformField::SRKUniformField(SRKFieldSettings inpFS) :
	SRKField(inpFS)
{
	fieldComponents[0] = scalingFactorWUnits * inpFS.direction.x();
	fieldComponents[1] = scalingFactorWUnits * inpFS.direction.y();
	fieldComponents[2] = scalingFactorWUnits * inpFS.direction.z();
}

SRKUniformField::~SRKUniformField()
{
}

// ------------------------------------------------------------------------

void SRKUniformField::addFieldValue(const double globalPoint[4], double fieldValue[9])
{
	fieldValue[g4FieldX] += fieldComponents[0];
	fieldValue[g4FieldY] += fieldComponents[1];
	fieldValue[g4FieldZ] += fieldComponents[2];
}
