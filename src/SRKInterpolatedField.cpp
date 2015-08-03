#include "SRKInterpolatedField.h"
using namespace std;
SRKInterpolatedField::SRKInterpolatedField(FieldSettings inpFS) :
	SRKField(inpFS)
{

	theROOTField = new SRKROOTField();
	offset.SetXYZ(fs.offset[0], fs.offset[1], fs.offset[2]); //Convert to ROOT vector format
	cout << "Loading field: " << fs.fieldFilePath << "... ";
	if(theROOTField->loadFieldFromFile(fs.fieldFilePath, fs.histName, scalingFactorWUnits, fs.spaceDim, fs.fieldDim)) cout << "error loading field!" << endl;
	cout << " loading complete." << endl;

	theROOTField->setRotation(fs.angleX, fs.angleY, fs.angleZ);
	theROOTField->setSymmetry(fs.symmetry[0], fs.symmetry[1], fs.symmetry[2]);

}

SRKInterpolatedField::~SRKInterpolatedField()
{
	// TODO Auto-generated destructor stub
}

void SRKInterpolatedField::addFieldValue(const double Point[4], double* BField) const
{
	if(TMath::Abs(Point[0]) > fs.extents[0] || TMath::Abs(Point[1]) > fs.extents[1] || TMath::Abs(Point[2]) > fs.extents[2]) //TODO: Need to decide whether to add this
	{

		return; // don't change if outside field boundary (currently centered)
	}
	TVector3 posIn, vecOut;
	posIn.SetXYZ(Point[0], Point[1], Point[2]);
	vecOut.SetXYZ(0, 0, 0);
	if(fs.useCubicInterpolation)
	{
		theROOTField->cubicInterp3D(posIn - offset, vecOut);
	}
	else
	{
		theROOTField->linearInterp3D(posIn - offset, vecOut);
	}
	BField[g4FieldX] += vecOut.x();
	BField[g4FieldY] += vecOut.y();
	BField[g4FieldZ] += vecOut.z();
}
