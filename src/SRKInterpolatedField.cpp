#include "SRKInterpolatedField.h"
using namespace std;
SRKInterpolatedField::SRKInterpolatedField(SRKFieldSettings inpFS) :
	SRKField(inpFS)
{

	offset.SetXYZ(fs.offset[0], fs.offset[1], fs.offset[2]); //Convert to ROOT vector format
	cout << "Loading field: " << fs.fieldFilePath << "... ";
	if(theROOTField.loadFieldFromFile(fs.fieldFilePath, fs.histName, scalingFactorWUnits, fs.spaceDim, fs.fieldDim))
	{
		cout << "error loading field!" << endl;
	}
	else
	{
		cout << " loading complete." << endl;
	}

	theROOTField.setRotation(fs.angleX, fs.angleY, fs.angleZ);
	theROOTField.setSymmetry(fs.symmetry[0], fs.symmetry[1], fs.symmetry[2]);

}

void SRKInterpolatedField::addFieldValue(const double globalPoint[4], double fieldValue[9])
{
	if(TMath::Abs(globalPoint[0]) > fs.extents[0] || TMath::Abs(globalPoint[1]) > fs.extents[1] || TMath::Abs(globalPoint[2]) > fs.extents[2])
	{

		return; // don't change if outside field boundary (currently centered)
	}
	TVector3 posIn, vecOut;
	posIn.SetXYZ(globalPoint[0], globalPoint[1], globalPoint[2]);
	vecOut.SetXYZ(0, 0, 0);
	if(fs.useCubicInterpolation)
	{
		theROOTField.cubicInterp3D(posIn - offset, vecOut);
	}
	else
	{
		theROOTField.linearInterp3D(posIn - offset, vecOut);
	}
	fieldValue[g4FieldX] += vecOut.x();
	fieldValue[g4FieldY] += vecOut.y();
	fieldValue[g4FieldZ] += vecOut.z();
}
