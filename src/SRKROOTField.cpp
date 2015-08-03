#include "SRKROOTField.h"

using namespace std;
SRKROOTField::SRKROOTField()
{
	reset();
}

SRKROOTField::SRKROOTField(string filePath, string histName, double scalingValue, int inpSpaceDim, int inpFieldDim)
{
	reset();
	loadFieldFromFile(filePath, histName, scalingValue, inpSpaceDim, inpFieldDim);
}

SRKROOTField::~SRKROOTField()
{
	reset();
}

void SRKROOTField::reset()
{
	for (int i = 0; i < 3; i++)
	{
		numBins[i] = 1;
		symmetry[i] = false;

	}
	fieldLoaded = false;
	setSize(0, 1, 1);
	startCorner.SetXYZ(0., 0., 0.);
	endCorner.SetXYZ(0., 0., 0.);
	spacing.SetXYZ(0., 0., 0.);
	rotateX = rotateY = rotateZ = 0;
	dimOfSpace = 0;
	dimOfField = 0;

}

void SRKROOTField::setSize(int inpRows, int inpColumns, int inpLayers)
{

	if(fieldLoaded && numBins[0] > 0 && numBins[1] > 0 && numBins[2] > 0)
	{
		for (int k = 0; k < dimOfField; k++)
		{
			for (int i = 0; i < numBins[0]; i++)
			{
				for (int j = 0; j < numBins[1]; j++)
				{
					theField[k][i][j].clear();
				}
				theField[k][i].clear();
			}
			theField[k].clear();
		}
	}

	numBins[0] = 0;
	numBins[1] = numBins[2] = 1;

	if(inpRows > 0 && inpColumns > 0 && inpLayers > 0)
	{
		numBins[0] = inpRows;
		numBins[1] = inpColumns;
		numBins[2] = inpLayers;

		for (int i = 0; i < dimOfField; i++)
		{
			theField[i].resize(numBins[0]);
			for (int j = 0; j < numBins[0]; j++)
			{
				theField[i][j].resize(numBins[1]);
				for (int k = 0; k < numBins[1]; k++)
				{
					theField[i][j][k].resize(numBins[2]);
				}
			}

		}
		fieldLoaded = true;

	}
}

bool SRKROOTField::isPositionInsideField(TVector3 pos)
{
	TVector3 posInField = getPosInField(pos);

	if(dimOfSpace == 3)
	{
		if(posInField.Z() >= startCorner.Z() && posInField.Z() <= endCorner.Z())
		{
			if(posInField.X() >= startCorner.X() && posInField.X() <= endCorner.X())
			{

				if(posInField.Y() >= startCorner.Y() && posInField.Y() <= endCorner.Y())
				{

					return true;
				}

			}

		}
	}
	else
	{
		double radius = sqrt(posInField.x() * posInField.x() + posInField.y() * posInField.y());
		if(posInField.Z() >= startCorner.Y() && posInField.Z() <= startCorner.Y())
		{
			if(radius >= startCorner.X() && radius <= endCorner.X())
			{
				return true;
			}

		}
	}

	return false;
}

int SRKROOTField::loadFieldFromFile(string filePath, string histName, double scalingValue, int inpSpaceDim, int inpFieldDim)
{
	if(fieldLoaded)
	{
		reset();
	}

	if(!fieldFileExists(filePath))
	{
		cout << "Error loading field file: " << filePath << endl;
		return -1;
	}

	if(inpSpaceDim < 2 || inpSpaceDim > 3)
	{
		cout << "Error in loading field, space dimension not 2 or 3" << endl;
		return -1;
	}

	cout << "Loading field from file: " << filePath << " with scale factor: " << scalingValue << endl;

	dimOfSpace = inpSpaceDim;
	dimOfField = inpFieldDim;

	//Determine file extension.
	string extension;
	size_t found;
	found = filePath.find_last_of(".");
	if(found == string::npos)
	{
		extension = "";
	}
	else
	{
		extension = filePath.substr(found + 1);
	}

	if(extension == "root")
	{
		return loadFieldROOTFile(filePath, histName, scalingValue);
	}

	return loadFieldTXTFile(filePath, scalingValue);

}

//Assumes file format is:
//Normal txt file: # comment lines
//X comment lines
//0-3 dimension lines
//X comment lines
//N field definition lines in the iteration order Z Y X
//Assumes field is stored as cm for text file
//Stored in class as meters
int SRKROOTField::loadFieldTXTFile(string fieldPath, double scalingValue)
{

	ifstream fieldFile;

	fieldFile.open(fieldPath.data());

	if(!fieldFile.is_open())
	{
		cout << "Error opening: " << fieldPath << endl;
		return 1;
	}

	bool isComsolFileFormat = (fieldFile.peek() == '%');

	//Skip comment lines
	string sBuffer;
	while (fieldFile.peek() == '#' || fieldFile.peek() == '%')
	{
		getline(fieldFile, sBuffer);
	}

	double tempDouble; //Temp Doubles

	int numLines = 1;
	//Read grid spacing and start location information from file.  For COMSOL file intepretation
	for (int i = 0; i < dimOfSpace; i++)
	{

		fieldFile >> tempDouble;
		startCorner[i] = tempDouble;
		fieldFile >> tempDouble >> numBins[i];
		spacing[i] = tempDouble;

		fieldFile.ignore(1024, '\n');

		length[i] = spacing[i] * (numBins[i] - 1);

		endCorner[i] = startCorner[i] + length[i];
		numLines *= numBins[i];

		if(numBins[i] <= 0)
		{
			cout << "Incorrect number of bins read from file: " << fieldPath << endl;
			fieldFile.close();
			reset();
			return -1;
		}
	}

	//Convert meters to cm
	startCorner *= 0.01;
	spacing *= 0.01;
	endCorner *= 0.01;
	length *= 0.01;

	if(isComsolFileFormat) //deal with reversed X direction in field model (RDK only)
	{
		tempDouble = -startCorner[0];
		startCorner[0] = -endCorner[0];
		endCorner[0] = tempDouble;
	}

	cout << "startCorner";
	startCorner.Print();
	cout << "spacing";
	spacing.Print();
	cout << "endCorner";
	endCorner.Print();
	cout << "length";
	length.Print();

	setSize(numBins[0], numBins[1], numBins[2]);

	//Skip comment lines
	while (fieldFile.peek() == '#' || fieldFile.peek() == '%')
	{
		getline(fieldFile, sBuffer);
	}

	int numDoublesToSkip;
	string tempString;
	for (int j = 0; j < numLines; j++)
	{

		numDoublesToSkip = dimOfSpace; //In my text files I stored the position
		if(isComsolFileFormat)
		{
			numDoublesToSkip = 4; //Comsol stores all three positions and voltage regardless of dimension
		}

		for (int i = 0; i < numDoublesToSkip; i++)
		{
			fieldFile >> tempDouble; //There are stored positions in text file....these do nothing but they're there so we skip them
		}

		int x = 0, y = 0, z = 0;
		if(isComsolFileFormat)
		{
			//x= j%numBins[0];
			x = (numBins[0] - 1) - j % numBins[0];  //TEMPORARY CHANGE SOON
			if(numBins[1])
			{
				y = (j / numBins[0]) % numBins[1];
			}
			if(numBins[2])
			{
				z = (j / (numBins[0] * numBins[1])) % numBins[2];
			}
		}
		else
		{
			if(dimOfSpace == 3)
			{
				z = j % numBins[2];
				y = (j / numBins[2]) % numBins[1];
				x = (j / (numBins[2] * numBins[1])) % numBins[0];
			}
			else
			{
				y = j % numBins[1];
				x = (j / numBins[1]) % numBins[0];
			}

		}

		for (int i = 0; i < dimOfField; i++)
		{
			fieldFile >> tempString;

			if(tempString == "NaN")
			{
				tempDouble = 0;
			}
			else
			{
				std::istringstream convStream(tempString);
				convStream >> tempDouble;
			}

			if(!fieldFile.is_open() || fieldFile.fail() || fieldFile.eof())
			{
				cout << "Data not found!" << endl;
				return -1;
			}
			if(i == 0 && isComsolFileFormat && dimOfField == 3) //reversed x direction accidently...need to figure out what to do for the 2D case later
			{
				tempDouble = -tempDouble;
			}
			theField[i][x][y][z] = tempDouble * scalingValue;
		}

	}
	fieldLoaded = true;
	fieldFile.close();
	return 0;
}

//Stored in meters
int SRKROOTField::loadFieldROOTFile(string filePath, string histName, double scalingValue)
{
	TFile inpFile(filePath.data(), "READ");
	TH1* histArray[3] =
	{ NULL, NULL, NULL };
	string specificHistName;
	for (int i = 0; i < dimOfField; i++)
	{
		std::ostringstream convStream;
		convStream << i;
		specificHistName = histName + convStream.str();
		if(dimOfField == 2)
		{
			histArray[i] = (TH2D*) inpFile.Get(specificHistName.data());
		}
		else
		{
			histArray[i] = (TH3D*) inpFile.Get(specificHistName.data());

		}

	}

	numBins[0] = histArray[0]->GetNbinsX();
	numBins[1] = histArray[0]->GetNbinsY();

	startCorner.SetX(histArray[0]->GetXaxis()->GetBinCenter(1));
	startCorner.SetY(histArray[0]->GetYaxis()->GetBinCenter(1));

	spacing.SetX(histArray[0]->GetXaxis()->GetBinWidth(1));
	spacing.SetY(histArray[0]->GetYaxis()->GetBinWidth(1));
	endCorner.SetX(startCorner.X() + numBins[0] * spacing.X());
	endCorner.SetY(startCorner.Y() + numBins[1] * spacing.Y());
	if(dimOfSpace == 3)
	{

		spacing.SetZ(histArray[0]->GetZaxis()->GetBinWidth(1));
		startCorner.SetZ(histArray[0]->GetZaxis()->GetBinCenter(1));
		numBins[2] = histArray[0]->GetNbinsZ();
		endCorner.SetZ(startCorner.Z() + numBins[2] * spacing.Z());
	}

	length = endCorner - startCorner;
	setSize(numBins[0], numBins[1], numBins[2]);

	for (int i = 0; i < dimOfField; i++)
	{
		for (int x = 0; x < numBins[0]; x++)
		{
			for (int y = 0; y < numBins[1]; y++)
			{
				for (int z = 0; z < numBins[2]; z++)
				{

					if(dimOfSpace == 2)
					{
						theField[i][x][y][z] = histArray[i]->GetBinContent(x + 1, y + 1) * scalingValue;
					}
					else
					{
						theField[i][x][y][z] = histArray[i]->GetBinContent(x + 1, y + 1, z + 1) * scalingValue;
					}

				}

			}

		};

	}

	fieldLoaded = true;
	inpFile.Close();
	return 0;
}

TVector3 SRKROOTField::getPosInField(TVector3 inp)
{
	TVector3 out = inp;

	out.RotateX(rotateX);
	out.RotateY(rotateY);
	out.RotateZ(rotateZ);

	for (int i = 0; i < 3; i++)
	{
		if(symmetry[i])
		{
			out[i] = abs(inp[i]);
		}

	}

	return out;
}

//Adds interpolated field to vecOut
void SRKROOTField::linearInterp3D(TVector3 posInField, TVector3& vecOut) const
{
	double radius;
	TVector3 resultVec;
	double modx, mody, modz, modxm, modym, modzm;
	double dx, dy, dz;
	int x, y, z, xp1, yp1, zp1;
	double a1, a2, a3, a4, a5, a6, a7, a8;
//	double rotateX,rotateY,rotateZ;  //Rotate angle theta aroudn each axis in order to transform coordinates.  Transformed in order X,Y,Z.

	if(!fieldLoaded)
	{
		cout << "Error field not loaded" << endl;
		return;
	}

	//Rotate coordinates of sampling location
	posInField.RotateX(rotateX);
	posInField.RotateY(rotateY);
	posInField.RotateZ(rotateZ);

//    posInField.Print();
//    startCorner.Print();
//    spacing.Print();

	if(dimOfSpace == 3)
	{
//        if(posInField.Z() < startCorner.Z() || posInField.Z() > endCorner.Z() || posInField.X() < startCorner.X() || posInField.X() > endCorner.X() || posInField.Y() < startCorner.Y() || posInField.Y() > endCorner.Y())
//        {
//            return;
//        }
		if(symmetry[0])
		{
			modx = modf((abs(posInField[0]) - startCorner[0]) / spacing[0], &(dx));
		}
		else
		{
			modx = modf((posInField[0] - startCorner[0]) / spacing[0], &(dx));
		}
		if(symmetry[1])
		{
			mody = modf((abs(posInField[1]) - startCorner[1]) / spacing[1], &(dy));
		}
		else
		{
			mody = modf((posInField[1] - startCorner[1]) / spacing[1], &(dy));
		}
		if(symmetry[2])
		{
			modz = modf((abs(posInField[2]) - startCorner[2]) / spacing[2], &(dz));
		}
		else
		{
			modz = modf((posInField[2] - startCorner[2]) / spacing[2], &(dz));
		}

		modxm = 1. - modx;
		modym = 1. - mody;
		modzm = 1. - modz;

		x = int(dx);
		y = int(dy);
		z = int(dz);
		xp1 = x + 1;
		yp1 = y + 1;
		zp1 = z + 1;

		a1 = (modxm) * (modym) * (modzm);
		a2 = (modxm) * (modym) * (modz);
		a3 = (modxm) * (mody) * (modzm);
		a4 = (modxm) * (mody) * (modz);
		a5 = (modx) * (modym) * (modzm);
		a6 = (modx) * (modym) * (modz);
		a7 = (modx) * (mody) * (modzm);
		a8 = (modx) * (mody) * (modz);

		if(x >= numBins[0] - 1 || x < 0 || y >= numBins[1] - 1 || y < 0 || z >= numBins[2] - 1 || z < 0)
		{
			return;
		}

		for (int i = 0; i < dimOfField; i++)
		{

			resultVec[i] = theField[i][x][y][z] * a1 + theField[i][x][y][zp1] * a2 + theField[i][x][yp1][z] * a3 + theField[i][x][yp1][zp1] * a4 + theField[i][xp1][y][z] * a5 + theField[i][xp1][y][zp1] * a6 + theField[i][xp1][yp1][z] * a7 + theField[i][xp1][yp1][zp1] * a8;

			if(dimOfField > 1)
			{
				if(symmetry[i] && posInField[i] < 0)
				{
					resultVec[i] *= -1.;
				}
			}
		}

	}

	else
	{

		radius = sqrt(posInField.x() * posInField.x() + posInField.y() * posInField.y());
//        cout << "New:" << endl;
//        cout << "PosInField " << posInField.x() << " "<< posInField.y() << " "<< posInField.z() << endl;
//        cout << "Radius : " << radius << " zPos: " <<  posInField[2] << endl ;//<< "Start Corner: ";
//        startCorner.Print();
//        spacing.Print();

		modx = modf((radius - startCorner[0]) / spacing[0], &(dx));
		mody = modf((posInField[2] - startCorner[1]) / spacing[1], &(dy));

		modxm = 1. - modx;
		modym = 1. - mody;

		x = int(dx);
		y = int(dy);
		xp1 = x + 1;
		yp1 = y + 1;

		if(x >= numBins[0] - 1 || x < 0 || y >= numBins[1] - 1 || y < 0)
		{
			return;
		}

		for (int i = 0; i < dimOfField; i++)
		{

			resultVec[i] = theField[i][x][y][0] * (modxm) * (modym) + theField[i][x][yp1][0] * (modxm) * (mody) + theField[i][xp1][y][0] * (modx) * (modym) + theField[i][xp1][yp1][0] * (modx) * (mody);

		}

		if(dimOfField == 2)
		{
			resultVec.SetZ(resultVec.Y());
			if(radius == 0)
			{
				resultVec.SetX(0);
				resultVec.SetY(0);
			}
			else
			{
				resultVec.SetY((resultVec.X() * posInField.Y()) / radius);
				resultVec.SetX((resultVec.X() * posInField.X()) / radius);

			}
		}

		resultVec.RotateZ(-rotateZ);
		resultVec.RotateY(-rotateY);
		resultVec.RotateX(-rotateX);

	}

	vecOut += resultVec;

	return;
}

TVector3 SRKROOTField::cubicInterpolate(TVector3 p[4], double xIn) const
{
	return p[1] + 0.5 * xIn * (p[2] - p[0] + xIn * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + xIn * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

TVector3 SRKROOTField::bicubicInterpolate(TVector3 p[4][4], double xIn, double yIn) const
{
	TVector3 arr[4];
	arr[0] = cubicInterpolate(p[0], yIn);
	arr[1] = cubicInterpolate(p[1], yIn);
	arr[2] = cubicInterpolate(p[2], yIn);
	arr[3] = cubicInterpolate(p[3], yIn);
	return cubicInterpolate(arr, xIn);
}

TVector3 SRKROOTField::tricubicInterpolate(TVector3 p[4][4][4], double xIn, double yIn, double zIn) const
{
	TVector3 arr[4];
	arr[0] = bicubicInterpolate(p[0], yIn, zIn);
	arr[1] = bicubicInterpolate(p[1], yIn, zIn);
	arr[2] = bicubicInterpolate(p[2], yIn, zIn);
	arr[3] = bicubicInterpolate(p[3], yIn, zIn);
	return cubicInterpolate(arr, xIn);
}

void SRKROOTField::cubicInterp3D(TVector3 posInField, TVector3& vecOut) const
{
	double modx, mody, modz;
	double dx, dy, dz;
	int x, y, z;
	TVector3 resultVec;

	if(!fieldLoaded)
	{
		cout << "Error field not loaded" << endl;
		return;
	}

	if(dimOfSpace != 3)
	{
		cout << "Error: 2D/cylindrical cubic interpolation not implemented" << endl;
	}

	//Rotate coordinates of sampling location
	posInField.RotateX(rotateX);
	posInField.RotateY(rotateY);
	posInField.RotateZ(rotateZ);

	if(symmetry[0])
	{
		modx = modf((abs(posInField[0]) - startCorner[0]) / spacing[0], &(dx));
	}
	else
	{
		modx = modf((posInField[0] - startCorner[0]) / spacing[0], &(dx));
	}
	if(symmetry[1])
	{
		mody = modf((abs(posInField[1]) - startCorner[1]) / spacing[1], &(dy));
	}
	else
	{
		mody = modf((posInField[1] - startCorner[1]) / spacing[1], &(dy));
	}
	if(symmetry[2])
	{
		modz = modf((abs(posInField[2]) - startCorner[2]) / spacing[2], &(dz));
	}
	else
	{
		modz = modf((posInField[2] - startCorner[2]) / spacing[2], &(dz));
	}

	x = int(dx);
	y = int(dy);
	z = int(dz);

	//Presume zero field if outside field boundaries.  I.E. Add nothing to result vec
	if(x >= numBins[0] - 1 || x < 0 || y >= numBins[1] - 1 || y < 0 || z >= numBins[2] - 1 || z < 0)
	{
		return;
	}

	TVector3 p[4][4][4]; //Data array around point.

	//Fill data array around point.  Assume that
	int xP, yP, zP;
	for (int i = 0; i < 4; i++)
	{
		xP = x + i - 1;
		if(xP < 0) xP = 0;
		if(xP > numBins[0] - 1) xP = numBins[0] - 1;

		for (int j = 0; j < 4; j++)
		{
			yP = y + j - 1;
			if(yP < 0) yP = 0;
			if(yP > numBins[1] - 1) yP = numBins[1] - 1;

			for (int k = 0; k < 4; k++)
			{
				zP = z + k - 1;
				if(zP < 0) zP = 0;
				if(zP > numBins[2] - 1) zP = numBins[2] - 1;
				p[i][j][k] = getVector(xP, yP, zP);
			}
		}

	}

	resultVec = tricubicInterpolate(p, modx, mody, modz);

	//If using symmetry, we must mirror the vectors appropriately.
	for (int i = 0; i < dimOfField; i++)
	{

		if(dimOfField > 1)
		{
			if(symmetry[i] && posInField[i] < 0)
			{
				resultVec[i] *= -1.;
			}
		}
	}

	vecOut += resultVec;

	return;
}

int SRKROOTField::saveFieldToFile(string filePath, string histName)
{
	if(!fieldLoaded)
	{
		cout << "Can't save, no field loaded." << endl;
		return -1;
	}

	TH1* histArray[3];
	TFile fieldFileOut(filePath.data(), "RECREATE");
	for (int i = 0; i < dimOfField; i++)
	{
		std::ostringstream convStream;
		convStream << i;
		string nameString = histName + convStream.str();
		if(dimOfSpace == 3)
		{
			histArray[i] = new TH3D(nameString.data(), filePath.data(), numBins[0], startCorner[0] - spacing[0] * .5, endCorner[0] + spacing[0] * .5, numBins[1], startCorner[1] - spacing[1] * .5, endCorner[1] + spacing[1] * .5, numBins[2], startCorner[2] - spacing[2] * .5, endCorner[2] + spacing[2] * .5);
		}
		else if(dimOfSpace == 2)
		{
			histArray[i] = new TH2D(nameString.data(), filePath.data(), numBins[0], startCorner[0] - spacing[0] * .5, endCorner[0] + spacing[0] * .5, numBins[1], startCorner[1] - spacing[1] * .5, endCorner[1] + spacing[1] * .5);
		}
		else
		{
			cout << "Fields only set for 2D and 3D" << endl;
			fieldFileOut.Close();
			reset();
			return -1;
		}

		for (int x = 0; x < numBins[0]; x++)
		{
			for (int y = 0; y < numBins[1]; y++)
			{
				for (int z = 0; z < numBins[2]; z++)
				{
					if(dimOfSpace == 2)
					{
						histArray[i]->SetBinContent(x + 1, y + 1, theField[i][x][y][z]);
					}
					else
					{
						histArray[i]->SetBinContent(x + 1, y + 1, z + 1, theField[i][x][y][z]);
					}

				}
			}
		}

		histArray[i]->Write("", TObject::kOverwrite);
	}

	fieldFileOut.Close();

	return 0;
}

bool SRKROOTField::fieldFileExists(string strFileName)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
#ifndef __CINT__
	// Attempt to get the file attributes
	intStat = stat(strFileName.c_str(), &stFileInfo);
	if(intStat == 0)
	{
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	}
	else
	{
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}

	return (blnReturn);
}
#else //__CINT__
Long_t *id,*size,*flags,*mt;
return !(gSystem->GetPathInfo(strFileName.data(),id,size,flags,mt));
}
#endif //__CINT__

