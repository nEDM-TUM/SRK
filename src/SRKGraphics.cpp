#include "SRKGraphics.h"

using namespace std;

void convertTH1ToTXT(TH1* inpHist, TString outputFileName)
{
	ofstream outFile;
	outFile.open(outputFileName);
	double binLow, data, error;
	outFile << "#Bin Center\tBin Content\tBin Error" << endl;
	for (int i = 1; i <= inpHist->GetNbinsX(); i++)
	{
		outFile.precision(15);
		data = inpHist->GetBinContent(i);
		binLow = inpHist->GetBinCenter(i);
		error = inpHist->GetBinError(i);
		outFile << scientific << binLow << "\t" << scientific << data << "\t" << scientific << error;
		if(i != inpHist->GetNbinsX()) outFile << endl;
	}
	outFile.close();
}

void convertXYDataWErrorToTXT(int numData, double* x, double* y, double* eX, double* eY, TString titleString, TString outputFileName)
{
	ofstream outFile;
	outFile.open(outputFileName);

	outFile << titleString << endl;
	for (int i = 0; i < numData; i++)
	{
		outFile.precision(15);
		outFile << scientific << x[i];
		if(eX != NULL)
		{
			outFile << "\t" << scientific << eX[i];
		}
		else
		{
			outFile << "\t" << scientific << 0;
		}

		outFile << "\t" << scientific << y[i];
		if(eY != NULL)
		{
			outFile << "\t" << scientific << eY[i];
		}
		else
		{
			outFile << "\t" << scientific << 0;
		}

		if(i != numData - 1)
		{
			outFile << endl;
		}

	}
	outFile.close();
}

void convertTGraphErrorsToTXT(TGraphErrors* inpGraph, TString outputFileName)
{
	ofstream outFile;
	outFile.open(outputFileName);

	double* x = inpGraph->GetX();
	double* y = inpGraph->GetY();
	double* eX = inpGraph->GetEX();
	double* eY = inpGraph->GetEY();
	TString titleString = TString("#") + inpGraph->GetTitle() + ";" + inpGraph->GetXaxis()->GetTitle() + ";" + inpGraph->GetYaxis()->GetTitle();
	convertXYDataWErrorToTXT(inpGraph->GetN(), x, y, eX, eY, titleString, outputFileName);

}

int getNonCommentLine(ifstream& inpFileStream, TString& outString, char delim)
{
	outString = "";
	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;

	TString sBuffer("#");

	while (sBuffer[0] == '#' || sBuffer[0] == '%')
	{
		sBuffer.ReadToDelim(inpFileStream, delim);
		if(inpFileStream.fail()) return -1;
	}

	outString = sBuffer;
	return 0;
}

bool fileExists(TString strFileName)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
#ifndef __CINT__
	// Attempt to get the file attributes
	intStat = stat(strFileName.Data(), &stFileInfo);
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
return !(gSystem->GetPathInfo(strFileName,id,size,flags,mt));
}
#endif //__CINT__

TString fileNameFromFullPath(TString fullPath)  //return fullPath if not found
{
	Ssiz_t found;
	found = fullPath.Last('/');
	if(found == kNPOS) return fullPath;
	return fullPath.Remove(0, found + 1);
}

int countTString(TString inpString, char inpChar)
{
	int numCharFound = 0;
	for (int i = 0; i < inpString.Sizeof(); i++)
	{
		if(inpString[i] == inpChar)
		{
			numCharFound++;
		}
	}
	return numCharFound;
}

TGraphErrors* getTabSeperatedTGraphErrors(TString filePath, char delim)
{
	std::vector<double> values[4];  //x,xE,y,yE depending on numColumns;

	if(!fileExists(filePath))
	{
		cout << "File, " << filePath << ", does not exist." << endl;
		return NULL;
	}

	ifstream inpFile;
	inpFile.open(filePath);

	//If comment line is present presume that it is the title
	TString titleString = fileNameFromFullPath(filePath);
	if(inpFile.peek() == '#')
	{
		titleString = "";
		titleString.ReadToDelim(inpFile, delim);
		titleString.Remove(0, 1);
	}

	TString temp, line;
	double value;

	//Get first line to check for column count and then return to beginning;
	getNonCommentLine(inpFile, line, delim);
	const int numColumns = countTString(line, '\t') + 1; //counts how many delims there are to determine number of columns
	inpFile.seekg(0);

	if(numColumns < 2 || numColumns > 4)
	{
		cout << "Data expected in X,Y or X,Y,Yerror   or  X,Xerror,Y,Yerror format!" << endl;
		cout << "Number of columns set, " << numColumns << ", is incorrected." << endl;
		return NULL;
	}

	//Read data
	for (int i = 0; i < 1000000; i++)  //Arbitrary 1 mill. data point max
	{

		if(getNonCommentLine(inpFile, line, delim))
		{
			break;  //Usually end of file
		}

		stringstream theSS(line.Data());

		theSS >> value;

		values[0].push_back(value); //X

		if(numColumns == 2)
		{
			theSS >> value;
			values[2].push_back(value); //Y;
			values[1].push_back(0); //X error is 0
			values[3].push_back(0); //Y error is 0
		}
		else if(numColumns == 3)
		{
			theSS >> value;
			values[2].push_back(value); //Y;
			theSS >> value;
			values[3].push_back(value); //Y error
			values[1].push_back(0); //X error is 0
		}
		else //4 columns
		{
			theSS >> value;
			values[1].push_back(value); //X error
			theSS >> value;
			values[2].push_back(value); //Y;
			theSS >> value;
			values[3].push_back(value); //Y error
		}
	}

	int numData = values[0].size();

	TGraphErrors* outGraph;

	if(numColumns == 2)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], NULL, NULL);
	}
	else if(numColumns == 3)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], NULL, &values[3][0]);
	}
	else if(numColumns == 4)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], &values[1][0], &values[3][0]);
	}

	outGraph->SetName(fileNameFromFullPath(filePath));
	outGraph->SetTitle(titleString);
	cout << "Title: " << titleString << endl;

	return outGraph;
}
