#include "SRKGraphics.h"
#include <algorithm>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include "TF1.h"
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
		if(eX != nullptr)
		{
			outFile << "\t" << scientific << eX[i];
		}
		else
		{
			outFile << "\t" << scientific << 0;
		}

		outFile << "\t" << scientific << y[i];
		if(eY != nullptr)
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
		return nullptr;
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
		return nullptr;
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
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], nullptr, nullptr);
	}
	else if(numColumns == 3)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], nullptr, &values[3][0]);
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

double makeMeanPhasePlot(TString filePath, TString imagePath, bool useWrapping, double& errorOut, double& stdevOut)
{
	TFile* rootFile = new TFile(filePath);

	double phi;
	TTree* theTree = (TTree*) rootFile->Get("hitTree");
	theTree->SetMakeClass(1);
	theTree->SetBranchAddress("phi", &(phi));

	gROOT->cd(0);

	const int numEntries = theTree->GetEntries();

	vector<double> phiVec;
	//theTree->GetEntry(0, 1);
	for (int i = 0; i < numEntries; i++)
	{

		theTree->GetEntry(i, 1);

		phiVec.push_back(phi);

	}
	rootFile->Close();
	delete rootFile;

	double mean;
	if(useWrapping)
	{
		mean = reducePeriodicToMeanInVector(phiVec);
	}
	else
	{
		mean = carefulMeanVector(phiVec);
	}

	stdevOut = carefullStDevVector(phiVec, true);
	errorOut = stdevOut / sqrt((double) phiVec.size());

	if(imagePath != "")
	{
		TH1* phaseHist = new TH1D("phaseHist", "Phase Hist; Final phase (radians);Prob (counts/total counts)", 2048, mean - 3 * stdevOut, mean + 3 * stdevOut);

		for (int i = 0; i < numEntries; ++i)
		{
			phaseHist->Fill(phiVec[i]);
		}

		TCanvas c2;
		gStyle->SetOptStat("eMKS");
		phaseHist->SetStats(true);
		phaseHist->Rebin(8);

		c2.SetGrid(1, 1);
		gPad->SetTickx(1);
		gPad->SetTicky(1);
		gPad->SetLogy(1);
		phaseHist->Draw();

		//TsallisFunc
		TF1 f1("TsallisFunc", "[0]/pow(1+((x-[3])/[1])*((x-[3])/[1]),[2])", mean - 3 * stdevOut, mean + 3 * stdevOut);
		f1.SetParNames("Amplitude","Sigma","Power","Mean");
		f1.SetParameters(phaseHist->GetMaximum()*5, stdevOut*0.005,.76,mean);
		f1.SetParLimits(3,mean,mean);

//		TF1 f1("Gaussian", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", mean - 3 * stdev, mean + 3 * stdev);
//		f1.SetParNames("Amplitude", "Sigma", "Mean");
//		f1.SetParameters(100, mean, stdev);
//		f1.SetParLimits(0, 80, 1000);
//		f1.SetParLimits(1, 0, 0.0001);
//		f1.SetParLimits(2, mean-stdev, mean+stdev);

		f1.SetLineColor(kRed);
		phaseHist->Fit(&f1, "VMR+");
		f1.Draw("same");

		c2.SaveAs(imagePath);
		gStyle->SetOptStat("");

		delete phaseHist;
	}
	return mean;
}

double meanVector(const vector<double>& theData)
{
	double mean = 0;
	for (double x : theData)
	{
		mean += x;
	}
	mean /= theData.size();
	return mean;
}

//Worried about precision
double carefulMeanVector(const vector<double>& theData)
{
	//First calc mean to use in subtraction
	double tempMean = meanVector(theData);

	double mean = 0;
	for (double x : theData)
	{
		mean += (x - tempMean);
	}
	mean /= theData.size();
	mean += tempMean;
	return mean;
}

double stDevVector(const vector<double>& theData, const bool useBesselCorrection)
{
	double mean = meanVector(theData);
	double stdev = 0;

	for (double x : theData)
	{
		stdev += pow(x - mean, 2);
	}

	if(useBesselCorrection)
	{
		stdev /= theData.size() - 1;
	}
	else
	{
		stdev /= theData.size();
	}

	stdev = sqrt(stdev);

	return stdev;

}

double carefullStDevVector(const vector<double>& theData, const bool useBesselCorrection)
{
	double mean = carefulMeanVector(theData);
	double stdev = 0;

	for (double x : theData)
	{
		stdev += pow(x - mean, 2);
	}

	if(useBesselCorrection)
	{
		stdev /= theData.size() - 1;
	}
	else
	{
		stdev /= theData.size();
	}

	stdev = sqrt(stdev);

	return stdev;

}

double minVector(const vector<double>& theData)
{
	return *(min_element(begin(theData), end(theData)));
}

double maxVector(const vector<double>& theData)
{
	return *(max_element(begin(theData), end(theData)));
}

double reducePeriodicToMeanInVector(vector<double>& theData)
{
	double tempMean = carefulMeanVector(theData);
	double mean = 0;
	for (double& x : theData)
	{
		x = reducePeriodicNumber(x, -TMath::Pi() + tempMean, TMath::Pi() + tempMean);
		mean += x;
	}
	mean /= theData.size();

	return mean;

}

double reducePeriodicNumber(const double inp, const double start, const double end)
{
	double period = end - start;
	double answer = inp;
	double fractPart, intPart;

	fractPart = modf((inp - start) / period, &intPart);

	if(inp > end)
	{
		answer = fractPart * period + start;
	}
	else if(inp < start)
	{
		answer = fractPart * period + end;
	}
	return answer;
}

void makeSteyerlPlot(TString titleString, vector<TGraphErrors*> theGraphs, vector<TString> legendList, TString imagePath)
{

	TCanvas* theCanvas = new TCanvas("theCanvas", "theCanvas", 1200, 1200);
	theCanvas->cd();
	theCanvas->SetLeftMargin(.09);
	theCanvas->SetRightMargin(.05);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);

	TLegend* theLegend;
	double legendSize = .1 * theGraphs.size();
	if(legendSize < .15) legendSize = .15;
	theLegend = new TLegend(0.55, .9 - legendSize, 0.8705314, 0.8381924);

	for (unsigned int i = 0; i < theGraphs.size(); ++i)
	{
		TGraphErrors* theGraph = theGraphs[i];
		if(i == 0)
		{
			theGraph->SetTitle(titleString);
			theGraph->Draw("AP");
			theGraph->GetYaxis()->SetTitle("#frac{#Delta#omega}{#zeta#eta#omega_{0}}");
			theGraph->GetXaxis()->SetTitle("#omega_{r} / #omega_{0}");
			theGraph->GetXaxis()->SetLimits(0.01, 100);
			theGraph->GetXaxis()->SetRangeUser(0.01, 100);
			theGraph->GetXaxis()->SetNoExponent();
			theGraph->GetYaxis()->SetLimits(-3.4999, 2);
			theGraph->GetYaxis()->SetRangeUser(-3.4999, 2);
			theGraph->SetLineColor(kRed);
			theGraph->SetLineStyle(2);

		}
		else
		{
			theGraph->Draw("PE same");
			theGraph->SetLineColor(kBlue);

		}
		theLegend->AddEntry(theGraph, legendList[i], "L");
	}

	theLegend->Draw();
	theCanvas->SetLogx();
	theCanvas->SaveAs(imagePath);
	delete theCanvas;
}

void saveGraphsToImage(vector<TGraphErrors*> theGraphs, vector<TString> legendList, TString imagePath)
{
	TCanvas* theCanvas = new TCanvas("theCanvas", "theCanvas", 1200, 1200);
//	theGraphs[0]->Draw("ALX");
	theCanvas->cd();
	theCanvas->SetLeftMargin(.13);
	theCanvas->SetRightMargin(.05);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);
	theCanvas->SetLogy();
	theCanvas->SetLogx();

	TLegend* theLegend;
	double legendSize = .1 * theGraphs.size();
	if(legendSize < .15) legendSize = .15;
	theLegend = new TLegend(0.2, .9 - legendSize, 0.4, 0.9);

	for (unsigned int i = 0; i < theGraphs.size(); ++i)
	{
		TGraphErrors* theGraph = theGraphs[i];

		if(i == 0)
		{

			theGraph->SetLineColor(kRed);
			theGraph->SetLineStyle(2);
			theGraph->SetMarkerStyle(4);
			theGraph->GetYaxis()->SetTitleOffset(1.6);
			theGraph->GetYaxis()->SetRangeUser(1e-34, 1e-27);
			theGraph->GetXaxis()->SetLimits(1.58E-011, 1.e-6);
			theGraph->Draw("ALX");
		}
		else
		{
			theGraph->Draw("LX same");
			theGraph->SetLineColor(kBlue);

		}
		theLegend->AddEntry(theGraph, legendList[i], "L");
	}

	theLegend->Draw();

	theCanvas->SaveAs(imagePath);
	delete theCanvas;
}
