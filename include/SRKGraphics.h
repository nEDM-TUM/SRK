#ifndef SRKGRAPHICS_H_
#define SRKGRAPHICS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <vector>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TH1.h"
#include "TROOT.h"
#include "TCanvas.h"


void convertTH1ToTXT(TH1* inpHist, TString outputFileName);
void convertXYDataWErrorToTXT(int numData, double* x, double* y, double* eX, double* eY, TString titleString, TString outputFileName);
void convertTGraphErrorsToTXT(TGraphErrors* inpGraph, TString outputFileName);
bool fileExists(TString strFileName);
TString fileNameFromFullPath(TString fullPath);
int countTString(TString inpString, char inpChar);
int getNonCommentLine(ifstream& inpFileStream, TString& outString, char delim);
TGraphErrors* getTabSeperatedTGraphErrors(TString filePath, char delim = '\n');
double makeMeanPhasePlot(TString filePath, TString imagePath, bool useWrapping,double& errorOut);
double meanVector(const std::vector<double>& theData);
double carefulMeanVector(const std::vector<double>& theData);
double stDevVector(const std::vector<double>& theData,const bool useBesselCorrection);
double carefullStDevVector(const std::vector<double>& theData, const bool useBesselCorrection);
double minVector(const std::vector<double>& theData);
double maxVector(const std::vector<double>& theData);
double reducePeriodicToMeanInVector(std::vector<double>& theData);
double reducePeriodicNumber(const double inp, const double start, double end);
void makeSteyerlPlot(TString titleString, std::vector<TGraphErrors*> theGraphs, std::vector<TString> legendList,TString imagePath);
void saveGraphsToImage(std::vector<TGraphErrors*> theGraphs, std::vector<TString> legendList, TString imagePath);
#endif /* SRKGRAPHICS_H_ */
