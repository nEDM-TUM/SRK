#ifndef SRKGRAPHICS_H_
#define SRKGRAPHICS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TH1.h"

void convertTH1ToTXT(TH1* inpHist, TString outputFileName);
void convertXYDataWErrorToTXT(int numData, double* x, double* y, double* eX, double* eY, TString titleString, TString outputFileName);
void convertTGraphErrorsToTXT(TGraphErrors* inpGraph, TString outputFileName);
bool fileExists(TString strFileName);
TString fileNameFromFullPath(TString fullPath);
int countTString(TString inpString, char inpChar);
int getNonCommentLine(ifstream& inpFileStream, TString& outString, char delim);
TGraphErrors* getTabSeperatedTGraphErrors(TString filePath, char delim = '\n');

#endif /* SRKGRAPHICS_H_ */
