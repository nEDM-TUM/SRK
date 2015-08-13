#include <iostream>
#include <time.h>
#include <functional>
#include <string>
#include <sstream>

#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TParameter.h"
#include "TFile.h"
#include "TROOT.h"
#include "TObjString.h"

#include "SRKManager.h"
#include "SRKMacroManager.h"
#include "SRKMotionTracker.h"
#include "SRKGraphics.h"
#include "TMap.h"

using namespace std;

int main(int argc, char* argv[])
{
	delete gRandom;
	gRandom = new TRandom3(0);

	/*Testing*/
//	delete gRandom;
//	gRandom=new TRandom3(0);
//
//	SRKMotionTracker theMotionTracker;
//	theMotionTracker.setReflectionLimit(1000);
//	//theMotionTracker.setTimeLimit(10);
//	theMotionTracker.setMeanVel(193);
//	theMotionTracker.setDiffuseReflectionProb(1);
//	theMotionTracker.makeTracks(1000);
//	theMotionTracker.writeTrackToFile("/home/mjbales/work/code/testproj/output/trackTest.root");
//
//	SRKManager theManager;
//	//theManager.setRecordAllSteps(true);
////	theManager.trackSpins(1000, "/home/mjbales/work/code/testproj/output/trackTest.root", "/home/mjbales/work/code/testproj/output/testRecordAll.root",true);
////	theManager.trackSpins(1000, "/home/mjbales/work/code/testproj/output/trackTest.root", "/home/mjbales/work/code/testproj/output/spinTestP.root",true);
////	theManager.trackSpins(1000, "/home/mjbales/work/code/testproj/output/trackTest.root", "/home/mjbales/work/code/testproj/output/spinTestA.root",false);
//
//	//Test Alternative tracking
//
//	theManager.trackSpins(1000, "/home/mjbales/work/code/testproj/output/trackTest.root", "/home/mjbales/work/code/testproj/output/spinTestP.root");
	/* Simpler Omega Styerl plot command =D */
//	SRKManager theManager;
//	theManager.setUseAltStepping(false);
//	theManager.setConstStepper(false);
//	theManager.trackSpins(100, "", SRKRESULTSDIR+"test.root");
	//theManager.setE0FieldStrength(0);
	//theManager.setBGradFieldStrength(0);
//	theManager.setTimeLimit(10);
//	theManager.setMeanVel(7.26);
//	vector< double > deltaMeanOmega;
//	vector< double > processingTime;
//	for(int i=0;i<50;i++)
//	{
//		clock_t t1,t2;
//		t1=clock();
//		double theEPS=1.0e-6*pow(10.,-i*.1);
//		theManager.setPerStepError(theEPS,theEPS);
//		cout << theEPS << endl;
//
////		double initStep=pow(10,-i*.1);
////		cout << "Init Step: " << initStep << endl;
////		theManager.setInitialStepSize(initStep);
//		theManager.trackSpins(1000, SRKRESULTSDIR+"diffuseHgOmega1.root", SRKRESULTSDIR+"diffuseHgResult.root");
//		deltaMeanOmega.push_back(theManager.getMeanOmega());
//		t2=clock();
//		double diff =((double)t2-(double)t1)/ (double)CLOCKS_PER_SEC;
//		processingTime.push_back(diff);
	//cout << scientific << "EPS: " << theEPS << "        Total Error: " << theManager.getStepsTaken()*theEPS/1000. << endl;
//	}
//	for(unsigned int i=0; i<deltaMeanOmega.size(); i++)
//	{
//		double theEPS=1.0e-6*pow(10.,-i*.1);
//		cout << setprecision(13) << theEPS << "\t" << deltaMeanOmega[deltaMeanOmega.size()-1] - deltaMeanOmega[i] << "\t" << processingTime[i] << endl;
//	}
	/* Steyerl plots*/
//	SRKManager theManager;
//	theManager.setUseAltStepping(false);
//	theManager.setConstStepper(false);
//	TString idString="Steyerl1";
//
//	TGraphErrors* diffuseGraph=theManager.trackSpinsDeltaOmegaSteyerlPlot(10000, idString+"_Diffuse", 100,  0.1,10,true,200);
//	theManager.setDiffuseReflectionProb(0);
//	TGraphErrors* specularGraph=theManager.trackSpinsDeltaOmegaSteyerlPlot(10000, idString+"_Specular", 100,  0.1,10,true,200);
//	vector<TGraphErrors*> theGraphs={specularGraph,diffuseGraph};
//	vector<TString> legendList={"Specular","Diffuse"};
//	makeSteyerlPlot( "Hg 193 m/s #approx200 reflections ILL chamber",  theGraphs, legendList, SRKGRAPHSDIR+idString+".png");
	/* Test dipole ILL setup.*/
	 //	const int numID = 8;
	 //	double dipoleValues[numID] = { 1e-15, 1e-15, 1e-14, 1e-14, 1e-16, 1e-16, 1e-15, 1e-15 };
	 //	double distances[numID] = { 0.05, 0.1, 0.2, 0.3, 0.05, 0.1, 0.2, 0.3 };
	 //	TString theArg = argv[1];
	 //	int idToRun = theArg.Atoi();
	 //	TString temp;
	 //	TString runName = "SRK_ID" + temp.Itoa(idToRun+8,10)+"B";
	 //
	 //	SRKManager theManager;
	 //	theManager.setConstStepper(false);
	 //	theManager.setUseAltStepping(false);
	 //	theManager.setPerStepError(1e-7, 1e-7);
	 //	//theManager.setInitialStepSize(3e-3);
	 //	theManager.setUse2D(false);
	 //	theManager.setChamberHeight(.12);
	 //	theManager.setChamberRadius(0.235);
	 //	theManager.setDiffuseReflectionProb(1);
	 //	theManager.setDipolePosition(TVector3(0, 0, -distances[idToRun] - 0.05 * 0.12));
	 //	theManager.setE0FieldStrength(1.e7);
	 //	theManager.setTimeLimit(100);
	 //	theManager.setMeanVel(193);
	 //	theManager.setDipoleFieldStrength(dipoleValues[idToRun]);
	 //	theManager.setBGradFieldStrength(0);
	 //
	 //	double falseEDMError;
	 //	double falseEDM = theManager.trackSpinsFalseEDM(1000, runName, falseEDMError);
	 //	TString resultPath = SRKHISTSDIR + runName + ".txt";
	 //	ofstream resultFile(resultPath.Data());
	 //	resultFile << scientific << setprecision(5) << falseEDM << "\t" << falseEDMError << endl;
	 //	resultFile.close();
	 //
	 //	double theError;
	 //	makeMeanPhasePlot(SRKRESULTSDIR + runName + "P.root", SRKGRAPHSDIR + runName + "P_Reduced.pdf",true, theError);
	 //	makeMeanPhasePlot(SRKRESULTSDIR + runName + "P.root", SRKGRAPHSDIR + runName + "P.pdf", false, theError);
	 /*Solve memory issue*/
//	SRKMotionTracker theMotionTracker;
//	theMotionTracker.setTimeLimit(1);
//	theMotionTracker.setUse2D(false);
//	cout << "2D" << endl;
//	theMotionTracker.makeTracks(10,SRKTRACKSDIR+"test.root");
	/* Remake plots for Steyerl Omega */

//	SRKManager theManager;
//	theManager.setDiffuseReflectionProb(1);
//	theManager.setTimeLimit(100);
//	theManager.trackSpinsDeltaOmegaSteyerlPlot(1000, "DiffuseConstantTime", 100,  0.1,1.3475,true,-1);
	/*Playing with lambda maps*/

//

//	std::unordered_map<std::string,std::function<void(std::string)>> commandMap;
//	std::string out="Matthew";
//	commandMap["testCommand"]=[&](std::string last){out += " " + last;};
//	commandMap["testCommand"]("Bales");
//	cout << out << endl;
	/*Make phase hists*/
//	TString idString="Dipole_1en11_0cm_Dynamic_EPS1en10_100s_193m_s_E1e7A";
//	makeMeanPhasePlot(SRKRESULTSDIR+idString+".root", SRKGRAPHSDIR+idString+"_Reduced.pdf", true);
//	makeMeanPhasePlot(SRKRESULTSDIR+idString+".root", SRKGRAPHSDIR+idString+".pdf", false);
	/*plot stepping changes */
//	vector<TGraphErrors*> theGraphs= {getTabSeperatedTGraphErrors(SRKHISTSDIR+"CustomDynamicSteppingSizeVsFalseEDM.txt")};
//	vector<TString> legendList={"CustomDynamic"};
//	saveGraphsToImage(theGraphs, legendList, SRKGRAPHSDIR+"CustomDynamicSteppingSizeVsFalseEDM.png");
	/* perform Peter's straight line path test 150807 */
	//Create track file
//	SRKMotionTracker theMotionTracker;
//	theMotionTracker.setTimeLimit(0.01);
//	theMotionTracker.setUse2D(false);
//	theMotionTracker.setManualTracking(true);
//	theMotionTracker.setChamberRadius(1000);
//	theMotionTracker.setChamberHeight(1000);
//
//
//
//	theMotionTracker.openTrackFile(SRKTRACKSDIR+"StraightLineTracksFromCenter150807.root");
//	for(int i=0;i<999;i++)
//	{
//		theMotionTracker.setMeanVel(200.);
//		theMotionTracker.setVel(TVector3(200.,0,0));
//		theMotionTracker.setPos(TVector3(0,0,1.-i*.001));
//		theMotionTracker.makeTracks(1);
//	}
//	theMotionTracker.writeTrackToFile();
//	theMotionTracker.closeTrackFile();
	//Create simple chambered trackFile
//	SRKMotionTracker theMotionTracker;
//	theMotionTracker.setTimeLimit(1.01);
//	theMotionTracker.setUse2D(true);
//	theMotionTracker.setDiffuseReflectionProb(0);
//	theMotionTracker.setManualTracking(true);
//	theMotionTracker.setChamberRadius(.235);
//	theMotionTracker.setChamberHeight(1000);
//	theMo
//
//
//
//	theMotionTracker.openTrackFile(SRKTRACKSDIR+"StraightInChamber150807.root");
//	for(int i=0;i<999;i++)
//	{
//		theMotionTracker.setMeanVel(200.);
//		theMotionTracker.setVel(TVector3(235.,0,0));
//		//theMotionTracker.setPos(TVector3(.235,0,1.-i*.001));
//		theMotionTracker.setPos(TVector3(.235,0,1.-i*.001));
//		theMotionTracker.makeTracks(1);
//	}
//	theMotionTracker.writeTrackToFile();
//	theMotionTracker.closeTrackFile();
	//Spin Tracking
//	double dipoleFieldStrength = 1e-11;
//	double gradFieldStrength = 0;
//	TString strengthString;
//	int number, exponent;
//	if(dipoleFieldStrength != 0)
//	{
//		exponent = abs(log10(dipoleFieldStrength)) + .999999;
//		number = round(dipoleFieldStrength / pow(10, -exponent));
//
//	}
//	else
//	{
//		exponent = round(abs(log10(gradFieldStrength)) + .999999);
//		number = round(gradFieldStrength / pow(10, -exponent));
//
//	}
//	strengthString = TString::Format("%in%i", number, exponent);
//	TString runString = "StraightInChamber150807";
//	SRKManager theManager;
//	theManager.setConstStepper(false);
//	theManager.setUseAltStepping(false);
//	theManager.setPerStepError(1e-7, 1e-7);
//
//	//theManager.setE0FieldStrength(1.e7);
//	theManager.setB0FieldStrength(1.e-6);
//
//	theManager.setDipolePosition(TVector3(0, 0, 0));
//	theManager.setDipoleFieldStrength(dipoleFieldStrength);
//	theManager.setBGradFieldStrength(gradFieldStrength);
//
////	theManager.setDipoleFieldStrength(0);
////	theManager.setBGradFieldStrength(1e-9);
////	theManager.setE0FieldStrength(0);
//
//	theManager.trackSpins(999, SRKTRACKSDIR + runString + ".root", SRKRESULTSDIR + "Results_" + runString + "_" + strengthString + ".root");
//
//	//Analyze
//	TFile resultsFile(SRKRESULTSDIR + "Results_" + runString + "_" + strengthString + ".root");
//	TTree* hitTree = (TTree*) resultsFile.Get("hitTree");
//
//	TCanvas theCanvas("theCanvas", "theCanvas", 1200, 1200);
//	theCanvas.SetLeftMargin(.2);
//	theCanvas.SetRightMargin(.05);
////	theCanvas.SetLogy(1);
//	gPad->SetTickx(1);
//	gPad->SetTicky(1);
//
//	hitTree->Draw("(phi-(48.4578839927*1.01)):pos.z()", "", "L");
//	TH2F *htemp = (TH2F*) gPad->GetPrimitive("htemp");
//	htemp->GetYaxis()->SetTitleOffset(2.2);
//	if(dipoleFieldStrength != 0)
//		htemp->SetTitle("Dipole: " + TString::Format("%e", dipoleFieldStrength) + " T/m^{3};Distance from Dipole (m);#Delta #phi (radians)");
//	else
//		htemp->SetTitle("Gradient: " + TString::Format("%e", gradFieldStrength) + " T/m;Distance from Center (m);#Delta #phi (radians)");
//
//	theCanvas.SaveAs(SRKGRAPHSDIR + "Plot_SRK_" + runString + "_" + strengthString + ".png");
//
//	hitTree->Draw("theta:pos.z()", "", "L");
//	htemp = (TH2F*) gPad->GetPrimitive("htemp");
//	htemp->GetYaxis()->SetTitleOffset(2.2);
//	if(dipoleFieldStrength != 0)
//		htemp->SetTitle("Dipole: " + TString::Format("%e", dipoleFieldStrength) + " T/m^{3};Distance from Dipole (m);#Delta #theta (radians)");
//	else
//		htemp->SetTitle("Gradient: " + TString::Format("%e", gradFieldStrength) + " T/m;Distance from Center (m);#Delta #theta (radians)");
//	theCanvas.SaveAs(SRKGRAPHSDIR + "Plot_SRK_" + runString + "_" + strengthString + "_theta.png");
//	resultsFile.Close();
	/* fit dipole distributions */
//	const int numID = 16;
//	double dipoleValues[numID] = { 1e-13, 1e-13, 1e-12, 1e-12, 1e-14, 1e-14, 1e-13, 1e-13, 1e-15, 1e-15, 1e-14, 1e-14, 1e-16, 1e-16, 1e-15, 1e-15 };
//	double distances[numID] = { 0.05, 0.1, 0.2, 0.3, 0.05, 0.1, 0.2, 0.3, 0.05, 0.1, 0.2, 0.3, 0.05, 0.1, 0.2, 0.3 };
//	for (int i = 0; i < numID; i++)
//	{
//		if(i==8)
//			continue;
//		TString idString = Form("SRK_ID%iP",i);
//		cout << idString << endl;
//		TString filePath = SRKRESULTSDIR + idString + ".root";
//
//		TFile* rootFile = new TFile(filePath);
//
//		double phi;
//		TTree* theTree = (TTree*) rootFile->Get("hitTree");
//		theTree->SetMakeClass(1);
//		theTree->SetBranchAddress("phi", &(phi));
//
//		gROOT->cd(0);
//
//		const int numEntries = theTree->GetEntries();
//
//		vector<double> phiVec;
//		//theTree->GetEntry(0, 1);
//		for (int i = 0; i < numEntries; i++)
//		{
//
//			theTree->GetEntry(i, 1);
//
//			phiVec.push_back(phi);
//
//		}
//		rootFile->Close();
//		delete rootFile;
//
//		double mean;
//		if(true)
//		{
//			mean = reducePeriodicToMeanInVector(phiVec);
//		}
//		else
//		{
//			mean = carefulMeanVector(phiVec);
//		}
//
//		double stdev = carefullStDevVector(phiVec, true);
//		double errorOut = stdev / sqrt((double) phiVec.size());
//
//		cout << "-----------------" << endl;
//
//		cout << "For file: " << filePath << "    Number of particles measured: " << numEntries << endl;
//		cout << "Mean Phase: " << setprecision(15) << mean << " +/- " << errorOut << endl;
//		cout << "Standard Deviation: " << stdev << " +/- " << errorOut << endl;
//		cout << "-----------------" << endl;
//
//
//		TH1* phaseHist = new TH1D("phaseHist", "Dipole: "+TString(Form("%4.2e",dipoleValues[i]))+" T/m  Distance: "+TString(Form("%4.2F",distances[i]))+" m; Final phase (radians);Prob (counts/total counts)", 2048, mean - 3 * stdev, mean + 3 * stdev);
//
//		for (int i = 0; i < numEntries; ++i)
//		{
//			phaseHist->Fill(phiVec[i]);
//		}
//
//		TCanvas c2;
//		gStyle->SetOptStat("eMKS");
//		phaseHist->SetStats(true);
//		phaseHist->Rebin(8);
//
//		c2.SetGrid(1, 1);
//		gPad->SetTickx(1);
//		gPad->SetTicky(1);
//		gPad->SetLogy(1);
//		phaseHist->Draw();
//
//		//TsallisFunc
//		TF1 f1("TsallisFunc", "[0]/pow(1+((x-[3])/[1])*((x-[3])/[1]),[2])", mean - 3 * stdev, mean + 3 * stdev);
//		f1.SetParNames("Amplitude", "Sigma", "Power", "Mean");
//		f1.SetParameters(phaseHist->GetMaximum() * 5, stdev * 0.005, .76, mean);
//		f1.SetParLimits(3, mean, mean);
//
//
//		//		TF1 f1("Gaussian", "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", mean - 3 * stdev, mean + 3 * stdev);
//		//		f1.SetParNames("Amplitude", "Sigma", "Mean");
//		//		f1.SetParameters(100, mean, stdev);
//		//		f1.SetParLimits(0, 80, 1000);
//		//		f1.SetParLimits(1, 0, 0.0001);
//		//		f1.SetParLimits(2, mean-stdev, mean+stdev);
//
//		f1.SetLineColor(kRed);
//		phaseHist->Fit(&f1, "VMR+");
//
//		double power=f1.GetParameter(2);
//		double powerError=f1.GetParError(2);
//		phaseHist->SetTitle(phaseHist->GetTitle()+TString(Form("  Tsallis N=%4.2F#pm%4.2F",power,powerError)));
//		f1.Draw("same");
//
//		c2.SaveAs(SRKGRAPHSDIR + idString + "_Reduced.png");
//		gStyle->SetOptStat("");
//
//		delete phaseHist;
//
//	}
	/*Double check dipoles*/
//	SRKManager theManager;
//	theManager.setDipoleFieldStrength(1e-11);
//	theManager.setB0FieldStrength(0);
//	theManager.setManualTracking(true);
//	theManager.setPos(TVector3(0,0,-0.05));
//	theManager.trackSpins(1,"",SRKRESULTSDIR+"test.root");

	/*Command Line Execution*/
	if(argc >2)
	{
		cout << "Error! SRK only accepts 0 or 1 arguments." << endl;
		return 1;
	}

	SRKManager theManager;
	SRKMacroManager theMacroManager(&theManager);

	if(argc ==2) //argument is the macro file path
	{
		TString macroFilePath = argv[1];
		//TString macroFilePath = SRKMACROSDIR + "test.mac";
		theMacroManager.openMacroFile(macroFilePath);
		theMacroManager.runMacroCommands();
	}
	else
	{
		theMacroManager.enterInteractiveMode();
	}


	return 0;
}
