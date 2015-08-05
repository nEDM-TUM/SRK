#include <iostream>
#include <time.h>
#include <functional>
#include <unordered_map>
#include <string>

#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"

#include "SRKManager.h"
#include "SRKMotionTracker.h"
#include "SRKGraphics.h"

using namespace std;

int main(int argc, char* argv[])
{
	delete gRandom;
	gRandom=new TRandom3(0);

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
//	//theManager.setE0FieldStrength(0);
//	//theManager.setBGradFieldStrength(0);
//
//	theManager.setTimeLimit(10);
//	theManager.setMeanVel(7.26);
//
//	vector< double > deltaMeanOmega;
//	vector< double > processingTime;
//
//
//
//
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
//
//		//cout << scientific << "EPS: " << theEPS << "        Total Error: " << theManager.getStepsTaken()*theEPS/1000. << endl;
//	}
//
//	for(unsigned int i=0; i<deltaMeanOmega.size(); i++)
//	{
//		double theEPS=1.0e-6*pow(10.,-i*.1);
//		cout << setprecision(13) << theEPS << "\t" << deltaMeanOmega[deltaMeanOmega.size()-1] - deltaMeanOmega[i] << "\t" << processingTime[i] << endl;
//	}



//	TGraphErrors* theGraph=theManager.trackSpinsDeltaOmegaSteyerlPlot(100, "DiffuseTest4", 10,  0.1,1.3475,true,100);
//	convertTGraphErrorsToTXT(theGraph,SRKHISTSDIR+"DiffuseTest4.txt");
//
////	TGraphErrors* theGraph=getTabSeperatedTGraphErrors(SRKHISTSDIR+"DiffuseTest2.txt");
////
//	theGraph->SetMarkerStyle(4);
//	TCanvas theCanvas("theCanvas","theCanvas",1600,1200);
//	theCanvas.SetLogx();
//	theGraph->Draw("APE1");
//	theCanvas.SaveAs(SRKGRAPHSDIR+"DiffuseTest4.png");

/* Test dipole ILL setup*/
	TString runName="Dipole_1en12_30cm_Dynamic_100s_193m_s_E1e7";

	SRKManager theManager;
	theManager.setConstStepper(false);
	theManager.setUseAltStepping(false);
	theManager.setPerStepError(1e-7,1e-7);
	//theManager.setInitialStepSize(3e-3);
	theManager.setUse2D(false);
	theManager.setChamberHeight(.12);
	theManager.setChamberRadius(0.235);
	theManager.setDiffuseReflectionProb(1);
	theManager.setDipolePosition(TVector3(0,0,-0.3-0.05*0.12));
	theManager.setE0FieldStrength(1.e7);
	theManager.setTimeLimit(100);
	theManager.setMeanVel(193);
	theManager.setDipoleFieldStrength(1e-12);
	theManager.setBGradFieldStrength(0);
	double falseEDMError;
	double falseEDM=theManager.trackSpinsFalseEDM(10000,runName, falseEDMError);

	TString resultPath=SRKHISTSDIR+runName+".txt";
	ofstream resultFile(resultPath.Data());
	resultFile << scientific << setprecision(5) << falseEDM << "\t" << falseEDMError <<  endl;
	resultFile.close();

	makeMeanPhasePlot(SRKRESULTSDIR+runName+"P.root", SRKGRAPHSDIR+runName+"P_Reduced.pdf", true);
	makeMeanPhasePlot(SRKRESULTSDIR+runName+"P.root", SRKGRAPHSDIR+runName+"P.pdf", false);

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



	return 0;
}
