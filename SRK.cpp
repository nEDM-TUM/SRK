#include <iostream>
#include <time.h>

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
	TString runName="Dipole_1en11_30cm_Dynamic_1s_193m_s_E1e8";

	SRKManager theManager;
	theManager.setConstStepper(false);
	//theManager.setPerStepError(1e-8,1e-8);
	theManager.setUse2D(false);
	theManager.setChamberHeight(.12);
	theManager.setChamberRadius(0.235);
	theManager.setDiffuseReflectionProb(1);
	theManager.setDipolePosition(TVector3(0,0,-0.3-0.5*0.12));//30 cm below bottom of chamber
//	theManager.setDipolePosition(TVector3(0,0,0.5*0.12));//Chamber top
	theManager.setE0FieldStrength(1.e8);
	//theManager.setE0FieldStrength(0);
	theManager.setTimeLimit(1);
	theManager.setMeanVel(193);
	theManager.setDipoleFieldStrength(1e-11);
	theManager.setBGradFieldStrength(0);
	double falseEDMError;
	double falseEDM=theManager.trackSpinsFalseEDM(2300000,runName, falseEDMError);

	TString resultPath=SRKHISTSDIR+runName+".txt";
	ofstream resultFile(resultPath.Data());
	resultFile << scientific << setprecision(5) << falseEDM << "\t" << falseEDMError <<  endl;
	resultFile.close();

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

	return 0;
}
