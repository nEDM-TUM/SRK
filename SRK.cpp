#include <iostream>
#include <time.h>

#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"

#include "SRKManager.h"
#include "SRKTrack.h"
#include "SRKGraphics.h"

using namespace std;

int main() {
	delete gRandom;
	gRandom=new TRandom3(0);

/*Testing*/
//	delete gRandom;
//	gRandom=new TRandom3(0);
//
//	SRKTrack theTrack;
//	theTrack.setReflectionLimit(1000);
//	//theTrack.setTimeLimit(10);
//	theTrack.setMeanVel(193);
//	theTrack.setDiffuseReflectionProb(1);
//	theTrack.makeTracksCylinder(1000);
//	theTrack.writeTrackToFile("/home/mjbales/work/code/testproj/output/trackTest.root");
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
	SRKManager theManager;
	theManager.setUseAltStepping(false);
	theManager.setConstStepper(true);
	theManager.setE0FieldStrength(0);
	theManager.setBGradFieldStrength(0);

	theManager.setTimeLimit(10);
	theManager.setMeanVel(7.26);

	vector< double > deltaMeanOmega;
	vector< double > processingTime;




	for(int i=0;i<4;i++)
	{
		clock_t t1,t2;
		t1=clock();
//		double theEPS=1.0e-1*pow(10.,-i);
//		theManager.setPerStepError(theEPS,theEPS);
//		cout << theEPS << endl;

		double initStep=pow(10,-i);
		cout << "Init Step: " << initStep << endl;
		theManager.setInitialStepSize(initStep);
		theManager.trackSpins(1000, SRKRESULTSDIR+"diffuseHgOmega1.root", SRKRESULTSDIR+"diffuseHgResult.root");
		deltaMeanOmega.push_back(theManager.getMeanOmega());
		t2=clock();
		double diff =((double)t2-(double)t1)/ (double)CLOCKS_PER_SEC;
		processingTime.push_back(diff);

		//cout << scientific << "EPS: " << theEPS << "        Total Error: " << theManager.getStepsTaken()*theEPS/1000. << endl;
	}

	for(unsigned int i=0; i<deltaMeanOmega.size(); i++)
	{
		cout << setprecision(13) << 1.0e-1*pow(10.,-i) << "\t" << deltaMeanOmega[0] - deltaMeanOmega[i] << "\t" << processingTime[i] << endl;
	}



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


	return 0;
}
