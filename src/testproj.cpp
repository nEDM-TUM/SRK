#include <iostream>

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
//	// Operations at 128-bit precision and full numeric_limits support:
//	float128 b = 2;
//	// There are 113-bits of precision:
//	std::cout << std::numeric_limits<float128>::digits << std::endl;
//	// Or 34 decimal places:
//	std::cout << std::numeric_limits<float128>::digits10 << std::endl;
//	// We can use any C++ std lib function, lets print all the digits as well:
//	std::cout << std::setprecision(std::numeric_limits<float128>::max_digits10) << log(b) << std::endl; // print log(2) = 0.693147180559945309417232121458176575

//	double x = 1;
//	double y = 1e-10;
//	double z = x;
//
//	for (long i = 0; i < 10000000000; ++i)//10000000000 = 7 sec
//	{
//		z=z+y;
//	}
//
//	cout << std::setprecision(20)<< std::scientific << z << endl;


//	float128 x = 1;
//	float128 y = 1e-10Q;
//	float128 z = x;
//
//	for (long i = 0; i < 1000000000; ++i)//1000000000 = 17 s-I/home/mjbales/work/software/root/root_v5.34.21/include
//	{
//		z=z+y;
//	}
//
//	cout << std::setprecision(std::numeric_limits<float128>::max_digits10) << std::scientific << z << endl;

	//Conclusion double is ~ 24x faster than float128

//	float128 x = 1;
//	double y = 1e-10Q;
//	float128 z = x;
//
//	for (long i = 0; i < 1000000000; ++i)//1000000000 = 17 s
//	{
//		z=z+y;
//	}

//	cout << std::setprecision(std::numeric_limits<float128>::max_digits10) << std::scientific << z << endl;


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


/*Omega loop*/
//	SRKManager theManager;
//
//	delete gRandom;
//	gRandom=new TRandom3(0);
//	int numTracksPerOmega=1000;
//	const int numOmega=114;
//	double begin=.1,end =1.347500000000000;
//	double logBegin= log10(begin);
//	double logEnd= log10(end);
//	double logIncrement=(logEnd-logBegin)/(numOmega-1);
//	double omegaR,velocity,timeLimit;
//	double radius=0.24;
//	double targetReflections=100;
//	double bGradientStrength=1e-9;
//	double bFieldStrength=1e-6;
//	double eFieldStrength=1e6;
//	double omega0=theManager.getGyromagneticRatio()*bFieldStrength;
//	TString idString="Diffuse";
//	double diffuseProb=1;
//	double x[numOmega];
//	double y[numOmega];
//	double eY[numOmega];
//	for (int i = 0; i < numOmega; ++i)
//	{
//		x[i]=pow(10,logBegin+logIncrement*i);
//		cout <<" ------------------------------" << endl;
//		cout << "Set#: " << i << "    Omega = " << x[i] << endl;
//		cout <<" ------------------------------" << endl;
//		omegaR =  x[i]* omega0;
//		velocity = omegaR * radius;
//		timeLimit = targetReflections * 0.75 / abs(velocity);
//		double Zeta=bGradientStrength*radius/(2*bFieldStrength);
//		double Eta=radius*omega0*eFieldStrength/(bFieldStrength*299792458*299792458);
//
//		//Make track
//		SRKTrack theTrack;
//		theTrack.setTimeLimit(timeLimit);
//		theTrack.setMeanVel(velocity);
//		theTrack.setDiffuseReflectionProb(diffuseProb);
//		theTrack.makeTracksCylinder(numTracksPerOmega);
//		theTrack.writeTrackToFile("/home/mjbales/work/code/testproj/output/trackTest.root");  //overwritten each time
//
//		//Run Parallel
//		theManager.setParallelFields(true);
//		theManager.trackSpins(numTracksPerOmega, "/home/mjbales/work/code/testproj/output/trackTest.root", "/home/mjbales/work/code/testproj/output/"+idString+"P"+TString::Format("%d",i)+".root");
//
//		y[i]=theManager.getMeanOmega();
//		double error1=theManager.getErrorOmega();
//
//		//Run Anti-Parallel
//		theManager.setParallelFields(false);
//		theManager.trackSpins(numTracksPerOmega, "/home/mjbales/work/code/testproj/output/trackTest.root", "/home/mjbales/work/code/testproj/output/"+idString+"A"+TString::Format("%d",i)+".root");
//
//		y[i]-=theManager.getMeanOmega();
//		y[i] /=  omega0*Zeta*Eta;
//		double error2=theManager.getErrorOmega();
//		eY[i]=sqrt(error1*error1+error2*error2)/(omega0*Zeta*Eta);
//
//	}
//
//	ofstream outFile;
//	outFile.open("/home/mjbales/work/nedm/hists/SRKTest"+ idString + ".txt");
//	double* eX=NULL;
//	outFile << "#SRKTest;#omega_{r} / #omega_{0};#frac{#Delta#omega}{#zeta#eta#omega_{0}}" << endl;
//	for (int i = 0; i < numOmega; i++)
//	{
//
//		outFile.precision(15);
//		outFile << scientific << x[i];
//		if(eX != NULL)
//		{
//			outFile << "\t" << scientific << eX[i];
//		}
//		else
//		{
//			outFile << "\t" << scientific << 0;
//		}
//
//		outFile << "\t" << scientific << y[i];
//		if(eY != NULL)
//		{
//			outFile << "\t" << scientific << eY[i];
//		}
//		else
//		{
//			outFile << "\t" << scientific << 0;
//		}
//
//		if(i != numOmega - 1)
//		{
//			outFile << endl;
//		}
//
//	}
//	outFile.close();

//	SRKTrack theTrack;
//	theTrack.setTimeLimit(10);
//	theTrack.setMeanVel(7.26);
//	theTrack.setDiffuseReflectionProb(1);
//	theTrack.makeTracksCylinder(1000);
//	theTrack.writeTrackToFile(SRKRESULTSDIR+"diffuseHgOmega1.root");
//	theTrack.closeTrackFile();


	/* Simpler Omega Styerl plot command =D */
	SRKManager theManager;
	theManager.setUseAltStepping(true);
	theManager.setConstStepper(false);
	theManager.setE0FieldStrength(8.33e5);

//	theManager.setBGradFieldStrength(1e-9);
//	theManager.setTimeLimit(10);
//	theManager.setMeanVel(7.26);
//
//	vector< double > deltaMeanOmega;
//	for(int i=0;i<9;i++)
//	{
//		double theEPS=1.0e-1*pow(10,-i);
//		theManager.setPerStepError(theEPS,theEPS);
//
//		//double initStep=pow(10,-i);
//		//cout << "Init Step: " << initStep << endl;
//		//theManager.setInitialStepSize(initStep);
//		theManager.trackSpins(1000, SRKRESULTSDIR+"diffuseHgOmega1.root", SRKRESULTSDIR+"diffuseHgResult.root");
//		deltaMeanOmega.push_back(theManager.getMeanOmega());
//
//		//cout << scientific << "EPS: " << theEPS << "        Total Error: " << theManager.getStepsTaken()*theEPS/1000. << endl;
//	}

//	for(unsigned int i=0;i<deltaMeanOmega.size();i++)
//	{
//		cout << "Delta Mean Omega:  " << setprecision(15)<< deltaMeanOmega[deltaMeanOmega.size()-1] - deltaMeanOmega[i] << endl;
//	}



	TGraphErrors* theGraph=theManager.trackSpinsDeltaOmegaSteyerlPlot(100, "DiffuseTest4", 10,  0.1,1.3475,true,100);
	convertTGraphErrorsToTXT(theGraph,SRKHISTSDIR+"DiffuseTest4.txt");

//	TGraphErrors* theGraph=getTabSeperatedTGraphErrors(SRKHISTSDIR+"DiffuseTest2.txt");
//
	theGraph->SetMarkerStyle(4);
	TCanvas theCanvas("theCanvas","theCanvas",1600,1200);
	theCanvas.SetLogx();
	theGraph->Draw("APE1");
	theCanvas.SaveAs(SRKGRAPHSDIR+"DiffuseTest4.png");


	return 0;
}
