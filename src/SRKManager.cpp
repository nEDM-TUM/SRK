#include "SRKManager.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>
#include <time.h>
#include <iomanip>
#include <string>

#include "SRKGraphics.h"
#include "SRKMacroManager.h"
#include "TF1.h"

using namespace std;

SRKManager::SRKManager()
{
	resultsFile = NULL;
	hitTree = NULL;
	recordAllSteps = false;
	useAltStepping = false;
	parallelFields = true;
	theGlobalField = new SRKGlobalField();
	theSpinTracker = new SRKSpinTracker(theGlobalField);
	theMotionTracker = new SRKMotionTracker();
	b0FieldStrength = 1e-6;
	e0FieldStrength = 1e6;
	bGradFieldStrength = 0;
	dipoleFieldStrength = 0;
	dipolePosition.SetXYZ(0, 0, 0);
	dipoleDirection.SetXYZ(0, 0, 1);
	deltaPhaseMean=deltaPhaseError=phaseMean = phaseError = phi = phi0 = theta = theta0 = time = time0 = 0.;
	trackID = 0;
	trackFilePath = "!dynamic";
	defaultResultsDir="/home/mjbales/work/nedm/results/"; //only when resultsFilePath is not appropriate
	resultsFilePath = defaultResultsDir + "test.root";
	runID="RIDX";
}

SRKManager::~SRKManager()
{

	delete theGlobalField;
	delete theSpinTracker;
	delete theMotionTracker;
}

void SRKManager::createResultsFile(TString resultsFilePath)
{
	if(resultsFile != NULL && resultsFile->IsOpen())
	{
		resultsFile->Close();
	}
	resultsFile = new TFile(resultsFilePath, "RECREATE");

	if(resultsFile->IsZombie() || !resultsFile->IsOpen())
	{
		cout << "Error opening results file: " << resultsFilePath << endl;
		return;

	}
	resultsFile->cd();

	hitTree = new TTree("hitTree", "Initial and final states after reflections and spin tracking");
	hitTree->Branch("trackID", &trackID, "trackID/I");
	hitTree->Branch("time0", &time0, "time0/D");
	hitTree->Branch("pos0", &pos0);
	hitTree->Branch("vel0", &vel0);
	hitTree->Branch("phi0", &phi0, "phi0/D");
	hitTree->Branch("theta0", &theta0, "theta0/D");
	hitTree->Branch("time", &time, "time/D");
	hitTree->Branch("pos", &pos);
	hitTree->Branch("vel", &vel);
	hitTree->Branch("phi", &phi, "phi/D");
	hitTree->Branch("theta", &theta, "theta/D");

}

void SRKManager::closeResultsFile()
{

	TList* userInfoList=hitTree->GetUserInfo(); //Every TTree has a list that you can add TObjects to

	userInfoList->Add(new TNamed("RecordAllSteps", Form("%i", (int) getRecordAllSteps())));
	userInfoList->Add(new TNamed("UseAltStepping", Form("%i", (int) getUseAltStepping())));
	userInfoList->Add(new TNamed("ParallelFields", Form("%i", (int) getParallelFields())));
	userInfoList->Add(new TNamed("ConstStepper", Form("%i", (int) getConstStepper())));
	userInfoList->Add(new TNamed("Use2D", Form("%i", (int) getUse2D())));
	userInfoList->Add(new TNamed("ManualTracking", Form("%i", (int) getManualTracking())));
	userInfoList->Add(new TNamed("B0FieldStrength", Form("%e", getB0FieldStrength())));
	userInfoList->Add(new TNamed("E0FieldStrength", Form("%e", getE0FieldStrength())));
	userInfoList->Add(new TNamed("BGradFieldStrength", Form("%e", getBGradFieldStrength())));
	userInfoList->Add(new TNamed("DipoleFieldStrength", Form("%e", getDipoleFieldStrength())));
	userInfoList->Add(new TNamed(TString("TrackFilePath"), getTrackFilePath()));
	userInfoList->Add(new TNamed(TString("ResultsFilePath"), getResultsFilePath()));
	userInfoList->Add(new TNamed(TString("RunID"), getRunID()));
	userInfoList->Add(new TNamed("GyromagneticRatio", Form("%e", getGyromagneticRatio())));
	userInfoList->Add(new TNamed("StepsTaken", Form("%i", getStepsTaken())));
	userInfoList->Add(new TNamed("TimeLimit", Form("%e", getTimeLimit())));
	userInfoList->Add(new TNamed("DiffuseReflectionProb", Form("%e", getDiffuseReflectionProb())));
	userInfoList->Add(new TNamed("ChamberRadius", Form("%e", getChamberRadius())));
	userInfoList->Add(new TNamed("ChamberHeight", Form("%e", getChamberHeight())));
	userInfoList->Add(new TNamed("MeanVel", Form("%e", getMeanVel())));
	userInfoList->Add(new TNamed("ReflectionLimit", Form("%i", getReflectionLimit())));
	userInfoList->Add(new TNamed("Pos", Form("%f %f %f", getPos().X(), getPos().Y(), getPos().Z())));
	userInfoList->Add(new TNamed("Vel", Form("%f %f %f", getVel().X(), getVel().Y(), getVel().Z())));
	userInfoList->Add(new TNamed("DipolePosition", Form("%f %f %f", getDipolePosition().X(), getDipolePosition().Y(), getDipolePosition().Z())));
	userInfoList->Add(new TNamed("DipoleDirection", Form("%f %f %f", getDipoleDirection().X(), getDipoleDirection().Y(), getDipoleDirection().Z())));

	double eps_abs, eps_rel;
	getPerStepError(eps_abs, eps_rel);
	userInfoList->Add(new TNamed("PerStepError", Form("%f %f", eps_abs, eps_rel)));

	resultsFile->Write("", TObject::kOverwrite);

	resultsFile->Close();
	delete resultsFile;
	resultsFile = NULL;
	hitTree = NULL;
}

void SRKManager::loadParametersFromResultsFile(TString filePath)
{
	TFile theFile(filePath,"READ");
	TTree* theTree = (TTree*) theFile.Get("hitTree");
	TList* userInfoList=theTree->GetUserInfo();
	SRKMacroManager tempManager(this);

	for(int i=0;i<userInfoList->GetEntries();i++)
	{
		TNamed* par=(TNamed*) userInfoList->At(i);
		string command=string("set")+par->GetName();
		string value = par->GetTitle();
		if(command != "setStepsTaken" && command != "setPerStepError")
		{
			cout << "SRK: " << command << " " << value << endl;
			tempManager.runMacroCommand(command,value);
		}
	}
	theFile.Close();
}

void SRKManager::writeEvent()
{
	hitTree->Fill();
}


void SRKManager::writeAllSteps(std::vector<SRKMotionState>* stepRecord, std::vector<double>* stepTimes)
{
	for (unsigned int i = 0; i < stepRecord->size(); i++)
	{
		(stepRecord->at(i))[8] = stepTimes->at(i); //we record the time in the last slot
		setFinalState(stepRecord->at(i));
		writeEvent();
	}

}


bool SRKManager::trackSpins(int numTracks)
{
	bool useDynamic = trackFilePath == "!dynamic";

	clock_t t1, t2;
	t1 = clock();

	if(!useDynamic)
	{
		theMotionTracker->loadTrackFile(trackFilePath);
	}

	TVector3 pos;
	TVector3 vel;
	TVector3 posOut;
	TVector3 velOut;
	double currentTime = 0;
	bool lastTrack = false;

	loadFields();

	SRKMotionState theState(9);
	SRKMotionState initialState(9);
	std::vector<SRKMotionState>* stepRecord = NULL; //Used for recording all steps
	std::vector<double>* stepTimes = NULL; //Used for recording all steps
	if(recordAllSteps)
	{
		stepRecord = new std::vector<SRKMotionState>;
		stepTimes = new std::vector<double>;
	}

	createResultsFile(resultsFilePath);

	for (int i = 0; i < numTracks; i++)  //Track loop
	{
		//Starting point
		if(useDynamic)
		{
			if(theMotionTracker->getManualTracking())
			{
				pos = theMotionTracker->getPos();
				vel = theMotionTracker->getVel();
			}
			else
			{
				theMotionTracker->getRandomDirectionAndPointInCylinder(pos, vel);
			}
//			pos.SetXYZ(0,0,0);
			trackID = i;
			lastTrack = false;

		}
		else
		{
			theMotionTracker->getNextTrackTreeEntry(pos, vel, currentTime, trackID, lastTrack);
		}
		updateMotionStatePosVel(theState, pos, vel, currentTime);
		theState[6] = 0; //Phi
		theState[7] = 0; //Theta
		setInitialState(theState);

		if(i % 100 == 0) cout << "Spinning track: " << trackID << endl;

		//Reflection point loop
		do
		{
			if(useDynamic)
			{
				lastTrack = theMotionTracker->getNextTrackingPoint(pos, vel, currentTime);
			}
			else
			{
				theMotionTracker->getNextTrackTreeEntry(pos, vel, currentTime, trackID, lastTrack);

			}

			if(useAltStepping)
			{
				theSpinTracker->trackSpinAltA(theState, currentTime - static_cast<double>(theState[8]), stepRecord, stepTimes); //Runge Kutta on Phi and Theta up to currentTime
			}
			else
			{
				theSpinTracker->trackSpin(theState, currentTime - static_cast<double>(theState[8]), stepRecord, stepTimes); //Runge Kutta on Phi and Theta up to currentTime
			}
			updateMotionStatePosVel(theState, pos, vel, currentTime); //Use the next reflection point for next step
			if(lastTrack) //Record at last point
			{

				if(stepRecord == NULL)
				{
					setFinalState(theState);
					writeEvent(); //Write the final state only
				}
				else
				{
					writeAllSteps(stepRecord, stepTimes);
				}
				currentTime = 0;
			}
		} while (!lastTrack);
	}
	theMotionTracker->closeTrackFile();
	closeResultsFile();
	phaseMean = makeMeanPhasePlot(resultsFilePath, "", true, phaseError); //Prints mean and stdev, no plot

	cout << "-----------------" << endl;
	cout << "For file: " << resultsFilePath << "    Number of particles measured: " << numTracks << endl;
	cout << "Mean Phase: " << setprecision(15) << phaseMean << " +/- " << phaseError << endl;
	cout << "-----------------" << endl;

	theMotionTracker->closeTrackFile();

	if(stepRecord != NULL)
	{
		delete stepRecord;
	}
	if(stepTimes != NULL)
	{
		delete stepTimes;
	}

	t2 = clock();

	double diff = ((double) t2 - (double) t1) / (double) CLOCKS_PER_SEC;
	cout << "Computation time: " << fixed << setprecision(3) << diff << "   Tracks per Second: " << (double) numTracks / diff << "   Tracks per hour: " << (float) numTracks * 3600 / diff << endl;
	cout.unsetf(ios_base::floatfield);

	return true;
}

void SRKManager::setInitialState(SRKMotionState& initialState)
{
	//cout << "-----------Initial------------" << endl;
	//printMotionState(initialState);

	pos0.SetXYZ(static_cast<double>(initialState[0]), static_cast<double>(initialState[1]), static_cast<double>(initialState[2]));
	vel0.SetXYZ(static_cast<double>(initialState[3]), static_cast<double>(initialState[4]), static_cast<double>(initialState[5]));
	phi0 = static_cast<double>(initialState[6]);
	theta0 = static_cast<double>(initialState[7]);
	time0 = static_cast<double>(initialState[8]);
}

void SRKManager::setFinalState(SRKMotionState& finalState)
{
	//cout << "-----------Final------------" << endl;
	//printMotionState(finalState);

	pos.SetXYZ(static_cast<double>(finalState[0]), static_cast<double>(finalState[1]), static_cast<double>(finalState[2]));
	vel.SetXYZ(static_cast<double>(finalState[3]), static_cast<double>(finalState[4]), static_cast<double>(finalState[5]));
	phi = static_cast<double>(finalState[6]);
	theta = static_cast<double>(finalState[7]);
	time = static_cast<double>(finalState[8]);

}

void SRKManager::loadFields()
{
	theGlobalField->setCurrentFieldSettingsToModify(0); //B0 Field
	theGlobalField->setFieldScalingValue(b0FieldStrength);
	theGlobalField->setFieldDirection(TVector3(0, 0, 1));

	theGlobalField->setCurrentFieldSettingsToModify(1); //E0 Field
	theGlobalField->setFieldType(FIELD_ELECTRIC);
	theGlobalField->setFieldDirection(TVector3(0, 0, 1));

	if(parallelFields)
		theGlobalField->setFieldScalingValue(e0FieldStrength);
	else
		theGlobalField->setFieldScalingValue(-e0FieldStrength);

	theGlobalField->setCurrentFieldSettingsToModify(2); //B Grad Field
	theGlobalField->setFieldClass(FIELDCLASS_GRADIENT);
	theGlobalField->setFieldScalingValue(bGradFieldStrength);
	theGlobalField->setFieldDirection(TVector3(0, 0, 1));

	theGlobalField->setCurrentFieldSettingsToModify(3); //B Dipole
	theGlobalField->setFieldClass(FIELDCLASS_DIPOLE);
	theGlobalField->setFieldScalingValue(dipoleFieldStrength);
	theGlobalField->setFieldMoment(dipoleDirection);
	theGlobalField->setFieldCenterPos(dipolePosition);

	theGlobalField->updateField();

}

void SRKManager::calcDeltaPhaseMean(TString inpRunID)
{
	double parMean,parError,antiMean,antiError;
	parMean=makeMeanPhasePlot(defaultResultsDir+"Results_"+inpRunID+"_P.root", "", true, parError); //Prints mean and stdev, no plot
	antiMean=makeMeanPhasePlot(defaultResultsDir+"Results_"+inpRunID+"_A.root", "", true, antiError);//Prints mean and stdev, no plot
	deltaPhaseMean = parMean - antiMean;
	deltaPhaseError = sqrt(parError * parError + antiError * antiError);

}

SRKRunStats SRKManager::calcResultsFileStats(TString filePath, bool useWrapping)
{
	SRKRunStats theStats;

	TFile* rootFile = new TFile(filePath);

	double phi, theta;
	TTree* theTree = (TTree*) rootFile->Get("hitTree");
	theTree->SetMakeClass(1);
	theTree->SetBranchAddress("phi", &(phi));
	theTree->SetBranchAddress("theta", &(theta));

	gROOT->cd(0);

	theStats.numEvents = theTree->GetEntries();

	vector<double> phiVec;
	vector<double> thetaVec;

	for (int i = 0; i < theStats.numEvents; i++)
	{

		theTree->GetEntry(i, 1);

		phiVec.push_back(phi);
		thetaVec.push_back(phi);

	}
	rootFile->Close();
	delete rootFile;

	if(useWrapping)
	{
		theStats.phiMean = reducePeriodicToMeanInVector(phiVec);
		theStats.thetaMean = reducePeriodicToMeanInVector(thetaVec);
	}
	else
	{
		theStats.phiMean = carefulMeanVector(phiVec);
		theStats.thetaMean = carefulMeanVector(thetaVec);
	}

	theStats.phiStDev = carefullStDevVector(phiVec, true);
	theStats.thetaStDev = carefullStDevVector(thetaVec, true);

	theStats.phiError = theStats.phiStDev / sqrt((double) theStats.numEvents);
	theStats.thetaError = theStats.thetaStDev / sqrt((double) theStats.numEvents);

	TH1D phiHist("phiHist", "phiHist", 10000, theStats.phiMean - TMath::Pi(), theStats.phiMean + TMath::Pi());
	TH1D thetaHist("phiHist", "phiHist", 10000, TMath::Pi(), TMath::Pi());

	for (int i = 0; i < theStats.numEvents; ++i)
	{
		phiHist.Fill(phiVec[i]);
		thetaHist.Fill(thetaVec[i]);
	}

	theStats.phiKurtosis = phiHist.GetKurtosis(1);
	theStats.phiKurtosisError = phiHist.GetKurtosis(11);
	theStats.phiSkewness = phiHist.GetSkewness(1);
	theStats.phiSkewnessError = phiHist.GetSkewness(11);

	theStats.thetaKurtosis = thetaHist.GetKurtosis(1);
	theStats.thetaKurtosisError = thetaHist.GetKurtosis(11);
	theStats.thetaSkewness = thetaHist.GetSkewness(1);
	theStats.thetaSkewnessError = thetaHist.GetSkewness(11);

	TF1 phiTsallisFunc("phiTsallisFunc", "[0]/pow(1+((x-[3])/[1])*((x-[3])/[1]),[2])", theStats.phiMean - TMath::Pi(), theStats.phiMean + TMath::Pi());
	phiTsallisFunc.SetParNames("Amplitude", "Sigma", "Power", "Mean");
	phiTsallisFunc.SetParameters(phiHist.GetMaximum() * 5, theStats.phiStDev * 0.005, .76, theStats.phiMean);
	phiTsallisFunc.SetParLimits(3, theStats.phiMean, theStats.phiMean);
	phiHist.Fit(&phiTsallisFunc, "VMR+");
	theStats.phiTsallisPower = phiTsallisFunc.GetParameter(2);
	theStats.phiTsallisPowerError = phiTsallisFunc.GetParError(2);

	TF1 thetaTsallisFunc("thetaTsallisFunc", "[0]/pow(1+((x-[3])/[1])*((x-[3])/[1]),[2])", theStats.thetaMean - TMath::Pi(), theStats.thetaMean + TMath::Pi());
	thetaTsallisFunc.SetParNames("Amplitude", "Sigma", "Power", "Mean");
	thetaTsallisFunc.SetParameters(thetaHist.GetMaximum() * 5, theStats.thetaStDev * 0.005, .76, theStats.thetaMean);
	thetaTsallisFunc.SetParLimits(3, theStats.thetaMean, theStats.thetaMean);
	thetaHist.Fit(&thetaTsallisFunc, "VMR+");
	theStats.thetaTsallisPower = thetaTsallisFunc.GetParameter(2);
	theStats.thetaTsallisPowerError = thetaTsallisFunc.GetParError(2);


	return theStats;
}

void SRKManager::outputDataForRIDs(TString rangeString)  //Format of int int
{
	stringstream theStringStream(rangeString.Data());
	int startRID, finalRID;
	theStringStream >> startRID >> finalRID;

	TString outFilePath = defaultResultsDir + Form("Data_RID%i_RID%i.txt", startRID, finalRID);
	ofstream outFile(outFilePath.Data());

	outFile << "#RunID\tNumEvents\tDeltaOmega\tDeltaOmegaError\tFalseEDM\tFalseEDMError";
	outFile << "\tphiStDev\tphiStDevError\tphiKurtosis\tphiKurtosisError\tphiSkewness\tphiSkewnessError\tphiTsallisPower\tPhiTsallisPowerError";
	outFile << "\tthetaStDev\tthetaStDevError\tthetaKurtosis\tthetaKurtosisError\tthetaSkewness\tthetaSkewnessError\tthetaTsallisPower\tthetaTsallisPowerError" << endl;

	for (int i = startRID; i < finalRID + 1; i++)
	{
		TString pFilePath = defaultResultsDir + Form("Results_RID%i_P.root", i);
		TString aFilePath = defaultResultsDir + Form("Results_RID%i_A.root", i);
		if(fileExists(pFilePath) && fileExists(aFilePath))
		{

			loadParametersFromResultsFile(defaultResultsDir + Form("Results_RID%i_P.root", i));

			SRKRunStats parStats = calcResultsFileStats(pFilePath, true);
			SRKRunStats antiStats = calcResultsFileStats(aFilePath, true);
			double deltaPhaseMean = parStats.phiMean - antiStats.phiMean;
			double deltaPhaseError = sqrt(parStats.phiError * parStats.phiError + antiStats.phiError * antiStats.phiError);
			double scaleFactor = 100. * 6.58211928E-016 / (4. * e0FieldStrength);
			double deltaOmega = deltaPhaseMean / getTimeLimit();
			double deltaOmegaError = deltaPhaseError / getTimeLimit();

			outFile << i << "\t" << parStats.numEvents << "\t" << deltaOmega << "\t" << deltaOmegaError << "\t" << deltaOmega * scaleFactor << "\t" << deltaOmegaError * scaleFactor;
			outFile << "\t" << parStats.phiStDev << "\t" << parStats.phiError << "\t" << parStats.phiKurtosis << "\t" << parStats.phiKurtosisError << "\t" << parStats.phiSkewness << "\t" << parStats.phiSkewnessError << "\t"
				<< parStats.phiTsallisPower << "\t" << parStats.phiTsallisPowerError;
			outFile << "\t" << parStats.thetaStDev << "\t" << parStats.thetaError << "\t" << parStats.thetaKurtosis << "\t" << parStats.thetaKurtosisError << "\t" << parStats.thetaSkewness << "\t" << parStats.thetaSkewnessError << "\t"
				<< parStats.thetaTsallisPower << "\t" << parStats.thetaTsallisPowerError;
		}
		else  //File does not exist
		{
			cout << "Cannot find one or both files for RID" << i <<"!" << endl;
			outFile << i;
		}
		if(i < finalRID)
		{
			outFile << endl;
		}
	}

	outFile.close();

}

void SRKManager::trackSpinsDeltaOmega(int numTracks)
{
	//Run Parallel
	setParallelFields(true);
	resultsFilePath = defaultResultsDir + "Results_" + runID + "_P.root";
	trackSpins(numTracks);


	//Run Anti-Parallel
	setParallelFields(false);
	resultsFilePath = defaultResultsDir + "Results_" + runID + "_A.root";
	trackSpins(numTracks);

	calcDeltaPhaseMean(runID);
	double deltaOmega=deltaPhaseMean/getTimeLimit();
	double deltaOmegaError=deltaPhaseError/getTimeLimit();


	cout << "Delta \\omega [rad // s]: " << scientific << setprecision(5) << deltaOmega << " +/- " << deltaOmegaError << endl;
	double scaleFactor = 100. * 6.58211928E-016 / (4. * e0FieldStrength);
	cout << "False EDM [e cm]: " << scientific << setprecision(5) << deltaOmega*scaleFactor << " +/- " << deltaOmegaError*scaleFactor << endl;
	double zetaEtaOmega0 = getZeta() * getEta() * getOmega0();

	if(bGradFieldStrength > 0)
	{
		cout << "Delta \\omega_Steyerl: " << scientific << setprecision(5) << deltaOmega/zetaEtaOmega0 << " +/- " << deltaOmegaError/zetaEtaOmega0 << endl;
		cout.unsetf(ios_base::floatfield);
	}

}

bool SRKManager::fileExists(TString strFileName)
{
	struct stat stFileInfo;
	return stat(strFileName.Data(), &stFileInfo) == 0;
}

//TGraphErrors* SRKManager::trackSpinsDeltaOmegaSteyerlPlot(int numTracksPerPoint, int numOmegaSteyerl, double OmegaSteyerlStart, double OmegaSteyerlEnd, bool useLog, int approximateReflectionsFixedTime)
//{
//	vector<double> x(numOmegaSteyerl);
//	vector<double> y(numOmegaSteyerl);
//	vector<double> eY(numOmegaSteyerl);
//
//	double oldTimeLimit = getTimeLimit(); //in case we are using approximateRefelectionsFixedTime
//
//	double increment;
//	if(useLog)
//	{
//		increment = (log10(OmegaSteyerlEnd) - log10(OmegaSteyerlStart)) / (numOmegaSteyerl - 1);
//	}
//	else
//	{
//		increment = (OmegaSteyerlEnd - OmegaSteyerlStart) / (numOmegaSteyerl - 1);
//	}
//
//	for (int i = 0; i < numOmegaSteyerl; i++)
//	{
//		TString subRunID=runID+"_" + Form("%i",i);
//		//If we want to emulate a fixed number of reflections but still not bias, we'll use a time limit to get us close to the same reflections
//
//		if(useLog)
//		{
//			x[i] = pow(10, log10(OmegaSteyerlStart) + increment * i);
//		}
//		else
//		{
//			x[i] = OmegaSteyerlStart + increment * i;
//		}
//		setVelByOmegaSteyerl(x[i]);
//		if(approximateReflectionsFixedTime > 0)
//		{
//			double timeLimit = approximateReflectionsFixedTime * .4 / abs(getMeanVel());
//			setTimeLimit(timeLimit);
//		}
//		cout << " ------------------------------" << endl;
//		cout << "Set#: " << i << "    Omega = " << x[i] << "    Vel = " << getMeanVel() << endl;
//		cout << " ------------------------------" << endl;
//		trackSpinsDeltaOmega(numTracksPerPoint);
//		double zetaEtaOmega0 = getZeta() * getEta() * getOmega0();
//		y[i] = deltaPhaseMean/zetaEtaOmega0;
//		eY[i] = deltaPhaseError/zetaEtaOmega0;
//	}
//
//	TCanvas theCanvas("theCanvas", "theCanvas", 1600, 1200);
//	theCanvas.SetLogx();
//
//	TGraphErrors* outGraph = new TGraphErrors(numOmegaSteyerl, x.data(), y.data(), NULL, eY.data());
//	outGraph->SetName(runNameString + "OmegaGraph");
//	outGraph->SetTitle(runNameString + " Geometric Phase Plot");
//	outGraph->GetYaxis()->SetTitle("#frac{#Delta#omega}{#zeta#eta#omega_{0}}");
//	outGraph->GetXaxis()->SetTitle("#omega_{r} / #omega_{0}");
//	outGraph->GetXaxis()->SetLimits(0.01, 100);
//	outGraph->GetXaxis()->SetRangeUser(0.01, 100);
//	outGraph->GetXaxis()->SetNoExponent();
//	//outGraph->GetYaxis()->SetLimits(-3.4999,2);
//	//outGraph->GetYaxis()->SetRangeUser(-3.4999,2);
//
//	outGraph->Draw("APE1");
//	theCanvas.SaveAs(SRKGRAPHSDIR + runNameString + ".png");
//
//	//Eventually I should figure out where to throw this function...
//	ofstream outFile;
//	outFile.open(SRKHISTSDIR + runNameString + ".txt");
//	outFile << TString("#") << outGraph->GetTitle() << ";" << outGraph->GetXaxis()->GetTitle() << ";" << outGraph->GetYaxis()->GetTitle() << endl;
//	for (int i = 0; i < numOmegaSteyerl; i++)
//	{
//
//		outFile.precision(15);
//		outFile << scientific << x[i];
//		outFile << "\t" << scientific << 0;
//		outFile << "\t" << scientific << y[i];
//		outFile << "\t" << scientific << eY[i];
//		if(i != numOmegaSteyerl - 1)
//		{
//			outFile << endl;
//		}
//
//	}
//	outFile.close();
//
//	setTimeLimit(oldTimeLimit);
//
//	return outGraph;
//}
