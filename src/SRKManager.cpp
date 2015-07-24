#include "SRKManager.h"

#include <iostream>
#include <time.h>
#include <iomanip>

using namespace std;

SRKManager::SRKManager()
{
	resultsFile=NULL;
	hitTree=NULL;
	recordAllSteps=false;
	useAltStepping=false;
	parallelFields=true;
	theGlobalField=new SRKGlobalField();
	theSpinTracker=new SRKSpinTracker(theGlobalField);
	theTrack=new SRKTrack();
	b0FieldStrength=1e-6;
	e0FieldStrength=1e6;
	bGradFieldStrength=1e-9;
	dipoleFieldStrength=0;
}

SRKManager::~SRKManager()
{

	delete theGlobalField;
	delete theSpinTracker;
	delete theTrack;
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
		cout << "Error opening track file: " << resultsFilePath << endl;
		return;

	}
	resultsFile->cd();

	hitTree= new TTree("hitTree", "Initial and final states after reflections and spin tracking");
	hitTree->Branch("trackID",&trackID,"trackID/I");
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
	hitTree->Write("", TObject::kOverwrite);

	resultsFile->Close();
	delete resultsFile;
	resultsFile = NULL;
	hitTree = NULL;
}

void SRKManager::writeEvent()
{
	hitTree->Fill();
}

void SRKManager::writeAllSteps(std::vector<SRKMotionState>* stepRecord, std::vector<double>* stepTimes)
{
	for(unsigned int i=0;i<stepRecord->size();i++)
	{
		(stepRecord->at(i))[8]=stepTimes->at(i); //we record the time in the last slot
		setFinalState(stepRecord->at(i));
		writeEvent();
	}

}

double SRKManager::trackSpins(int numTracks, TString trackFilePath, TString resultsFilePath)
{
	clock_t t1,t2;
	t1=clock();

	vector<double> omegaArray; // For calculating mean and error
	omegaArray.reserve(numTracks);


	if(trackFilePath != "")
	{
		theTrack->loadTrackFromFile(trackFilePath);
	}


	TVector3 pos;
	TVector3 vel;
	double currentTime;
	bool lastTrack;


	loadFields();


	SRKMotionState theState(9);
	SRKMotionState initialState(9);
	std::vector<SRKMotionState>* stepRecord=NULL;//Used for recording all steps
	std::vector<double>* stepTimes=NULL; //Used for recording all steps
	if(recordAllSteps)
	{
		stepRecord=new std::vector<SRKMotionState>;
		stepTimes=new std::vector<double>;
	}


	int previousTrackID = -1;

	createResultsFile(resultsFilePath);
	for (int i = 0; i < theTrack->getTrackTreeEntries(); i++)
	{
		theTrack->getTrackTreeEntry(i,pos,vel,currentTime,trackID,lastTrack);
		if(trackID >= numTracks)
		{
			break;
		}
		if(previousTrackID != trackID) //New Track
		{
			if(trackID%100==0) cout << "Spinning track: " << trackID << endl;
//			pos.Print();
//			vel.Print();
			updateMotionStatePosVel(theState, pos, vel, currentTime);
			theState[6] = 0; //Phi
			theState[7] = 0; //Theta
			setInitialState(theState);
			previousTrackID = trackID;
		}
		else
		{
			if(useAltStepping)
			{
				theSpinTracker->trackSpinAltA(theState, currentTime - static_cast<double>(theState[8]),stepRecord,stepTimes); //Runge Kutta on Phi and Theta
			}
			else
			{
				theSpinTracker->trackSpin(theState, currentTime - static_cast<double>(theState[8]),stepRecord,stepTimes); //Runge Kutta on Phi and Theta
			}

			updateMotionStatePosVel(theState, pos, vel, currentTime); //Use the next reflection point for next step
//			pos.Print();
//			vel.Print();
			if(lastTrack) //Record at last point
			{

				if(stepRecord == NULL)
				{
					setFinalState(theState);
					writeEvent(); //Write the final state only
					double deltaOmega=phi/currentTime;
					omegaArray.push_back(deltaOmega);

				}
				else
				{
					writeAllSteps(stepRecord,stepTimes);
				}
//				cout << "Steps taken: " << theSpinTracker.getStepsTaken() << endl;
			//	theSpinTracker.resetStepsTaken();
			}

			previousTrackID = trackID;
		}

	}

	//Calc mean and Error and print
	errorOmega=0;
	meanOmega=0;
	for(unsigned int i=0;i<omegaArray.size();i++)
	{
		meanOmega+=omegaArray[i];
	}
	meanOmega/=numTracks;
	for(int i=0;i<numTracks;i++)
	{
		errorOmega += pow(omegaArray[i]-meanOmega,2);
	}

	errorOmega=sqrt(errorOmega);
	errorOmega /= numTracks-1;

	cout << resultsFilePath << "-- Mean omega: " << setprecision(15) << fixed << meanOmega << " +/- " << scientific << errorOmega << endl;
	closeResultsFile();

	if(stepRecord != NULL)
	{
		delete stepRecord;
	}
	if(stepTimes != NULL)
	{
		delete stepTimes;
	}

	t2=clock();

	double diff =((double)t2-(double)t1)/ (double)CLOCKS_PER_SEC;
	cout << "Computation time: " << fixed << setprecision(3) << diff << "   Tracks per Second: " << (double) numTracks / diff   << "   Tracks per hour: " << (float) numTracks*3600 / diff << endl;;
	cout.unsetf(ios_base::floatfield);
	theTrack->closeTrackFile();
	return static_cast<double>(theState[6]); //Phi
}

void SRKManager::setInitialState(SRKMotionState& initialState)
{
//	cout << "-----------Initial------------" << endl;
//	printMotionState(initialState);

	pos0.SetXYZ(static_cast<double>(initialState[0]),static_cast<double>(initialState[1]),static_cast<double>(initialState[2]));
	vel0.SetXYZ(static_cast<double>(initialState[3]),static_cast<double>(initialState[4]),static_cast<double>(initialState[5]));
	phi0=static_cast<double>(initialState[6]);
	theta0=static_cast<double>(initialState[7]);
	time0=static_cast<double>(initialState[8]);
}

void SRKManager::setFinalState(SRKMotionState& finalState)
{
//	cout << "-----------Final------------" << endl;
//	printMotionState(finalState);

	pos.SetXYZ(static_cast<double>(finalState[0]),static_cast<double>(finalState[1]),static_cast<double>(finalState[2]));
	vel.SetXYZ(static_cast<double>(finalState[3]),static_cast<double>(finalState[4]),static_cast<double>(finalState[5]));
	phi=static_cast<double>(finalState[6]);
	theta=static_cast<double>(finalState[7]);
	time=static_cast<double>(finalState[8]);

}

void SRKManager::loadFields()
{
	theGlobalField->setCurrentFieldSettingsToModify(0); //B0 Field
	theGlobalField->setFieldScalingValue(b0FieldStrength);
	theGlobalField->setFieldDirection(TVector3(0,0,1));

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
	theGlobalField->setFieldDirection(TVector3(0,0,1));

	theGlobalField->setCurrentFieldSettingsToModify(3); //B Dipole
	theGlobalField->setFieldClass(FIELDCLASS_DIPOLE);
	theGlobalField->setFieldScalingValue(dipoleFieldStrength);
	theGlobalField->setFieldMoment(TVector3(0,0,1));
	theGlobalField->setFieldCenterPos(TVector3(0,0,theTrack->getChamberHeight()*.5));

	theGlobalField->updateField();

}

double SRKManager::trackSpinsDeltaOmega(int numTracks, TString runNameString, double& deltaOmegaError)
{
	//Run Parallel
	theTrack->closeTrackFile();
	theTrack->makeTracksCylinder(numTracks);
	setParallelFields(true);
	trackSpins(numTracks, "", SRKRESULTSDIR+runNameString+"P.root");
	double parMean=getMeanOmega();
	double parError=getErrorOmega();



	//Run Anti-Parallel
	theTrack->closeTrackFile();
	theTrack->makeTracksCylinder(numTracks);
	setParallelFields(false);
	trackSpins(numTracks, "", SRKRESULTSDIR+runNameString+"P.root");
	double antiMean=getMeanOmega();
	double antiError=getErrorOmega();

	deltaOmegaError=sqrt(parError*parError+antiError*antiError);
	double deltaOmega=parMean-antiMean;
	cout << "Delta \\omega: " << scientific << setprecision(5) << deltaOmega << " +/- " << deltaOmegaError <<  endl;
	cout.unsetf(ios_base::floatfield);
	return deltaOmega;



}

double SRKManager::trackSpinsDeltaOmegaSteyerl(int numTracks, TString runNameString, double& deltaOmegaSteyerlError)
{
	double meanDeltaOmegaSteyerl = trackSpinsDeltaOmega(numTracks, runNameString, deltaOmegaSteyerlError);
	double zetaEtaOmega0 = getZeta() * getEta()*getOmega0();
	deltaOmegaSteyerlError /= zetaEtaOmega0;
	meanDeltaOmegaSteyerl /= zetaEtaOmega0;

	cout << "Delta \\omega_Steyerl: " << scientific << setprecision(5) << meanDeltaOmegaSteyerl << " +/- " << deltaOmegaSteyerlError << endl;
	cout.unsetf(ios_base::floatfield);
	return meanDeltaOmegaSteyerl;

}

TGraphErrors* SRKManager::trackSpinsDeltaOmegaSteyerlPlot(int numTracksPerPoint, TString runNameString, int numOmegaSteyerl, double OmegaSteyerlStart, double OmegaSteyerlEnd,bool useLog,int approximateReflectionsFixedTime)
{
	vector<double> x(numOmegaSteyerl);
	vector<double> y(numOmegaSteyerl);
	vector<double> eY(numOmegaSteyerl);

	double oldTimeLimit=getTimeLimit(); //in case we are using approximateRefelectionsFixedTime

	double increment;
	if(useLog)
	{
		increment = (log10(OmegaSteyerlEnd) - log10(OmegaSteyerlStart)) / (numOmegaSteyerl - 1);
	}
	else
	{
		increment = (OmegaSteyerlEnd - OmegaSteyerlStart) / (numOmegaSteyerl - 1);
	}

	for (int i = 0; i < numOmegaSteyerl; i++)
	{
		//If we want to emulate a fixed number of reflections but still not bias, we'll use a time limit to get us close to the same reflections

		if(useLog)
		{
			x[i] = pow(10, log10(OmegaSteyerlStart) + increment * i);
		}
		else
		{
			x[i] = OmegaSteyerlStart + increment * i;
		}
		setVelByOmegaSteyerl(x[i]);
		if(approximateReflectionsFixedTime > 0)
		{
			double timeLimit=approximateReflectionsFixedTime * .4 / abs(getMeanVel());
			setTimeLimit( timeLimit);
		}
		cout << " ------------------------------" << endl;
		cout << "Set#: " << i << "    Omega = " << x[i] << "    Vel = "<< getMeanVel() <<endl;
		cout << " ------------------------------" << endl;

		y[i] = trackSpinsDeltaOmegaSteyerl(numTracksPerPoint, runNameString + "_" + TString(i), eY[i]);
	}

	TGraphErrors* outGraph = new TGraphErrors(numOmegaSteyerl, x.data(), y.data(), NULL, eY.data());
	outGraph->SetName(runNameString + "OmegaGraph");
	outGraph->SetTitle(runNameString + " Geometric Phase Plot");
	outGraph->GetYaxis()->SetTitle("#frac{#Delta#omega}{#zeta#eta#omega_{0}}");
	outGraph->GetXaxis()->SetTitle("#omega_{r} / #omega_{0}");
	outGraph->GetXaxis()->SetLimits(0.01,100);
	outGraph->GetXaxis()->SetRangeUser(0.01,100);
	outGraph->GetXaxis()->SetNoExponent();
	//outGraph->GetYaxis()->SetLimits(-3.4999,2);
	//outGraph->GetYaxis()->SetRangeUser(-3.4999,2);

	//Eventually I should figure out where to throw this function...
	ofstream outFile;
	outFile.open(SRKHISTSDIR + runNameString + ".txt");
	outFile << TString("#") << outGraph->GetTitle() << ";" << outGraph->GetXaxis()->GetTitle() << ";" << outGraph->GetYaxis()->GetTitle() << endl;
	for (int i = 0; i < numOmegaSteyerl; i++)
	{

		outFile.precision(15);
		outFile << scientific << x[i];
		outFile << "\t" << scientific << 0;
		outFile << "\t" << scientific << y[i];
		outFile << "\t" << scientific << eY[i];
		if(i != numOmegaSteyerl - 1)
		{
			outFile << endl;
		}

	}
	outFile.close();

	setTimeLimit(oldTimeLimit);



	return outGraph;

}
