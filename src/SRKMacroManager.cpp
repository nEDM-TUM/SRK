#include "SRKMacroManager.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "TRandom.h"

using namespace std;

SRKMacroManager::SRKMacroManager(SRKManager* inpManager)
{
	theManager = inpManager;
	defineCommands();

}

SRKMacroManager::~SRKMacroManager()
{

}

bool SRKMacroManager::openMacroFile(const TString filePath)
{
	ifstream theFile(filePath);
	TString sBuffer("#");
	if(!theFile.is_open() || theFile.fail() || theFile.eof())
	{
		cout << "Macro file: " << filePath << " not recognized" << endl;
		return false;
	}

	cout << "Opening macro file: " << filePath << endl;
	list<TString> fileCommandStringList;
	list<TString> fileCommandValueStringList;

	while (!theFile.fail() || !theFile.eof())
	{
		skipCommentLines(theFile);
		if(theFile.eof()) break;
		TString commandString, commandValueString;
		theFile >> commandString;
		theFile.ignore(1);
		commandValueString.ReadLine(theFile); //rest is value
		//cout << "Grabbed: " << commandString << "  and value: " << commandValueString << endl;

		fileCommandStringList.push_back(commandString);
		fileCommandValueStringList.push_back(commandValueString);
	}

	theFile.close();

	commandStringList.splice(commandStringList.begin(), fileCommandStringList);
	cout << "Number of commands: " << commandStringList.size() << endl;
	commandValueStringList.splice(commandValueStringList.begin(), fileCommandValueStringList);

	return true;
}

void SRKMacroManager::runMacroCommands()
{
	while (!commandStringList.empty())
	{

		TString command = commandStringList.front();
		TString value = commandValueStringList.front();
		commandStringList.pop_front();
		commandValueStringList.pop_front();

		cout << "SRK: " << command << " " << value << endl;
		if(!runMacroCommand(command.Data(), value.Data()))
		{
			break;  //Break if there is a command issue;
		}

	}
}

void SRKMacroManager::enterInteractiveMode()
{
	cout << endl << endl;
	cout << "*******************************************************" << endl;
	cout << "*  SRK Interactive Mode" << endl;
	cout << "*  Type \"q\" to exit." << endl;
	cout << "*******************************************************" << endl;
	bool endInteractiveMode = false;
	stringstream lineStream;
	while (!endInteractiveMode)
	{
		cout << "SRK: ";
		string theLine;
		getline(cin, theLine);
		if(theLine == "end" || theLine == "quit" || theLine == "exit" || theLine == "q")
		{
			cout << "Exiting Interactive Mode..." << endl;
			endInteractiveMode = true;
			break;
		}
		lineStream.clear();
		lineStream << theLine;
		string commandString, commandValueString;
		lineStream >> commandString;
		lineStream.ignore(1);
		getline(lineStream, commandValueString);
		//cout << "Command: " << commandString << "  Value: " << commandValueString << endl;
		runMacroCommand(commandString, commandValueString);

		if(commandString == "runMacroFile") //Need to follow up and run the commands that got added but still return to interactive mode.
		{
			runMacroCommands();
		}
	}
}

bool SRKMacroManager::runMacroCommand(const string command, const string value)
{
	if(commandMap.count(command) == 0) //Command doesn't exist?
	{
		cout << "Command: " << command << " doesn't exist!" << endl;
		return false;
	}

	commandMap[command](value);
	return true;
}

bool SRKMacroManager::skipCommentLines(ifstream& inpFileStream)
{

	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;
	char peek = inpFileStream.peek();
	while (peek == '#' || peek == '%' || peek == '\n')
	{
		inpFileStream.ignore(numeric_limits<streamsize>::max(), '\n'); //discards everything in that line
		peek = inpFileStream.peek();
	}

	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return false;

	return true;
}

bool SRKMacroManager::getNonCommentLine(ifstream& inpFileStream, TString& outString)
{
	outString = "";
	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;

	TString sBuffer("#");

	while (sBuffer[0] == '#' || sBuffer[0] == '%')
	{
		sBuffer.ReadLine(inpFileStream);
		if(inpFileStream.fail()) return false;
	}

	outString = sBuffer;
	return true;
}

bool SRKMacroManager::stobool(const string inp)
{
	if(inp == "0" || inp == "false" || inp == "FALSE" || inp == "False")
	{
		return false;
	}
	else if(inp == "1" || inp == "true" || inp == "TRUE" || inp == "True")
	{
		return true;
	}
	else
	{
		cout << "Truth value: " << inp << " not recognized" << endl;
	}
	return false;
}

void SRKMacroManager::sto3double(const string inp, double& x, double& y, double& z)
{
	stringstream aStringStream(inp);
	aStringStream >> x >> y >> z;
	return;

}

TVector3 SRKMacroManager::stoTVector3(const string inp)
{
	double x, y, z;
	sto3double(inp, x, y, z);
	return TVector3(x, y, z);

}

void SRKMacroManager::defineCommands()
{
	commandMap["setRecordAllSteps"] = [&](string inp){	theManager->setRecordAllSteps(stobool(inp));};
	commandMap["setUseAltStepping"] = [&](string inp){	theManager->setUseAltStepping(stobool(inp));};
	commandMap["setParallelFields"] = [&](string inp){	theManager->setParallelFields(stobool(inp));};
	commandMap["setUse2D"] = [&](string inp){	theManager->setUse2D(stobool(inp));};
	commandMap["setConstStepper"] = [&](string inp){	theManager->setConstStepper(stobool(inp));};
	commandMap["setManualTracking"] = [&](string inp){	theManager->setManualTracking(stobool(inp));};
	commandMap["setGyromagneticRatio"] = [&](string inp){	theManager->setGyromagneticRatio(stod(inp));};
	commandMap["setMass"] = [&](string inp){	theManager->setMass(stod(inp));};
	commandMap["setMeanFreePath"] = [&](string inp){	theManager->setMeanFreePath(stod(inp));};
	commandMap["setPhiStart"] = [&](string inp){	theManager->setPhiStart(stod(inp));};
	commandMap["setThetaStart"] = [&](string inp){	theManager->setThetaStart(stod(inp));};
	commandMap["setTimeLimit"] = [&](string inp){	theManager->setTimeLimit(stod(inp));};
	commandMap["setDiffuseReflectionProb"] = [&](string inp){	theManager->setDiffuseReflectionProb(stod(inp));};
	commandMap["setMeanVel"] = [&](string inp){	theManager->setMeanVel(stod(inp));};
	commandMap["setReflectionLimit"] = [&](string inp){	theManager->setReflectionLimit(stod(inp));};
	commandMap["setChamberRadius"] = [&](string inp){	theManager->setChamberRadius(stod(inp));};
	commandMap["setChamberHeight"] = [&](string inp){	theManager->setChamberHeight(stod(inp));};
	commandMap["setChamberRotation"] = [&](string inp){	double x, y, z; sto3double(inp,x, y, z);theManager->setChamberRotation(x,y,z);};
	commandMap["setB0FieldStrength"] = [&](string inp){	theManager->setB0FieldStrength(stod(inp));};
	commandMap["setE0FieldStrength"] = [&](string inp){	theManager->setE0FieldStrength(stod(inp));};
	commandMap["setBGradFieldStrength"] = [&](string inp){	theManager->setBGradFieldStrength(stod(inp));};
	commandMap["setEGradFieldStrength"] = [&](string inp){	theManager->setEGradFieldStrength(stod(inp));};
	commandMap["setDipoleFieldStrength"] = [&](string inp){	theManager->setDipoleFieldStrength(stod(inp));};
	commandMap["setEPSAbs"] = [&](string inp){	theManager->setEPSAbs(stod(inp));};
	commandMap["setEPSRel"] = [&](string inp){	theManager->setEPSRel(stod(inp));};
	commandMap["setInitialStepSize"] = [&](string inp){	theManager->setInitialStepSize(stod(inp));};
	commandMap["setDipolePosition"] = [&](string inp){	theManager->setDipolePosition(stoTVector3(inp));};
	commandMap["setDipoleDirection"] = [&](string inp){	theManager->setDipoleDirection(stoTVector3(inp));};
	commandMap["setE0FieldDirection"] = [&](string inp){	theManager->setE0FieldDirection(stoTVector3(inp));};
	commandMap["setB0FieldDirection"] = [&](string inp){	theManager->setB0FieldDirection(stoTVector3(inp));};
	commandMap["setPos"] = [&](string inp){	theManager->setPos(stoTVector3(inp));};
	commandMap["setVel"] = [&](string inp){	theManager->setVel(stoTVector3(inp));};
	commandMap["setResultsFilePath"] = [&](string inp){	theManager->setResultsFilePath(inp);};
	commandMap["setTrackFilePath"] = [&](string inp){	theManager->setTrackFilePath(inp);};
	commandMap["setRandomSeed"] = [&](string inp){	theManager->setRandomSeed(stoi(inp));};
	commandMap["setDefaultResultsDir"] = [&](string inp){	theManager->setDefaultResultsDir(inp);};
	commandMap["setVelProfHistPath"] = [&](string inp){	theManager->setVelProfHistPath(inp);};
	commandMap["setRunID"] = [&](string inp){	theManager->setRunID(inp);};
	commandMap["makeTracks"] = [&](string inp){	theManager->makeTracks(stoi(inp));};
	commandMap["trackSpins"] = [&](string inp) {	theManager->precessSpinsAlongTracks(stoi(inp));};          // Depricated
	commandMap["precessSpinsAlongTracks"] = [&](string inp){	theManager->precessSpinsAlongTracks(stoi(inp));};
	commandMap["trackSpinsDeltaOmega"] = [&](string inp)  {	theManager->precessSpinsAlongTracksParAndAnti(stoi(inp));};// Depricated
	commandMap["precessSpinsAlongTracksParAndAnti"] = [&](string inp){	theManager->precessSpinsAlongTracksParAndAnti(stoi(inp));};
	commandMap["openMacroFile"] = [&](string inp){	openMacroFile(inp);};
	commandMap["loadParametersFromResultsFile"] = [&](string inp){	theManager->loadParametersFromResultsFile(inp);};

}
