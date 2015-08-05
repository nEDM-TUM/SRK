#include "SRKMacroManager.h"
#include <iostream>
#include <fstream>
using namespace std;

SRKMacroManager::SRKMacroManager(SRKManager* inpManager)
{
	theManager = inpManager;

}

SRKMacroManager::~SRKMacroManager()
{

}

bool SRKMacroManager::openMacroFile(TString filePath)
{
	ifstream theFile(filePath);
	TString sBuffer("#");
	if(!theFile.is_open() || theFile.fail() || theFile.eof()) return false;

	while (!theFile.fail() || !theFile.eof())
	{
		skipCommentLines(theFile);
		if(theFile.eof()) return true;
		TString commandString, commandValueString;
		theFile >> commandString >> commandValueString;
		commandStringList.push(commandString);
		commandValueStringList.push(commandValueString);
	}

	return true;
}

void SRKMacroManager::doMacroCommands()
{
	while (!commandStringList.empty())
	{
		cout << "SRK: " << commandStringList.front() << " " << commandValueStringList.front() << endl;

		commandStringList.pop();
		commandValueStringList.pop();

	}
}

bool SRKMacroManager::skipCommentLines(ifstream& inpFileStream)
{
	TString sBuffer;
	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;
	while (inpFileStream.peek() == '#' || inpFileStream.peek() == '%')
	{
		sBuffer.ReadLine(inpFileStream);
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

bool SRKMacroManager::stobool(string inp)
{
	if(inp=="0" || inp=="false" || inp=="FALSE" || inp=="False")
	{
		return false;
	}
	else if(inp=="1" || inp=="true" || inp=="TRUE" || inp=="True")
	{
		return true;
	}
	else
	{
		cout << "Truth value: " << inp << " not recognized" << endl;
	}
	return false;
}

void SRKMacroManager::defineCommands()
{
	commandMap["recordAllSteps"]=[this](string inp){theManager->setRecordAllSteps(stobool(inp));};
	commandMap["recordAllSteps"]=[this](string inp){theManager->setRecordAllSteps(stobool(inp));};
}
