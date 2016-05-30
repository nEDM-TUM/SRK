#ifndef SRKMACROMANAGER_H_
#define SRKMACROMANAGER_H_

#include <list>
#include <functional>
#include <unordered_map>
#include <string>
#include <iostream>

#include "TString.h"

#include "SRKManager.h"

////////////////////////////////////////////////////////////////
/// SRKMacroManager
/// Given macro commands (either directly or from a file)
/// will tell SRKManager to execute the commands
///
/// Commands are implemented as lambda functions in defineCommands()
/// and stored in an unordered_map
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

class SRKMacroManager
{
public:
	SRKMacroManager(SRKManager* inpManager);
	~SRKMacroManager();
	bool openMacroFile(const TString filePath); /// Open macro file and stores commands and values in lists
	void runMacroCommands(); /// Runs commands from list
	bool runMacroCommand(const std::string command, const std::string value); /// Runs a single command
	void enterInteractiveMode();  /// Using standard output, enter commands
private:

	void defineCommands(); /// Defines all commands

	bool getNonCommentLine(std::ifstream& inpFileStream, TString& outString); ///  Gets the next noncomment line from a string.  Comments begin with "#".
	bool skipCommentLines(std::ifstream& inpFileStream); ///  Skips next comment line from a string.  Comments begin with "#".
	bool stobool(const std::string inp);  /// Converts a string to a bool (several options given)
	void sto3double(const std::string inp, double& x, double& y, double& z); /// Converts a string containing three numbers separated by whitespace to doubles
	TVector3 stoTVector3(const std::string inp);  /// Converts a string containing three numbers separated by whitespace to a TVector3

	std::list<TString> commandStringList;  //List of command strings
	std::list<TString> commandValueStringList; //List of command values

	std::unordered_map<std::string, std::function<void(std::string)>> commandMap;  //Map of command strings to the commands (using lambda functions)
	SRKManager* theManager;  //Pointer to SRKManager which runs the commands
};

#endif /* SRKMACROMANAGER_H_ */
