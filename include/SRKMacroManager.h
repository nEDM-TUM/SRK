#ifndef SRKMACROMANAGER_H_
#define SRKMACROMANAGER_H_

#include <list>
#include <functional>
#include <unordered_map>
#include <string>

#include "TString.h"

#include "SRKManager.h"

class SRKMacroManager
{
public:
	SRKMacroManager(SRKManager* inpManager);
	virtual ~SRKMacroManager();
	bool openMacroFile(TString filePath);
	void runMacroCommands(); //From commmandString/ValueList
	bool runMacroCommand(std::string command, std::string value);
	void enterInteractiveMode();
private:

	void defineCommands();

	bool getNonCommentLine(ifstream& inpFileStream, TString& outString);
	bool skipCommentLines(ifstream& inpFileStream);
	bool stobool(std::string inp);
	TVector3 stoTVector3(std::string inp);

	std::list<TString> commandStringList;
	std::list<TString> commandValueStringList;

	std::unordered_map<std::string,std::function<void(std::string)>> commandMap;  //Map of command strings to the commands (using lambda functions)
	SRKManager* theManager;
};

#endif /* SRKMACROMANAGER_H_ */
