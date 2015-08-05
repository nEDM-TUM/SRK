#ifndef SRKMACROMANAGER_H_
#define SRKMACROMANAGER_H_

#include <queue>
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
	void doMacroCommands();
private:

	void defineCommands();

	bool getNonCommentLine(ifstream& inpFileStream, TString& outString);
	bool skipCommentLines(ifstream& inpFileStream);
	bool stobool(std::string inp);

	std::queue<TString> commandStringList;
	std::queue<TString> commandValueStringList;

	std::unordered_map<std::string,std::function<void(std::string)>> commandMap;  //Map of command strings to the commands (using lambda functions)
	SRKManager* theManager;
};

#endif /* SRKMACROMANAGER_H_ */
