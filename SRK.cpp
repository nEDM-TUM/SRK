#include <iostream>

#include "TRandom3.h"

#include "SRKManager.h"
#include "SRKMacroManager.h"

using namespace std;

int main(int argc, char* argv[])
{
	delete gRandom;
	gRandom = new TRandom3(0);

	/*Command Line Execution*/
	if(argc >2)
	{
		cout << "Error! SRK only accepts 0 or 1 arguments." << endl;
		return 1;
	}

	SRKManager theManager;
	SRKMacroManager theMacroManager(&theManager);

	if(argc ==2) //argument is the macro file path
	{
		TString macroFilePath = argv[1];
		theMacroManager.openMacroFile(macroFilePath);
		theMacroManager.runMacroCommands();
	}
	else
	{
		theMacroManager.enterInteractiveMode();
	}

	return 0;
}
