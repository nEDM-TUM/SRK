////////////////////////////////////////////////////////////////
/// SRK - Spin Runge Kutta - MC simulation of spin precession
/// Created to study spin precession in chambered experiments
/// in arbitrary electric and magnetic fields.  Specifically
/// intended to study chambered experiments
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

#include <iostream>

#include "TRandom3.h"

#include "SRKManager.h"
#include "SRKMacroManager.h"

int main(int argc, char* argv[])
{
	//Generate random numbers automatically from lowest 4 bytes of TUUUID
	//If other instances are started within ~100 ns they might have the same seed!
	delete gRandom;
	gRandom = new TRandom3(0);

	std::cout << std::scientific << std::endl;


	//Command Line Execution
	if(argc >2)
	{
		std::cout << "Error! SRK only accepts 0 or 1 arguments." << std::endl;
		return 1;
	}

	SRKManager theManager;  //Monolithic class that runs everything
	SRKMacroManager theMacroManager(&theManager);  //Reads macro files and sends manager to execute commands

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
