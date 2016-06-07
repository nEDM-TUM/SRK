#include "SRKMotionState.h"

#include <iostream>

using namespace std;

SRKMotionState::SRKMotionState()
{
	// TODO Auto-generated constructor stub

}

SRKMotionState::~SRKMotionState()
{
	// TODO Auto-generated destructor stub
}

void SRKMotionState::print() const
{
	cout << scientific << endl;
	cout << "-------------------------------------------" << endl;
	cout << "Pos: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << endl;
	cout << "Vel: " << vel.X() << ", " << vel.Y() << ", " << vel.Z() << endl;
	cout << "Time: " << time << ", " << SRKStepPointTypeToTString(type) << endl;
	cout << "-------------------------------------------" << endl;

}

TString SRKStepPointTypeToTString(SRKStepPointType inp)
{
	switch( inp )
	{
		case SRKStepPointType::UNDEFINED:
			return "UNDEFINED";
			break;

		case SRKStepPointType::START:
			return "START";
			break;

		case SRKStepPointType::REFLECTION:
			return "REFLECTION";
			break;

		case SRKStepPointType::PERIODICSTOP:
			return "PERIODICSTOP";
			break;

		case SRKStepPointType::TIMELIMIT:
			return "TIMELIMIT";
			break;

		case SRKStepPointType::ABSORBED:
			return "ABSORBED";
			break;

		case SRKStepPointType::DEPOLARIZED:
			return "DEPOLARIZED";
			break;

		case SRKStepPointType::GASSCATTER:
			return "GASSCATTER";
			break;

	}
	return "";
}
