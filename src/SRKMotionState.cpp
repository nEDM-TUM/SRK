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

	cout << "-------------------------------------------" << endl;
	cout << "Pos: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << endl;
	cout << "Vel: " << vel.X() << ", " << vel.Y() << ", " << vel.Z() << endl;
	cout << "Time: " << time << endl;
	cout << "-------------------------------------------" << endl;

}
