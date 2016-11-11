#ifndef INCLUDE_SRKODESTATE_H_
#define INCLUDE_SRKODESTATE_H_

//////////////////////////////////////////////////////////
/// SRKMotionState
/// Setting up a state of particle motion/precession which
/// allows for SRKSpinFloat to be a higher precision number
///////////////////////////////////////////////////////////

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

//Float128 Option
//#include <boost/multiprecision/float128.hpp>
//extern "C"
//{
//#include <quadmath.h>
//}
//typedef boost::multiprecision::float128 SRKSpinFloat;


//Higher floating point option
//#include <boost/multiprecision/cpp_dec_float.hpp>
//typedef boost::multiprecision::cpp_dec_float_50 SRKSpinFloat;

using SRKSpinFloat=double;

#include <array>
#include <vector>

#include "TVector3.h"

#include "SRKMotionState.h"

using SRKODEState= std::vector<SRKSpinFloat>;  ///< Presumes 9 entries: X,Y,Z,Vx,Vy,Vz,spin_phi,spin_theta,time

void setMotionState(SRKODEState& outState, const TVector3* pos, const TVector3* vel, const double phi, const double theta);
void updateMotionStatePosVel(SRKODEState& outState, const SRKMotionState& inpMotionState);
void printMotionState(const SRKODEState& theState);

////////////////////////////////////////////////////////////////
/// struct push_back_state_and_time
/// Observer for integration (records all steps)
/// see http://www.boost.org/doc/libs/1_58_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/short_example.html
////////////////////////////////////////////////////////////////
struct pushBackStateAndTime
{
	std::vector<SRKODEState>* m_states;
	std::vector<double>* m_times;

	pushBackStateAndTime(std::vector<SRKODEState>* states, std::vector<double>* times) :
		m_states(states), m_times(times)
	{
	}

	void operator()(const SRKODEState &x, SRKSpinFloat t)
	{
		m_states->push_back(x);
		m_times->push_back(static_cast<double>(t));
	}
};

#endif /* INCLUDE_SRKODESTATE_H_ */
