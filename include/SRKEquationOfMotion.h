#ifndef SRKEQUATIONOFMOTION_H_
#define SRKEQUATIONOFMOTION_H_

#include <vector>

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <boost/multiprecision/float128.hpp>
//extern "C"
//{
//#include <quadmath.h>
//}
//typedef boost::multiprecision::float128 SRKSpinFloat;

//#include <boost/multiprecision/cpp_dec_float.hpp>
//typedef boost::multiprecision::cpp_dec_float_50 SRKSpinFloat;

typedef double SRKSpinFloat;
//typedef boost::multiprecision::float128 SRKSpinFloat;

#include "TVector3.h"

#include "SRKGlobalField.h"

typedef std::vector<SRKSpinFloat> SRKMotionState;  //Presumes 9 entries: X,Y,Z,Vx,Vy,Vz,spin_phi,spin_theta,time   time not used for runge kutta but for recording state
//typedef std::vector<SRKSpinFloat> SRKInternalMotionState;  //Uses higher floating point numbers.  Presumes 9 entries: X,Y,Z,Vx,Vy,Vz,spin_phi,spin_theta,time   time not used for runge kutta but for recording state

class SRKEquationOfMotion
{
public:
	SRKEquationOfMotion(SRKGlobalField* inpGlobalField);
	virtual ~SRKEquationOfMotion();
	void operator()(const SRKMotionState& x, SRKMotionState& dxdt, const SRKSpinFloat /* t */); //THE equation of motion for use with BOOST odeint
	void loadFields(bool EBParallel);

	inline void setGyromagneticRatio(double inp){ gyromagneticRatio = SRKSpinFloat(inp); }
	inline double getGyromagneticRatio(){ return static_cast<double> (gyromagneticRatio); }

protected:
	SRKGlobalField* theGlobalField;

	//For retrieval of field
	double fieldDouble[9];
	double posDouble[3];

	SRKSpinFloat gyromagneticRatio;
	SRKSpinFloat theField[9];

};

void setMotionState(SRKMotionState& outState, TVector3* pos, TVector3* vel, double phi, double theta);
void updateMotionStatePosVel(SRKMotionState& outState, TVector3 pos, TVector3 vel, double currentTime);
void printMotionState(const SRKMotionState& theState);
void printMotionState(const SRKMotionState& theState);

//Observer for integration (records all steps)
//see http://www.boost.org/doc/libs/1_58_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/short_example.html
struct push_back_state_and_time
{
	std::vector<SRKMotionState>* m_states;
	std::vector<double>* m_times;

	push_back_state_and_time(std::vector<SRKMotionState>* states, std::vector<double>* times) :
		m_states(states), m_times(times)
	{
	}

	void operator()(const SRKMotionState &x, SRKSpinFloat t)
	{
		m_states->push_back(x);
		m_times->push_back(static_cast<double>(t));
	}
};
#endif /* SRKEQUATIONOFMOTION_H_ */
