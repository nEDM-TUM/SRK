#ifndef SRKMOTIONSTATE_H_
#define SRKMOTIONSTATE_H_

#include "TVector3.h"

/// Defines why the tracker ended the step
enum class SRKStepPointType : char { UNDEFINED = 0, START=1, REFLECTION = 2, MAXTRACKTIME = 3, TIMELIMIT = 4, ABSORBED=5, DEPOLARIZED=6, GASSCATTER=7};

class SRKMotionState
{
public:
	SRKMotionState();
	virtual ~SRKMotionState();

	TVector3 pos;
	TVector3 vel;
	double time=0.;
	SRKStepPointType type=SRKStepPointType::UNDEFINED;

	void print() const;

};

#endif /* SRKMOTIONSTATE_H_ */
