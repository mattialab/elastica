/*
 * RollingFrictionInclinedPlane.h
 *
 *  Created on: May 25, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_ROLLINGFRICTIONINCLINEDPLANE_H_
#define SOURCE_ROLLINGFRICTIONINCLINEDPLANE_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class RollingFrictionInclinedPlane : public Test
{
protected:
	static bool _test(REAL& tEnergy, REAL& rEnergy, const REAL angle);
	void _sweepAngle(const REAL start, const REAL stop, const REAL step);

public:
	RollingFrictionInclinedPlane(const int argc, const char ** argv);
	virtual ~RollingFrictionInclinedPlane();

	void run();
	void paint(){};
};

#endif /* SOURCE_ROLLINGFRICTIONINCLINEDPLANE_H_ */
