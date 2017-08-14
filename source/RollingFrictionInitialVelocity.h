/*
 * RollingFriction.h
 *
 *  Created on: May 19, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_ROLLINGFRICTIONINITIALVELOCITY_H_
#define SOURCE_ROLLINGFRICTIONINITIALVELOCITY_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class RollingFrictionInitialVelocity : public Test
{
protected:
	static bool _test(REAL& tEnergy, REAL& rEnergy, REAL IFactor);
	void _sweepInertia(const REAL start, const REAL stop, const REAL step);

public:
	RollingFrictionInitialVelocity(const int argc, const char ** argv);
	virtual ~RollingFrictionInitialVelocity();

	void run();
	void paint(){};
};

#endif /* SOURCE_ROLLINGFRICTIONINITIALVELOCITY_H_ */
