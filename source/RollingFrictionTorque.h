/*
 * RollingFrictionTorque.h
 *
 *  Created on: May 28, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_ROLLINGFRICTIONTORQUE_H_
#define SOURCE_ROLLINGFRICTIONTORQUE_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class RollingFrictionTorque : public Test
{
protected:
	static bool _test(REAL& tEnergy, REAL& rEnergy, const REAL torque);
	void _sweepTorque(const REAL start, const REAL stop, const REAL step);

public:
	RollingFrictionTorque(const int argc, const char ** argv);
	virtual ~RollingFrictionTorque();

	void run();
	void paint(){};
};

#endif /* SOURCE_ROLLINGFRICTIONTORQUE_H_ */
