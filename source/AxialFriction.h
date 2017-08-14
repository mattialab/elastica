/*
 * AxialFriction.h
 *
 *  Created on: Sep 17, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_AXIALFRICTION_H_
#define SOURCE_AXIALFRICTION_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class AxialFriction : public Test
{
protected:
	static bool _test(REAL& tEnergy, REAL& rEnergy, const REAL angle);
	void _sweepAngle(const REAL start, const REAL stop, const REAL step);

public:
	AxialFriction(const int argc, const char ** argv);
	virtual ~AxialFriction();

	void run();
	void paint(){};
};

#endif /* SOURCE_AXIALFRICTION_H_ */
