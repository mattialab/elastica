/*
 * ValidationEulerBuckling.h
 *
 *  Created on: Jul 28, 2014
 *      Author: mgazzola
 */

#ifndef VALIDATIONEULERBUCKLING_H_
#define VALIDATIONEULERBUCKLING_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class EulerBuckling : public Test
{
protected:
	static bool _testEulerBuckling(const REAL force, const REAL lengthRod, const REAL alpha, const REAL timeSimulation, const REAL bendLimit, REAL &bendingEnergy);
	void _sweepEulerBucklingForceLength();
	void _sweepEulerBucklingForceAlpha();

public:
	EulerBuckling(const int argc, const char ** argv);
	virtual ~EulerBuckling();

	void run();
	void paint(){};
};

#endif /* VALIDATIONEULERBUCKLING_H_ */
