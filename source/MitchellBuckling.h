/*
 * ValidationMitchellBuckling.h
 *
 *  Created on: Aug 4, 2014
 *      Author: mgazzola
 */

#ifndef VALIDATIONMITCHELLBUCKLING_H_
#define VALIDATIONMITCHELLBUCKLING_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "VelocityVerlet2nd.h"
#include "PositionVerlet2nd.h"

using namespace std;

class MitchellBuckling : public Test
{
protected:
	static bool _testMitchellBuckling(const REAL alpha, const REAL beta, const REAL twist, const REAL L, const REAL timeSimulation, const REAL translationalLimit, REAL &translationalEnergy);
	void _sweepMitchellBuckling(const REAL alpha, const REAL L, string outfileName);

public:
	MitchellBuckling(const int argc, const char ** argv);
	virtual ~MitchellBuckling();

	void run();
	void paint(){};
};

#endif /* VALIDATIONMITCHELLBUCKLING_H_ */
