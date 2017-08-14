/*
 * Knot.h
 *
 *  Created on: Nov 28, 2015
 *      Author: ldudte
 */

#ifndef SOURCE_KNOT_H_
#define SOURCE_KNOT_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "VelocityVerlet2nd.h"
#include "PositionVerlet2nd.h"
#include "PositionVerletGeneral.h"
#include "AdaptiveVelocityVerlet.h"
#include "MRAGProfiler.h"
#include "Tolerance.h"

using namespace std;

class Knot : public Test
{
protected:
	static bool _testLocalizedHelicalBuckling(const int nEdges, const REAL _timeTwist, const REAL _timeRelax, const REAL _nu, string outfileName);
	static void _test(const unsigned int nEdges, string outputfilename);

public:
	Knot(const int argc, const char ** argv);
	virtual ~Knot();

	void run();
	void paint(){};
};


#endif /* SOURCE_KNOT_H_ */
