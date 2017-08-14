/*
 * LongitudinalWaves.h
 *
 *  Created on: Mar 30, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_LONGITUDINALWAVES_H_
#define SOURCE_LONGITUDINALWAVES_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodBoundaryConditions.h"
#include "MRAGProfiler.h"
#include "Tolerance.h"

using namespace std;

class LongitudinalWaves : public Test
{
protected:
	static bool _test(	const int nEdges, const REAL _dt, const REAL _L, const REAL _f, const REAL _amplitude, const REAL _timeSimulation, const REAL _E, const REAL _rho,
						const REAL _nu, const REAL _relaxationNu, const string outfileName);
	static void _longWaveTest(const int nEdges, const string outputdata);

public:
	LongitudinalWaves(const int argc, const char ** argv);
	virtual ~LongitudinalWaves();

	void run();
	void paint(){};
};

#endif /* SOURCE_LONGITUDINALWAVES_H_ */
