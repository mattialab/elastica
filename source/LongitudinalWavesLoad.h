/*
 * LongitudinalWavesLoad.h
 *
 *  Created on: May 14, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_LONGITUDINALWAVESLOAD_H_
#define SOURCE_LONGITUDINALWAVESLOAD_H_

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

class LongitudinalWavesLoad : public Test
{
protected:
	static bool _test(	const int nEdges, const REAL _dt, const REAL _L, const REAL _r, const REAL _f, const REAL _F, const REAL _timeSimulation, const REAL _E, const REAL _rho,
						const REAL _nu, const REAL _relaxationNu, const string outfileName);
	static void _longWaveTest(const int nEdges, const string outputdata);

public:
	LongitudinalWavesLoad(const int argc, const char ** argv);
	virtual ~LongitudinalWavesLoad();

	void run();
	void paint(){};
};

#endif /* SOURCE_LONGITUDINALWAVESLOAD_H_ */
