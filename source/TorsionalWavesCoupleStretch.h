/*
 * TorsionalWavesCoupleStretch.h
 *
 *  Created on: May 15, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_TORSIONALWAVESCOUPLESTRETCH_H_
#define SOURCE_TORSIONALWAVESCOUPLESTRETCH_H_

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

class TorsionalWavesCoupleStretch : public Test
{
protected:
	static bool _test(	const int nEdges, const REAL _dt, const REAL _L, const REAL _r, const REAL _f, const REAL _M, const REAL _timeSimulation, const REAL _G, const REAL _rho,
						const REAL _nu, const REAL _relaxationNu, const string outfileName);
	static void _longWaveTest(const int nEdges, const string outputdata);

public:
	TorsionalWavesCoupleStretch(const int argc, const char ** argv);
	virtual ~TorsionalWavesCoupleStretch();

	void run();
	void paint(){};
};

#endif /* SOURCE_TORSIONALWAVESCOUPLESTRETCH_H_ */
