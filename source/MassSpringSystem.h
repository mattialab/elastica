/*
 * StretchRelease.h
 *
 *  Created on: May 7, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_MASSSPRINGSYSTEM_H_
#define SOURCE_MASSSPRINGSYSTEM_H_

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

class MassSpringSystem : public Test
{
protected:
	static bool _test(	const int nEdges, const REAL _dt, const REAL _M, const REAL _L0, const REAL _A0, const REAL _timeSimulation,
						const REAL _E, const REAL _rho, const string outfileName);
	static void _largeAttachedMassTest(const int nEdges, const REAL E, const string outputdata);
	static void _smallAttachedMassTest(const int nEdges, const REAL E, const string outputdata);

public:
	MassSpringSystem(const int argc, const char ** argv);
	virtual ~MassSpringSystem();

	void run();
	void paint(){};
};




#endif /* SOURCE_MASSSPRINGSYSTEM_H_ */
