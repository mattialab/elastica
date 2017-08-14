/*
 * Solenoids.h
 *
 *  Created on: Nov 23, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_SOLENOIDSJCP_H_
#define SOURCE_SOLENOIDSJCP_H_

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

class SolenoidsJCP : public Test
{
protected:
	static bool _test(const int nEdges, const REAL _E, const REAL _R, const REAL _F, const REAL _timeTwist, const REAL _timeRelax, const REAL _nu, string outfileName);

public:
	SolenoidsJCP(const int argc, const char ** argv);
	virtual ~SolenoidsJCP();

	void run();
	void paint(){};
};


#endif /* SOURCE_SOLENOIDSJCP_H_ */
