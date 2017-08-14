/*
 * BendingWavesCouple.h
 *
 *  Created on: Jun 21, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_EULERBEAM_H_
#define SOURCE_EULERBEAM_H_

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

class EulerBeam : public Test
{
protected:
	static bool _cantileverPointLoad( const string outfileName );
	static bool _cantileverPointTorque( const string outfileName );
	static bool _cantileverTorqueAtFreeEnds( const string outfileName );
	static bool _cantileverConstantMuscularActivity( const string outfileName );

public:
	EulerBeam(const int argc, const char ** argv);
	virtual ~EulerBeam();

	void run();
	void paint(){};
};

#endif /* SOURCE_QUASISTATICTIMOSHENKOBEAM_H_ */
