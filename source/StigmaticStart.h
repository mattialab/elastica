/*
 * SelfAssembly.h
 *
 *  Created on: Dec 16, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_STIGMATICSTART_H_
#define SOURCE_STIGMATICSTART_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class StigmaticStart: public Test
{
protected:
	MRAG::ArgumentParser parser;
	InterfaceCma interfaceCma;
	vector<double> paramtersCMA;

	REAL _snakeRun(const REAL _W, const REAL _M, const REAL _Froude, const REAL _L, const REAL _T, const REAL _d, const REAL _E, const REAL _rho);

#ifdef SNAKE_VIZ
	void _paint(){};
#endif

public:

	StigmaticStart(const int argc, const char ** argv);
	~StigmaticStart(){};

	void run();
	void paint(){};
};


#endif /* SOURCE_STIGMATICSTART_H_ */
