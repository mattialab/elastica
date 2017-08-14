/*
 * AnisotropicFriction.h
 *
 *  Created on: Jul 8, 2015
 *      Author: mgazzola
 */

#ifndef SOURCE_EELSLENDERBODY_H_
#define SOURCE_EELSLENDERBODY_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class SlenderBodyStokes: public Test
{
protected:
	vector <double> amp;
	REAL w;
	REAL v;

	REAL ncycles;
	unsigned int framesPerUnitTime;

	MRAG::ArgumentParser parser;
	InterfaceCma interfaceCma;

	REAL _snakeRun();

#ifdef SNAKE_VIZ
	void _paint(){};
#endif

public:

	SlenderBodyStokes(const int argc, const char ** argv);
	~SlenderBodyStokes(){};

	void run();
	void paint(){};
};


#endif /* SOURCE_ANISOTROPICFRICTION_H_ */
