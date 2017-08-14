/*
 * Snake.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef SNAKE_H_
#define SNAKE_H_

#include "UsualHeaders.h"
#include "Test.h"
#include "ArgumentParser.h"
#include "InterfaceCma.h"
#include "RodInitialConfigurations.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"

using namespace std;

class Snake: public Test
{
protected:
	vector <double> amp;
	double w, v;
	double ncycles;
	unsigned int framesPerUnitTime;

	MRAG::ArgumentParser parser;
	InterfaceCma interfaceCma;

	REAL _snakeRun();

#ifdef SNAKE_VIZ
	void _paint(){};
#endif

public:

	Snake(const int argc, const char ** argv);
	~Snake(){};

	void run();
	void paint(){};
};



#endif /* SNAKE_H_ */
