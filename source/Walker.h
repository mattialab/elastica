/*
 * Walker.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef WALKER_H_
#define WALKER_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class Walker : public Test {
 protected:
  double YoungM;
  double LegL;
  double Location;
  int NOR;

  REAL _walkerRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  Walker(const int argc, const char **argv);
  ~Walker(){};

  void run();
  void paint(){};
};

#endif /* WALKER_H_ */
