/*
 * MuscularSnake.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef WING_H_
#define WING_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class Wing : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _wingRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  Wing(const int argc, const char **argv);
  ~Wing(){};

  void run();
  void paint(){};
};

#endif /* WING_H_ */
