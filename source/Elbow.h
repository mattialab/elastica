/*
 * MuscularSnake.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef ELBOW_H_
#define ELBOW_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class Elbow : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _elbowRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  Elbow(const int argc, const char **argv);
  ~Elbow(){};

  void run();
  void paint(){};
};

#endif /* ELBOW_H_ */
