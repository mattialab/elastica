/*
 * MuscularSnake.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef MUSCULARSNAKE_H_
#define MUSCULARSNAKE_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class MuscularSnake : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _muscularsnakeRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  MuscularSnake(const int argc, const char **argv);
  ~MuscularSnake(){};

  void run();
  void paint(){};
};

#endif /* MUSCULARSNAKE_H_ */
