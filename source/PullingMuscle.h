/*
 * PullingMuscle.h
 *
 */

#ifndef PULLINGMUSCLE_H_
#define PULLINGMUSCLE_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class PullingMuscle : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _pullingmuscleRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  PullingMuscle(const int argc, const char **argv);
  ~PullingMuscle(){};

  void run();
  void paint(){};
};

#endif /* PULLINGMUSCLE_H_ */
