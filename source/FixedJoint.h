/*
 * FixedJoint.h
 *
 */

#ifndef FIXEDJOINT_H_
#define FIXEDJOINT_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class FixedJoint : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _fixedjointRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  FixedJoint(const int argc, const char **argv);
  ~FixedJoint(){};

  void run();
  void paint(){};
};

#endif /* FIXEDJOINT_H_ */
