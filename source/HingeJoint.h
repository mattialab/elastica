/*
 * SphericalJoint.h
 *
 */

#ifndef HINGEJOINT_H_
#define HINGEJOINT_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class HingeJoint : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _hingejointRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  HingeJoint(const int argc, const char **argv);
  ~HingeJoint(){};

  void run();
  void paint(){};
};

#endif /* HINGEJOINT_H_ */
