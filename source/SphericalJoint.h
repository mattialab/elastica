/*
 * SphericalJoint.h
 *
 */

#ifndef SPHERICALJOINT_H_
#define SPHERICALJOINT_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class SphericalJoint : public Test {
 protected:
  vector<double> amp;
  double w, v;
  int NOR;

  vector<REAL> _sphericaljointRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  SphericalJoint(const int argc, const char **argv);
  ~SphericalJoint(){};

  void run();
  void paint(){};
};

#endif /* SPHERICALJOINT_H_ */
