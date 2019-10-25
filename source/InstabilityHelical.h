/*
 * InstabilityHelical.h
 *
 *  Created on: Sep 8, 2014
 *      Author: mgazzola
 */

#ifndef INSTABILITYHELICAL_H_
#define INSTABILITYHELICAL_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodBoundaryConditions.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "Tolerance.h"
#include "UsualHeaders.h"

using namespace std;

class InstabilityHelical : public Test {
 protected:
  static bool _testLocalizedHelicalBuckling(const int nEdges,
                                            const REAL _timeTwist,
                                            const REAL _timeRelax,
                                            const REAL _nu, string outfileName);
  static void _test(const unsigned int nEdges, string outputfilename);

 public:
  InstabilityHelical(const int argc, const char** argv);
  virtual ~InstabilityHelical();

  void run();
  void paint(){};
};

#endif /* INSTABILITYHELICAL_H_ */
