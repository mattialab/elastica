/*
 * Flagella.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef FLAGELLA_H_
#define FLAGELLA_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class Flagella : public Test {
 protected:
  vector<double> amp;
  double rsmall, rlarge;
  double head;
  double cell1;
  double cell2;
  double ncycles;
  int NOR;

  REAL _flagellaRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  Flagella(const int argc, const char **argv);
  ~Flagella(){};

  void run();
  void paint(){};
};

#endif /* FLAGELLA_H_ */
