/*
 * Snake.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef SNAKE_H_
#define SNAKE_H_

#include "ArgumentParser.h"
#include "Polymer.h"
#include "PositionVerlet2nd.h"
#include "RodInitialConfigurations.h"
#include "Test.h"
#include "UsualHeaders.h"

using namespace std;

class Snake : public Test {
 protected:
  vector<double> amp;
  double w, v;

  REAL _snakeRun();

#ifdef SNAKE_VIZ
  void _paint(){};
#endif

 public:
  Snake(const int argc, const char **argv);
  ~Snake(){};

  void run();
  void paint(){};
};

#endif /* SNAKE_H_ */