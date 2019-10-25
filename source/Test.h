/*
 * Test.h
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#ifndef TEST_H_
#define TEST_H_

#include "ExternalContact.h"
#include "Interaction.h"
#include "MRAGProfiler.h"
#include "Matrix3.h"
#include "Rod.h"
#include "RodBoundaryConditions.h"
#include "RodExternalForces.h"
#include "UsualHeaders.h"
#include "Vector3.h"

class Test {
 protected:
  typedef std::vector<Vector3> VV3;
  typedef std::vector<Matrix3> VM3;
  typedef std::vector<REAL> VREAL;
  typedef std::vector<bool> VBOOL;
  typedef std::vector<int> VINT;
  typedef std::vector<RodBC *> Vbcptr;
  typedef std::vector<ExternalForces *> Vefptr;
  typedef std::vector<Rod *> Vrodptr;
  typedef std::vector<Interaction *> Vinterptr;

 public:
  Test(){};
  virtual ~Test(){};

  virtual void run() = 0;
  virtual void paint() = 0;
};

#endif /* TEST_H_ */
