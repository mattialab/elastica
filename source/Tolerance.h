/*
 * Tolerance.h
 *
 *  Created on: Jul 28, 2014
 *      Author: mgazzola
 */

#ifndef TOLERANCE_H_
#define TOLERANCE_H_

#include "ArithmeticPrecision.h"
#include "UsualHeaders.h"

class Tolerance {
 public:
  Tolerance();
  virtual ~Tolerance();

  inline static REAL tol() {
    return 10.0 * std::numeric_limits<REAL>::epsilon();
  }
};

#endif /* TOLERANCE_H_ */
