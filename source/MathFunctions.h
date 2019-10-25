#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include "ArithmeticPrecision.h"
#include "UsualHeaders.h"

// Returns fmod shifted into [0,r]
REAL posMod(const REAL a, const REAL r);

// Calls the normal asin func, but if x>1 returns PI/2, and if x<-1 returns
// -PI/2, so as to get around precision errors making asin(sin(PI/2)) =
// nan
REAL arcSin(const REAL x);

// Calls the normal acos func, but if x>1 returns 0, and if x<-1 returns PI,
// so as to get around precision errors making acos(cos(0)) = nan
REAL arcCos(const REAL x);

// Gaussian random normal generator without trigonometric calls using Box-Muller
// transform
double randn_notrig(const double mu = 0.0, const double sigma = 1.0);

#endif
