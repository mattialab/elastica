#ifndef GEOMETRYFUNCTIONS_H_
#define GEOMETRYFUNCTIONS_H_

#include <limits>
#include "Matrix3.h"
#include "Rod.h"
#include "Tolerance.h"
#include "UsualHeaders.h"
#include "Vector3.h"
#include "VectorFunctions.h"

std::vector<Matrix3> alignFrames(const std::vector<Vector3> &x,
                                 const Vector3 &u0);

void applyTwists(std::vector<Matrix3> &Q, const std::vector<REAL> &twists);

std::vector<Vector3> findMinDistVectors(const Vector3 &x1, const Vector3 &e1,
                                        const Vector3 &x2, const Vector3 &e2,
                                        const REAL tol = 1e-6);

#endif
