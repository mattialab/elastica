#ifndef SPEEDFUNCTIONS_H_
#define SPEEDFUNCTIONS_H_

#include <algorithm>
#include <typeinfo>
#include "ArithmeticPrecision.h"
#include "Matrix3.h"
#include "UsualHeaders.h"
#include "Vector3.h"

using namespace std;

//////////////////// Mattia's new vector functions for speed (no memory
/// allocation) ////////////////////////////

// No need for inline functions here (I checked - no gain)

void vLength(const vector<Vector3> &A, vector<REAL> &B);
void vUnitize(const vector<Vector3> &A, vector<Vector3> &B);
void vDelta(const vector<Vector3> &A, vector<Vector3> &B);
void vSqrt(const vector<REAL> &A, vector<REAL> &B);
void vDiff(const vector<Vector3> &A, vector<Vector3> &B);
void vRotDiff(const vector<Matrix3> &A, vector<Matrix3> &B);
void vLog(const vector<Matrix3> &A, vector<Vector3> &B);
void vExp(const vector<Vector3> &A, vector<Matrix3> &B,
          const REAL tol = Tolerance::tol());
void vMidAvg(const vector<Vector3> &A, vector<Vector3> &B);
void vMidAvgInterior(const vector<Vector3> &A, vector<Vector3> &B);
void vMidAvgInterior(const vector<REAL> &A, vector<REAL> &B);
void vFromPointsToElements(const vector<Vector3> &A, vector<Vector3> &B);
void vFromPointsToElements(const vector<REAL> &A, vector<REAL> &B);
void vT(const vector<Matrix3> &A, vector<Matrix3> &B);
void vDiagI(const vector<Matrix3> &A, vector<Matrix3> &B);
void v_timesequal(const vector<Matrix3> &A, vector<Matrix3> &B);  // a = b * a
void v_a_plus_b_equal_c(const vector<Vector3> &A, const vector<Vector3> &B,
                        vector<Vector3> &C);  // c = a + b
void v_plusequal_a_times_b(const REAL A, const vector<Vector3> &B,
                           vector<Vector3> &C);  // c += a * b
void v_a_plus_b_plus_c_equal_d(const vector<Vector3> &A,
                               const vector<Vector3> &B,
                               const vector<Vector3> &C,
                               vector<Vector3> &D);  // d = a + b + c
void v_a_minus_b_equal_c(const vector<Vector3> &A, const vector<Vector3> &B,
                         vector<Vector3> &C);  // c = a - b
void v_a_times_b_equal_c(const vector<Vector3> &A, const REAL B,
                         vector<Vector3> &C);  // c = a * b
void v_a_times_b_equal_c(const vector<REAL> &A, const REAL B,
                         vector<REAL> &C);  // c = a * b
void v_a_times_b_equal_c(const vector<Matrix3> &A, const vector<Vector3> &B,
                         vector<Vector3> &C);  // c = a * b
void v_a_cross_b_equal_c(const vector<Vector3> &A, const vector<Vector3> &B,
                         vector<Vector3> &C);  // c = a x b
void v_a_times_b_times_c_equal_d(const REAL A, const vector<REAL> &B,
                                 const vector<Vector3> &C,
                                 vector<Vector3> &D);  // d = a * b * c
void v_a_times_b_cross_c_equal_d(const vector<REAL> &A,
                                 const vector<Vector3> &B,
                                 const vector<Vector3> &C,
                                 vector<Vector3> &D);  // d = a * b x c
void v_a_times_b_cross_c_equal_d(const REAL A, const vector<Vector3> &B,
                                 const vector<Vector3> &C,
                                 vector<Vector3> &D);  // d = a * b x c
void v_a_divide_b_equal_c(const vector<Vector3> &A, const vector<REAL> &B,
                          vector<Vector3> &C);  // c = a / b
void v_a_dividepar_b_times_c_equal_d(const vector<REAL> &A,
                                     const vector<REAL> &B, const REAL C,
                                     vector<REAL> &D);  // d = a / (b*c)
void v_sqrt_a_dividepar_b_times_c_equal_d(
    const vector<REAL> &A, const vector<REAL> &B, const REAL C,
    vector<REAL> &D);  // d = sqrt( a / (b*c) )
void v_a_times_b_divide_c_minus_d_equal_c(
    const vector<Matrix3> &A, const vector<Vector3> &B, const vector<REAL> &C,
    const Vector3 &D, vector<Vector3> &E);  // e = (a * b / c) - d

#endif
