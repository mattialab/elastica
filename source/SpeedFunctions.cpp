/*
 * VectorFunctions.cpp
 *
 *  Created on: Feb 27, 2015
 *      Author: mgazzola
 */

#include "SpeedFunctions.h"

void vLength(const vector<Vector3> &A, vector<REAL> &B) {
  assert(A.size() == B.size());
  assert(B.size() != 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<REAL>::iterator iterB = B.begin();
  ;
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterA++)
    (*iterB) =
        (*iterA)
            .length();  // length returns a real, no need for any in-place thing
}

void vUnitize(const vector<Vector3> &A, vector<Vector3> &B) {
  assert(A.size() == B.size());
  assert(B.size() != 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Vector3>::iterator iterB = B.begin();
  ;
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterA++)
    (*iterB) = (*iterA).unitize();  // length returns a real, no need for any
                                    // in-place thing
}

void vSqrt(const vector<REAL> &A, vector<REAL> &B) {
  assert(A.size() == B.size());
  assert(B.size() != 0);
  vector<REAL>::const_iterator iterA = A.begin();
  vector<REAL>::iterator iterB = B.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterA++)
    (*iterB) =
        sqrt((*iterA));  // operation with reals, no need for any in-place thing
}

void vDiff(const vector<Vector3> &A, vector<Vector3> &B) {
  assert(B.size() == A.size() - 1);
  assert(A.size() >= 1);
  vector<Vector3>::iterator iterB;
  vector<Vector3>::const_iterator iterA = A.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++)
    Vector3::a_minus_b_equal_c(
        (*++iterA), (*iterA),
        (*iterB));  // in-place operator, no memory allocation
}

void vDelta(const vector<Vector3> &A, vector<Vector3> &B) {
  assert(B.size() == A.size() + 1);
  assert(A.size() >= 1);
  const unsigned int n = A.size();
  B[0] = A[0];
  for (unsigned int i = 1; i < n; i++)
    Vector3::a_minus_b_equal_c(
        A[i], A[i - 1],
        B[i]);  // in-place operator, no memory allocation
  B.back() = -A[n - 1];
}

void vRotDiff(const vector<Matrix3> &A, vector<Matrix3> &B) {
  assert(A.size() - 1 == B.size());
  assert(A.size() > 1);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Matrix3>::iterator iterB = B.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++)
    Matrix3::matrix_times_transposed(
        (*++iterA), (*iterA),
        (*iterB));  // in-place operator, no memory allocation
}

void vLog(const vector<Matrix3> &A, vector<Vector3> &B) {
  assert(A.size() == B.size());
  assert(B.size() >= 1);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Vector3>::iterator iterB = B.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++)
    (*iterA).log((*iterB));  // in-place
}

void vExp(const vector<Vector3> &A, vector<Matrix3> &B, const REAL tol) {
  assert(A.size() == B.size());
  assert(B.size() >= 1);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Matrix3>::iterator iterB = B.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++)
    exp((*iterA), (*iterB));  // in-place
}

void vMidAvg(const vector<Vector3> &A, vector<Vector3> &B) {
  assert(A.size() + 1 == B.size());
  assert(B.size() >= 1);
  const unsigned int n = A.size();
  B[0] = A[0] / 2.0;
  for (unsigned int i = 1; i < n; i++)
    Vector3::a_plus_b_times_c_equal_d(A[i], A[i - 1], 0.5, B[i]);
  B.back() = A[n - 1] / 2.0;
}

void vMidAvgInterior(const vector<Vector3> &A, vector<Vector3> &B) {
  assert(A.size() == B.size() + 1);
  assert(B.size() >= 1);
  const unsigned int n = A.size();
  for (unsigned int i = 1; i < n; i++)
    Vector3::a_plus_b_times_c_equal_d(A[i], A[i - 1], 0.5, B[i - 1]);
}

void vMidAvgInterior(const vector<REAL> &A, vector<REAL> &B) {
  assert(A.size() == B.size() + 1);
  assert(B.size() >= 1);
  const unsigned int n = A.size();
  for (unsigned int i = 1; i < n; i++) B[i - 1] = (A[i] + A[i - 1]) / 2.0;
}

void vFromPointsToElements(const vector<Vector3> &A, vector<Vector3> &B) {
  assert(A.size() == B.size() + 1);
  assert(B.size() >= 1);
  const unsigned int n = A.size();
  B.front() = A.front();
  for (unsigned int i = 1; i < n - 1; i++) {
    B[i] = A[i] / 2.0;
    B[i - 1] += A[i] / 2.0;
  }
  B.back() += A.back();
}

void vFromPointsToElements(const vector<REAL> &A, vector<REAL> &B) {
  assert(A.size() == B.size() + 1);
  assert(B.size() >= 1);
  const unsigned int n = A.size();
  B.front() = A.front();
  for (unsigned int i = 1; i < n - 1; i++) {
    B[i] = A[i] / 2.0;
    B[i - 1] += A[i] / 2.0;
  }
  B.back() += A.back();
}

void vT(const vector<Matrix3> &A, vector<Matrix3> &B) {
  assert(A.size() == B.size());
  assert(B.size() >= 1);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Matrix3>::iterator iterB = B.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++)
    (*iterA).T((*iterB));  // replace with in-place T()
}

void vDiagI(const vector<Matrix3> &A, vector<Matrix3> &B) {
  assert(A.size() == B.size());
  assert(B.size() >= 1);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Matrix3>::iterator iterB = B.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++)
    (*iterA).diagI((*iterB));  // replace with in-place I()
}

// a = b * a
void v_timesequal(const vector<Matrix3> &A, vector<Matrix3> &B) {
  assert(A.size() == B.size());
  assert(B.size() >= 1);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Matrix3>::iterator iterB = B.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++)
    (*iterB) = (*iterA) * (*iterB);
}

// c = a + b
void v_a_plus_b_equal_c(const vector<Vector3> &A, const vector<Vector3> &B,
                        vector<Vector3> &C) {
  assert(A.size() == B.size());
  assert(B.size() == C.size());
  assert(C.size() > 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++, iterC++)
    Vector3::a_plus_b_equal_c(
        (*iterA), (*iterB),
        (*iterC));  // in-place operator, no memory allocation
}

// c += a * b
void v_plusequal_a_times_b(const REAL A, const vector<Vector3> &B,
                           vector<Vector3> &C) {
  assert(B.size() == C.size());
  assert(C.size() > 0);
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterC++)
    Vector3::plusequal_a_times_b(A, (*iterB), (*iterC));
}

// d = a + b + c
void v_a_plus_b_plus_c_equal_d(const vector<Vector3> &A,
                               const vector<Vector3> &B,
                               const vector<Vector3> &C, vector<Vector3> &D) {
  assert(A.size() == B.size());
  assert(B.size() == C.size());
  assert(C.size() == D.size());
  assert(C.size() > 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::const_iterator iterC = C.begin();
  vector<Vector3>::iterator iterD = D.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++, iterC++, iterD++)
    Vector3::a_plus_b_plus_c_equal_d(
        (*iterA), (*iterB), (*iterC),
        (*iterD));  // in-place operator, no memory allocation
}

// c = a - b
void v_a_minus_b_equal_c(const vector<Vector3> &A, const vector<Vector3> &B,
                         vector<Vector3> &C)  // c = a - b
{
  assert(A.size() == C.size());
  assert(A.size() == B.size());
  assert(C.size() != 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterC = C.begin(); iterC != C.end(); iterC++, iterA++, iterB++)
    Vector3::a_minus_b_equal_c(
        (*iterA), (*iterB),
        (*iterC));  // in-place operator, no memory allocation
}

// c = a * b
void v_a_times_b_equal_c(const vector<Vector3> &A, const REAL B,
                         vector<Vector3> &C) {
  assert(A.size() == C.size());
  assert(C.size() != 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterC++)
    Vector3::a_times_b_equal_c((*iterA), B, (*iterC));
}

// c = a * b
void v_a_times_b_equal_c(const vector<REAL> &A, const REAL B, vector<REAL> &C) {
  assert(A.size() == C.size());
  assert(C.size() != 0);
  vector<REAL>::const_iterator iterA = A.begin();
  vector<REAL>::iterator iterC = C.begin();
  for (iterC = C.begin(); iterC != C.end(); iterC++, iterA++)
    (*iterC) =
        (*iterA) * B;  // operation with reals, no need for any in-place thing
}

// c = a * b
void v_a_times_b_equal_c(const vector<Matrix3> &A, const vector<Vector3> &B,
                         vector<Vector3> &C) {
  assert(A.size() == C.size());
  assert(A.size() == B.size());
  assert(C.size() != 0);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++, iterC++)
    Matrix3::matrix_times_vector(
        (*iterA), (*iterB),
        (*iterC));  // in-place operator, no memory allocation
}

// d = a * b * c
void v_a_times_b_times_c_equal_d(const REAL A, const vector<REAL> &B,
                                 const vector<Vector3> &C, vector<Vector3> &D) {
  assert(B.size() == C.size());
  assert(C.size() == D.size());
  assert(C.size() > 0);
  vector<REAL>::const_iterator iterB = B.begin();
  vector<Vector3>::const_iterator iterC = C.begin();
  vector<Vector3>::iterator iterD = D.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterC++, iterD++)
    Vector3::a_times_b_equal_c(
        (*iterC), A * (*iterB),
        (*iterD));  // in-place operator, no memory allocation
}

// c = a x b
void v_a_cross_b_equal_c(const vector<Vector3> &A, const vector<Vector3> &B,
                         vector<Vector3> &C) {
  assert(A.size() == B.size());
  assert(B.size() == C.size());
  assert(C.size() > 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++, iterC++)
    Vector3::a_cross_b_equal_c(
        (*iterA), (*iterB),
        (*iterC));  // in-place operator, no memory allocation
}

// d = a * b x c
void v_a_times_b_cross_c_equal_d(const vector<REAL> &A,
                                 const vector<Vector3> &B,
                                 const vector<Vector3> &C, vector<Vector3> &D) {
  assert(A.size() == B.size());
  assert(B.size() == C.size());
  assert(C.size() == D.size());
  assert(C.size() > 0);
  vector<REAL>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::const_iterator iterC = C.begin();
  vector<Vector3>::iterator iterD = D.begin();
  for (iterA = A.begin(); iterA != A.end();
       iterA++, iterB++, iterC++, iterD++) {
    Vector3::a_cross_b_equal_c(
        (*iterB), (*iterC),
        (*iterD));  // in-place operator, no memory allocation
    (*iterD) *= (*iterA);
  }
}

// d = a * b x c
void v_a_times_b_cross_c_equal_d(const REAL A, const vector<Vector3> &B,
                                 const vector<Vector3> &C, vector<Vector3> &D) {
  assert(B.size() == C.size());
  assert(C.size() == D.size());
  assert(C.size() > 0);
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<Vector3>::const_iterator iterC = C.begin();
  vector<Vector3>::iterator iterD = D.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterC++, iterD++) {
    Vector3::a_cross_b_equal_c(
        (*iterB), (*iterC),
        (*iterD));  // in-place operator, no memory allocation
    (*iterD) *= A;
  }
}

// c = a / b
void v_a_divide_b_equal_c(const vector<Vector3> &A, const vector<REAL> &B,
                          vector<Vector3> &C) {
  assert(A.size() == C.size());
  assert(A.size() == B.size());
  assert(C.size() != 0);
  vector<Vector3>::const_iterator iterA = A.begin();
  vector<REAL>::const_iterator iterB = B.begin();
  vector<Vector3>::iterator iterC = C.begin();
  for (iterA = A.begin(); iterA != A.end(); iterA++, iterB++, iterC++)
    Vector3::a_divide_b_equal_c(
        (*iterA), (*iterB),
        (*iterC));  // in-place operator, no memory allocation
}

// e = a * b / c - d
void v_a_times_b_divide_c_minus_d_equal_c(const vector<Matrix3> &A,
                                          const vector<Vector3> &B,
                                          const vector<REAL> &C,
                                          const Vector3 &D,
                                          vector<Vector3> &E) {
  assert(A.size() == C.size());
  assert(A.size() == B.size());
  assert(A.size() == E.size());
  assert(C.size() != 0);
  vector<Matrix3>::const_iterator iterA = A.begin();
  vector<Vector3>::const_iterator iterB = B.begin();
  vector<REAL>::const_iterator iterC = C.begin();
  vector<Vector3>::iterator iterE = E.begin();

  for (iterA = A.begin(); iterA != A.end();
       iterA++, iterB++, iterC++, iterE++) {
    Matrix3::matrix_times_vector((*iterA), (*iterB), (*iterE));
    (*iterE) /= (*iterC);  // in-place operator, no memory allocation
    (*iterE) -= D;         // in-place operator, no memory allocation
  }
}

// d = a / (b*c)
void v_a_dividepar_b_times_c_equal_d(const vector<REAL> &A,
                                     const vector<REAL> &B, const REAL C,
                                     vector<REAL> &D) {
  assert(A.size() == B.size());
  assert(D.size() == B.size());
  assert(B.size() != 0);
  vector<REAL>::const_iterator iterA = A.begin();
  vector<REAL>::const_iterator iterB = B.begin();
  vector<REAL>::iterator iterD = D.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterA++, iterD++)
    (*iterD) =
        (*iterA) /
        ((*iterB) * C);  // operation with reals, no need for any in-place thing
}

// d = sqrt( a / (b*c) )
void v_sqrt_a_dividepar_b_times_c_equal_d(const vector<REAL> &A,
                                          const vector<REAL> &B, const REAL C,
                                          vector<REAL> &D) {
  assert(A.size() == B.size());
  assert(D.size() == B.size());
  assert(B.size() != 0);
  vector<REAL>::const_iterator iterA = A.begin();
  vector<REAL>::const_iterator iterB = B.begin();
  vector<REAL>::iterator iterD = D.begin();
  for (iterB = B.begin(); iterB != B.end(); iterB++, iterA++, iterD++)
    (*iterD) =
        sqrt((*iterA) /
             ((*iterB) *
              C));  // operation with reals, no need for any in-place thing
}
