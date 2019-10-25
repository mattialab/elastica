/*
 * Matrix3.cpp
 *
 *  Created on: Feb 26, 2015
 *      Author: mgazzola
 */

#include "Matrix3.h"

// constructors
Matrix3::Matrix3()
    : r1c1(1.),
      r1c2(0.),
      r1c3(0.),
      r2c1(0.),
      r2c2(1.),
      r2c3(0.),
      r3c1(0.),
      r3c2(0.),
      r3c3(1.) {}

Matrix3::Matrix3(REAL r1c1, REAL r1c2, REAL r1c3, REAL r2c1, REAL r2c2,
                 REAL r2c3, REAL r3c1, REAL r3c2, REAL r3c3)
    : r1c1(r1c1),
      r1c2(r1c2),
      r1c3(r1c3),
      r2c1(r2c1),
      r2c2(r2c2),
      r2c3(r2c3),
      r3c1(r3c1),
      r3c2(r3c2),
      r3c3(r3c3) {}

Matrix3::Matrix3(const Vector3 &S)
    :  // skew matrix <-> vector mapping
      r1c1(0.),
      r1c2(S.z),
      r1c3(-S.y),
      r2c1(-S.z),
      r2c2(0.),
      r2c3(S.x),
      r3c1(S.y),
      r3c2(-S.x),
      r3c3(0.) {}

Matrix3::Matrix3(const Vector3 &r1, const Vector3 &r2, const Vector3 &r3)
    : r1c1(r1.x),
      r1c2(r1.y),
      r1c3(r1.z),
      r2c1(r2.x),
      r2c2(r2.y),
      r2c3(r2.z),
      r3c1(r3.x),
      r3c2(r3.y),
      r3c3(r3.z) {}

Matrix3::Matrix3(const Matrix3 &M)
    : r1c1(M.r1c1),
      r1c2(M.r1c2),
      r1c3(M.r1c3),
      r2c1(M.r2c1),
      r2c2(M.r2c2),
      r2c3(M.r2c3),
      r3c1(M.r3c1),
      r3c2(M.r3c2),
      r3c3(M.r3c3) {}

// I/O
std::ostream &operator<<(std::ostream &out, const Matrix3 &A) {
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r1c1 << ", ";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r1c2 << ", ";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r1c3 << "\n";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r2c1 << ", ";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r2c2 << ", ";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r2c3 << "\n";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r3c1 << ", ";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r3c2 << ", ";
  out.setf(std::ios::fixed | std::ios::right);
  out.width(A.PRINTPRECISION + 3);
  out.precision(A.PRINTPRECISION);
  out << A.r3c3;

  return out;
}

////////////////// operators not declarable inside Matrix3 //////////

Matrix3 operator*(REAL a, const Matrix3 &M) { return M * a; }

// treating Vector3's as row vectors here: we need to be careful in using it!
Vector3 operator*(const Vector3 &v, const Matrix3 &M) {
  return Vector3(v.x * M.r1c1 + v.y * M.r2c1 + v.z * M.r3c1,
                 v.x * M.r1c2 + v.y * M.r2c2 + v.z * M.r3c2,
                 v.x * M.r1c3 + v.y * M.r2c3 + v.z * M.r3c3);
}

// exponential of a vector (\in so(3)) maps to an orthogonal matrix (\in SO(3))
// using rodiguez formula. The size (angle) of the rotation is just the
// magnitude of the vector.
Matrix3 exp(const Vector3 &S, const REAL tol) {
  const REAL theta = S.length();
  if (theta < tol) return Matrix3();
  Matrix3 Sx = Matrix3(S.unitize());
  return Matrix3() + Sx * sin(theta) + (1 - cos(theta)) * Sx * Sx;
}

void exp(const Vector3 &S, Matrix3 &B) {
  // Clear matrix container and initialize to identity matrix
  B.r1c1 = 1.0;
  B.r1c2 = 0.0;
  B.r1c3 = 0.0;
  B.r2c1 = 0.0;
  B.r2c2 = 1.0;
  B.r2c3 = 0.0;
  B.r3c1 = 0.0;
  B.r3c2 = 0.0;
  B.r3c3 = 1.0;

  // Fast in-place evaluation of the rodriguez formula
  const REAL theta = S.length();
  const Matrix3 Sx = Matrix3(S.unitize());
  const REAL sintheta = sin(theta);
  const REAL oneminuscostheta = 1.0 - cos(theta);

  B.r1c1 += oneminuscostheta * (Sx.r1c2 * Sx.r2c1 + Sx.r1c3 * Sx.r3c1);
  B.r1c2 += Sx.r1c2 * sintheta + oneminuscostheta * (Sx.r1c3 * Sx.r3c2);
  B.r1c3 += Sx.r1c3 * sintheta + oneminuscostheta * (Sx.r1c2 * Sx.r2c3);
  B.r2c1 += Sx.r2c1 * sintheta + oneminuscostheta * (Sx.r2c3 * Sx.r3c1);
  B.r2c2 += oneminuscostheta * (Sx.r2c1 * Sx.r1c2 + Sx.r2c3 * Sx.r3c2);
  B.r2c3 += Sx.r2c3 * sintheta + oneminuscostheta * (Sx.r2c1 * Sx.r1c3);
  B.r3c1 += Sx.r3c1 * sintheta + oneminuscostheta * (Sx.r3c2 * Sx.r2c1);
  B.r3c2 += Sx.r3c2 * sintheta + oneminuscostheta * (Sx.r3c1 * Sx.r1c2);
  B.r3c3 += oneminuscostheta * (Sx.r3c1 * Sx.r1c3 + Sx.r3c2 * Sx.r2c3);
}

// generates a uniform gaussian random variable in each element
Matrix3 randMatrix3() {
  return Matrix3(randn_notrig(0, 1), randn_notrig(0, 1), randn_notrig(0, 1),
                 randn_notrig(0, 1), randn_notrig(0, 1), randn_notrig(0, 1),
                 randn_notrig(0, 1), randn_notrig(0, 1), randn_notrig(0, 1));
}
