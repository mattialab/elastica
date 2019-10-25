#ifndef VECTOR3_H_
#define VECTOR3_H_

#include "ArithmeticPrecision.h"
#include "MathFunctions.h"
#include "Tolerance.h"
#include "UsualHeaders.h"

class Vector3 {
 private:
  static const unsigned int PRINTPRECISION = 6;

 public:
  // This should be completely recoded in a HPC code
  REAL x, y, z;

  // constructors (never inline)
  Vector3();
  Vector3(REAL x, REAL y, REAL z);
  Vector3(const Vector3 &V);

  // accessor
  inline REAL operator[](int element) const {
    switch (element) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        throw "Bad [] access attempt on Vector3.\n";
    }
  }

  // vector operators
  inline Vector3 operator+(const Vector3 &V) const {
    return Vector3(x + V.x, y + V.y, z + V.z);
  }
  inline Vector3 operator-(const Vector3 &V) const {
    return Vector3(x - V.x, y - V.y, z - V.z);
  }
  inline Vector3 operator*(const Vector3 &V) const {
    return Vector3(y * V.z - z * V.y, z * V.x - x * V.z, x * V.y - y * V.x);
  }  // cross product
  inline REAL operator%(const Vector3 &V) const {
    return x * V.x + y * V.y + z * V.z;
  }  // dot product

  // scalar operators
  inline Vector3 operator*(REAL a) const {
    return Vector3(x * a, y * a, z * a);
  }
  inline Vector3 operator/(REAL a) const {
    return Vector3(x / a, y / a, z / a);
  }
  inline Vector3 operator-() const { return Vector3(-x, -y, -z); }

  // in-place modifiers
  inline Vector3 &operator+=(const Vector3 &V) {
    x += V.x;
    y += V.y;
    z += V.z;
    return *this;
  }
  inline Vector3 &operator-=(const Vector3 &V) {
    x -= V.x;
    y -= V.y;
    z -= V.z;
    return *this;
  }
  inline Vector3 &operator*=(const Vector3 &V) {
    *this = *this * V;
    return *this;
  }

  inline Vector3 &operator*=(REAL a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }
  inline Vector3 &operator/=(REAL a) {
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }

  // utility functions
  inline REAL length() const { return sqrt(x * x + y * y + z * z); }
  inline bool isUnit(const REAL tol = Tolerance::tol()) const {
    return fabs(length() - 1.0) < tol;
  }

  inline Vector3 unitize() const {
    const REAL l = length();
    if (l < Tolerance::tol()) return Vector3();
    return *this / l;
  }

  // speed functions
  static inline void plusequal_a_times_b(const REAL A, const Vector3 &B,
                                         Vector3 &C) {
    C.x += A * B.x;
    C.y += A * B.y;
    C.z += A * B.z;
  }

  static inline void a_plus_b_equal_c(const Vector3 &A, const Vector3 &B,
                                      Vector3 &C) {
    C.x = A.x + B.x;
    C.y = A.y + B.y;
    C.z = A.z + B.z;
  }

  static inline void a_plus_b_plus_c_equal_d(const Vector3 &A, const Vector3 &B,
                                             const Vector3 &C, Vector3 &D) {
    D.x = A.x + B.x + C.x;
    D.y = A.y + B.y + C.y;
    D.z = A.z + B.z + C.z;
  }

  static inline void a_minus_b_equal_c(const Vector3 &A, const Vector3 &B,
                                       Vector3 &C) {
    C.x = A.x - B.x;
    C.y = A.y - B.y;
    C.z = A.z - B.z;
  }

  static inline void a_times_b_equal_c(const Vector3 &A, const REAL B,
                                       Vector3 &C) {
    C.x = A.x * B;
    C.y = A.y * B;
    C.z = A.z * B;
  }

  static inline void a_cross_b_equal_c(const Vector3 &A, const Vector3 &B,
                                       Vector3 &C)  // cross product
  {
    C.x = A.y * B.z - A.z * B.y;
    C.y = A.z * B.x - A.x * B.z;
    C.z = A.x * B.y - A.y * B.x;
  }

  static inline void a_divide_b_equal_c(const Vector3 &A, const REAL B,
                                        Vector3 &C) {
    C.x = A.x / B;
    C.y = A.y / B;
    C.z = A.z / B;
  }

  static inline void a_plus_b_times_c_equal_d(const Vector3 &A,
                                              const Vector3 &B, const REAL C,
                                              Vector3 &D) {
    D.x = (A.x + B.x) * C;
    D.y = (A.y + B.y) * C;
    D.z = (A.z + B.z) * C;
  }

  // It performs the following operation a = v*(I-0.5*t'*t), where v is the
  // velocity vector, t is the tangent and I the unit matrix
  static inline void projectionSBT(const Vector3 &v, const Vector3 &t,
                                   Vector3 &f) {
    const REAL tx = t.x;
    const REAL ty = t.y;
    const REAL tz = t.z;

    const REAL vx = v.x;
    const REAL vy = v.y;
    const REAL vz = v.z;

    // t'*t
    const REAL a11 = tx * tx;
    const REAL a12 = tx * ty;
    const REAL a13 = tx * tz;

    const REAL a21 = ty * tx;
    const REAL a22 = ty * ty;
    const REAL a23 = ty * tz;

    const REAL a31 = tz * tx;
    const REAL a32 = tz * ty;
    const REAL a33 = tz * tz;

    // v*(I-0.5*t'*t)
    f.x = (1 - 0.5 * a11) * vx + (0 - 0.5 * a12) * vy + (0 - 0.5 * a13) * vz;
    f.y = (0 - 0.5 * a21) * vx + (1 - 0.5 * a22) * vy + (0 - 0.5 * a23) * vz;
    f.z = (0 - 0.5 * a31) * vx + (0 - 0.5 * a32) * vy + (1 - 0.5 * a33) * vz;
  }

  // IO (never inline)
  friend std::ostream &operator<<(std::ostream &out, const Vector3 &A);
};

///// Additional vector functions not defined as member functions //////

Vector3 operator*(REAL a, const Vector3 &V);

// gets the angle between two vector in range(0, PI)
REAL angle(const Vector3 &A, const Vector3 &B);

// generates a uniform gaussian random variable in 3 dimensions
Vector3 randVector3();
Vector3 randVector3(REAL mu, REAL sigma);

// random only in x direction
Vector3 randXOnly();

#endif
