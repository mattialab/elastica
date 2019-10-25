#ifndef RODINITIALCONFIGURATIONS_H_
#define RODINITIALCONFIGURATIONS_H_

#include "GeometryFunctions.h"
#include "Matrix3.h"
#include "Rod.h"
#include "Tolerance.h"
#include "UsualHeaders.h"
#include "Vector3.h"

using namespace std;

class RodInitialConfigurations {
 public:
  static Rod *straightRod(const int n, const REAL totalMass, const REAL r0,
                          const Matrix3 _J0, const Matrix3 _B0,
                          const Matrix3 _S0, const REAL L0, const REAL totTwist,
                          const Vector3 origin, const Vector3 direction,
                          const Vector3 normal, const REAL nu,
                          const REAL relaxationNu, const bool useSelfContact);

  static Rod compressedRod(const int n, const REAL totalMass, const REAL r,
                           const Matrix3 J, const Matrix3 B, const Matrix3 S,
                           const REAL L, const REAL F, const REAL totTwist,
                           const REAL nu, const Vector3 origin,
                           const Vector3 direction, const Vector3 normal,
                           const bool useSelfContact);

  static Rod circleRod(const int n, const REAL totalMass, const REAL r,
                       const Matrix3 I, const Matrix3 B, const Matrix3 S,
                       const REAL L, const REAL totTwist, const REAL nu,
                       const bool useSelfContact);

  static Rod *curvedRod(const int n, const REAL density, const vector<REAL> r0,
                        const vector<Matrix3> J0, const vector<Matrix3> B0,
                        const vector<Matrix3> S0, const REAL totTwist,
                        const vector<Vector3> points, const Vector3 normal,
                        const REAL nu, const REAL relaxationNu,
                        const bool useSelfContact);

  static Rod *straightRod_v(const int n, const REAL density,
                            const vector<REAL> r0, const vector<Matrix3> J0,
                            const vector<Matrix3> B0, const vector<Matrix3> S0,
                            const REAL L0, const REAL totTwist,
                            const Vector3 origin, const Vector3 direction,
                            const Vector3 normal, const REAL nu,
                            const REAL relaxationNu, const bool useSelfContact);
  static Rod *straightRod_vscale(
      const int n, const REAL density, const vector<REAL> r0,
      const vector<Matrix3> J0, const vector<Matrix3> B0,
      const vector<Matrix3> S0, const REAL L0, const REAL totTwist,
      const Vector3 origin, const Vector3 direction, const Vector3 normal,
      const REAL nu, const REAL relaxationNu, const bool useSelfContact);
  // Twisted, helical artificial muscle
  static Rod *helicalRod(const int n, const REAL density, const REAL r0,
                         const REAL pitch_factor, const REAL helix_radius,
                         const REAL n_turns, const REAL E, const REAL G,
                         const REAL totTwist, const Vector3 origin,
                         const Vector3 direction, const REAL nu,
                         const REAL relaxationNu, const bool useSelfContact);
};

#endif
