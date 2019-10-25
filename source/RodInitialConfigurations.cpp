/*
 * RodInitialConfiguations.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: mgazzola
 */

#include "RodInitialConfigurations.h"

Rod *RodInitialConfigurations::straightRod(
    const int n, const REAL totalMass, const REAL r0, const Matrix3 _J0,
    const Matrix3 _B0, const Matrix3 _S0, const REAL L0, const REAL totTwist,
    const Vector3 origin, const Vector3 direction, const Vector3 normal,
    const REAL nu, const REAL relaxationNu, const bool useSelfContact) {
  // Bunch of sanity checks
  assert(n > 1);
  assert(totalMass > Tolerance::tol());
  assert(r0 > Tolerance::tol());
  assert(L0 > Tolerance::tol());
  assert(nu >= 0.0);
  assert(relaxationNu >= 0.0);
  assert(direction.length() > Tolerance::tol());
  assert(normal.length() > Tolerance::tol());

  assert(_B0.r2c1 == 0.0);
  assert(_B0.r3c1 == 0.0);
  assert(_B0.r1c2 == 0.0);
  assert(_B0.r3c2 == 0.0);
  assert(_B0.r1c3 == 0.0);
  assert(_B0.r2c3 == 0.0);
  assert(_B0.r1c1 >= Tolerance::tol());
  assert(_B0.r2c2 >= Tolerance::tol());
  assert(_B0.r3c3 >= Tolerance::tol());

  assert(_S0.r2c1 == 0.0);
  assert(_S0.r3c1 == 0.0);
  assert(_S0.r1c2 == 0.0);
  assert(_S0.r3c2 == 0.0);
  assert(_S0.r1c3 == 0.0);
  assert(_S0.r2c3 == 0.0);
  assert(_S0.r1c1 >= Tolerance::tol());
  assert(_S0.r2c2 >= Tolerance::tol());
  assert(_S0.r3c3 >= Tolerance::tol());

  // Density
  const REAL density = totalMass / (M_PI * r0 * r0 * L0);

  // Initialize discretization points
  vector<Vector3> x0 = vRange(origin, origin + direction * L0, n + 1);

  // Set velocities to zero
  vector<Vector3> v0 = vector<Vector3>(n + 1);

  // u0 is a vector orthogonal to the tangent, .i.e. the direction.
  // If not provided I just compute one, using a random non-zero vector.
  // Note that to be able to compare codes, I fix the random vector and
  // I do not generate one, because that is BAD!
#ifndef NDEBUG
  const Vector3 t0 = (x0[1] - x0[0]).unitize();
  assert((t0 * normal).length() > Tolerance::tol());
#endif

  // Now I can align frames using the orthogonal vector provided/generated
  vector<Matrix3> Q0 = alignFrames(x0, normal.unitize());

  // Apply twist about the orthonormal vector d3
  applyTwists(Q0, vRange((REAL)0.0, totTwist, n));

  // Set angular velocity to zero
  vector<Vector3> w0 = vector<Vector3>(n);

  // Set rest edge lengths
  vector<REAL> l0 = vLength(vDiff(x0));
  const REAL dl0 = l0[0];

  // Set volume discretization elements
  const vector<REAL> V0 = vector<REAL>(n, M_PI * r0 * r0 * dl0);

  // Set shear vector to zero
  vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n);

  // Mass of vertex point wise element.
  // VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
  // This is why m is obtained dividing the total mass by n, and then to
  // conserve the total mass the masses of the first and last verteces are
  // divideb by 2
  const REAL m = totalMass / (double)(n);
  vector<REAL> masses = vector<REAL>(n + 1, m);
  masses.front() /= 2.0;
  masses.back() /= 2.0;

  // Intrinsic curvature in rest configuration
  const vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1);

  // Mass second moment of inertia matrix in rest configuration
  vector<Matrix3> J0 = vector<Matrix3>(n, _J0);

  // Bending matrix in reference configuration
  vector<Matrix3> B0 = vector<Matrix3>(n - 1, _B0);

  // Shear matrix in reference configuration
  vector<Matrix3> S0 = vector<Matrix3>(n, _S0);

  return new Rod(n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0,
                 masses, V0, density, J0, B0, S0, nu, relaxationNu,
                 useSelfContact);
}

// A rod on the z-axis compressed so as to push back with forces F
Rod RodInitialConfigurations::compressedRod(
    const int n, const REAL totalMass, const REAL r, const Matrix3 J,
    const Matrix3 B, const Matrix3 S, const REAL L, const REAL F,
    const REAL totTwist, const REAL nu, const Vector3 origin,
    const Vector3 direction, const Vector3 normal, const bool useSelfContact) {
  // Bunch of sanity checks
  assert(n > 1);
  assert(totalMass > Tolerance::tol());
  assert(r > Tolerance::tol());
  assert(L > Tolerance::tol());
  assert(nu >= 0.0);
  assert(direction.length() > Tolerance::tol());
  assert(normal.length() > Tolerance::tol());

  assert(B.r2c1 == 0.0);
  assert(B.r3c1 == 0.0);
  assert(B.r1c2 == 0.0);
  assert(B.r3c2 == 0.0);
  assert(B.r1c3 == 0.0);
  assert(B.r2c3 == 0.0);
  assert(B.r1c1 >= Tolerance::tol());
  assert(B.r2c2 >= Tolerance::tol());
  assert(B.r3c3 >= Tolerance::tol());

  assert(S.r2c1 == 0.0);
  assert(S.r3c1 == 0.0);
  assert(S.r1c2 == 0.0);
  assert(S.r3c2 == 0.0);
  assert(S.r1c3 == 0.0);
  assert(S.r2c3 == 0.0);
  assert(S.r1c1 >= Tolerance::tol());
  assert(S.r2c2 >= Tolerance::tol());
  assert(S.r3c3 >= Tolerance::tol());

  // Density
  const REAL density = totalMass / (M_PI * r * r * L);

  // fractional compression of rod needed to resist force F
  const REAL c = F / S[2][2];

  // Initialize discretization points
  vector<Vector3> x = vRange(origin, origin + direction * L * (1 - c), n + 1);

  // Set velocities to zero
  vector<Vector3> v = vector<Vector3>(n + 1);

  // u0 is a vector orthogonal to the tangent, .i.e. the direction.
  // If not provided I just compute one, using a random non-zero vector.
  // Note that to be able to compare codes, I fix the random vector and
  // I do not generate one, because that is BAD!
#ifndef NDEBUG
  const Vector3 t0 = (x[1] - x[0]).unitize();
  assert((t0 * normal).length() > Tolerance::tol());
#endif

  // Now I can align frames using the orthogonal vector provided/generated
  vector<Matrix3> Q = alignFrames(x, normal.unitize());

  // Apply twist about the orthonormal vector d3
  applyTwists(Q, vRange((REAL)0.0, totTwist, n));

  // Set angular velocity to zero
  vector<Vector3> w = vector<Vector3>(n);

  // Set rest edge lengths
  vector<REAL> l = vector<REAL>(n, L / n);

  // Set volume discretization elements
  vector<REAL> V = vector<REAL>(n, l[0] * r * r * M_PI);

  // Set shear vector to zero
  vector<Vector3> sigma0 = vector<Vector3>(n);

  const REAL relaxationNu = 0.0;

  // Mass of vertex point wise element.
  // VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
  // This is why m is obtained dividing the total mass by n, and then to
  // conserve the total mass the masses of the first and last verteces are
  // divideb by 2
  const REAL m = totalMass / (double)(n);
  vector<REAL> masses = vector<REAL>(n + 1, m);
  masses.front() /= 2.0;
  masses.back() /= 2.0;

  // Intrinsic curvature in rest configuration
  const vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1);

  // Mass second moment of inertia matrix in rest configuration
  vector<Matrix3> J0 = vector<Matrix3>(n, J);

  // Bending matrix in reference configuration
  vector<Matrix3> B0 = vector<Matrix3>(n - 1, B);

  // Shear matrix in reference configuration
  vector<Matrix3> S0 = vector<Matrix3>(n, S);

  return Rod(n, x, v, Q, w, l, intrinsic_k0, sigma0, masses, V, density, J0, B0,
             S0, nu, relaxationNu, useSelfContact);
}

Rod RodInitialConfigurations::circleRod(const int n, const REAL totalMass,
                                        const REAL r, const Matrix3 I,
                                        const Matrix3 B, const Matrix3 S,
                                        const REAL L, const REAL totTwist,
                                        const REAL nu,
                                        const bool useSelfContact) {
  // Density
  const REAL density = totalMass / (M_PI * r * r * L);

  // Masses
  const REAL m = totalMass / (double)(n);
  vector<REAL> masses = vector<REAL>(n + 1, m);
  masses.front() /= 2.0;
  masses.back() /= 2.0;

  vector<Vector3> x;

  for (int i = 0; i < n + 1; i++)
    x.push_back(L / (2 * M_PI) *
                Vector3(cos(2 * M_PI / n * i), sin(2 * M_PI / n * i), 0));

  vector<Vector3> v = vector<Vector3>(n + 1);
  vector<Matrix3> Q = alignFrames(x, Vector3(0, 0, 1));
  applyTwists(Q, vRange((REAL)0.0, totTwist, n));
  vector<Vector3> w = vector<Vector3>(n);
  vector<REAL> l = vLength(vDiff(x));
  vector<REAL> V = vector<REAL>(n, l[0] * r * r * M_PI);
  vector<Vector3> k0 = vector<Vector3>(n - 1);
  vector<Vector3> sigma0 = vector<Vector3>(n);

  const REAL relaxationNu = 0.0;

  return Rod(n, x, v, Q, w, l, k0, sigma0, masses, V, density,
             vector<Matrix3>(n, I), vector<Matrix3>(n - 1, B),
             vector<Matrix3>(n, S), nu, relaxationNu, useSelfContact);
}

// This function allows rod to be initialized with arbitrary 3D curvature.
Rod *RodInitialConfigurations::curvedRod(
    const int n, const REAL density, const vector<REAL> r0,
    const vector<Matrix3> J0, const vector<Matrix3> B0,
    const vector<Matrix3> S0, const REAL totTwist, const vector<Vector3> points,
    const Vector3 normal, const REAL nu, const REAL relaxationNu,
    const bool useSelfContact) {
  assert(n == (points.size() - 1));

  // Initialize discretization points
  vector<Vector3> x0 = points;

  // Set velocities to zero
  vector<Vector3> v0 = vector<Vector3>(n + 1);

  // Now I can align frames using the orthogonal vector provided/generated
  vector<Matrix3> Q0 = alignFrames(x0, normal.unitize());

  // Apply twist about the orthonormal vector d3
  applyTwists(Q0, vRange((REAL)0.0, totTwist, n));

  // Set angular velocity to zero
  vector<Vector3> w0 = vector<Vector3>(n);

  // Set rest edge lengths
  vector<REAL> l0 = vLength(vDiff(x0));

  // Set shear vector to zero
  vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n);

  // Mass of vertex point wise element.
  // VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
  // This is why m is obtained dividing the total mass by n, and then to
  // conserve the total mass the masses of the first and last verteces are
  // divideb by 2

  // Set volume discretization elements
  vector<REAL> V0 = vector<REAL>(n);
  vector<REAL> masses = vector<REAL>(n + 1);
  REAL totalmass = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    V0[i] = M_PI * r0[i] * r0[i] * l0[i];
    REAL mass = V0[i] * density;
    masses[i] += mass / 2.0;
    masses[i + 1] += mass / 2.0;
    totalmass += mass;
  }

  // Voronaoi Domian
  vector<REAL> d0 = vector<REAL>(n - 1);
  vMidAvgInterior(l0, d0);

  // Intrinsic curvature in rest configuration
  vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1);
  vector<Matrix3> tempVM3_nminus = vector<Matrix3>(n - 1);
  vector<Vector3> tempVV3_nminus = vector<Vector3>(n - 1);

  vRotDiff(Q0, tempVM3_nminus);
  vLog(tempVM3_nminus, tempVV3_nminus);
  v_a_divide_b_equal_c(tempVV3_nminus, d0, intrinsic_k0);

  return new Rod(n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0,
                 masses, V0, density, J0, B0, S0, nu, relaxationNu,
                 useSelfContact);
}

// A vector version of straight rod initialization allows varying radii.
Rod *RodInitialConfigurations::straightRod_v(
    const int n, const REAL density, const vector<REAL> r0,
    const vector<Matrix3> J0, const vector<Matrix3> B0,
    const vector<Matrix3> S0, const REAL L0, const REAL totTwist,
    const Vector3 origin, const Vector3 direction, const Vector3 normal,
    const REAL nu, const REAL relaxationNu, const bool useSelfContact) {
  // Bunch of sanity checks
  assert(n > 1);
  assert(density > Tolerance::tol());
  assert(r0[0] > Tolerance::tol());
  assert(L0 > Tolerance::tol());
  assert(nu >= 0.0);
  assert(relaxationNu >= 0.0);
  assert(direction.length() > Tolerance::tol());
  assert(normal.length() > Tolerance::tol());

  assert(B0[0].r2c1 == 0.0);
  assert(B0[0].r3c1 == 0.0);
  assert(B0[0].r1c2 == 0.0);
  assert(B0[0].r3c2 == 0.0);
  assert(B0[0].r1c3 == 0.0);
  assert(B0[0].r2c3 == 0.0);
  assert(B0[0].r1c1 >= Tolerance::tol());
  assert(B0[0].r2c2 >= Tolerance::tol());
  assert(B0[0].r3c3 >= Tolerance::tol());

  assert(S0[0].r2c1 == 0.0);
  assert(S0[0].r3c1 == 0.0);
  assert(S0[0].r1c2 == 0.0);
  assert(S0[0].r3c2 == 0.0);
  assert(S0[0].r1c3 == 0.0);
  assert(S0[0].r2c3 == 0.0);
  assert(S0[0].r1c1 >= Tolerance::tol());
  assert(S0[0].r2c2 >= Tolerance::tol());
  assert(S0[0].r3c3 >= Tolerance::tol());

  // Initialize discretization points
  vector<Vector3> x0 = vRange(origin, origin + direction * L0, n + 1);

  // Set velocities to zero
  vector<Vector3> v0 = vector<Vector3>(n + 1);

  // u0 is a vector orthogonal to the tangent, .i.e. the direction.
  // If not provided I just compute one, using a random non-zero vector.
  // Note that to be able to compare codes, I fix the random vector and
  // I do not generate one, because that is BAD!
#ifndef NDEBUG
  const Vector3 t0 = (x0[1] - x0[0]).unitize();
  assert((t0 * normal).length() > Tolerance::tol());
#endif

  // Now I can align frames using the orthogonal vector provided/generated
  vector<Matrix3> Q0 = alignFrames(x0, normal.unitize());

  // Apply twist about the orthonormal vector d3
  applyTwists(Q0, vRange((REAL)0.0, totTwist, n));

  // Set angular velocity to zero
  vector<Vector3> w0 = vector<Vector3>(n);

  // Set rest edge lengths
  vector<REAL> l0 = vLength(vDiff(x0));
  const REAL dl0 = l0[0];

  // Set shear vector to zero
  vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n);

  // Mass of vertex point wise element.
  // VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
  // This is why m is obtained dividing the total mass by n, and then to
  // conserve the total mass the masses of the first and last verteces are
  // divideb by 2

  // const REAL m = totalMass/(double)(n);
  // Set volume discretization elements
  vector<REAL> V0 = vector<REAL>(n);
  vector<REAL> masses = vector<REAL>(n + 1);
  REAL totalmass = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    V0[i] = M_PI * r0[i] * r0[i] * dl0;
    REAL mass = V0[i] * density;
    masses[i] += mass / 2.0;
    masses[i + 1] += mass / 2.0;
  }
  totalmass = vSum(masses);
  // Intrinsic curvature in rest configuration
  const vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1);

  return new Rod(n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0,
                 masses, V0, density, J0, B0, S0, nu, relaxationNu,
                 useSelfContact);
}

// The function is based on "straightRod_v", but does not scale the first and
// last element by 0.5. Masses are computed exactly corresponding to the local
// radii. So that the first and last elements are not lightweight. This will
// match better the actual physics when rod is tapered and is very important for
// numerical stability. Only used in the wing case.
Rod *RodInitialConfigurations::straightRod_vscale(
    const int n, const REAL density, const vector<REAL> r0,
    const vector<Matrix3> J0, const vector<Matrix3> B0,
    const vector<Matrix3> S0, const REAL L0, const REAL totTwist,
    const Vector3 origin, const Vector3 direction, const Vector3 normal,
    const REAL nu, const REAL relaxationNu, const bool useSelfContact) {
  // Bunch of sanity checks
  assert(n > 1);
  assert(density > Tolerance::tol());
  assert(r0[0] > Tolerance::tol());
  assert(L0 > Tolerance::tol());
  assert(nu >= 0.0);
  assert(relaxationNu >= 0.0);
  assert(direction.length() > Tolerance::tol());
  assert(normal.length() > Tolerance::tol());

  assert(B0[0].r2c1 == 0.0);
  assert(B0[0].r3c1 == 0.0);
  assert(B0[0].r1c2 == 0.0);
  assert(B0[0].r3c2 == 0.0);
  assert(B0[0].r1c3 == 0.0);
  assert(B0[0].r2c3 == 0.0);
  assert(B0[0].r1c1 >= Tolerance::tol());
  assert(B0[0].r2c2 >= Tolerance::tol());
  assert(B0[0].r3c3 >= Tolerance::tol());

  assert(S0[0].r2c1 == 0.0);
  assert(S0[0].r3c1 == 0.0);
  assert(S0[0].r1c2 == 0.0);
  assert(S0[0].r3c2 == 0.0);
  assert(S0[0].r1c3 == 0.0);
  assert(S0[0].r2c3 == 0.0);
  assert(S0[0].r1c1 >= Tolerance::tol());
  assert(S0[0].r2c2 >= Tolerance::tol());
  assert(S0[0].r3c3 >= Tolerance::tol());

  // Initialize discretization points
  vector<Vector3> x0 = vRange(origin, origin + direction * L0, n + 1);

  // Set velocities to zero
  vector<Vector3> v0 = vector<Vector3>(n + 1);

#ifndef NDEBUG
  const Vector3 t0 = (x0[1] - x0[0]).unitize();
  assert((t0 * normal).length() > Tolerance::tol());
#endif

  // Now I can align frames using the orthogonal vector provided/generated
  vector<Matrix3> Q0 = alignFrames(x0, normal.unitize());

  // Apply twist about the orthonormal vector d3
  applyTwists(Q0, vRange((REAL)0.0, totTwist, n));

  // Set angular velocity to zero
  vector<Vector3> w0 = vector<Vector3>(n);

  // Set rest edge lengths
  vector<REAL> l0 = vLength(vDiff(x0));
  const REAL dl0 = l0[0];

  // Set shear vector to zero
  vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n);

  vector<REAL> V0 = vector<REAL>(n);
  vector<REAL> masses = vector<REAL>(n + 1);
  REAL totalmass = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    V0[i] = M_PI * r0[i] * r0[i] * dl0;
    REAL mass = V0[i] * density;
    masses[i] = mass;
  }
  REAL mass = V0[n - 1] * density;
  masses[n] = mass;

  totalmass = vSum(masses);
  // Intrinsic curvature in rest configuration
  const vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1);

  return new Rod(n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0,
                 masses, V0, density, J0, B0, S0, nu, relaxationNu,
                 useSelfContact);
}

// Artificial Muscle
Rod *RodInitialConfigurations::helicalRod(
    const int n, const REAL density, const REAL r0, const REAL pitch_factor,
    const REAL helix_radius, const REAL n_turns, const REAL E, const REAL G,
    const REAL totTwist, const Vector3 origin, const Vector3 direction,
    const REAL nu, const REAL relaxationNu, const bool useSelfContact) {
  // Bunch of sanity checks
  assert(n > 1);
  assert(density > Tolerance::tol());
  assert(r0 > Tolerance::tol());
  assert(pitch_factor > Tolerance::tol());
  assert(nu >= 0.0);
  assert(relaxationNu >= 0.0);
  assert(helix_radius > 0.0);
  assert(n_turns > 0.0);
  assert(direction.length() > Tolerance::tol());

  // Initialize discretization at nodes to zero
  // Node positions
  vector<Vector3> x0 = vector<Vector3>(n + 1);
  // Node velocities
  vector<Vector3> v0 = vector<Vector3>(n + 1);

  // Initialize discretization at elements to zero
  // Angular velocities
  vector<Vector3> w0 = vector<Vector3>(n);

  // Geometrical parameters
  REAL rod_cross_area = M_PI * r0 * r0;

  // Max of curve parameter s, min is 0.0
  REAL max_s = 2 * M_PI * n_turns;
  auto s = vRange(0.0, max_s, n + 1);

  REAL s_offset = 0.5 * M_PI;

  // Fill in points wrt parametrized curve
  for (size_t i = 0; i < n + 1; ++i) {
    auto temp =
        // Vector3(helix_radius * cos(s[i]), helix_radius * sin(s[i]),
        //                     pitch_factor * s[i]);
        Vector3(helix_radius * cos(s[i] + s_offset),
                helix_radius * sin(s[i] + s_offset), pitch_factor * s[i]);
    // auto temp = Vector3(s[i] / (M_PI * n_turns) - 1.,
    //                     std::tanh(s[i] - M_PI * n_turns), 0.0);
    // Todo: Fill in direction change
    x0[i] = temp + origin;
  }

  // std::ofstream my_out;
  // my_out.open("test.csv");
  // for (size_t i = 0; i < n + 1; ++i) {
  //   my_out << i << "," << x0[i].x << "," << x0[i].y << "," << x0[i].z <<
  //   "\n";
  // }
  // my_out.close();

  // Now I can align frames using the orthogonal vector provided/generated
  // vector<Matrix3> Q1 = alignFrames(x0, direction.unitize());

  std::vector<Matrix3> Q0;
  {
    // Done this way because no one knows what basis coordinate system we
    // are using
    Q0.reserve(n);
    // Index i starts from 1, because the first frame has been filled already
    REAL d3_prefac =
        1. / sqrt(helix_radius * helix_radius + pitch_factor * pitch_factor);
    for (size_t i = 0; i < n; ++i) {
      // Create body frames like so
      //    d2  d3
      //     | /
      // d1 __|/
      // where the helix is (cos(t),sin(t),t)
      auto new_s = 0.5 * (s[i] + s[i + 1]);
      auto d1 = Vector3(-cos(new_s + s_offset), -sin(new_s + s_offset), 0.0);
      auto d3 = d3_prefac * Vector3(-helix_radius * sin(new_s + s_offset),
                                    helix_radius * cos(new_s + s_offset),
                                    pitch_factor);
      // auto d1 = Vector3(-cos(new_s), -sin(new_s), 0.0);
      // auto d3 = d3_prefac * Vector3(-helix_radius * sin(new_s),
      //                               helix_radius * cos(new_s), pitch_factor);

      Q0.push_back(Matrix3(d1, d3 * d1, d3));
    }
  }

  // my_out.open("test_mat.csv");
  // for (size_t i = 0; i < n; ++i) {
  //   my_out << Q0[i] << "\n";
  // }
  // my_out.close();

  // my_out.open("test_mat_native.csv");
  // for (size_t i = 0; i < n; ++i) {
  //   my_out << Q1[i] << "\n";
  // }
  // my_out.close();

  // Set rest edge lengths
  vector<REAL> l0 = vLength(vDiff(x0));
  const REAL dl0 = l0[0];

  // Voronoi, because I need to code up intrinsic curvature
  // without touching the rod class
  auto d0 = vector<REAL>(n - 1);  // --> (n-1)
  vMidAvgInterior(l0, d0);        // (n) --> (n-1)

  auto tempVM3_nminus = vector<Matrix3>(n - 1);
  auto tempVV3_nminus = vector<Vector3>(n - 1);

  // Set shear vector to zero
  vector<Vector3> intrinsic_k0 = vector<Vector3>(n - 1);
  vRotDiff(Q0, tempVM3_nminus);
  vLog(tempVM3_nminus, tempVV3_nminus);
  v_a_divide_b_equal_c(tempVV3_nminus, d0, intrinsic_k0);

  // my_out.open("curvature_test.csv");
  // for (size_t i = 0; i < n - 1; ++i) {
  //   // auto temp_two = Q0[i].T() * intrinsic_k0[i];
  //   auto temp_two = intrinsic_k0[i];
  //   my_out << i << "," << temp_two.x << "," << temp_two.y << "," <<
  //   temp_two.z
  //          << "\n";
  // }
  // my_out.close();

  REAL twist_per_meter = 2. * 5. * M_PI;
  // T sqrt(a^2+b^2)
  REAL helix_length =
      max_s * sqrt(helix_radius * helix_radius + pitch_factor * pitch_factor);
  // REAL tot_twist = twist_per_meter * helix_length;
  // REAL tot_twist = 2. * 10. * M_PI;
  // REAL tot_twist = 0. * 10. * M_PI;
  REAL tot_twist = 0.5 * M_PI;
  auto twist_per_step = tot_twist / (n - 1);
  auto twist_values = vector<REAL>(n);
  for (size_t i = 0; i < n; ++i) {
    twist_values[i] = (double)(i)*twist_per_step;
  }

  // // Apply twists now
  // // Apply twists function not working properly too
  // {
  //   // Done this way because no one knows what basis coordinate system we
  //   // are using
  //   // Index i starts from 1, because the first frame has been filled already
  //   REAL d_prefac =
  //       1. / sqrt(helix_radius * helix_radius + pitch_factor * pitch_factor);
  //   for (size_t i = 0; i < n; ++i) {
  //     // Create body frames like so
  //     //    d2  d3
  //     //     | /
  //     // d1 __|/
  //     // where the helix is (cos(t),sin(t),t)
  //     auto new_s = 0.5 * (s[i] + s[i + 1]);
  //     auto d1 = Vector3(-cos(new_s), -sin(new_s), 0.0);
  //     auto d2 = d_prefac * Vector3(pitch_factor * sin(new_s),
  //                                  -pitch_factor * cos(new_s), helix_radius);
  //     auto new_d1 = d1 * cos(twist_values[i]) + d2 * sin(twist_values[i]);
  //     auto new_d2 = -d1 * sin(twist_values[i]) + d2 * cos(twist_values[i]);
  //     auto d3 = d_prefac * Vector3(-helix_radius * sin(new_s),
  //                                  helix_radius * cos(new_s), pitch_factor);
  //     Q0[i] = (Matrix3(new_d1, new_d2, d3));
  //   }
  // }

  // my_out.open("test_mat_native.csv");
  // for (size_t i = 0; i < n; ++i) {
  //   my_out << Q0[i] << "\n";
  // }
  // my_out.close();

  // Set shear vector to zero
  vector<Vector3> new_k0 = vector<Vector3>(n - 1);
  vRotDiff(Q0, tempVM3_nminus);
  vLog(tempVM3_nminus, tempVV3_nminus);
  v_a_divide_b_equal_c(tempVV3_nminus, d0, new_k0);

  // my_out.open("curvature_test_native.csv");
  // for (size_t i = 0; i < n - 1; ++i) {
  //   // auto temp_two = Q0[i].T() * new_k0[i];
  //   auto temp_two = new_k0[i];
  //   my_out << i << "," << temp_two.x << "," << temp_two.y << "," <<
  //   temp_two.z
  //          << "\n";
  // }
  // my_out.close();

  vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n);

  // Mass of vertex point wise element.
  // VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
  // This is why m is obtained dividing the total mass by n, and then to
  // conserve the total mass the masses of the first and last verteces are
  // divideb by 2

  // const REAL m = totalMass/(double)(n);
  // Set volume discretization elements
  vector<REAL> V0 = vector<REAL>(n);
  vector<REAL> masses = vector<REAL>(n + 1);
  REAL totalmass = 0.0;
  for (unsigned int i = 0; i < n; i++) {
    V0[i] = M_PI * r0 * r0 * dl0;
    REAL mass = V0[i] * density;
    masses[i] += mass / 2.0;
    masses[i + 1] += mass / 2.0;
    totalmass += mass;
  }

  // Mass second moment of inertia matrix in rest configuration
  REAL I0_1 = rod_cross_area * rod_cross_area / (4.0 * M_PI);
  REAL I0_2 = I0_1;
  REAL I0_3 = 2.0 * I0_1;

  // clang-format off
  auto J =
    density * dl0 * Matrix3(I0_1, 0.0, 0.0,
    							0.0, I0_2, 0.0,
    							0.0, 0.0, I0_3);

	// Bending matrix
	auto B = Matrix3(	E*I0_1,	0.0,	0.0,
							0.0,	E*I0_2,	0.0,
							0.0,	0.0,	G*I0_3);
 // 0.039270,  0.000000,  0.000000
 // 0.000000,  0.039270,  0.000000
 // 0.000000,  0.000000,  0.052360

	// Shear matrix
	auto S = Matrix3(	(4.0/3.0)*G*rod_cross_area,	0.0,			0.0,
							0.0,			(4.0/3.0)*G*rod_cross_area,	0.0,
							0.0,			0.0,			E*rod_cross_area);
  // clang-format on

  vector<Matrix3> J0 = vector<Matrix3>(n, J);

  // Bending matrix in reference configuration
  vector<Matrix3> B0 = vector<Matrix3>(n - 1, B);

  // Shear matrix in reference configuration
  vector<Matrix3> S0 = vector<Matrix3>(n, S);

  return new Rod(n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0,
                 masses, V0, density, J0, B0, S0, nu, relaxationNu,
                 useSelfContact);
}
