/*
 * RodInitialConfiguations.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: mgazzola
 */

#include "RodInitialConfigurations.h"

Rod RodInitialConfigurations::straightRod(	const int n, const REAL totalMass, const REAL r0, const Matrix3 _J0, const Matrix3 _B0, const Matrix3 _S0, const REAL L0, const REAL totTwist,
		const Vector3 origin, const Vector3 direction, const Vector3 normal, const REAL nu, const REAL relaxationNu, const bool useSelfContact)
{
	// Bunch of sanity checks
	assert(n>1);
	assert(totalMass>Tolerance::tol());
	assert(r0>Tolerance::tol());
	assert(L0>Tolerance::tol());
	assert(nu>=0.0);
	assert(relaxationNu>=0.0);
	assert(direction.length()>Tolerance::tol());
	assert(normal.length()>Tolerance::tol());

	assert(_B0.r2c1==0.0);
	assert(_B0.r3c1==0.0);
	assert(_B0.r1c2==0.0);
	assert(_B0.r3c2==0.0);
	assert(_B0.r1c3==0.0);
	assert(_B0.r2c3==0.0);
	assert(_B0.r1c1>=Tolerance::tol());
	assert(_B0.r2c2>=Tolerance::tol());
	assert(_B0.r3c3>=Tolerance::tol());

	assert(_S0.r2c1==0.0);
	assert(_S0.r3c1==0.0);
	assert(_S0.r1c2==0.0);
	assert(_S0.r3c2==0.0);
	assert(_S0.r1c3==0.0);
	assert(_S0.r2c3==0.0);
	assert(_S0.r1c1>=Tolerance::tol());
	assert(_S0.r2c2>=Tolerance::tol());
	assert(_S0.r3c3>=Tolerance::tol());

	// Density
	const REAL density = totalMass / (M_PI*r0*r0*L0);

	// Initialize discretization points
	vector<Vector3> x0 = vRange(origin, origin + direction * L0, n+1);

	// Set velocities to zero
	vector<Vector3> v0 = vector<Vector3>(n+1);

	// u0 is a vector orthogonal to the tangent, .i.e. the direction.
	// If not provided I just compute one, using a random non-zero vector.
	// Note that to be able to compare codes, I fix the random vector and
	// I do not generate one, because that is BAD!
#ifndef NDEBUG
	const Vector3 t0 = (x0[1]-x0[0]).unitize();
	assert((t0*normal).length()>Tolerance::tol());
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
	const vector<REAL> V0 = vector<REAL>(n, M_PI*r0*r0*dl0);

	// Set shear vector to zero
	vector<Vector3> intrinsicShearStrain0 = vector<Vector3>(n);

	// Mass of vertex point wise element.
	// VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
	// This is why m is obtained dividing the total mass by n, and then to conserve the total mass
	// the masses of the first and last verteces are divideb by 2
	const REAL m = totalMass/(double)(n);
	vector<REAL> masses = vector<REAL>(n+1,m);
	masses.front() /= 2.0;
	masses.back() /= 2.0;

	// Intrinsic curvature in rest configuration
	const vector<Vector3> intrinsic_k0 = vector<Vector3>(n-1);

	// Mass second moment of inertia matrix in rest configuration
	vector<Matrix3> J0 = vector<Matrix3>(n, _J0);

	// Bending matrix in reference configuration
	vector<Matrix3> B0 = vector<Matrix3>(n-1, _B0);

	// Shear matrix in reference configuration
	vector<Matrix3> S0 = vector<Matrix3>(n, _S0);

	return Rod(	n, x0, v0, Q0, w0, l0, intrinsic_k0, intrinsicShearStrain0, masses, V0, density, J0, B0, S0, nu, relaxationNu, useSelfContact);
}

// A rod on the z-axis compressed so as to push back with forces F
Rod RodInitialConfigurations::compressedRod(const int n, const REAL totalMass, const REAL r, const Matrix3 J, const Matrix3 B, const Matrix3 S, const REAL L, const REAL F,
											const REAL totTwist, const REAL nu, const Vector3 origin, const Vector3 direction, const Vector3 normal, const bool useSelfContact)
{
	// Bunch of sanity checks
	assert(n>1);
	assert(totalMass>Tolerance::tol());
	assert(r>Tolerance::tol());
	assert(L>Tolerance::tol());
	assert(nu>=0.0);
	assert(direction.length()>Tolerance::tol());
	assert(normal.length()>Tolerance::tol());

	assert(B.r2c1==0.0);
	assert(B.r3c1==0.0);
	assert(B.r1c2==0.0);
	assert(B.r3c2==0.0);
	assert(B.r1c3==0.0);
	assert(B.r2c3==0.0);
	assert(B.r1c1>=Tolerance::tol());
	assert(B.r2c2>=Tolerance::tol());
	assert(B.r3c3>=Tolerance::tol());

	assert(S.r2c1==0.0);
	assert(S.r3c1==0.0);
	assert(S.r1c2==0.0);
	assert(S.r3c2==0.0);
	assert(S.r1c3==0.0);
	assert(S.r2c3==0.0);
	assert(S.r1c1>=Tolerance::tol());
	assert(S.r2c2>=Tolerance::tol());
	assert(S.r3c3>=Tolerance::tol());

	// Density
	const REAL density = totalMass / (M_PI*r*r*L);

	//fractional compression of rod needed to resist force F
	const REAL c = F/S[2][2];

	// Initialize discretization points
	vector<Vector3> x = vRange(origin, origin + direction * L * (1 - c), n + 1);

	// Set velocities to zero
	vector<Vector3> v = vector<Vector3>(n+1);

	// u0 is a vector orthogonal to the tangent, .i.e. the direction.
	// If not provided I just compute one, using a random non-zero vector.
	// Note that to be able to compare codes, I fix the random vector and
	// I do not generate one, because that is BAD!
#ifndef NDEBUG
	const Vector3 t0 = (x[1]-x[0]).unitize();
	assert((t0*normal).length()>Tolerance::tol());
#endif

	// Now I can align frames using the orthogonal vector provided/generated
	vector<Matrix3> Q = alignFrames(x, normal.unitize());

	// Apply twist about the orthonormal vector d3
	applyTwists(Q, vRange((REAL)0.0, totTwist, n));

	// Set angular velocity to zero
	vector<Vector3> w = vector<Vector3>(n);

	// Set rest edge lengths
	vector<REAL> l = vector<REAL>(n, L/n);

	// Set volume discretization elements
	vector<REAL> V = vector<REAL>(n, l[0]*r*r*M_PI);

	// Set shear vector to zero
	vector<Vector3> sigma0 = vector<Vector3>(n);

	const REAL relaxationNu = 0.0;

	// Mass of vertex point wise element.
	// VERY IMPORTANT TO BE CONSISTENT WITH CYLINDER MASS, BETTER CONVERGENCE!
	// This is why m is obtained dividing the total mass by n, and then to conserve the total mass
	// the masses of the first and last verteces are divideb by 2
	const REAL m = totalMass/(double)(n);
	vector<REAL> masses = vector<REAL>(n+1,m);
	masses.front() /= 2.0;
	masses.back() /= 2.0;

	// Intrinsic curvature in rest configuration
	const vector<Vector3> intrinsic_k0 = vector<Vector3>(n-1);

	// Mass second moment of inertia matrix in rest configuration
	vector<Matrix3> J0 = vector<Matrix3>(n, J);

	// Bending matrix in reference configuration
	vector<Matrix3> B0 = vector<Matrix3>(n-1, B);

	// Shear matrix in reference configuration
	vector<Matrix3> S0 = vector<Matrix3>(n, S);

	return Rod(	n, x, v, Q, w, l,intrinsic_k0, sigma0, masses, V, density, J0, B0, S0, nu, relaxationNu, useSelfContact);
}

Rod RodInitialConfigurations::circleRod(const int n, const REAL totalMass, const REAL r, const Matrix3 I, const Matrix3 B, const Matrix3 S, const REAL L, const REAL totTwist, const REAL nu, const bool useSelfContact)
{
	// Density
	const REAL density = totalMass / (M_PI*r*r*L);

	// Masses
	const REAL m = totalMass/(double)(n);
	vector<REAL> masses = vector<REAL>(n+1,m);
	masses.front() /= 2.0;
	masses.back() /= 2.0;

	vector<Vector3> x;

	for (int i=0; i<n+1; i++)
		x.push_back(L/(2*M_PI)*Vector3(cos(2*M_PI/n*i), sin(2*M_PI/n*i), 0));

	vector<Vector3> v = vector<Vector3>(n+1);
	vector<Matrix3> Q = alignFrames(x, Vector3(0,0,1));
	applyTwists(Q, vRange((REAL)0.0, totTwist, n));
	vector<Vector3> w = vector<Vector3>(n);
	vector<REAL> l = vLength(vDiff(x));
	vector<REAL> V = vector<REAL>(n, l[0]*r*r*M_PI);
	vector<Vector3> k0 = vector<Vector3>(n-1);
	vector<Vector3> sigma0 = vector<Vector3>(n);

	const REAL relaxationNu = 0.0;

	return Rod(	n, x, v, Q, w, l, k0, sigma0, masses, V, density, vector<Matrix3>(n, I), vector<Matrix3>(n-1, B), vector<Matrix3>(n, S), nu, relaxationNu, useSelfContact);
}


