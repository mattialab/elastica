/*
 * StretchRelease.cpp
 *
 *  Created on: May 7, 2015
 *      Author: mgazzola
 */

#include "MassSpringSystem.h"

MassSpringSystem::MassSpringSystem(const int argc, const char ** argv)
{
}

MassSpringSystem::~MassSpringSystem()
{
}

bool MassSpringSystem::_test(	const int nEdges, const REAL _dt, const REAL _M, const REAL _L0, const REAL _A0, const REAL _timeSimulation,
								const REAL _E, const REAL _rho, const string outfileName)
{
	// Input parameters
	const int n = nEdges;									// number of discretization edges (i.e. n+1 points) along the entire rod
	const REAL timeSimulation = _timeSimulation;			// total simulation time
	const REAL dt = _dt;									// time step
	const REAL L0 = _L0;									// total length of rod [m]
	const REAL density = _rho;								// [kg/m^3]
	const REAL A0 = _A0;
	const REAL E = _E;										// GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
	const REAL M = _M;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = ceil(1000.0/timeSimulation);
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL totalMass = L0*A0*density;
	const REAL dL0 = L0/(double)n;							// length of cross-section element
	const REAL r0 = sqrt(A0/M_PI);							// radius [m]
	const REAL nu = 0.0;
	const REAL relaxationNu = 0.0;
	const REAL g = 9.81;
	const REAL poissonRatio = 0.50;							// Incompressible
	const REAL G = E / (poissonRatio+1.0);					// Shear modulus
	const REAL initialTotalTwist = 0.0;
	const Vector3 originRod = Vector3(0.0,0.0,0.0);
	const Vector3 directionRod = Vector3(1.0,0.0,0.0);
	const Vector3 normalRod = Vector3(0.0,0.0,1.0);

	// Second moment of area for disk cross section
	const REAL I0_1 = A0*A0/(4.0*M_PI);
	const REAL I0_2 = I0_1;
	const REAL I0_3 = 2.0*I0_1;
	const Matrix3 I0 = Matrix3(	I0_1,	 0.0,	 0.0,
								 0.0,	I0_2,	 0.0,
								 0.0,	 0.0,	I0_3);

	// Mass inertia matrix for disk cross section
	const Matrix3 J0 = density*dL0*I0;

	// Bending matrix (TOD: change this is wrong!!)
	Matrix3 B0 = Matrix3(	E*I0_1,	0.0,	0.0,
							0.0,	E*I0_2,	0.0,
							0.0,	0.0,	G*I0_3);

	// Shear matrix
	Matrix3 S0 = Matrix3(	G*A0*4.0/3.0,	0.0,			0.0,
							0.0,			G*A0*4.0/3.0,	0.0,
							0.0,			0.0,			E*A0);

	// Initialize straight rod and pack it into a vector of pointers to rod.
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.m.back() += M;
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	StretchReleaseBC endBC = StretchReleaseBC(rodPtrs);
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	GravityForce gravity = GravityForce(g*directionRod);
	MultipleForces multipleForces;
	multipleForces.add(&gravity);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Empty interaction forces (no substrate in this case)
	vector<Interaction*> substrateInteractionsPtrs;

	// Set up time integrator
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs);

	// Simulate
	Polymer poly = Polymer(integrator);

	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

	return false;
}

void MassSpringSystem::_largeAttachedMassTest(const int nEdges, const REAL E, const string outputdata)
{
	const REAL rho = 1000.0;
	const REAL L0 = 1.0;
	const REAL m = 1.0;
	const REAL M = 100.0;
	const REAL A0 = m/(L0*rho);
	const REAL D = 3.0;
	const REAL f = sqrt(E*A0/(M + m/D))/(2.0*M_PI);
	const REAL simTime = 1.0/f;
	const REAL dt = simTime/1000000;
	_test(nEdges, dt, M, L0, A0, simTime, E, rho, outputdata);
}

void MassSpringSystem::_smallAttachedMassTest(const int nEdges, const REAL E, const string outputdata)
{
	const REAL rho = 1000.0;
	const REAL L0 = 1.0;
	const REAL m = 1;
	const REAL M = 0.0;
	const REAL A0 = m/(L0*rho);
	const REAL D = (M_PI*M_PI)/4.0;
	const REAL f = sqrt(E*A0/(M + m/D))/(2.0*M_PI);
	const REAL simTime = 1.0/f;
	const REAL dt = simTime/1000000;
	_test(nEdges, dt, M, L0, A0, simTime, E, rho, outputdata);
}

void MassSpringSystem::run()
{
	/*
	{
		// Large mass case
		_largeAttachedMassTest( 100,  1e7, "largemass_1E07");
		_largeAttachedMassTest( 100,  2e7, "largemass_2E07");
		_largeAttachedMassTest( 100,  3e7, "largemass_3E07");
		//_largeAttachedMassTest( 100,  4e7, "largemass_4E07");
		_largeAttachedMassTest( 100,  5e7, "largemass_5E07");
		//_largeAttachedMassTest( 100,  6e7, "largemass_6E07");
		//_largeAttachedMassTest( 100,  7e7, "largemass_7E07");
		//_largeAttachedMassTest( 100,  8e7, "largemass_8E07");
		//_largeAttachedMassTest( 100,  9e7, "largemass_9E07");

		_largeAttachedMassTest( 100,  1e8, "largemass_1E08");
		//_largeAttachedMassTest( 100,  2e8, "largemass_2E08");
		//_largeAttachedMassTest( 100,  3e8, "largemass_3E08");
		//_largeAttachedMassTest( 100,  4e8, "largemass_4E08");
		//_largeAttachedMassTest( 100,  5e8, "largemass_5E08");
		//_largeAttachedMassTest( 100,  6e8, "largemass_6E08");
		//_largeAttachedMassTest( 100,  7e8, "largemass_7E08");
		//_largeAttachedMassTest( 100,  8e8, "largemass_8E08");
		//_largeAttachedMassTest( 100,  9e8, "largemass_9E08");

		//_largeAttachedMassTest( 100,  1e9, "largemass_1E09");
		//_largeAttachedMassTest( 100,  2e9, "largemass_2E09");
		//_largeAttachedMassTest( 100,  3e9, "largemass_3E09");
		//_largeAttachedMassTest( 100,  4e9, "largemass_4E09");
		//_largeAttachedMassTest( 100,  5e9, "largemass_5E09");
		//_largeAttachedMassTest( 100,  6e9, "largemass_6E09");
		//_largeAttachedMassTest( 100,  7e9, "largemass_7E09");
		//_largeAttachedMassTest( 100,  8e9, "largemass_8E09");
		//_largeAttachedMassTest( 100,  9e9, "largemass_9E09");

		_largeAttachedMassTest( 100, 1e10, "largemass_1E10");
	}
	*/

	{
		// Small mass case
		_smallAttachedMassTest( 100,   1e4, "smallmass_1E04");
		_smallAttachedMassTest( 100,   2e4, "smallmass_2E04");
		_smallAttachedMassTest( 100,   3e4, "smallmass_3E04");
		//_smallAttachedMassTest( 100,   4e4, "smallmass_4E04");
		_smallAttachedMassTest( 100,   5e4, "smallmass_5E04");
		//_smallAttachedMassTest( 100,   6e4, "smallmass_6E04");
		//_smallAttachedMassTest( 100,   7e4, "smallmass_7E04");
		//_smallAttachedMassTest( 100,   8e4, "smallmass_8E04");
		//_smallAttachedMassTest( 100,   9e4, "smallmass_9E04");

		_smallAttachedMassTest( 100,   1e5, "smallmass_1E05");
		_smallAttachedMassTest( 100,   2e5, "smallmass_2E05");
		//_smallAttachedMassTest( 100,   3e5, "smallmass_3E05");
		//_smallAttachedMassTest( 100,   4e5, "smallmass_4E05");
		//_smallAttachedMassTest( 100,   5e5, "smallmass_5E05");
		//_smallAttachedMassTest( 100,   6e5, "smallmass_6E05");
		//_smallAttachedMassTest( 100,   7e5, "smallmass_7E05");
		//_smallAttachedMassTest( 100,   8e5, "smallmass_8E05");
		//_smallAttachedMassTest( 100,   9e5, "smallmass_9E05");

		_smallAttachedMassTest( 100,   1e6, "smallmass_1E06");
		//_smallAttachedMassTest( 100,   2e6, "smallmass_2E06");
		//_smallAttachedMassTest( 100,   3e6, "smallmass_3E06");
		//_smallAttachedMassTest( 100,   4e6, "smallmass_4E06");
		//_smallAttachedMassTest( 100,   5e6, "smallmass_5E06");
		//_smallAttachedMassTest( 100,   6e6, "smallmass_6E06");
		//_smallAttachedMassTest( 100,   7e6, "smallmass_7E06");
		//_smallAttachedMassTest( 100,   8e6, "smallmass_8E06");
		//_smallAttachedMassTest( 100,   9e6, "smallmass_9E06");

		_smallAttachedMassTest( 100,   1e7, "smallmass_1E07");
		//_smallAttachedMassTest( 100,   1e8, "smallmass_1E08");
		//_smallAttachedMassTest( 100,   1e9, "smallmass_1E09");
	}

	cout << "Mass spring test case concluded!" << endl;

	exit(0);
}
