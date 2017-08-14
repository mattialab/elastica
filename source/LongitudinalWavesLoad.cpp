/*
 * LongitudinalWavesLoad.cpp
 *
 *  Created on: May 14, 2015
 *      Author: mgazzola
 */

#include "LongitudinalWavesLoad.h"

LongitudinalWavesLoad::LongitudinalWavesLoad(const int argc, const char ** argv)
{
}

LongitudinalWavesLoad::~LongitudinalWavesLoad()
{
}

bool LongitudinalWavesLoad::_test(	const int nEdges, const REAL _dt, const REAL _L, const REAL _r, const REAL _f, const REAL _F, const REAL _timeSimulation, const REAL _E, const REAL _rho,
									const REAL _nu, const REAL _relaxationNu, const string outfileName)
{
	// Initialize profiler
	MRAG::Profiler profiler;

	// Input parameters
	const int n = nEdges;									// number of discretization edges (i.e. n+1 points) along the entire rod
	const REAL timeSimulation = _timeSimulation;			// total simulation time
	const REAL dt = _dt;									// time step
	const REAL f = _f;										// frequency [1/s]
	const REAL L0 = _L;										// total length of rod [m]
	const REAL density = _rho;								// [kg/m^3]
	const REAL E = _E;										// GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
	const REAL r0 = _r;
	const REAL nu = _nu;									// Numerical damping viscosity [m^2/s]
	const REAL relaxationNu = _relaxationNu;				// relaxation time for exponential decay of nu
	const REAL F = _F;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = ceil(11.0*f);
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL omega = 2.0*M_PI*f;							// angular frequency
	const REAL dL0 = L0/(double)n;							// length of cross-section element
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL totalMass = Vol * density;
	const REAL poissonRatio = 0.50;							// Incompressible
	const REAL G = E / (poissonRatio+1.0);
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

	// Initialize straight rod and pack it into a vector of pointers to rod --> Use linear load-strain (hence the true flag at the end)!!!
	const REAL attachedMass = 0.0;
	const bool useSelfContact = false;
	const bool useLinearLoadStrain = true;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, attachedMass, useSelfContact, useLinearLoadStrain);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	LongitudinalWavesLoadBC endBC = LongitudinalWavesLoadBC(rodPtrs);
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	PeriodicLoad endpointsForce = PeriodicLoad(F*directionRod, omega);
	MultipleForces multipleForces;
	multipleForces.add(&endpointsForce);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Empty interaction forces (no substrate in this case)
	vector<Interaction*> substrateInteractionsPtrs;

	// Set up time integrator
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs, &profiler);

	// Simulate
	Polymer poly = Polymer(integrator, &profiler);
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

	profiler.printSummary();

	return false;
}

void LongitudinalWavesLoad::_longWaveTest(const int nEdges, const string outputdata)
{
	// Very important! The viscosity should be zero at the end of the simulations and should start rather high
	// to get rid of initial spurious oscillations due to the fact that the simulation is started impulsively.
	// That is why nu has an initial high non realistic value, it doesnt matter since it will decay exponentially to zero.
	// In order to get good convergence one needs to play with the decaying period and with the initial viscosity
	// to make sure that for all simulations initial oscillations are wiped out and viscosity is zero at the end
	const REAL f = 1.0;
	const REAL rho = 10;
	const REAL E = 1e6;
	const REAL F = 1e7;
	const REAL L = 2.0*sqrt(E/rho)/(2.0*f);
	const REAL r = 0.01*L;
	const REAL T = 1.0/f;
	const REAL dL = L / nEdges;
	const REAL dt = 0.0001*dL;
	const REAL simTime = 2000*T;
	const REAL nu = 1000.0;
	const REAL halfLife = simTime / 10.0;
	const REAL relaxationNu = halfLife/log(2.0);
	_test(nEdges, dt, L, r, f, F, simTime, E, rho, nu, relaxationNu, outputdata);
}

void LongitudinalWavesLoad::run()
{
	testProfiler.push_start("EXEC");

	_longWaveTest( 10, "longwaves_0010");
	_longWaveTest( 25, "longwaves_0025");
	_longWaveTest( 50, "longwaves_0050");
	_longWaveTest( 100, "longwaves_0100");
	_longWaveTest( 200, "longwaves_0200");
	_longWaveTest( 400, "longwaves_0400");
	_longWaveTest( 800, "longwaves_0800");
	_longWaveTest( 1600, "longwaves_1600");

	testProfiler.pop_stop();
	testProfiler.printSummary();
	exit(0);
}
