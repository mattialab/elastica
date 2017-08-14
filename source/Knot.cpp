/*
 * Knot.cpp
 *
 *  Created on: Nov 28, 2015
 *      Author: ldudte
 */

#include "Knot.h"

Knot::Knot(const int argc, const char ** argv)
{
}

Knot::~Knot()
{
}

bool Knot::_testLocalizedHelicalBuckling(const int nEdges, const REAL _timeTwist, const REAL _timeRelax, const REAL _nu, string outfileName)
{

	// Initialize profiler
	MRAG::Profiler profiler;

	// Discretization parameters
	const int n = nEdges;									// number of discretization edges (i.e. n+1 points) along the entire rod
	const REAL timeTwist = _timeTwist;						// twisting time
	const REAL timeRelax = _timeRelax;						// relaxing time
	const REAL timeSimulation = timeTwist + timeRelax;		// total simulation time

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 1;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL massPerUnitLength = .05;									// mass per unit length
	const REAL L0 = 1.0;												// total length of rod
	//const REAL alpha = 1.345;											// flexural rigidity about d1 and d2 (often referred to as B - note that here we use B for the bending matrix instead)
	//const REAL beta = 0.789;											// torsional rigidity about d3 (often referred to as C)
	const REAL A = 3.0;													// Shake magnitude
	const REAL N = 200;													// Number of shakes
	const REAL nu = _nu;												// Numerical dumping viscosity
	const REAL nuRelaxation = 0.0;
	const REAL dL0 = L0/(double)n;										// length of cross-section element
	const REAL dt = 0.0005*dL0;											// time step
	const REAL totalMass = massPerUnitLength*L0;						// total mass of the rod
	const REAL poissonRatio = 0.50;							// Incompressible
	const REAL E = 1e5;
	const REAL G = E / (poissonRatio+1.0);					// Shear modulus
	const REAL r0 = 0.005;//5*dL0;
	const REAL A0 = M_PI*r0*r0;
	const REAL density  = totalMass/(L0*A0);
	const REAL totalInitialTwist = 0.0;
	const Vector3 originRod = Vector3(0.0,0.0,L0/2.0);
	const Vector3 directionRod = Vector3(0.0,0.0,-1.0);
	const Vector3 normalRod = Vector3(1.0,0.0,0.0);

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
	const REAL attachedMass = 0.0;
	const bool useSelfContact = true;
	const bool useLinearLoadStrain = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, totalInitialTwist, originRod, directionRod, normalRod, nu, nuRelaxation, attachedMass, useSelfContact, useLinearLoadStrain);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);

	// Perturb rod and update it
	rod.v[floor(n/2)].x+=1e-6;
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	KnotBC endBC = KnotBC(rodPtrs, timeTwist, A, N);
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	GravityForce gravity = GravityForce(9.8*Vector3(0.0,0.0,-1.0));
	MultipleForces multipleForces;
	multipleForces.add(&gravity);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Empty interaction forces (no substrate in this case)
	vector<Interaction*> substrateInteractionsPtrs;

	// Set up integrator (define integration order)
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs, &profiler);

	// Simulate
	Polymer poly = Polymer(integrator, &profiler);
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);
	//const bool goodRun = poly.simulate(220.0, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in knot formation, what is going on?";

	// Compute maximum envelope to be compared to analytical solution (See MATLAB script "vanHejiden.m")
	VV3 tangents = vUnitize(vDiff(rod.x));
	assert(tangents.size()==rod.l0.size());
	REAL maxPhi = numeric_limits<REAL>::min();
	for(unsigned int i=0; i<tangents.size(); i++)
	{
		const REAL argumentACos = tangents[i]%directionRod;
		assert( argumentACos>=-1.0-Tolerance::tol() && argumentACos<=1.0+Tolerance::tol() );
		const REAL argumentACosClamped = max((REAL)-1.0,(REAL)min(argumentACos,(REAL)1.0));
		const REAL phi = acos(argumentACosClamped);
		maxPhi = max(maxPhi,phi);
	}

	// Compute max velocity to make sure it is at equilibrium
	REAL maxVel = 0.0;
	for(unsigned int i=0; i<rod.v.size(); i++)
		maxVel = max(maxVel,fabs(rod.v[i].length()));

	// Dump post buckling helical shape
	{
		string helixShape = outfileName + "_shape.txt";
		FILE * fitnessFile = fopen(helixShape.c_str(),"w");
		assert(fitnessFile!=NULL);

		for(unsigned int i=0; i<rod.x.size(); i++)
			fprintf(fitnessFile,"%1.15e %1.15e %1.15e\n", rod.x[i].x, rod.x[i].y, rod.x[i].z);

		fclose(fitnessFile);
	}

	cout << "TENSION AT THE ENDS" << endl;
	cout << rod.totalInternalForces.front() <<endl;
	cout << rod.totalInternalForces.back() <<endl;

	cout << "TORQUE AT THE ENDS" << endl;
	cout << rod.totalInternalTorques.front() <<endl;
	cout << rod.totalInternalTorques.back() <<endl;

	profiler.printSummary();

	return false;
}

void Knot::_test(const unsigned int nEdges, string outputfilename)
{
	//const REAL timeTwist = 250.0;
	//const REAL timeRelax = 10000.0;
	const REAL timeShake = 100.0;
	const REAL timeRelax = 20.0;
	const REAL nu = 1e-2;
	_testLocalizedHelicalBuckling(nEdges, timeShake, timeRelax, nu, outputfilename);
}

void Knot::run()
{
	testProfiler.push_start("EXEC");

	{
		cout << "Knot formation" << endl << endl;
		_test(40, "knot_040");
		//_test(200, "helix_0200");
		//_test(400, "helix_0400");
		//_test(800, "helix_0800");
		//_test(1600, "helix_1600");
		//_test(3200, "helix_3200");
		cout << "study completed!" << endl << endl;
	}

	testProfiler.pop_stop();
	testProfiler.printSummary();
	exit(0);
}

