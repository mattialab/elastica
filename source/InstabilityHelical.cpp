/*
 * ValidationLocalizedHelicalBuckling_SpaceConvergence.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: mgazzola
 */

#include "InstabilityHelical.h"

InstabilityHelical::InstabilityHelical(const int argc, const char ** argv)
{
}

InstabilityHelical::~InstabilityHelical()
{
}

bool InstabilityHelical::_testLocalizedHelicalBuckling(const int nEdges, const REAL _timeTwist, const REAL _timeRelax, const REAL _nu, string outfileName)
{
	// Discretization parameters
	const int n = nEdges;									// number of discretization edges (i.e. n+1 points) along the entire rod
	const REAL timeTwist = _timeTwist;						// twisting time
	const REAL timeRelax = _timeRelax;						// relaxing time
	const REAL timeSimulation = timeTwist + timeRelax;		// total simulation time

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 1;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL massPerUnitLength = 1.0;									// mass per unit length
	const REAL L0 = 100.0;												// total length of rod
	const REAL alpha = 1.345;											// flexural rigidity about d1 and d2 (often referred to as B - note that here we use B for the bending matrix instead)
	const REAL beta = 0.789;											// torsional rigidity about d3 (often referred to as C)
	const REAL R = 27.0;												// Number of turns applied at the first extrema (see "Writhing instabilities of twisted rods: from infinite to finite length", 2001)
	const REAL D = 3.0;													// Imposed displacement or slack at the first extrema (see above reference)
	const REAL nu = _nu;												// Numerical dumping viscosity
	const REAL nuRelaxation = 0.0;
	const REAL dL0 = L0/(double)n;										// length of cross-section element
	const REAL dt = 0.001*dL0;											// time step
	const REAL totalMass = massPerUnitLength*L0;							// total mass of the rod
	const REAL r0 = 0.35;//5*dL0;
	const REAL A0 = M_PI*r0*r0;
	const REAL density  = totalMass/(L0*A0);
	const REAL totalInitialTwist = 0.0;
	const Vector3 originRod = Vector3(0.0,0.0,-L0/2.0);
	const Vector3 directionRod = Vector3(0.0,0.0,1.0);
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

	// Bending matrix
	Matrix3 B0 = Matrix3(	alpha,		0.0,	0.0,
							0.0,		alpha,	0.0,
							0.0,		0.0,	beta);

	// Shear matrix --> the high value numerically enforce both unshearability and unstreatchability!
	Matrix3 S0 = 1e5 * Matrix3(	1.0, 0.0, 0.0,
								0.0, 1.0, 0.0,
								0.0, 0.0, 1.0);

	// Initialize straight rod and pack it into a vector of pointers to rod.
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, totalInitialTwist, originRod, directionRod, normalRod, nu, nuRelaxation, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);

	// Perturb rod and update it
	rod.v[floor(n/2)].x+=1e-6;
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	HelicalBucklingBC endBC = HelicalBucklingBC(rodPtrs,timeTwist, D, R);
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	NoForces endpointsForce = NoForces();
	MultipleForces multipleForces;
	multipleForces.add(&endpointsForce);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Empty interaction forces (no substrate in this case)
	vector<Interaction*> substrateInteractionsPtrs;

	// Set up integrator (define integration order)
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs);

	// Simulate
	Polymer poly = Polymer(integrator);
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);
	//const bool goodRun = poly.simulate(220.0, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

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

	// Compute max velocity to make sure it is at quilibrium
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

	return false;
}

void InstabilityHelical::_test(const unsigned int nEdges, string outputfilename)
{
	//const REAL timeTwist = 250.0;
	//const REAL timeRelax = 10000.0;
	const REAL timeTwist = 500.0;
	const REAL timeRelax = 10000.0;
	const REAL nu = 1e-2;
	_testLocalizedHelicalBuckling(nEdges, timeTwist, timeRelax, nu, outputfilename);
}

void InstabilityHelical::run()
{
	{
		cout << "Localized helical buckling: time-space convergnece" << endl << endl;
		_test(100, "helix_0100");
		_test(200, "helix_0200");
		_test(400, "helix_0400");
		_test(800, "helix_0800");
		//_test(1600, "helix_1600");
		//_test(3200, "helix_3200");
		cout << "study completed!" << endl << endl;
	}

	exit(0);
}
