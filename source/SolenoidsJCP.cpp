/*
 * Solenoids.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: mgazzola
 */

#include "SolenoidsJCP.h"

SolenoidsJCP::SolenoidsJCP(const int argc, const char ** argv)
{
}

SolenoidsJCP::~SolenoidsJCP()
{
}

bool SolenoidsJCP::_test(const int nEdges, const REAL _E, const REAL _R, const REAL _F, const REAL _timeTwist, const REAL _timeRelax, const REAL _nu, string outfileName)
{
	// Discretization parameters
	const int n = nEdges;									// number of discretization edges (i.e. n+1 points) along the entire rod
	const REAL timeTwist = _timeTwist;						// twisting time
	const REAL timeRelax = _timeRelax;						// relaxing time
	const REAL timeSimulation = timeTwist + timeRelax;		// total simulation time

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 1;
	const unsigned int povrayPerUnitTime = 25;

	// Physical parameters
	const REAL massPerUnitLength = 1.0;									// mass per unit length
	const REAL L0 = 1.0;												// total length of rod
	const REAL F = _F;
	const REAL R = _R;												// Number of turns applied at the first extrema (see "Writhing instabilities of twisted rods: from infinite to finite length", 2001)
	const REAL nu = _nu;												// Numerical dumping viscosity
	const REAL dL0 = L0/(double)n;										// length of cross-section element
	const REAL dt = 0.01*dL0;											// time step
	const REAL totalMass = massPerUnitLength*L0;						// total mass of the rod
	const REAL poissonRatio = 0.50;							// Incompressible
	const REAL E = _E;
	const REAL G = E / (poissonRatio+1.0);					// Shear modulus
	const REAL r0 = 0.025;
	const REAL A0 = M_PI*r0*r0;
	const REAL density  = totalMass/(L0*A0);
	const REAL twistFraction = 0.0;
	const REAL totalInitialTwist = twistFraction*R*2*M_PI;
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
	const bool useSelfContact = true;
	Rod rod = RodInitialConfigurations::compressedRod(n, totalMass, r0, J0, B0, S0, L0, -F, totalInitialTwist, nu, originRod, directionRod, normalRod, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);

	// Perturb rod and update it
	rod.v[floor(n/2)].x+=1e-6;
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	SolenoidsBC_JCP endBC = SolenoidsBC_JCP(rodPtrs,timeTwist, R*(1-twistFraction));
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	EndpointForces endpointForce = EndpointForces(Vector3(0,0,0), F*directionRod);
	RandomForces rf = RandomForces(0.01);
	MultipleForces multipleForces;
	multipleForces.add(&endpointForce);
	//multipleForces.add(&rf);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Empty interaction forces (no substrate in this case)
	vector<Interaction*> substrateInteractionsPtrs;

	// Set up integrator (define integration order)
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs);

	// Simulate
	Polymer poly = Polymer(integrator);
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

	// Dump post buckling helical shape
	{
		string helixShape = outfileName + "_shape.txt";
		FILE * fitnessFile = fopen(helixShape.c_str(),"w");
		assert(fitnessFile!=NULL);

		for(unsigned int i=0; i<rod.x.size(); i++)
			fprintf(fitnessFile,"%1.15e %1.15e %1.15e\n", rod.x[i].x, rod.x[i].y, rod.x[i].z);

		fclose(fitnessFile);
	}

	return false;
}

void SolenoidsJCP::run()
{
	// Plectoneme
	{
		const int nEdges = 100;
		const REAL E = 1e6;
		const REAL R = 4;
		const REAL F = 0;
		const REAL timeTwist = 15*R;
		const REAL nu = 2;//1.0/nEdges;
		const REAL timeRelax = 5;
		_test(nEdges, E, R, F, timeTwist, timeRelax, nu, "solenoid");
	}

	/*
	// Solenoid
	{
		const int nEdges = 100;
		const REAL E = 1e6;
		const REAL R = 13;
		const REAL F = 300;
		const REAL timeTwist = 5*R;
		const REAL nu = 2;
		const REAL timeRelax = 50.0;
		_test(nEdges, E, R, F, timeTwist, timeRelax, nu, "solenoid");
	}
	*/

	exit(0);
}

