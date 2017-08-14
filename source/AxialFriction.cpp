/*
 * AxialFriction.cpp
 *
 *  Created on: Sep 17, 2015
 *      Author: mgazzola
 */

#include "AxialFriction.h"

AxialFriction::AxialFriction(const int argc, const char ** argv)
{
}

AxialFriction::~AxialFriction()
{
}

bool AxialFriction::_test(REAL& tEnergy, REAL& rEnergy, const REAL F)
{
	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 1;
	const unsigned int povrayPerUnitTime = 0;

	// Driving parameters
	const int n = 50;
	const REAL totalMass = 1.0;
	const REAL L0 = 1.0;
	const REAL r0 = 0.025;
	const REAL timeSimulation = 0.25;
	const REAL dt = 0.00001;
	const REAL E = 1e5;

	const REAL dL0 = L0/(double)n;							// length of cross-section element
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL density = totalMass/Vol;
	const REAL poissonRatio = 0.50;							// Incompressible
	const REAL G = E / (poissonRatio+1.0);

	// Define plane
	const Vector3 originPlane = Vector3(0.0,0.0,0.0);
	const Vector3 normalPlane = Vector3(0.0,sin(0.0),cos(0.0)).unitize();
	const Vector3 tangentPlane = Vector3(1.0,0.0,0.0);
	const Vector3 bitangentPlane = (normalPlane * tangentPlane).unitize();

	// Define rod
	const Vector3 directionRod = Vector3(1.0,0.0,0.0);
	const Vector3 normalRod = Vector3(0.0,0.0,1.0);
	const Vector3 originRod = (originPlane - L0/2.0*directionRod) - bitangentPlane + r0*normalPlane;

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
	Matrix3 B0 = Matrix3(	E*I0_1,	0.0,	0.0,
							0.0,	E*I0_2,	0.0,
							0.0,	0.0,	G*I0_3);

	// Shear matrix
	Matrix3 S0 = 10000*Matrix3(	1.0,	0.0,	0.0,
									0.0,	1.0,	0.0,
									0.0,	0.0,	1.0);

	// Initialize straight rod and pack it into a vector of pointers to rod
	const REAL initialTotalTwist = 0.0;
	const REAL nu = 1e-6;
	const REAL relaxationNu = 0.0;
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	FreeBC freeBC = FreeBC();
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&freeBC);

	// Pack all forces together (no forces applied)
	UniformForces endpointsForce = UniformForces( -(F/n)*directionRod );
	GravityForce gravity = GravityForce(Vector3(0.0,0.0,-9.81));
	MultipleForces multipleForces;
	multipleForces.add(&endpointsForce);
	multipleForces.add(&gravity);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Set up substrate properties and snake-plane interaction object
	const REAL kPlane = 10; // stiffness of ground
	const REAL nuPlane = 1e-4; // viscous damping of ground
	const REAL muKineticForward = 0.4;
	const REAL muKineticBackward = 0.2;
	const REAL muKineticSideways = 0.2;
	const REAL muStaticForward = 2.0*muKineticForward;
	const REAL muStaticBackward = 2.0*muKineticBackward;
	const REAL muStaticSideways = 2.0*muKineticSideways;
	const REAL vStatic = 1e-4;
	AnisotropicFrictionPlaneInteraction frictionPlane = AnisotropicFrictionPlaneInteraction(	rodPtrs, normalPlane, originPlane, kPlane, nuPlane,
																								muKineticForward, muKineticBackward, muKineticSideways,
																								muStaticForward, muStaticBackward, muStaticSideways, vStatic);
	vector<Interaction*> substrateInteractionsPtrs;
	substrateInteractionsPtrs.push_back(&frictionPlane);

	// Set up integrator (define integration order)
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs);

	// Simulate
	Polymer poly = Polymer(integrator);
	string outfileName = string("prova");
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

	rod.computeEnergies();
	tEnergy = rod.translationalEnergy;
	rEnergy = rod.rotationalEnergy;

	return false;
}

void AxialFriction::_sweepAngle(const REAL start, const REAL stop, const REAL step)
{
	REAL tEnergy, rEnergy;

	char buffer[1000];
	sprintf(buffer,"friction.txt");
	FILE * outfile = fopen( buffer, "w" );
	fclose( outfile );

	for (REAL angle = start; angle <= stop; angle += step)
	{
		_test(tEnergy, rEnergy, angle);
		printf("angle = %f, kinetic energy = %10.10e, rotation energy = %10.10e, tot = %10.10e\n", angle, tEnergy, rEnergy, tEnergy+rEnergy);

		sprintf(buffer,"friction.txt");
		FILE * outfile = fopen( buffer, "a" );
		fprintf( outfile, "%1.10e %1.10e %1.10e\n", angle, tEnergy, rEnergy );
		fclose( outfile );
	}

}

void AxialFriction::run()
{
	{
		// Quick test spanning across the jump
		const REAL startAngle = 0.0;
		const REAL stopAngle = 10.0;
		const REAL step = (stopAngle-startAngle)/100.0;
		_sweepAngle(startAngle, stopAngle, step);
	}

	/*
	{
		// Extensive scan
		const REAL startAngle = 0.0;
		const REAL stopAngle = M_PI/2.0;
		const REAL step = (stopAngle-startAngle)/300.0;
		_sweepAngle(startAngle, stopAngle, step);
	}
	*/

	exit(0);
}


