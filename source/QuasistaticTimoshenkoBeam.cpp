/*
 * BendingWavesCouple.cpp
 *
 *  Created on: Jun 21, 2015
 *      Author: mgazzola
 */

#include "QuasistaticTimoshenkoBeam.h"

QuasistaticTimoshenkoBeam::QuasistaticTimoshenkoBeam(const int argc, const char ** argv)
{
}

QuasistaticTimoshenkoBeam::~QuasistaticTimoshenkoBeam()
{
}

bool QuasistaticTimoshenkoBeam::_test(	const int nEdges, const REAL _dt, const REAL _L, const REAL _r, const REAL _P, const REAL _timeSimulation, const REAL _E, const REAL _G, const REAL _rho,
										const REAL _nu, const REAL _relaxationNu, const string outfileName)
{
	// Input parameters
	const int n = nEdges;									// number of discretization edges (i.e. n+1 points) along the entire rod
	const REAL timeSimulation = _timeSimulation;			// total simulation time
	const REAL dt = _dt;									// time step
	const REAL P = _P;
	const REAL L0 = _L;										// total length of rod [m]
	const REAL r0 = _r;										// radius [m]
	const REAL density = _rho;								// [kg/m^3]
	const REAL E = _E;										// GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
	const REAL G = _G;										// GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
	const REAL nu = _nu;									// Numerical damping viscosity [m^2/s]
	const REAL relaxationNu = _relaxationNu;				// relaxation time for exponential decay of nu

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 5;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL dL0 = L0/(double)n;							// length of cross-section element
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL totalMass = Vol * density;
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
	Matrix3 S0 = Matrix3(	(4.0/3.0)*G*A0,	0.0,			0.0,
							0.0,			(4.0/3.0)*G*A0,	0.0,
							0.0,			0.0,			E*A0);

	// Initialize straight rod and pack it into a vector of pointers to rod --> Use linear load-strain (hence the true flag at the end)!!!
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	TimoshenkoBeamBC endBC = TimoshenkoBeamBC(rodPtrs);
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	GradualEndpointForces endpointsForce = GradualEndpointForces(Vector3(), Vector3(0,-P,0), timeSimulation/2);
	MultipleForces multipleForces;
	multipleForces.add(&endpointsForce);
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

	// Dump post buckling helical shape
	{
		string rodShape = outfileName + "_shape.txt";
		FILE * fitnessFile = fopen(rodShape.c_str(),"w");
		assert(fitnessFile!=NULL);

		for(unsigned int i=0; i<rod.x.size(); i++)
			fprintf(fitnessFile,"%1.15e %1.15e %1.15e\n", rod.x[i].x, rod.x[i].y, rod.x[i].z);

		fclose(fitnessFile);
	}

	cout <<  "total internal energy = " << poly.getTotalEnergy() << endl;
	cout <<  "total translational energy = " << poly.getTotalTranslationalEnergy() << endl;
	cout <<  "total rotational energy = " << poly.getTotalRotationalEnergy() << endl;

	return false;
}

void QuasistaticTimoshenkoBeam::_longWaveTest(const int nEdges, const string outputdata)
{
	/*
	const REAL rho = 1000;
	const REAL L = 1.0;
	const REAL poissonRatio = 0.5; // incompressible material
	const REAL E = 1e9;
	const REAL G = E/(2.0*(1.0+poissonRatio));
	const REAL nu = 5e-2;
	const REAL P = 0.0005;
	const REAL dL = L / nEdges;
	const REAL r = 0.01;
	const REAL dt = 0.001*dL;
	const REAL simTime = 10.0;
	const REAL relaxationNu = 0.0;
	_test(nEdges, dt, L, r, P, simTime, E, G, rho, nu, relaxationNu, outputdata);
	*/

	/*
	const REAL rho = 5000;
	const REAL L = 2.0;
	const REAL E = 1e6;
	const REAL G = 1e3;
	const REAL nu = 1e-1;
	const REAL minNu = 0.0;
	const REAL P = 1.0;
	const REAL dL = L / nEdges;
	const REAL r = 0.2;
	const REAL dt = 0.03*dL;
	const REAL simTime = 20000;
	const REAL halfLife = 10.0;
	const REAL relaxationNu = 0.0;
	_test(nEdges, dt, L, r, P, simTime, E, G, rho, nu, relaxationNu, minNu, outputdata);
	*/

	const REAL rho = 5000;
	const REAL L = 3.0;
	const REAL E = 1e6;
	const REAL G = 1e4;
	const REAL nu = 1e-1;
	const REAL P = 15;
	const REAL dL = L / nEdges;
	const REAL r = 0.25;
	const REAL dt = 0.01*dL;
	const REAL simTime = 5000;
	const REAL relaxationNu = 0.0;
	_test(nEdges, dt, L, r, P, simTime, E, G, rho, nu, relaxationNu, outputdata);
}

void QuasistaticTimoshenkoBeam::run()
{
	_longWaveTest( 5, "longwaves_0005");
	_longWaveTest( 6, "longwaves_0006");
	_longWaveTest( 7, "longwaves_0007");
	_longWaveTest( 8, "longwaves_0008");
	_longWaveTest( 9, "longwaves_0009");
	_longWaveTest( 10, "longwaves_0010");
	_longWaveTest( 20, "longwaves_0020");
	_longWaveTest( 30, "longwaves_0030");
	_longWaveTest( 40, "longwaves_0040");
	_longWaveTest( 50, "longwaves_0050");
	_longWaveTest( 60, "longwaves_0060");
	_longWaveTest( 70, "longwaves_0070");
	_longWaveTest( 80, "longwaves_0080");
	_longWaveTest( 90, "longwaves_0090");
	_longWaveTest( 100, "longwaves_0100");
	_longWaveTest( 200, "longwaves_0200");
	_longWaveTest( 400, "longwaves_0400");

	exit(0);
}

