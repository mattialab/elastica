#include "Snake.h"

Snake::Snake(const int argc, const char ** argv) :
	amp(0.0), w(0.0), v(0.0), parser(argc, argv)
{
	const string ctrl = parser("-cmaes").asString();
	ncycles = parser("-ncycles").asDouble();
	framesPerUnitTime = parser("-framesPerUnitTime").asInt();

	vector<double> paramters = interfaceCma.parse(ctrl);
	interfaceCma.printParameters();

	amp.clear();
	amp.push_back( paramters[0] );
	amp.push_back( paramters[1] );
	amp.push_back( paramters[2] );
	amp.push_back( paramters[3] );
	v = paramters[4];
}

REAL Snake::_snakeRun()
{
	// Driving parameters
	const int n = 100;
	const REAL density = 1000.0;
	const REAL L0 = 1.0;
	const REAL r0 = 0.025*L0;
	const REAL totalMass = density*M_PI*r0*r0*L0;
	const REAL E = 1e7;
	const REAL g = 9.81;
	const REAL Froude = 0.1;
	const REAL T = 1.0;
	const REAL mu = L0/(T*T*g*Froude);
	const REAL dt = 1e-5*T;
	const REAL timeSimulation = (11.0+0.01)*T;
	w = 2.0*M_PI/T;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = (unsigned int)ceil(120/T);
	const unsigned int povrayPerUnitTime = (unsigned int)ceil(60/T);

	const REAL dL0 = L0/(double)n;							// length of cross-section element
	const REAL A0 = M_PI*r0*r0;
	const REAL poissonRatio = 0.5;							// Incompressible
	const REAL G = E / (poissonRatio+1.0);

	// Define plane
	const REAL angle = 0.0;
	const Vector3 originPlane = Vector3(0.0,0.0,0.0);
	const Vector3 normalPlane = Vector3(0.0,sin(angle),cos(angle)).unitize();

	// Define rod
	const Vector3 directionRod = Vector3(1.0,0.0,0.0);
	const Vector3 normalRod = Vector3(0.0,0.0,1.0);
	const Vector3 originRod = (originPlane - L0/2.0*directionRod) + r0*normalPlane;

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
	Matrix3 S0 = Matrix3(	(4.0/3.0)*G*A0,	0.0,			0.0,
							0.0,			(4.0/3.0)*G*A0,	0.0,
							0.0,			0.0,			E*A0);

	// Initialize straight rod and pack it into a vector of pointers to rod
	const REAL initialTotalTwist = 0.0;
	const REAL nu = 5;
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
	SplineMuscleTorques muscleLateral  = SplineMuscleTorques(amp, w, v, normalPlane, T);
	GravityForce gravity = GravityForce(Vector3(0.0,0.0,-g));
	MultipleForces multipleForces;
	multipleForces.add(&muscleLateral);
	multipleForces.add(&gravity);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Set up substrate properties and snake-plane interaction object
	const REAL kPlane = 1.0; // stiffness of ground
	const REAL nuPlane = 1e-6; // viscous damping of ground
	const REAL muKineticForward = mu;
	const REAL muKineticSideways = 2.0*muKineticForward;
	const REAL muKineticBackward = 1.5*muKineticForward;
	const REAL muStaticForward = 2.0*mu;
	const REAL muStaticSideways = 2.0*muStaticForward;
	const REAL muStaticBackward = 1.5*muStaticForward;
	const REAL vStatic = 1e-8;
	AnisotropicFrictionPlaneInteraction frictionPlane = AnisotropicFrictionPlaneInteraction(	rodPtrs, normalPlane, originPlane, kPlane, nuPlane,
																								muKineticForward, muKineticBackward, muKineticSideways,
																								muStaticForward, muStaticBackward, muStaticSideways, vStatic);
	vector<Interaction*> substrateInteractionsPtrs;
	substrateInteractionsPtrs.push_back(&frictionPlane);

	// Set up integrator (define integration order)
	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs);

	// Instantiate simulator
	Polymer poly = Polymer(integrator);

	// I am goint go collect data over this time window
	poly.setWindowStats(4.0*T, 5.0*T);

	// Run simulation
	string outfileName = string("prova");
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

	const Vector3 avgVel= poly.getAverageVelocity();
	const REAL fwdAvgVel = avgVel % directionRod;
	//const REAL latAvgVel = fabs(avgVel % Vector3(0.0,1.0,0.0));
	const REAL fitness = fwdAvgVel;

	return (fitness);
}

void Snake::run()
{
	const REAL fitness = _snakeRun();
	cout << "Mean velocity: " << fitness << endl;
	interfaceCma.dumpFitness( fitness );
	exit(0);
}

