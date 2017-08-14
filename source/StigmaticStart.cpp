/*
 * SelfAssembly.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: mgazzola
 */

#include "StigmaticStart.h"

StigmaticStart::StigmaticStart(const int argc, const char ** argv) : parser(argc, argv)
{
	const string ctrl = parser("-cmaes").asString();

	paramtersCMA = interfaceCma.parse(ctrl);
	interfaceCma.printParameters();
}

REAL StigmaticStart::_snakeRun(const REAL _W, const REAL _M, const REAL _Froude, const REAL _L, const REAL _T, const REAL _d, const REAL _E, const REAL _rho)
{
	// Driving parameters
	const int n = 100;
	const REAL W = _W;
	const REAL M = _M;
	const REAL Froude = _Froude;
	const REAL L0 = _L;
	const REAL T = _T;
	const REAL d0 = _d;
	const REAL E = _E;
	const REAL density = _rho;

	// Derived quantities
	const REAL r0 = 0.5*d0;
	const REAL A0 = M_PI*r0*r0;
	const REAL totalMass = L0*A0*density;
	const REAL characteristicSpeed = L0/T;

	const REAL dt = 1e-6*T;
	const REAL timeSimulation = T;

	const REAL g = W*E*d0*d0/(density*L0*L0*L0);
	const REAL mu = L0/(T*T*g*Froude);
	const REAL Gamma = M*d0*E/L0;
	const REAL torque = Gamma*d0*d0*d0;

	cout << "g = " << g << endl;
	cout << "mu = " << mu << endl;
	cout << "Gamma = " << Gamma << endl;
	cout << "torque = " << torque << endl;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = (unsigned int)ceil(0/T);
	const unsigned int povrayPerUnitTime = (unsigned int)ceil(0/T);

	const REAL poissonRatio = 0.5;							// Incompressible
	const REAL G = E / (poissonRatio+1.0);

	// Define plane
	const Vector3 originPlane = Vector3(0.0,0.0,0.0);
	const Vector3 normalPlane = Vector3(0.0,0.0,1.0);

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
	const REAL dL0 = L0/(REAL)n;
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
	const REAL nu = 1.0;
	const REAL relaxationNu = 0.0;
	const bool useSelfContact = true;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	FreeBC freeBC = FreeBC();
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&freeBC);

	// Set gravity
	GravityForce gravity = GravityForce(Vector3(0.0,0.0,-g));

	// Set bending localized traveling muscle activity
	const REAL amp_firstBend = -torque;
	const REAL contractSpeed_firstBend = characteristicSpeed;
	const REAL widthFactor_firstBend = L0/10.0;
	const REAL position_firstBend = 0.1*L0;
	const REAL rumpUpTime_firstBend = 0.1*T;
	LocalizedMuscleForce muscle_firstBend = LocalizedMuscleForce(amp_firstBend, widthFactor_firstBend, contractSpeed_firstBend, position_firstBend, normalRod, rumpUpTime_firstBend);

	const REAL amp_secondBend = -amp_firstBend;
	const REAL contractSpeed_secondBend = contractSpeed_firstBend;
	const REAL widthFactor_secondBend = widthFactor_firstBend;
	const REAL position_secondBend = 0.3*L0;
	const REAL rumpUpTime_secondBend = rumpUpTime_firstBend;
	LocalizedMuscleForce muscle_secondBend = LocalizedMuscleForce(amp_secondBend, widthFactor_secondBend, contractSpeed_secondBend, position_secondBend, normalPlane, rumpUpTime_secondBend);

	const REAL amp_middle = 0.1*fabs(amp_firstBend);
	const REAL contractSpeed_middle = contractSpeed_firstBend;
	const REAL widthFactor_middle = widthFactor_firstBend/2.0;
	const REAL position_middle = 0.2*L0;
	const REAL rumpUpTime_middle = rumpUpTime_firstBend;
	const Vector3 d2 = Vector3(0.0,1.0,0.0);
	LocalizedMuscleForceLagrangian_ForcedToBeInPlane muscle_middle = LocalizedMuscleForceLagrangian_ForcedToBeInPlane(amp_middle, widthFactor_middle, contractSpeed_middle, position_middle, d2, normalPlane, rumpUpTime_middle);

	// Pack all forces together
	MultipleForces multipleForces;
	multipleForces.add(&muscle_firstBend);
	multipleForces.add(&muscle_secondBend);
	multipleForces.add(&muscle_middle);
	multipleForces.add(&gravity);
	vector<ExternalForces*> externalForcesPtrs = multipleForces.get();

	// Set up substrate properties and snake-plane interaction object
	const REAL kPlane = 1e5; // stiffness of ground
	const REAL nuPlane = 1.0; // viscous damping of ground
	const REAL muKineticForward = mu/2.0; // That is because the mu in the Froude number is relative to the side friction coeafficient
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
	poly.setWindowStats(0.1*T, 0.7*T);

	// Run simulation
	string outfileName = string("prova");
	const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime, povrayPerUnitTime, outfileName);

	// Throw exception if something went wrong
	if (!goodRun)
		throw "not good run in localized helical buckling, what is going on?";

	const Vector3 avgVel= poly.getAverageVelocity();
	const REAL fitness = avgVel.length();

	return (fitness);
}

void StigmaticStart::run()
{
	/*
	{
		// Individual testcase
		const REAL L = 1.0; // Length [m]
		const REAL T = 1.0; // Characterisitc time [s]
		const REAL d = 0.03*L; // Diameter [m]
		const REAL rho = 1000.0; // Density [kg/m^3]
		const REAL E = 1e7; // Young modulus [Pa] --> Bib tendons+muscles: Shinohara:2010, Ogneva:2010, Kot:2012 --> 5e6
		const REAL g = 9.81;
		const REAL Gamma = 2.5e5;

		const REAL W = rho*g*L*L*L/(E*d*d);
		const REAL M = Gamma*L/(d*E);
		const REAL Froude = 0.1; // Froude number Fr = L/(T^2*g*mu);

		const REAL fitness = _snakeRun(W, M, Froude, L, T, d, E, rho);
		cout << "Mean velocity: " << fitness << endl;
	}
	*/

	{
		// CMA simulation
		const REAL L = 1.0; // Length [m]
		const REAL T = 1.0; // Characterisitc time [s]
		const REAL d = 0.03*L; // Diameter [m]
		const REAL rho = 1000.0; // Density [kg/m^3]
		const REAL E = 1e7; // Young modulus [Pa] --> Bib tendons+muscles: Shinohara:2010, Ogneva:2010, Kot:2012 --> 5e6
		const REAL g = 9.81;
		const REAL Gamma = paramtersCMA[0];

		const REAL W = rho*g*L*L*L/(E*d*d);
		const REAL M = Gamma*L/(d*E);
		const REAL Froude = 0.1; // Froude number Fr = L/(T^2*g*mu);

		const REAL fitness = _snakeRun(W, M, Froude, L, T, d, E, rho);
		cout << "Mean velocity: " << fitness << endl;

		interfaceCma.dumpFitness( -fitness );
	}

	exit(0);
}
