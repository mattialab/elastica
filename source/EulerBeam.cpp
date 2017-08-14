/*
 * BendingWavesCouple.cpp
 *
 *  Created on: Jun 21, 2015
 *      Author: mgazzola
 */

#include "EulerBeam.h"

EulerBeam::EulerBeam(const int argc, const char ** argv)
{
}

EulerBeam::~EulerBeam()
{
}

bool EulerBeam::_cantileverPointLoad( const string outfileName )
{
	// Input parameters
	const int n = 100;
	const REAL timeSimulation = 200.0; // you need a long time and a small nu to fully converge to the exact steady state solution
	const REAL L0 = 1.0;
	const REAL r0 = 0.02*L0;
	const REAL dL0 = L0/(REAL)n;
	const REAL dt = 1e-3*dL0;
	const REAL P = 0.5;
	const REAL density = 10000;
	const REAL E = 1e7;
	const REAL poissonRatio = 0.25;
	const REAL G = E / (poissonRatio+1.0);
	const REAL nu = 0.1;
	const REAL relaxationNu = 0.0;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 5;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL totalMass = Vol * density;
	const REAL initialTotalTwist = 0.0;
	const Vector3 originRod = Vector3(-L0/2.0,0.0,0.0);
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
	Matrix3 S0 = 1e6*Matrix3(	1.0,	0.0,	0.0,
								0.0,	1.0,	0.0,
								0.0,	0.0,	1.0);

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

	poly.print_s_internalTorques("internalTorque");
	poly.print_s_coordinates("shape");
	poly.print_s_internalShears("internalShear");

	cout <<  "total internal energy = " << poly.getTotalEnergy() << endl;
	cout <<  "total translational energy = " << poly.getTotalTranslationalEnergy() << endl;
	cout <<  "total rotational energy = " << poly.getTotalRotationalEnergy() << endl;

	return false;
}

bool EulerBeam::_cantileverPointTorque( const string outfileName )
{
	// Input parameters
	const int n = 100;
	const REAL timeSimulation = 50.0; // you need a long time and a small nu to fully converge to the exact steady state solution
	const REAL L0 = 1.0;
	const REAL r0 = 0.02*L0;
	const REAL dL0 = L0/(REAL)n;
	const REAL dt = 1e-3*dL0;
	const REAL M = 0.5;
	const REAL density = 10000;
	const REAL E = 1e7;
	const REAL poissonRatio = 0.25;
	const REAL G = E / (poissonRatio+1.0);
	const REAL nu = 0.1;
	const REAL relaxationNu = 0.0;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 5;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL totalMass = Vol * density;
	const REAL initialTotalTwist = 0.0;
	const Vector3 originRod = Vector3(-L0/2.0,0.0,0.0);
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
	Matrix3 S0 = 1e6*Matrix3(	1.0,	0.0,	0.0,
								0.0,	1.0,	0.0,
								0.0,	0.0,	1.0);

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
	GradualEndpointTorques endpointsTorque = GradualEndpointTorques(Vector3(), Vector3(0,0,-M), timeSimulation/2);
	MultipleForces multipleForces;
	multipleForces.add(&endpointsTorque);
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

	poly.print_s_internalTorques("internalTorque");
	poly.print_s_coordinates("shape");
	poly.print_s_internalShears("internalShear");
	poly.print_s_curvatures("curvature");

	cout <<  "total internal energy = " << poly.getTotalEnergy() << endl;
	cout <<  "total translational energy = " << poly.getTotalTranslationalEnergy() << endl;
	cout <<  "total rotational energy = " << poly.getTotalRotationalEnergy() << endl;

	return false;
}

bool EulerBeam::_cantileverTorqueAtFreeEnds( const string outfileName )
{
	// Input parameters
	const int n = 100;
	const REAL timeSimulation = 20.0; // you need a long time and a small nu to fully converge to the exact steady state solution
	const REAL L0 = 1.0;
	const REAL r0 = 0.02*L0;
	const REAL dL0 = L0/(REAL)n;
	const REAL dt = 1e-3*dL0;
	const REAL M = 0.5;
	const REAL density = 10000;
	const REAL E = 1e7;
	const REAL poissonRatio = 0.25;
	const REAL G = E / (poissonRatio+1.0);
	const REAL nu = 0.1;
	const REAL relaxationNu = 0.0;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 5;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL totalMass = Vol * density;
	const REAL initialTotalTwist = 0.0;
	const Vector3 originRod = Vector3(-L0/2.0,0.0,0.0);
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
	Matrix3 S0 = 1e6*Matrix3(	1.0,	0.0,	0.0,
								0.0,	1.0,	0.0,
								0.0,	0.0,	1.0);

	// Initialize straight rod and pack it into a vector of pointers to rod --> Use linear load-strain (hence the true flag at the end)!!!
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	FreeBC endBC = FreeBC();
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	GradualEndpointTorques endpointsTorque = GradualEndpointTorques(Vector3(0,0,M), Vector3(0,0,-M), timeSimulation/2);
	MultipleForces multipleForces;
	multipleForces.add(&endpointsTorque);
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

	poly.print_s_internalTorques("internalTorque");
	poly.print_s_coordinates("shape");
	poly.print_s_internalShears("internalShear");
	poly.print_s_curvatures("curvature");

	cout <<  "total internal energy = " << poly.getTotalEnergy() << endl;
	cout <<  "total translational energy = " << poly.getTotalTranslationalEnergy() << endl;
	cout <<  "total rotational energy = " << poly.getTotalRotationalEnergy() << endl;

	return false;
}

bool EulerBeam::_cantileverConstantMuscularActivity( const string outfileName )
{
	// Input parameters
	const int n = 100;
	const REAL timeSimulation = 50; // you need a long time and a small nu to fully converge to the exact steady state solution
	const REAL L0 = 1.0;
	const REAL r0 = 0.02*L0;
	const REAL dL0 = L0/(REAL)n;
	const REAL dt = 1e-3*dL0;
	const REAL density = 10000;
	const REAL E = 1e7;
	const REAL poissonRatio = 0.25;
	const REAL G = E / (poissonRatio+1.0);
	const REAL nu = 0.1;
	const REAL relaxationNu = 0.0;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 5;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL A0 = M_PI*r0*r0;
	const REAL Vol = A0 * L0;
	const REAL totalMass = Vol * density;
	const REAL initialTotalTwist = 0.0;
	const Vector3 originRod = Vector3(-L0/2.0,0.0,0.0);
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
	Matrix3 S0 = 1e6*Matrix3(	1.0,	0.0,	0.0,
								0.0,	1.0,	0.0,
								0.0,	0.0,	1.0);

	const REAL M = (2.0*M_PI/L0)*E*I0_1;

	// Initialize straight rod and pack it into a vector of pointers to rod --> Use linear load-strain (hence the true flag at the end)!!!
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::straightRod(n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod, directionRod, normalRod, nu, relaxationNu, useSelfContact);
	vector<Rod*> rodPtrs;
	rodPtrs.push_back(&rod);
	rod.update();
	rod.computeEnergies();

	// Pack boundary conditions
	FreeBC endBC = FreeBC();
	vector<RodBC*> boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&endBC);

	// Pack all forces together (no forces applied)
	ConstantMuscleTorques constTorque = ConstantMuscleTorques(-M, Vector3(1.0,0.0,0.0), timeSimulation/2.0);
	MultipleForces multipleForces;
	multipleForces.add(&constTorque);
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

	poly.print_s_internalTorques("internalTorque");
	poly.print_s_coordinates("shape");
	poly.print_s_internalShears("internalShear");
	poly.print_s_curvatures("curvature");

	cout <<  "total internal energy = " << poly.getTotalEnergy() << endl;
	cout <<  "total translational energy = " << poly.getTotalTranslationalEnergy() << endl;
	cout <<  "total rotational energy = " << poly.getTotalRotationalEnergy() << endl;

	return false;
}

void EulerBeam::run()
{
	//_cantileverPointLoad("longwaves_0005");
	//_cantileverPointTorque("longwaves_0005");
	//_cantileverTorqueAtFreeEnds("longwaves_0005");
	_cantileverConstantMuscularActivity("longwaves_0005");

	exit(0);
}

