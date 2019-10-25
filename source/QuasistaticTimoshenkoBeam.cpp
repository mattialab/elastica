
#include "QuasistaticTimoshenkoBeam.h"

/*
This code simulates a cantilevered slender beam (clamped at the wall at one end,
free at another) bending under the influence of a downward point force like so:

// clang-format off

Wall
/|							      | Force
/|							      |
/|							      v
/|================================== Rod/Beam
/|
/|
/|

// clang-format on

        Any other effects like gravity are ignored. In our case, where the shear
modulus (G) is comparable to the Young's modulus (E), the static deflection of
the beam centerline is given by the Timoshenko bending solution (see Gazzola et
al, 2018).

        In the code below, we consider a Cosserat rod and clamp it at the end
using FixedBC (applies forces and torques). We time-step the code until it
reaches a steady-state (or in this case a quasi-steady state). A comparison of
the obtained solution with the theoretical results show a good match. For more
details the reader is referred to "Forward and inverse problems in the mechanics
of soft filaments" by Gazzola et.al, RSoS, 2018. This example  is instructive in
showing how a rod, along with simple boundary conditions can be coded up in the
framework of Elastica.

*/

QuasistaticTimoshenkoBeam::QuasistaticTimoshenkoBeam(const int argc,
                                                     const char **argv) {}

QuasistaticTimoshenkoBeam::~QuasistaticTimoshenkoBeam() {}

bool QuasistaticTimoshenkoBeam::_test(const int nEdges, const REAL _dt,
                                      const REAL _L, const REAL _r,
                                      const REAL _P, const REAL _timeSimulation,
                                      const REAL _E, const REAL _G,
                                      const REAL _rho, const REAL _nu,
                                      const REAL _relaxationNu,
                                      const string outfileName) {
  // Input parameters
  const int n = nEdges;  // number of discretization edges (i.e. n+1 points)
                         // along the entire rod
  const REAL timeSimulation = _timeSimulation;  // total simulation time
  const REAL dt = _dt;                          // time step
  const REAL P = _P;
  const REAL L0 = _L;         // total length of rod [m]
  const REAL r0 = _r;         // radius [m]
  const REAL density = _rho;  // [kg/m^3]
  const REAL E = _E;          // GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
  const REAL G = _G;          // GPa --> rubber~0.01-0.1Gpa, iron~200Gpa
  const REAL nu = _nu;        // Numerical damping viscosity [m^2/s]
  const REAL relaxationNu =
      _relaxationNu;  // relaxation time for exponential decay of nu

  // Dumping frequencies (number of frames/dumps per unit time)
  const unsigned int diagPerUnitTime = 5;
  const unsigned int povrayPerUnitTime = 0;

  // Physical parameters
  const REAL dL0 = L0 / (double)n;  // length of cross-section element
  const REAL A0 = M_PI * r0 * r0;
  const REAL Vol = A0 * L0;
  const REAL totalMass = Vol * density;
  const REAL initialTotalTwist = 0.0;
  const Vector3 originRod = Vector3(0.0, 0.0, 0.0);
  const Vector3 directionRod = Vector3(1.0, 0.0, 0.0);
  const Vector3 normalRod = Vector3(0.0, 0.0, 1.0);

  // Second moment of area for disk cross section
  const REAL I0_1 = A0 * A0 / (4.0 * M_PI);
  const REAL I0_2 = I0_1;
  const REAL I0_3 = 2.0 * I0_1;
  const Matrix3 I0 = Matrix3(I0_1, 0.0, 0.0, 0.0, I0_2, 0.0, 0.0, 0.0, I0_3);

  // Mass inertia matrix for disk cross section
  const Matrix3 J0 = density * dL0 * I0;

  // Bending matrix (TOD: change this is wrong!!)
  Matrix3 B0 =
      Matrix3(E * I0_1, 0.0, 0.0, 0.0, E * I0_2, 0.0, 0.0, 0.0, G * I0_3);

  // Shear matrix
  Matrix3 S0 = Matrix3((4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0,
                       (4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0, E * A0);

  // Initialize straight rod and pack it into a vector of pointers to rod -->
  // Use linear load-strain (hence the true flag at the end)!!!
  const bool useSelfContact = false;
  Rod *rod = RodInitialConfigurations::straightRod(
      n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod,
      directionRod, normalRod, nu, relaxationNu, useSelfContact);
  vector<Rod *> rodPtrs;
  rodPtrs.push_back(rod);
  rod->update(0.0);
  rod->computeEnergies();

  // Pack boundary conditions
  TimoshenkoBeamBC endBC = TimoshenkoBeamBC(rodPtrs);
  vector<RodBC *> boundaryConditionsPtrs;
  boundaryConditionsPtrs.push_back(&endBC);

  // Pack all forces together
  vector<ExternalForces *> externalForcesPtrs;
  GradualEndpointForces endpointsForce =
      GradualEndpointForces(Vector3(), Vector3(0, -P, 0), timeSimulation / 2);
  MultipleForces multipleForces;
  multipleForces.add(&endpointsForce);
  MultipleForces *multipleForcesPtr = multipleForces.get();
  externalForcesPtrs.push_back(multipleForcesPtr);

  // Empty interaction forces (no substrate in this case)
  vector<Interaction *> substrateInteractionsPtrs;

  // No external contact function needed
  vector<pair<int, int>> attachpoint;
  vector<ExternalContact *> externalcontactPtrs;
  ExternalContact externalcontact =
      ExternalContact(rodPtrs, 0.0, 0.0, attachpoint);
  externalcontactPtrs.push_back(&externalcontact);

  // No Simple Connection needed
  vector<SimpleConnection *> simpleconnectionPtrs;
  SimpleConnection simpleconnection = SimpleConnection(rodPtrs);
  simpleconnectionPtrs.push_back(&simpleconnection);
  //-----------------------------------------------------------------------------------------------------------------
  // Set up time integrator
  PolymerIntegrator *integrator = new PositionVerlet2nd(
      rodPtrs, externalForcesPtrs, boundaryConditionsPtrs,
      substrateInteractionsPtrs, externalcontactPtrs, simpleconnectionPtrs);

  // Simulate
  Polymer poly = Polymer(integrator);
  const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime,
                                     povrayPerUnitTime, outfileName);

  // Throw exception if something went wrong
  if (!goodRun)
    throw "not good run in localized helical buckling, what is going on?";

  // Dump post buckling helical shape
  {
    string rodShape = outfileName + "_shape.txt";
    FILE *fitnessFile = fopen(rodShape.c_str(), "w");
    assert(fitnessFile != NULL);

    for (unsigned int i = 0; i < rod->x.size(); i++) {
      const REAL s_hat = rod->x[i].x;
      /* Timoshenko theory analytical solution for shear hardened rod
       (set in _longWaveTest below)
      */
      const REAL yposition = -P / (4.0 / 3.0 * A0 * G) * s_hat -
                             P * L0 / (2 * E * I0_1) * s_hat * s_hat +
                             P / (6 * E * I0_1) * s_hat * s_hat * s_hat;
      /* Euler-Bernoulli analytical solution for shear hardened rod
       (set in _longWaveTest below)
      */
      // const REAL yposition = -P * L0 / (2 * E * I0_1) * s_hat * s_hat +
      //                        P / (6 * E * I0_1) * s_hat * s_hat * s_hat;
      fprintf(fitnessFile, "%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n",
              rod->x[i].x, rod->x[i].y, rod->x[i].z, rod->x[i].x, yposition,
              rod->x[i].z);
    }

    fclose(fitnessFile);
  }

  cout << "total internal energy = " << poly.getTotalEnergy() << endl;
  cout << "total translational energy = " << poly.getTotalTranslationalEnergy()
       << endl;
  cout << "total rotational energy = " << poly.getTotalRotationalEnergy()
       << endl;

  return false;
}

void QuasistaticTimoshenkoBeam::_longWaveTest(const int nEdges,
                                              const string outputdata) {
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
  _test(nEdges, dt, L, r, P, simTime, E, G, rho, nu, relaxationNu, minNu,
  outputdata);
  */

  const REAL rho = 5000;
  const REAL L = 3.0;
  const REAL E = 1e6;

  /*
     Here we make a choice for the variable G, the shear modulus.

     If G is set low such that (E \cdot I)/\alpha \cdot L^2 \cdot A \cdot G not
     << 1, then the deflection of the beam under gravity is governed by the
     Timoshenko theorem of beams, which include the effect of shear.

     However, if we harden G or equivalently if (E \cdot I)/\alpha \cdot L^2
     \cdot A \cdot G << 1, then shear is unimportant in the deflection. Then we
     can calculate the final position of the beam using the (simplified) Euler-
     Bernoulli beam theory.

     Two values of G are specified by default. The choice governs which model
     is used to compare with Elastica's rod solution.
  */

  const REAL G = 1e4;  // Calculate using Timoshenko beam theory
  // const REAL G = 1e7; // Calculate using Euler-Bernoulli theory
  const REAL nu = 1e-1;
  const REAL P = 15;
  const REAL dL = L / nEdges;
  const REAL r = 0.25;
  const REAL dt = 0.01 * dL;
  const REAL simTime = 5000;
  const REAL relaxationNu = 0.0;
  _test(nEdges, dt, L, r, P, simTime, E, G, rho, nu, relaxationNu, outputdata);
}

void QuasistaticTimoshenkoBeam::run() {
  // Teja : replaced names here for ease in plotting script
  /*
  _longWaveTest( 5, "timoshenko_final");
  _longWaveTest( 6, "timoshenko_final");
  _longWaveTest( 7, "timoshenko_final");
  _longWaveTest( 8, "timoshenko_final");
  _longWaveTest( 9, "timoshenko_final");
  _longWaveTest( 10, "timoshenko_final");
  _longWaveTest( 20, "timoshenko_final");
  _longWaveTest( 30, "timoshenko_final");
  _longWaveTest( 40, "timoshenko_final");
  _longWaveTest(50, "timoshenko_final");
  */
  /*
  _longWaveTest( 60, "timoshenko_final");
  _longWaveTest( 70, "timoshenko_final");
  _longWaveTest( 80, "timoshenko_final");
  _longWaveTest( 90, "timoshenko_final");
  */
  _longWaveTest(100, "timoshenko_final");
  /*
  _longWaveTest( 200, "timoshenko_final");
  _longWaveTest( 400, "timoshenko_final");
  */
  exit(0);
}
