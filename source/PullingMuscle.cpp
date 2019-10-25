#include "PullingMuscle.h"

/*
This example demonstrates how we connect two Cosserat rods together using an
active muscle along their centers like so:

// clang-format off

     Hinge joint	Hinge joint
         (x)           (x)
         ||            ||
         ||   Muscle   ||
         ||  ++++++++  ||
     Rod1||++++++++++++|| Rod2
         ||  ++++++++  ||
         ||            ||
         ||            ||


   ^ z
   |
   |
   o — — —> y
  /
 /
x

// clang-format on

No other forces (like gravity) are considered.

The rods are attached to the ground via a hinge joint, which only allows
relative orientation changes in the x-z plane (no displacement of the rod
connected point is allowed). The muscle connecting them is initially in a
relaxed configuration. The muscle is then activated and so it starts
contracting, exerting forces on the rod COM. This muscle force causes a torque
about the hinge joint, which causes the rods to rotate in the y-z plane. After
reaching maximal contraction, the muscle starts relaxing and reaches its
original configuration, as a result of which the rods are returned to their
initial configuration as well. This contraction-relaxation cycles is repeated
sinusoidally. This example is instructive in showing how muscles can be
implemented using our Cosserat rod model.


*/
PullingMuscle::PullingMuscle(const int argc, const char **argv)
    : amp(0.0), w(0.0), v(0.0) {}

// Units in this case are mm/g/s

vector<REAL> PullingMuscle::_pullingmuscleRun() {
  vector<Rod *> rodPtrs;

  // Dumping frequencies (number of frames/dumps per unit time)
  const REAL diagPerUnitTime = 120;
  const REAL povrayPerUnitTime = 50;
  const REAL dt = 1e-6;
  const REAL timeSimulation = (10.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Links
  // Driving parameters
  // Link One shape
  const int n = 10;
  const REAL density = 1.75e-3;  // 1.75g/cm^3
  const REAL L0 = 200.0;
  const REAL r0 = 7.0;
  const REAL E = 3e7;
  const REAL totalMass = density * M_PI * r0 * r0 * L0;

  const REAL dL0 = L0 / (double)n;  // length of cross-section element
  const REAL poissonRatio = 0.5;    // Incompressible
  const REAL G = E / (poissonRatio + 1.0);

  // Define rod
  const Vector3 Linkdirection = Vector3(0.0, 0.0, -1.0);
  const Vector3 Linknormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 Linkoneorigin = Vector3(0.0, L0, 1.5 * L0);
  const Vector3 Linktwoorigin = Linkoneorigin + L0 * Vector3(0.0, -1.0, 0.0);

  // Second moment of area for disk cross section
  const REAL A0 = M_PI * r0 * r0;
  const REAL I0_1 = A0 * A0 / (4.0 * M_PI);
  const REAL I0_2 = I0_1;
  const REAL I0_3 = 2.0 * I0_1;
  const Matrix3 I0 = Matrix3(I0_1, 0.0, 0.0, 0.0, I0_2, 0.0, 0.0, 0.0, I0_3);
  // Mass inertia matrix for disk cross section
  const Matrix3 J0 = density * dL0 * I0;

  // Bending matrix
  Matrix3 B0 =
      Matrix3(E * I0_1, 0.0, 0.0, 0.0, E * I0_2, 0.0, 0.0, 0.0, G * I0_3);
  // Shear matrix
  Matrix3 S0 = Matrix3((4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0,
                       (4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0, E * A0);

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist = 0.0;
  const REAL nu = 0.1;
  const REAL relaxationNu = 0.0;
  const bool useSelfContact = false;

  Rod *rod1 = RodInitialConfigurations::straightRod(
      n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, Linkoneorigin,
      Linkdirection, Linknormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod1);
  rod1->update(0.0);
  Rod *rod2 = RodInitialConfigurations::straightRod(
      n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, Linktwoorigin,
      Linkdirection, Linknormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod2);
  rod2->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Soft Contracting Element
  // two ends are stiff --- representing tendons
  // middle elements are soft --- representing muscle
  // Driving parameters
  const int nm = 10;
  const REAL densitym = 1.06e-3;
  const REAL Em = 1e4;
  const REAL Et = 5e8;

  const REAL poissonRatiom = 0.5;  // Incompressible
  const REAL Gm = Em / (poissonRatiom + 1.0);
  const REAL Gt = Et / (poissonRatiom + 1.0);

  Vector3 originmuscle = rodPtrs[0]->x[6];
  Vector3 directionmuscle = (rodPtrs[1]->x[6] - rodPtrs[0]->x[6]).unitize();
  Vector3 normalmuscle = Vector3(1.0, 0.0, 0.0);

  REAL Lm = (rodPtrs[1]->x[6] - rodPtrs[0]->x[6]).length();
  REAL dLm = Lm / (double)nm;

  // geometrical properties of tendons
  // Second moment of area for disk cross section for tendons
  const REAL rm_t = 6.0;
  const REAL Am_t = M_PI * rm_t * rm_t;
  const REAL I0_1m_t = Am_t * Am_t / (4.0 * M_PI);
  const REAL I0_2m_t = I0_1m_t;
  const REAL I0_3m_t = 2.0 * I0_1m_t;
  const Matrix3 I0m_t =
      Matrix3(I0_1m_t, 0.0, 0.0, 0.0, I0_2m_t, 0.0, 0.0, 0.0, I0_3m_t);

  // Mass inertia matrix for disk cross section
  const Matrix3 Jm_t = densitym * dLm * I0m_t;

  // Bending matrix
  const Matrix3 Bm_t = Matrix3(Et * I0_1m_t, 0.0, 0.0, 0.0, Et * I0_2m_t, 0.0,
                               0.0, 0.0, Gt * I0_3m_t);

  // Shear matrix
  const Matrix3 Sm_t =
      Matrix3((4.0 / 3.0) * Gt * Am_t, 0.0, 0.0, 0.0, (4.0 / 3.0) * Gt * Am_t,
              0.0, 0.0, 0.0, Et * Am_t);

  // geometrical properties of muscle
  // Second moment of area for disk cross section for muscle
  const REAL rm = 12.0;
  const REAL Am = M_PI * rm * rm;
  const REAL I0_1m = Am * Am / (4.0 * M_PI);
  const REAL I0_2m = I0_1m;
  const REAL I0_3m = 2.0 * I0_1m;
  const Matrix3 I0m =
      Matrix3(I0_1m, 0.0, 0.0, 0.0, I0_2m, 0.0, 0.0, 0.0, I0_3m);
  // Mass inertia matrix for disk cross section
  const Matrix3 Jm = densitym * dLm * I0m;

  // Bending matrix
  const Matrix3 Bm =
      Matrix3(Em * I0_1m, 0.0, 0.0, 0.0, Em * I0_2m, 0.0, 0.0, 0.0, Gm * I0_3m);
  // Shear matrix
  const Matrix3 Sm = Matrix3((4.0 / 3.0) * Gm * Am, 0.0, 0.0, 0.0,
                             (4.0 / 3.0) * Gm * Am, 0.0, 0.0, 0.0, Em * Am);

  vector<REAL> rm_v = vector<REAL>(nm);  //(mm)
  vector<Matrix3> Jm_v = vector<Matrix3>(nm);
  vector<Matrix3> Bm_v = vector<Matrix3>(nm - 1);
  vector<Matrix3> Sm_v = vector<Matrix3>(nm);

  // tendon
  for (unsigned int i = 0; i < 3; i++) {
    rm_v[i] = rm_t;
    Jm_v[i] = Jm_t;
    Bm_v[i] = Bm_t;
    Sm_v[i] = Sm_t;
  }
  Bm_v[2] = (Bm_t + Bm) / 2;

  for (unsigned int i = 3; i < 7; i++) {
    rm_v[i] = rm;
    Jm_v[i] = Jm;
    Bm_v[i] = Bm;
    Sm_v[i] = Sm;
  }
  Bm_v[6] = (Bm_t + Bm) / 2;
  // muscle SP-SP
  for (unsigned int i = 7; i < 10; i++) {
    rm_v[i] = rm_t;
    Jm_v[i] = Jm_t;
    if (i < 9) Bm_v[i] = Bm_t;
    Sm_v[i] = Sm_t;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist_m = 0.0;
  const REAL nu_m = 10.0;
  const REAL relaxationNu_m = 0.0;
  const bool useSelfContact_m = false;

  Rod *rod_m = RodInitialConfigurations::straightRod_v(
      nm, densitym, rm_v, Jm_v, Bm_v, Sm_v, Lm, initialTotalTwist_m,
      originmuscle, directionmuscle, normalmuscle, nu_m, relaxationNu_m,
      useSelfContact_m);
  rodPtrs.push_back(rod_m);
  rod_m->update(0.0);
  //-----------------------------------------------------------------------------------------------------------------
  // Pack boundary conditions
  vector<RodBC *> boundaryConditionsPtrs;
  // Link One -- first node pinned
  HingeBC hinge1 = HingeBC(rodPtrs[0]);
  boundaryConditionsPtrs.push_back(&hinge1);
  // Link Two -- first node pinned
  HingeBC hinge2 = HingeBC(rodPtrs[1]);
  boundaryConditionsPtrs.push_back(&hinge2);
  // Contracting element
  FreeBC freeBC = FreeBC();
  boundaryConditionsPtrs.push_back(&freeBC);

  // Pack all forces together (no forces applied)
  vector<ExternalForces *> externalForcesPtrs;

  // Gravity
  MultipleForces multipleForces1;
  GravityForce gravity = GravityForce(Vector3(0.0, 0.0, 0.0));
  multipleForces1.add(&gravity);
  MultipleForces *multipleForcesPtr1 = multipleForces1.get();
  for (unsigned int i = 0; i < 3; i++) {
    externalForcesPtrs.push_back(multipleForcesPtr1);
  }

  vector<Interaction *> substrateInteractionsPtrs;

  // Set up External Contact -- This is for the five cases in the paper, not
  // used in this case
  vector<pair<int, int>> attachpoint;
  vector<ExternalContact *> externalcontactPtrs;
  /* The second and third argument are unimportant, but
         are preserved here for legacy purposes. Hence we simply
         set it to 0.0
  */
  ExternalContact externalcontact =
      ExternalContact(rodPtrs, 0.0, 0.0, attachpoint);
  externalcontactPtrs.push_back(&externalcontact);

  // Set up Simple Connection
  vector<SimpleConnection *> simpleconnectionPtrs;
  SimpleConnection simpleconnection = SimpleConnection(rodPtrs);
  simpleconnectionPtrs.push_back(&simpleconnection);
  //-----------------------------------------------------------------------------------------------------------------
  // Set up integrator (define integration order)
  PolymerIntegrator *integrator = new PositionVerlet2nd(
      rodPtrs, externalForcesPtrs, boundaryConditionsPtrs,
      substrateInteractionsPtrs, externalcontactPtrs, simpleconnectionPtrs);

  // Instantiate simulator
  Polymer poly = Polymer(integrator);

  // I am goint go collect data over this time window
  poly.setWindowStats(1.0, 2.0);

  // Run simulation
  string outfileName = string("prova");
  const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime,
                                     povrayPerUnitTime, outfileName);

  // Throw exception if something went wrong
  if (!goodRun)
    throw "not good run in localized helical buckling, what is going on?";

  const vector<Vector3> avgVel = poly.getAverageVelocity();

  vector<REAL> fwdAvgVel;
  for (unsigned int i = 0; i < 1; i++) {
    fwdAvgVel.push_back(avgVel[i] % Linkdirection);
  }
  const vector<REAL> fitness = fwdAvgVel;

  return (fitness);
}

void PullingMuscle::run() {
  const vector<REAL> fitness = _pullingmuscleRun();
  exit(0);
}
