#include "Flagella.h"

Flagella::Flagella(const int argc, const char **argv)
    : rsmall(0.0), rlarge(0.0), head(0.0), cell1(0.0), cell2(0.0) {}

// Units in this case are mm/g/s

REAL Flagella::_flagellaRun() {
  vector<Rod *> rodPtrs;

  // Dumping frequencies (number of frames/dumps per unit time)
  const REAL diagPerUnitTime = 120;
  const REAL povrayPerUnitTime = 50;

  REAL dt = 0.5e-7;                   // ms
  const REAL timeSimulation = (6.5);  // ms

  //-----------------------------------------------------------------------------------------------------------------
  // Define PDMS Body

  // Driving parameters
  const int n = 18;
  const REAL density = 0.965e-3;
  // const REAL density = 0.965*1e-3;          //(g/mm^-3)
  const REAL L0 = 1.927;  //(mm)
  // const REAL totalMass = density*M_PI*r0*r0*L0;     //(g)
  const REAL E = 3.86e6;          //(g*mm^-1*ms^-2)
  const REAL poissonRatio = 0.5;  // Incompressible
  const REAL G = E / (poissonRatio + 1.0);

  // fluid property
  const REAL T = 1.0 / 1.8;           //(s)
  const REAL fluidDensity = 1.15e-3;  //(g/mm^3)
  const REAL RE = 1.8e-2;
  const REAL dynamicViscosity = 1.2e-3;

  const REAL w = 2.0 * M_PI / T;

  const REAL dL0 = L0 / (double)n;  // length of cross-section element

  // Define rod
  const Vector3 spermdirection = Vector3(1.0, 0.0, 0.0);
  const Vector3 spermnormal = Vector3(0.0, 0.0, 1.0);
  const Vector3 spermorigin = Vector3(0.0, 0.0, 0.1);

  // vector values
  vector<REAL> r0_v2 = vector<REAL>(n);  //(mm)
  vector<Matrix3> J0_v2 = vector<Matrix3>(n);
  vector<Matrix3> B0_v2 = vector<Matrix3>(n - 1);
  vector<Matrix3> S0_v2 = vector<Matrix3>(n);

  // The actual radius for the head and tail
  rsmall = 0.007;  //(mm)
  rlarge = 0.02;   //(mm)

  // This is the virtual radius for matching the bending stiffness
  // of the tail.
  // See
  // Neuromuscular actuation of biohybrid motile bots.
  // by Aydin, O., Zhang, X., Nuethong, S., Pagan-Diaz, G. J., et al.
  // PNAS (2019).for more details.
  REAL requal = 0.0053;
  REAL Asmall = M_PI * requal * requal;
  REAL Alarge = M_PI * rlarge * rlarge;

  // Second moment of area for disk cross section
  REAL I0_1_s = Asmall * Asmall / (4.0 * M_PI);
  REAL I0_2_s = I0_1_s;
  REAL I0_3_s = 2.0 * I0_1_s;
  Matrix3 I0_s = Matrix3(I0_1_s, 0.0, 0.0, 0.0, I0_2_s, 0.0, 0.0, 0.0, I0_3_s);
  // Mass inertia matrix for disk cross section
  Matrix3 J0_s = density * dL0 * I0_s;

  // Bending matrix
  Matrix3 B0_s =
      Matrix3(E * I0_1_s, 0.0, 0.0, 0.0, E * I0_2_s, 0.0, 0.0, 0.0, G * I0_3_s);
  // Shear matrix
  Matrix3 S0_s = Matrix3((4.0 / 3.0) * G * Asmall, 0.0, 0.0, 0.0,
                         (4.0 / 3.0) * G * Asmall, 0.0, 0.0, 0.0, E * Asmall);

  REAL I0_1_l = Alarge * Alarge / (4.0 * M_PI);
  REAL I0_2_l = I0_1_l;
  REAL I0_3_l = 2.0 * I0_1_l;
  Matrix3 I0_l = Matrix3(I0_1_l, 0.0, 0.0, 0.0, I0_2_l, 0.0, 0.0, 0.0, I0_3_l);
  // Mass inertia matrix for disk cross section
  Matrix3 J0_l = density * dL0 * I0_l;

  // Bending matrix
  Matrix3 B0_l =
      Matrix3(E * I0_1_l, 0.0, 0.0, 0.0, E * I0_2_l, 0.0, 0.0, 0.0, G * I0_3_l);
  // Shear matrix
  Matrix3 S0_l = Matrix3((4.0 / 3.0) * G * Alarge, 0.0, 0.0, 0.0,
                         (4.0 / 3.0) * G * Alarge, 0.0, 0.0, 0.0, E * Alarge);

  int headindex = 4;
  // Set properties for the head section
  if (headindex != 0) {
    for (unsigned int i = 0; i < headindex; i++) {
      r0_v2[i] = rlarge;
      J0_v2[i] = J0_l;
      B0_v2[i] = B0_l;
      S0_v2[i] = S0_l;
    }
  }
  // Set properties for the tail section
  for (unsigned int i = headindex; i < n; i++) {
    r0_v2[i] = rsmall;
    J0_v2[i] = J0_s;
    if (i < (n - 1)) {
      B0_v2[i] = B0_s;
    }
    S0_v2[i] = S0_s;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist = 0.0;
  const REAL nu = 0.0;
  const REAL relaxationNu = 0.0;
  const bool useSelfContact = false;

  Rod *rod1 = RodInitialConfigurations::straightRod_v(
      n, density, r0_v2, J0_v2, B0_v2, S0_v2, L0, initialTotalTwist,
      spermorigin, spermdirection, spermnormal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod1);
  rod1->update(0.0);
  rod1->computeEnergies();

  // This enables the computation of hydrodynamic forces on the
  // rod using Slender Body Theory (SBT). This vector of Rod* is
  // later on passed into the SlenderBodyTheoryEnvironment Class
  vector<Rod *> spermPtr;
  spermPtr.push_back(rodPtrs[0]);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Cell

  // number of cells
  NOR = 1;
  vector<pair<int, int>> attachpoint;

  // Driving parameters
  const int nm = 2;
  const REAL densitym = 2.6e-4;
  const REAL Lm = 0.107056;
  const REAL rm = 0.01;
  const REAL totalMassm = densitym * M_PI * rm * rm * Lm;
  const REAL Em = 0.3e5;  // 30Kpa

  const REAL dLm = Lm / (double)nm;  // length of cross-section element
  const REAL Am = M_PI * rm * rm;
  const REAL poissonRatiom = 0.5;  // Incompressible
  const REAL Gm = Em / (poissonRatiom + 1.0);

  vector<Vector3> directionmuscle;
  vector<Vector3> normalmuscle;
  vector<Vector3> originmuscle;

  pair<int, int> attach;

  directionmuscle.push_back(Vector3(1.0, 0.0, 0.0));
  normalmuscle.push_back(Vector3(0.0, 0.0, 1.0));
  originmuscle.push_back(Vector3(4.5 * Lm, -0.0053, 0.1));

  // The attaching index and direction for the cell
  // The first element is the index on the body
  // The second index represents the 'side' at which
  // the cell is attached to the body
  attach.first = 4;
  attach.second = -1;
  attachpoint.push_back(attach);

  // Second moment of area for disk cross section
  const REAL I0_1m = Am * Am / (4.0 * M_PI);
  const REAL I0_2m = I0_1m;
  const REAL I0_3m = 2.0 * I0_1m;
  const Matrix3 I0m =
      Matrix3(I0_1m, 0.0, 0.0, 0.0, I0_2m, 0.0, 0.0, 0.0, I0_3m);

  // Mass inertia matrix for disk cross section
  const Matrix3 Jm = densitym * dLm * I0m;

  // Bending matrix
  Matrix3 Bm =
      Matrix3(Em * I0_1m, 0.0, 0.0, 0.0, Em * I0_2m, 0.0, 0.0, 0.0, Gm * I0_3m);

  // Shear matrix
  Matrix3 Sm = Matrix3((4.0 / 3.0) * Gm * Am, 0.0, 0.0, 0.0,
                       (4.0 / 3.0) * Gm * Am, 0.0, 0.0, 0.0, Em * Am);

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwistm = 0.0;
  const REAL num = 1e-6;
  const REAL relaxationNum = 0.0;
  const bool useSelfContactm = false;

  // Just one rod
  for (unsigned int i = 0; i < NOR; i++) {
    Rod *rod = RodInitialConfigurations::straightRod(
        nm, totalMassm, rm, Jm, Bm, Sm, Lm, initialTotalTwistm, originmuscle[i],
        directionmuscle[i], normalmuscle[i], num, relaxationNum,
        useSelfContactm);
    rodPtrs.push_back(rod);
    rod->update(0.0);
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Pack boundary conditions
  FreeBC freeBC = FreeBC();
  vector<RodBC *> boundaryConditionsPtrs;
  // Both body rod and cell rod are free
  for (unsigned int i = 0; i < 2; i++) {
    boundaryConditionsPtrs.push_back(&freeBC);
  }

  // Pack all forces together
  vector<ExternalForces *> externalForcesPtrs;

  // No force and torques on the body, all actuation comes from the cell
  // Hence all these coefficients are set to 0, and are required only for
  // legacy purposes
  int cellindex1 = 0;
  cellindex1 = cellindex1 + headindex;
  LocalTorque local1 = LocalTorque(Vector3(0.0, 0.0, 0.0e-3), w, 0);
  LocalTorque local2 = LocalTorque(Vector3(0.0, 0.0, 0.0e-3), w, 0);

  MultipleForces multipleForces1;
  multipleForces1.add(&local1);
  multipleForces1.add(&local2);

  MultipleForces *multipleForcesPtr1 = multipleForces1.get();
  externalForcesPtrs.push_back(multipleForcesPtr1);

  // muscular activity of cell
  // define contracting groups
  // This represents a 12 microN force amplitude beating at w=3.6Hz
  LocalForce muscle0 = LocalForce(12, w);
  MultipleForces multipleForces2;
  multipleForces2.add(&muscle0);

  MultipleForces *multipleForcesPtr2 = multipleForces2.get();
  externalForcesPtrs.push_back(multipleForcesPtr2);

  // Set up substrate properties and slenderbodytheory environment
  SlenderBodyTheoryEnvironment sbt =
      SlenderBodyTheoryEnvironment(spermPtr, dynamicViscosity);
  vector<Interaction *> substrateInteractionsPtrs;
  substrateInteractionsPtrs.push_back(&sbt);

  // Set up External Contact
  vector<ExternalContact *> externalcontactPtrs;
  /* The second and third argument are unimportant, but
         are preserved here for legacy purposes. Hence we simply
         set it to 0.0
  */
  /* This function takes care of assembly of body and cell
   */
  ExternalContact externalcontact =
      ExternalContact(rodPtrs, 0.0, 0.0, attachpoint);
  externalcontactPtrs.push_back(&externalcontact);

  // Set up Simple Connection -- Not used in this case, only in simple toy
  // examples
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

  // This call is done for legacy purposes (for data collection). In this
  // case we don't collect any data, but it is dumped in the flagella.txt
  // file below.
  poly.setWindowStats(0.5, 1.5);

  // Run simulation
  string outfileName = string("prova");
  const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime,
                                     povrayPerUnitTime, outfileName);

  // Throw exception if something went wrong
  if (!goodRun)
    throw "not good run in localized helical buckling, what is going on?";

  const vector<Vector3> avgVel = poly.getAverageVelocity();

  vector<REAL> fwdAvgVel;
  for (unsigned int i = 0; i < NOR; i++) {
    fwdAvgVel.push_back(avgVel[i] % spermdirection);
  }
  const REAL fitness = fwdAvgVel[0];

  return (fitness);
}

void Flagella::run() {
  const REAL fitness = _flagellaRun();
  exit(0);
}
