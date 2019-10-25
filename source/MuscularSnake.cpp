#include "MuscularSnake.h"

MuscularSnake::MuscularSnake(const int argc, const char **argv)
    : amp(0.0), w(0.0), v(0.0) {}

// Units in this case are m/kg/s

vector<REAL> MuscularSnake::_muscularsnakeRun() {
  vector<Rod *> rodPtrs;

  // Dumping frequencies (number of frames/dumps per unit time)
  const REAL diagPerUnitTime = 120;
  const REAL povrayPerUnitTime = 50;
  const REAL dt = 0.2e-7;
  const REAL timeSimulation = (5.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Snake Body -- same with the continuum snake
  // Driving parameters
  const int n = 100;
  const REAL density = 1000.0;
  const REAL L0 = 1.0;
  const REAL r0 = 0.025;
  // const REAL totalMass = density*M_PI*r0*r0*L0;
  const REAL E = 1e7;
  const REAL g = 9.81;
  const REAL Froude = 0.1;
  const REAL T = 1.0;
  const REAL mu = L0 / (T * T * g * Froude);
  const REAL totalMass = density * M_PI * r0 * r0 * L0;
  w = 2.0 * M_PI / T;

  const REAL dL0 = L0 / (double)n;  // length of cross-section element
  const REAL poissonRatio = 0.5;    // Incompressible
  const REAL G = E / (poissonRatio + 1.0);

  // Define rod
  const Vector3 snakedirection = Vector3(1.0, 0.0, 0.0);
  const Vector3 snakenormal = Vector3(0.0, 0.0, 1.0);
  const Vector3 snakeorigin = Vector3(0.0, 0.0, r0);

  // vector values for passing into straightrod_v

  /*
    We use these vector (*_v) values only when we want to
    include effects like tapering (where the radii of the snake body
    changes with the position). In this case, for comparison with a
    continuum snake, we implement the snake body as a straight rod
    with uniform cross-section (no-tapering.). Hence these vector *_v
    values are not used.
  */
  vector<REAL> r0_v = vector<REAL>(n);  //(mm)
  vector<Matrix3> J0_v = vector<Matrix3>(n);
  vector<Matrix3> B0_v = vector<Matrix3>(n - 1);
  vector<Matrix3> S0_v = vector<Matrix3>(n);

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

  /*
  We use these vector (*_v) values only when we want to
  include effects like tapering (where the radii of the snake body
  changes with the position). In this case, for comparison with a
  continuum snake, we implement the snake body as a straight rod
  with uniform cross-section (NO-tapering.). Hence these vector *_v
  values are not used.
*/
  // Head
  for (unsigned int i = 0; i < 5; i++) {
    r0_v[i] = r0 * (0.5 + i * 0.15);
    J0_v[i] = 10 * J0;
    B0_v[i] = 10 * B0;
    S0_v[i] = 10 * S0;
  }

  // Body
  for (unsigned int i = 5; i < (n - 30); i++) {
    r0_v[i] = r0;
    J0_v[i] = J0;
    B0_v[i] = B0;
    S0_v[i] = S0;
  }

  // Tail
  REAL At = 0.0;
  REAL It_1 = 0.0;
  REAL It_2 = 0.0;
  REAL It_3 = 0.0;
  Matrix3 It = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  Matrix3 Jt = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  Matrix3 Bt = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  Matrix3 St = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  REAL rt = 0.0;

  for (unsigned int i = (n - 30); i < n; i++) {
    rt = r0 * (1 - 0.5 / 29 * (i - (n - 30)));
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * dL0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    r0_v[i] = r0;
    J0_v[i] = J0;
    if (i < (n - 1)) {
      B0_v[i] = B0;
    }
    S0_v[i] = S0;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist = 0.0;
  const REAL nu = 5;
  const REAL relaxationNu = 0.0;
  const bool useSelfContact = false;

  /*
  Please use straightRod_v for tapering effects where you have to
  pass all vector (*_v) values for rod initialization.
  e.g. RodInitialConfigurations::straightRod_v(
      n, totalMass, r0_v, J0_v, B0_v, .....);

  */
  Rod *rod = RodInitialConfigurations::straightRod(
      n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, snakeorigin,
      snakedirection, snakenormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod);
  rod->update(0.0);
  rod->computeEnergies();

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers --- Group 1 and 3
  // number of muscle fibers
  NOR = 8;
  vector<pair<int, int>> attachpoint;

  // Driving parameters
  const int nm = 13;
  /*
    The muscle density is higher than the physiological one, since
    we lump many muscles (SSP-SP, LD and IC) into one actuator. These rods
    also represent the two tendons on the sides of the muscle which biologically
    have a higher density than the muscle itself. For these reasons,we set the
    muscle density to approximately twice the biological value.
  */
  const REAL densitym = 2000.0;
  const REAL Lm = 0.39;
  // const REAL totalMassm = densitym*M_PI*rm*rm*Lm;
  const REAL Em = 1e4;
  const REAL dLm = Lm / (double)nm;  // length of cross-section element

  const REAL poissonRatiom = 0.5;  // Incompressible
  const REAL Gm = Em / (poissonRatiom + 1.0);

  const Vector3 directionmuscle = Vector3(1.0, 0.0, 0.0);
  const Vector3 normalmuscle = Vector3(0.0, 0.0, 1.0);
  vector<Vector3> originmuscle;
  // Sets up starting positions of two antagonistic muscle pairs using the same
  // vector
  originmuscle.push_back(Vector3(4 * dL0, -r0 / 2, r0));
  originmuscle.push_back(Vector3(4 * dL0, r0 / 2, r0));
  originmuscle.push_back(Vector3(33 * dL0, -r0 / 2, r0));
  originmuscle.push_back(Vector3(33 * dL0, r0 / 2, r0));

  // Second moment of area for disk cross section for tendon
  /*
  In our simulation, we lump many biological tendons into one computational
  tendon. As a result, our computational tendon is bigger in size, set as rm_t
  below
  */
  const REAL rm_t = 0.003;
  const REAL Am_t = M_PI * rm_t * rm_t;
  const REAL I0_1m_t = Am_t * Am_t / (4.0 * M_PI);
  const REAL I0_2m_t = I0_1m_t;
  const REAL I0_3m_t = 2.0 * I0_1m_t;
  const Matrix3 I0m_t =
      Matrix3(I0_1m_t, 0.0, 0.0, 0.0, I0_2m_t, 0.0, 0.0, 0.0, I0_3m_t);

  // Mass inertia matrix for disk cross section
  const Matrix3 Jm_t = densitym * dLm * I0m_t;

  /*
  The biological tendons have a high Young's modulus E.,but are very slender.
  As a result, they resist extension (stretch) but can bend easily.

  Due to our decision to lump tendons and in order to mimic the above behavior
  of the biological tendons, we use a lower Young's
  Modulus and harden the stiffness of the shear and stretch modes only.
  Numerically, this is done by putting a pre-factor of 50000 before the
  shear/stretch matrix below. The actual value of the prefactor does not matter,
  what is important is that it is a high value to high stretch/shear stiffness.
  */
  // Bending matrix
  const Matrix3 Bm_t = Matrix3(Em * I0_1m_t, 0.0, 0.0, 0.0, Em * I0_2m_t, 0.0,
                               0.0, 0.0, Gm * I0_3m_t);

  // Shear matrix
  const Matrix3 Sm_t =
      50000 * Matrix3((4.0 / 3.0) * Gm * Am_t, 0.0, 0.0, 0.0,
                      (4.0 / 3.0) * Gm * Am_t, 0.0, 0.0, 0.0, Em * Am_t);

  // Second moment of area for disk cross section for muscle
  const REAL rm = 0.006;
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
  for (unsigned int i = 0; i < 4; i++) {
    rm_v[i] = rm_t;
    Jm_v[i] = Jm_t;
    Bm_v[i] = Bm_t;
    Sm_v[i] = Sm_t;
  }
  Bm_v[3] = (Bm_t + Bm) / 2;  // Continuity between thin tendon and thick muscle
  // other tendon
  for (unsigned int i = 9; i < 13; i++) {
    rm_v[i] = rm_t;
    Jm_v[i] = Jm_t;
    if (i < 12) Bm_v[i] = Bm_t;
    Sm_v[i] = Sm_t;
  }
  // muscle
  for (unsigned int i = 4; i < 9; i++) {
    rm_v[i] = rm;
    Jm_v[i] = Jm;
    Bm_v[i] = Bm;
    Sm_v[i] = Sm;
  }
  Bm_v[8] = (Bm_t + Bm) / 2;

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist_m = 0.0;
  const REAL nu_m = 5.0;  // 5
  const REAL relaxationNu_m = 0.0;
  const bool useSelfContact_m = false;

  // Half of the muscles (first and third) initialized here
  for (unsigned int i = 0; i < (NOR / 2); i++) {
    Rod *rod_1 = RodInitialConfigurations::straightRod_v(
        nm, densitym, rm_v, Jm_v, Bm_v, Sm_v, Lm, initialTotalTwist_m,
        originmuscle[i], directionmuscle, normalmuscle, nu_m, relaxationNu_m,
        useSelfContact_m);
    rodPtrs.push_back(rod_1);
    rod_1->update(0.0);
  }
  // mindex is muscle index to regulate muscle activation and coordination while
  // firing. See usage in External contact.
  rodPtrs[1]->mindex = 4;
  rodPtrs[2]->mindex = 4;
  rodPtrs[3]->mindex = 33;
  rodPtrs[4]->mindex = 33;

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers --- Group 2 and 4 (same procedure as 1 and 3 above)

  // Driving parameters
  const int nm_2 = 11;
  const REAL Lm_2 = 0.33;

  originmuscle.push_back(Vector3(23 * dL0, -r0 / 2, r0));
  originmuscle.push_back(Vector3(23 * dL0, r0 / 2, r0));
  originmuscle.push_back(Vector3(61 * dL0, -r0 / 2, r0));
  originmuscle.push_back(Vector3(61 * dL0, r0 / 2, r0));

  vector<REAL> rm_v2 = vector<REAL>(nm_2);
  vector<Matrix3> Jm_v2 = vector<Matrix3>(nm_2);
  vector<Matrix3> Bm_v2 = vector<Matrix3>(nm_2 - 1);
  vector<Matrix3> Sm_v2 = vector<Matrix3>(nm_2);

  // tendon
  for (unsigned int i = 0; i < 4; i++) {
    rm_v2[i] = rm_t;
    Jm_v2[i] = Jm_t;
    Bm_v2[i] = Bm_t;
    Sm_v2[i] = Sm_t;
  }
  Bm_v2[3] = (Bm_t + Bm) / 2;
  for (unsigned int i = 9; i < 11; i++) {
    rm_v2[i] = rm_t;
    Jm_v2[i] = Jm_t;
    if (i < 10) Bm_v2[i] = Bm_t;
    Sm_v2[i] = Sm_t;
  }
  // muscle
  for (unsigned int i = 4; i < 9; i++) {
    rm_v2[i] = rm;
    Jm_v2[i] = Jm;
    Bm_v2[i] = Bm;
    Sm_v2[i] = Sm;
  }
  Bm_v2[8] = (Bm_t + Bm) / 2;

  for (unsigned int i = 4; i < 8; i++) {
    Rod *rod_2 = RodInitialConfigurations::straightRod_v(
        nm_2, densitym, rm_v2, Jm_v2, Bm_v2, Sm_v2, Lm_2, initialTotalTwist_m,
        originmuscle[i], directionmuscle, normalmuscle, nu_m, relaxationNu_m,
        useSelfContact_m);
    rodPtrs.push_back(rod_2);
    rod_2->update(0.0);
  }
  rodPtrs[5]->mindex = 23;
  rodPtrs[6]->mindex = 23;
  rodPtrs[7]->mindex = 61;
  rodPtrs[8]->mindex = 61;

  //-----------------------------------------------------------------------------------------------------------------
  // Pack boundary conditions
  FreeBC freeBC = FreeBC();
  vector<RodBC *> boundaryConditionsPtrs;
  // NOR + 1 , the extra 1 rod is the body
  for (unsigned int i = 0; i < (NOR + 1); i++) {
    boundaryConditionsPtrs.push_back(&freeBC);
  }

  // Pack all forces together (no forces applied)
  vector<ExternalForces *> externalForcesPtrs;

  // Gravity only for snake body -- to compare with the continuum case, we
  // assume snakes have the same overall weight
  GravityForce gravity = GravityForce(Vector3(0.0, 0.0, -g));
  MultipleForces multipleForces1;
  multipleForces1.add(&gravity);
  MultipleForces *multipleForcesPtr1 = multipleForces1.get();
  externalForcesPtrs.push_back(multipleForcesPtr1);

  // muscular activity of muscles are applied in ExternalContact, not here
  MuscleContraction muscle0 = MuscleContraction(0.0, w);
  MultipleForces multipleForces2;
  multipleForces2.add(&muscle0);
  MultipleForces *multipleForcesPtr2 = multipleForces2.get();

  // Kept for legacy purposes
  for (unsigned int i = 0; i < NOR; i++) {
    externalForcesPtrs.push_back(multipleForcesPtr2);
  }

  // Define plane
  const REAL angle = 0.0;
  const Vector3 originPlane = Vector3(0.0, 0.0, 0.0);
  const Vector3 normalPlane = Vector3(0.0, sin(angle), cos(angle)).unitize();

  // Set up substrate properties and snake-plane interaction object
  const REAL kPlane = 1.0;    // stiffness of ground
  const REAL nuPlane = 1e-6;  // viscous damping of ground
  const REAL muKineticForward = mu;
  const REAL muKineticSideways = 2.0 * muKineticForward;
  const REAL muKineticBackward = 1.5 * muKineticForward;
  const REAL muStaticForward = 2.0 * mu;
  const REAL muStaticSideways = 2.0 * muStaticForward;
  const REAL muStaticBackward = 1.5 * muStaticForward;
  const REAL vStatic = 1e-8;

  vector<Rod *> snakeptr;
  snakeptr.push_back(rodPtrs[0]);
  AnisotropicFrictionPlaneInteraction frictionPlane =
      AnisotropicFrictionPlaneInteraction(
          snakeptr, normalPlane, originPlane, kPlane, nuPlane, muKineticForward,
          muKineticBackward, muKineticSideways, muStaticForward,
          muStaticBackward, muStaticSideways, vStatic);
  vector<Interaction *> substrateInteractionsPtrs;
  substrateInteractionsPtrs.push_back(&frictionPlane);

  // Set up External Contact
  vector<ExternalContact *> externalcontactPtrs;
  /* The second and third argument are unimportant, but
         are preserved here for legacy purposes. Hence we simply
         set it to 0.0
  This function takes care of the various connections and muscle actuations
  within the muscularsnake case.
  */
  ExternalContact externalcontact =
      ExternalContact(rodPtrs, 0.0, 0.0, attachpoint);
  externalcontactPtrs.push_back(&externalcontact);

  // Set up Simple Connection -- Not used in this case, only in simple toy
  // examples
  vector<SimpleConnection *> simpleconnectionPtrs;
  SimpleConnection simpleconnection = SimpleConnection(rodPtrs);
  simpleconnectionPtrs.push_back(&simpleconnection);
  //******************************************************************************

  // Set up integrator (define integration order)
  PolymerIntegrator *integrator = new PositionVerlet2nd(
      rodPtrs, externalForcesPtrs, boundaryConditionsPtrs,
      substrateInteractionsPtrs, externalcontactPtrs, simpleconnectionPtrs);

  // Instantiate simulator
  Polymer poly = Polymer(integrator);

  // This call is done for legacy purposes (for data collection). In this
  // case we don't collect any data, but it is dumped in the flagella.txt
  // file below.
  poly.setWindowStats(3.0 * T, 4.0 * T);

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
    fwdAvgVel.push_back(avgVel[i] % directionmuscle);
  }
  const vector<REAL> fitness = fwdAvgVel;

  return (fitness);
}

void MuscularSnake::run() {
  const vector<REAL> fitness = _muscularsnakeRun();
  exit(0);
}
