#include "Walker.h"

Walker::Walker(const int argc, const char **argv)
    : YoungM(0.0), LegL(0.0), Location(0.0) {}

// Units in this case are mm/g/s

REAL Walker::_walkerRun() {
  vector<Rod *> rodPtrs;

  // Dumping frequencies (number of frames/dumps per unit time)
  const REAL diagPerUnitTime = 120;
  const REAL povrayPerUnitTime = 50;
  REAL dt = 0.75e-7;  // ms

  const REAL timeSimulation = 10.0;  // ms

  /*
  In the case of the bio-hybrid walker, we consider bots displacing in a shallow
  solution that reaches the muscle tissue. As a consequence, the tissue is
  immersed and suspended in the fluid (due to buoyant forces – density of
  muscles close to density of liquid). However, the scaffold is mostly exposed
  to the atmosphere so that we neglect buoyancy and hydrodynamic loads acting on
  it.
  */

  //-----------------------------------------------------------------------------------------------------------------
  // Define Beam
  // Driving parameters
  const int n_B = 16;
  const REAL L0_B = 14.0;                 //(mm)
  const REAL dL0_B = L0_B / (double)n_B;  // length of cross-section element
  const REAL density = 1.12e-3;           //(g/mm^-3)
  const REAL g = -9810.0;
  YoungM = 3.194;                 // Class member
  const REAL E = YoungM * 1e5;    //(g*mm^-1*ms^-2) --- original 3.194
  const REAL poissonRatio = 0.5;  // Incompressible
  const REAL G = E / (poissonRatio + 1.0);

  // fluid property -- no fluid when considering shallow water
  const REAL T = 1.0 / 1.8;           //(s)
  const REAL fluidDensity = 1.15e-3;  //(g/mm^3)
  const REAL RE = 1.8e-2;
  const REAL dynamicViscosity =
      1.1e-3;  // dynamicViscosity is not used in this simulation
  const REAL w = 2.0 * M_PI / T;

  // Define rod
  const Vector3 beamdirection = Vector3(0.0, -1.0, 0.0);
  const Vector3 beamnormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 beamorigin = Vector3(0.0, 0.0, 3.5);

  // vector values
  vector<REAL> rB_v = vector<REAL>(n_B);  //(mm)
  vector<Matrix3> JB_v = vector<Matrix3>(n_B);
  vector<Matrix3> BB_v = vector<Matrix3>(n_B - 1);
  vector<Matrix3> SB_v = vector<Matrix3>(n_B);

  const REAL rsmall = 0.8786;  //(mm)0.8786
  const REAL rlarge = 0.8786;  //(mm)0.8786

  /* The bending stiffness of the beam is set using the equivalent radius
  "requal" so as to match the experimental bending stiffness of the beam.
  // See Simulation and fabrication of stronger, larger, and faster walking
  // biohybrid machines. Pagan‐Diaz, G. J., Zhang, X., Grant, L., et al.
  // Advanced Functional Materials (2018) for more details
*/
  REAL requal = 0.5842;

  REAL Asmall = M_PI * requal * requal;
  REAL Alarge = M_PI * rlarge * rlarge;

  // Second moment of area for disk cross section
  REAL I0_1_s = Asmall * Asmall / (4.0 * M_PI);
  REAL I0_2_s = I0_1_s;
  REAL I0_3_s = 2.0 * I0_1_s;
  Matrix3 I0_s = Matrix3(I0_1_s, 0.0, 0.0, 0.0, I0_2_s, 0.0, 0.0, 0.0, I0_3_s);
  // Mass inertia matrix for disk cross section
  Matrix3 J0_s = density * dL0_B * I0_s;

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
  Matrix3 J0_l = density * dL0_B * I0_l;

  // Bending matrix
  Matrix3 B0_l =
      Matrix3(E * I0_1_l, 0.0, 0.0, 0.0, E * I0_2_l, 0.0, 0.0, 0.0, G * I0_3_l);
  // Shear matrix
  Matrix3 S0_l = Matrix3((4.0 / 3.0) * G * Alarge, 0.0, 0.0, 0.0,
                         (4.0 / 3.0) * G * Alarge, 0.0, 0.0, 0.0, E * Alarge);

  // 4 is the limit of segments connected to the leg
  for (unsigned int i = 0; i < 4; i++) {
    rB_v[i] = rlarge;
    JB_v[i] = J0_l;
    BB_v[i] = B0_l;
    if (i == 3) BB_v[i] = (B0_s + B0_l) / 2.0;
    SB_v[i] = S0_l;
  }
  // These represent the center segment of the beam
  for (unsigned int i = 4; i < 12; i++) {
    rB_v[i] = rsmall;
    JB_v[i] = J0_s;
    if (i < 11) BB_v[i] = B0_s;
    if (i == 11) BB_v[i] = (B0_s + B0_l) / 2.0;
    SB_v[i] = S0_s;
  }
  // These 12-16 elements are connected to the other leg
  for (unsigned int i = 12; i < 16; i++) {
    rB_v[i] = rlarge;
    JB_v[i] = J0_l;
    if (i < 15) BB_v[i] = B0_l;
    SB_v[i] = S0_l;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist = 0.0;
  const REAL nu = 0.3;
  const REAL relaxationNu = 0.1;
  const bool useSelfContact = false;
  const REAL factorW = 0.5;  // adjustment to match the overall weight, account
                             // for irregular shape of the actual walker

  // This call is kept for legacy purposes and is useful to include hydrodynamic
  // effects
  // We DO NOT use it in this case.
  vector<Rod *> spermPtr;

  for (unsigned int i = 0; i < 2; i++) {
    Rod *rod = RodInitialConfigurations::straightRod_v(
        n_B, density * factorW, rB_v, JB_v, BB_v, SB_v, L0_B, initialTotalTwist,
        beamorigin + i * 3.5 * Vector3(1.0, 0.0, 0.0), beamdirection,
        beamnormal, nu, relaxationNu, useSelfContact);
    rodPtrs.push_back(rod);
    rod->update(0.0);
    rod->computeEnergies();
    spermPtr.push_back(rod);
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Define Leg
  // Driving parameters
  const int n_L = 6;
  const REAL L0_L = 3.4;
  const REAL dL0_L = L0_L / (double)n_L;  // length of cross-section element
  const REAL r_L = 1.75;
  const REAL totalMass_L = density * M_PI * r_L * r_L * L0_L;

  vector<Vector3> Legdirection;
  Legdirection.push_back(-rodPtrs[0]->Q[2][1]);
  Legdirection.push_back(-rodPtrs[0]->Q[14][1]);
  Legdirection.push_back(-rodPtrs[1]->Q[2][1]);
  Legdirection.push_back(-rodPtrs[1]->Q[14][1]);

  const Vector3 Legnormal = Vector3(1.0, 0.0, 0.0);
  vector<Vector3> Legorigin;

  Legorigin.push_back(rodPtrs[0]->x[2]);
  Legorigin.push_back(rodPtrs[0]->x[14]);
  Legorigin.push_back(rodPtrs[1]->x[2]);
  Legorigin.push_back(rodPtrs[1]->x[14]);

  // Second moment of area for disk cross section for tendon
  REAL A_L = M_PI * r_L * r_L;
  REAL IL_1 = A_L * A_L / (4.0 * M_PI);
  REAL IL_2 = IL_1;
  REAL IL_3 = 2.0 * IL_1;
  Matrix3 IL = Matrix3(IL_1, 0.0, 0.0, 0.0, IL_2, 0.0, 0.0, 0.0, IL_3);

  // Mass inertia matrix for disk cross section
  Matrix3 J_L = density * dL0_L * IL;

  // Bending matrix
  Matrix3 B_L =
      Matrix3(E * IL_1, 0.0, 0.0, 0.0, E * IL_2, 0.0, 0.0, 0.0, G * IL_3);

  // Shear matrix
  Matrix3 S_L = Matrix3((4.0 / 3.0) * G * A_L, 0.0, 0.0, 0.0,
                        (4.0 / 3.0) * G * A_L, 0.0, 0.0, 0.0, E * A_L);

  REAL factorL = 1.0;
  LegL = 2.8;  // Class member
  // 4 here indicates the total number of legs
  for (unsigned int i = 0; i < 4; i++) {
    if ((i % 2) == 0) factorL = LegL / L0_L;
    if ((i % 2) == 1) factorL = 1.0;
    Rod *rodL = RodInitialConfigurations::straightRod(
        n_L, factorL * totalMass_L, r_L, factorL * J_L, B_L, S_L,
        factorL * L0_L, initialTotalTwist, Legorigin[i], Legdirection[i],
        Legnormal, nu, relaxationNu, useSelfContact);
    rodPtrs.push_back(rodL);
    rodL->update(0.0);
    spermPtr.push_back(rodL);
  }

  // MaxHeight represents a proxy variable that can be used to pass
  // information into the ExternalContact module
  rodPtrs[2]->MaxHeight = LegL;

  //-----------------------------------------------------------------------------------------------------------------
  // Define Cross -- connecting beams to beams and legs to legs
  // Driving parameters

  const int n_C = 5;
  const REAL L0_C = 3.5;
  const REAL dL0_C = L0_C / (double)n_C;  // length of cross-section element
  const REAL r_C = 0.72;
  const REAL totalMass_C = density * M_PI * r_C * r_C * L0_C;

  const Vector3 Crossdirection = Vector3(1.0, 0.0, 0.0);
  const Vector3 Crossnormal = Vector3(0.0, 0.0, 1.0);
  vector<Vector3> Crossorigin;

  // This takes input of the leg positions
  Crossorigin.push_back(rodPtrs[2]->x[0]);
  Crossorigin.push_back(rodPtrs[3]->x[0]);
  Crossorigin.push_back(rodPtrs[2]->x[4]);
  Crossorigin.push_back(rodPtrs[3]->x[4]);

  // Second moment of area for disk cross section for cross
  REAL A_C = M_PI * r_C * r_C;
  REAL IC_1 = A_C * A_C / (4.0 * M_PI);
  REAL IC_2 = IC_1;
  REAL IC_3 = 2.0 * IC_1;
  Matrix3 IC = Matrix3(IC_1, 0.0, 0.0, 0.0, IC_2, 0.0, 0.0, 0.0, IC_3);

  // Mass inertia matrix for disk cross section
  Matrix3 J_C = density * dL0_C * IC;

  // Bending matrix
  Matrix3 B_C =
      Matrix3(E * IC_1, 0.0, 0.0, 0.0, E * IC_2, 0.0, 0.0, 0.0, G * IC_3);

  // Shear matrix
  Matrix3 S_C = Matrix3((4.0 / 3.0) * G * A_C, 0.0, 0.0, 0.0,
                        (4.0 / 3.0) * G * A_C, 0.0, 0.0, 0.0, E * A_C);

  for (unsigned int i = 0; i < 4; i++) {
    Rod *rodC = RodInitialConfigurations::straightRod(
        n_C, totalMass_C, r_C, J_C, B_C, S_C, L0_C, initialTotalTwist,
        Crossorigin[i], Crossdirection, Crossnormal, nu, relaxationNu,
        useSelfContact);
    rodPtrs.push_back(rodC);
    rodC->update(0.0);
    spermPtr.push_back(rodC);  // Not used
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers

  // Driving parameters
  const int n_M = 24;
  const REAL density_M = 1.06e-3;  //(g/mm^-3)
  const REAL E_M = 1.0e4;          //(g*mm^-1*ms^-2)
  const REAL G_M = E_M / (poissonRatio + 1.0);

  // Define rod
  const Vector3 musclenormal = Vector3(0.0, 0.0, 1.0);
  vector<Vector3> Positions_M1 = vector<Vector3>(n_M + 1);
  vector<Vector3> Positions_M2 = vector<Vector3>(n_M + 1);

  const REAL rmuscle = 0.29;  // mm
  // Corner of the muscle ring structure
  const Vector3 corner1 = Vector3(-2.05, -14.35, 1.5);
  const Vector3 corner2 = Vector3(-2.05, -14.35, 0.9);

  // Define Shape of muscle ring as a rectangular path around the leg
  for (unsigned int i = 0; i < 8; i++) {
    Positions_M1[i] = corner1 + i * Vector3((7.6 / 7.0), 0.0, 0.0);
    Positions_M2[i] = corner2 + i * Vector3((7.6 / 7.0), 0.0, 0.0);
  }
  for (unsigned int i = 8; i < 12; i++) {
    Positions_M1[i] =
        Positions_M1[7] + (i - 7) * Vector3(0.0, (4.6 / 4.0), 0.0);
    Positions_M2[i] =
        Positions_M2[7] + (i - 7) * Vector3(0.0, (4.6 / 4.0), 0.0);
  }
  for (unsigned int i = 12; i < 19; i++) {
    Positions_M1[i] =
        Positions_M1[11] + (i - 11) * Vector3(-(7.6 / 7.0), 0.0, 0.0);
    Positions_M2[i] =
        Positions_M2[11] + (i - 11) * Vector3(-(7.6 / 7.0), 0.0, 0.0);
  }
  for (unsigned int i = 19; i < 23; i++) {
    Positions_M1[i] =
        Positions_M1[18] + (i - 18) * Vector3(0.0, -(4.6 / 4.0), 0.0);
    Positions_M2[i] =
        Positions_M2[18] + (i - 18) * Vector3(0.0, -(4.6 / 4.0), 0.0);
  }

  // Loop the muscle ring back to the origin and glue it
  Positions_M1[23] = Positions_M1[1];
  Positions_M2[23] = Positions_M2[1];
  Positions_M1[24] = Positions_M1[2];
  Positions_M2[24] = Positions_M2[2];

  REAL Amuscle = M_PI * rmuscle * rmuscle;

  // Second moment of area for disk cross section
  REAL I0_1_m = Amuscle * Amuscle / (4.0 * M_PI);
  REAL I0_2_m = I0_1_m;
  REAL I0_3_m = 2.0 * I0_1_m;
  Matrix3 I0_m = Matrix3(I0_1_m, 0.0, 0.0, 0.0, I0_2_m, 0.0, 0.0, 0.0, I0_3_m);
  // Mass inertia matrix for disk cross section
  Matrix3 J0_m =
      density_M * 0.1 *
      I0_m;  // Here J matrix is not exactly computed with varying
             // dl. We approximate dl=0.1. This doesn't effect the
             // dynamics at all since muscles do not rotate. Even with varying
             // dl, the dynamics was found to be exactly the same

  // Bending matrix
  Matrix3 B0_m = Matrix3(E_M * I0_1_m, 0.0, 0.0, 0.0, E_M * I0_2_m, 0.0, 0.0,
                         0.0, G_M * I0_3_m);
  // Shear matrix
  Matrix3 S0_m =
      Matrix3((4.0 / 3.0) * G_M * Amuscle, 0.0, 0.0, 0.0,
              (4.0 / 3.0) * G_M * Amuscle, 0.0, 0.0, 0.0, E_M * Amuscle);

  // Params for initializig muscle as curved rod
  vector<REAL> rM_v = vector<REAL>(n_M, rmuscle);
  // Rigidities
  vector<Matrix3> JM_v = vector<Matrix3>(n_M, J0_m);
  vector<Matrix3> BM_v = vector<Matrix3>(n_M - 1, B0_m);
  vector<Matrix3> SM_v = vector<Matrix3>(n_M, S0_m);

  const REAL nu_M = 0.05;
  const REAL relaxationNu_M = 0.0;

  // The first muscle ring layer
  Rod *rodm = RodInitialConfigurations::curvedRod(
      n_M, density_M, rM_v, JM_v, BM_v, SM_v, initialTotalTwist, Positions_M1,
      musclenormal, nu_M, relaxationNu_M, useSelfContact);
  rodPtrs.push_back(rodm);
  rodm->update(0.0);
  rodm->computeEnergies();
  // The second muscle ring layer
  Rod *rodm1 = RodInitialConfigurations::curvedRod(
      n_M, density_M, rM_v, JM_v, BM_v, SM_v, initialTotalTwist, Positions_M2,
      musclenormal, nu_M, relaxationNu_M, useSelfContact);
  rodPtrs.push_back(rodm1);
  rodm1->update(0.0);
  rodm1->computeEnergies();

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers

  // The ring on the other side, the process of assembly is similar
  // to the one described above
  vector<Vector3> Positions_M3 = vector<Vector3>(n_M + 1);
  vector<Vector3> Positions_M4 = vector<Vector3>(n_M + 1);
  // Define Shape
  const Vector3 corner3 = Vector3(-2.05, 0.33, 1.5);
  const Vector3 corner4 = Vector3(-2.05, 0.33, 0.9);

  for (unsigned int i = 0; i < 8; i++) {
    Positions_M3[i] = corner3 + i * Vector3((7.6 / 7.0), 0.0, 0.0);
    Positions_M4[i] = corner4 + i * Vector3((7.6 / 7.0), 0.0, 0.0);
  }
  for (unsigned int i = 8; i < 12; i++) {
    Positions_M3[i] =
        Positions_M3[7] + (i - 7) * Vector3(0.0, (-4.6 / 4.0), 0.0);
    Positions_M4[i] =
        Positions_M4[7] + (i - 7) * Vector3(0.0, (-4.6 / 4.0), 0.0);
  }
  for (unsigned int i = 12; i < 19; i++) {
    Positions_M3[i] =
        Positions_M3[11] + (i - 11) * Vector3(-(7.6 / 7.0), 0.0, 0.0);
    Positions_M4[i] =
        Positions_M4[11] + (i - 11) * Vector3(-(7.6 / 7.0), 0.0, 0.0);
  }
  for (unsigned int i = 19; i < 23; i++) {
    Positions_M3[i] =
        Positions_M3[18] + (i - 18) * Vector3(0.0, (4.6 / 4.0), 0.0);
    Positions_M4[i] =
        Positions_M4[18] + (i - 18) * Vector3(0.0, (4.6 / 4.0), 0.0);
  }
  Positions_M3[23] = Positions_M3[1];
  Positions_M4[23] = Positions_M4[1];
  Positions_M3[24] = Positions_M3[2];
  Positions_M4[24] = Positions_M4[2];

  Rod *rodm2 = RodInitialConfigurations::curvedRod(
      n_M, density_M, rM_v, JM_v, BM_v, SM_v, initialTotalTwist, Positions_M3,
      musclenormal, nu_M, relaxationNu_M, useSelfContact);
  rodPtrs.push_back(rodm2);
  rodm2->update(0.0);
  rodm2->computeEnergies();
  Rod *rodm3 = RodInitialConfigurations::curvedRod(
      n_M, density_M, rM_v, JM_v, BM_v, SM_v, initialTotalTwist, Positions_M4,
      musclenormal, nu_M, relaxationNu_M, useSelfContact);
  rodPtrs.push_back(rodm3);
  rodm3->update(0.0);
  rodm3->computeEnergies();

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers
  // Straight muscles
  const int n_M2 = 4;
  const REAL L0_M = 5.48;
  const REAL dL0_M = L0_M / (double)n_M2;  // length of cross-section element
  const REAL r_M = 0.227;
  const REAL totalMass_M = density_M * M_PI * r_M * r_M * L0_M;

  const Vector3 Muscledirection = Vector3(0.0, -1.0, 0.0);
  const Vector3 Musclenormal = Vector3(0.0, 0.0, 1.0);
  vector<Vector3> Muscleorigin;

  // Takes the position of the center of the ring
  Muscleorigin.push_back(Positions_M3[13]);
  Muscleorigin.push_back(Positions_M3[14]);
  Muscleorigin.push_back(Positions_M3[15]);
  Muscleorigin.push_back(Positions_M3[16]);

  for (unsigned int i = 0; i < 4; i++) {
    Rod *rodm4 = RodInitialConfigurations::straightRod(
        n_M2, totalMass_M, r_M, J0_m, B0_m, S0_m, L0_M, initialTotalTwist,
        Muscleorigin[i], Muscledirection, Musclenormal, nu_M, relaxationNu_M,
        useSelfContact);
    rodPtrs.push_back(rodm4);
    rodm4->update(0.0);
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Pack boundary conditions
  FreeBC freeBC = FreeBC();
  vector<RodBC *> boundaryConditionsPtrs;
  // FreeBCs for all (18) rods in the simulation
  for (unsigned int i = 0; i < 18; i++) {
    boundaryConditionsPtrs.push_back(&freeBC);
  }

  // Pack all forces together
  vector<ExternalForces *> externalForcesPtrs;

  MultipleForces multipleForces1;
  GravityForce gravity = GravityForce(Vector3(0.0, 0.0, g));
  multipleForces1.add(&gravity);
  MultipleForces *multipleForcesPtr1 = multipleForces1.get();
  for (unsigned int i = 0; i < 10; i++) {
    externalForcesPtrs.push_back(multipleForcesPtr1);
  }
  // Here we don't apply gravity to the muscles by considering shallow
  // water/solution will make the muscle tissue naturally buoyant Muscle force
  // is applied in ExternalContact function, not here.
  MultipleForces multipleForces2;
  GravityForce gravity2 = GravityForce(Vector3(0.0, 0.0, 0.0));
  multipleForces2.add(&gravity2);
  MultipleForces *multipleForcesPtr2 = multipleForces2.get();
  for (unsigned int i = 0; i < 10; i++) {
    externalForcesPtrs.push_back(multipleForcesPtr2);
  }

  // Set up substrate properties and slenderbodytheory environment
  // In this case, although we set it up here, we manually switch off
  // the effects of SlenderBodyTheoryEnvironment in the PolymerIntergrator.cpp
  // file, line 25:28, under the assumption of shallow water.
  SlenderBodyTheoryEnvironment sbt =
      SlenderBodyTheoryEnvironment(spermPtr, dynamicViscosity);
  vector<Interaction *> substrateInteractionsPtrs;
  substrateInteractionsPtrs.push_back(&sbt);

  // Set up External Contact
  vector<ExternalContact *> externalcontactPtrs;
  vector<pair<int, int>> attachpoint;
  /* The second and third argument are unimportant, but
         are preserved here for legacy purposes. Hence we simply
         set it to 0.0

    This function takes care of the various connections, ground setup (friction)
    and muscle actuations within the walker case.
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
  poly.setWindowStats(0.5, 2.0);

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
    fwdAvgVel.push_back(avgVel[i] % Vector3(0.0, 1.0, 0.0));
  }

  REAL fitness = fwdAvgVel[0];
  if (fitness < (-1.0)) fitness = 0.0;

  return (fitness);
}

void Walker::run() {
  const REAL fitness = _walkerRun();
  exit(0);
}
