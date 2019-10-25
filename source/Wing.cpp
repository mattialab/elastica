#include "Wing.h"

Wing::Wing(const int argc, const char **argv) : amp(0.0), w(0.0), v(0.0) {}

/*
This is an example study demonstrating how our numerical solver can be used to
qualitatively simulate a complex biological assembly consisting of O(1000) rods
with non-linear interactions. Due to numerical constraints, an exact replica of
a flapping pigeon wing is computationally expensive. To relax these constraints,
we rescale the system while retaining similar dynamics as follows:

Notation: hatted denotes real biological quantities of pigeon flapping.
          No hat denotes the corresponding computational quantities.

1. density (rho) is kept constant for all components --> rho = rho_hat.
2. Dimension (L) is scaled up by 2.5 times --> L = 2.5*L_hat.
3. Time (t) is scaled up by 4 times --> t = 4*t_hat.

Thus, all quantities are scaled according to these rules. We expect the dynamics
of our system to predominantly depend on input muscle force F and deformation of
the feathers (characterized by EI). Upon scaling of these critical parameters,
we see:
1. Muscle Force F = rho*l^4/t^2 = 2.44*F_hat.
2. Bending Stiffness EI = rho*l^3*l*l^2/t^2 = 15.26*EI_hat.
3. Young's Modulus E = rho*l^3*l/t^2/l^2 = 0.39*E_hat.

These mean that the input computational muscle force and feather bending
stiffness are 2.44 and 15.26 times the real biological values respectively. We
thus scale these two parameters accordingly to reflect accurate dynamics.

Here, the scaling of E for bones and tendons are less important because they are
almost rigid (very high E). The E scaling for muscles are also less important
since they are very soft (very low E). Nevertheless, we still scale the values
according to point 3 above.
*/

// Units in this case are mm/g/s

vector<REAL> Wing::_wingRun() {
  vector<Rod *> rodPtrs;

  // Dumping frequencies (number of frames/dumps per unit time)
  const REAL diagPerUnitTime = 120;
  const REAL povrayPerUnitTime = 200;
  const REAL dt = 0.5e-7;              // 0.5e-7
  const REAL timeSimulation = (0.51);  // Scaled Time

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- HUMERUS
  // Driving parameters
  const int n = 10;
  /*
  We assume the same density with human bone.
  Noticing that the bird bone is generally hollow, we divide the density by 2
  when initializing the rod by assuming 50% hollowness. This assumption is based
  on rough measurements from paper "Extreme lightweight structures: avian
  feathers and bones" Sullivan et al
  */
  const REAL density = 1.75e-3;
  const REAL L0 = 150;  // Scaled Length - same for all dimensions below
  const REAL E = 15e9;  // does not change with unit

  const REAL dL0 = L0 / (double)n;  // length of cross-section element
  const REAL poissonRatio = 0.5;    // Incompressible
  const REAL G = E / (poissonRatio + 1.0);

  // Define rod
  const Vector3 humerusdirection = Vector3(0.8, -1.0, 0.0).unitize();
  const Vector3 humerusnormal = Vector3(0.0, 0.0, 1.0);
  const Vector3 humerusorigin = Vector3(0.0, 0.0, 150);

  // Humerus shape
  // vector values
  vector<REAL> r0_v = vector<REAL>(n);  //(mm)
  vector<Matrix3> J0_v = vector<Matrix3>(n);
  vector<Matrix3> B0_v = vector<Matrix3>(n - 1);
  vector<Matrix3> S0_v = vector<Matrix3>(n);

  REAL r0 = 8.4;  // middle radius
  // Second moment of area for disk cross section
  const REAL A0 = M_PI * r0 * r0;
  const REAL I0_1 = A0 * A0 / (4.0 * M_PI);
  const REAL I0_2 = I0_1;
  const REAL I0_3 = 2.0 * I0_1;
  const Matrix3 I0 = Matrix3(I0_1, 0.0, 0.0, 0.0, I0_2, 0.0, 0.0, 0.0, I0_3);

  /*
  Given the fact that bones are rigid, we know that they resist bending and
  shearing modes. Hence we set their B and S matrices to be high and do not
  consider their internal structure (for e.g. hollowness). Nevertheless, we
  tested different B,S, and J (Inertia) and observed that the resulting dynamics
  remain essentially the same (<1% difference). More detailed bio-mechanical
  assessments of the bone structure is beyond the scope of this paper, but can
  be included by the user.
  */

  // Mass inertia matrix for disk cross section
  const Matrix3 J0 = density * 15.0 * I0;

  // Bending matrix
  Matrix3 B0 =
      Matrix3(E * I0_1, 0.0, 0.0, 0.0, E * I0_2, 0.0, 0.0, 0.0, G * I0_3);
  // Shear matrix
  Matrix3 S0 = Matrix3((4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0,
                       (4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0, E * A0);

  // middle shaft
  for (unsigned int i = 2; i < 8; i++) {
    r0_v[i] = r0;
    J0_v[i] = J0;
    if (i < 7) B0_v[i] = B0;
    S0_v[i] = S0;
  }

  // Heads
  REAL At = 0.0;
  REAL It_1 = 0.0;
  REAL It_2 = 0.0;
  REAL It_3 = 0.0;
  Matrix3 It = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  Matrix3 Jt = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  Matrix3 Bt = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  Matrix3 St = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  REAL rt = 0.0;

  // 2 indicates the first two segments of Humerus, for morphology
  for (unsigned int i = 0; i < 2; i++) {
    rt = (18 - 9.6 * (i * i) / 4.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * 15.0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    r0_v[i] = rt;
    r0_v[9 - i] = 0.8 * rt;
    J0_v[i] = Jt;
    J0_v[9 - i] = 0.64 * 0.64 * Jt;
    B0_v[i] = Bt;
    B0_v[8 - i] = 0.64 * 0.64 * Bt;
    S0_v[i] = St;
    S0_v[9 - i] = 0.64 * St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist = 0.0;
  const REAL nu = 1.0;  // 0.5
  const REAL relaxationNu = 0.0;
  const bool useSelfContact = false;

  Rod *rod0 = RodInitialConfigurations::straightRod_vscale(
      n, density / 2.0, r0_v, J0_v, B0_v, S0_v, L0, initialTotalTwist,
      humerusorigin, humerusdirection, humerusnormal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod0);
  rod0->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- ULNA
  // Driving parameters
  const int n_R = 12;
  const REAL L0_R = 180;

  const REAL Rdis = 8.4;

  // Define rod
  const Vector3 radiusdirection =
      Vector3(0.2, -1.0, 0.0).unitize();  //-0.8,-1.0,0.0
  const Vector3 radiusnormal = Vector3(0.0, 0.0, 1.0);
  // 10 is the last index of Humerus
  const Vector3 radiusorigin = rodPtrs[0]->x[10] + Vector3(Rdis, 0.0, 0.0);

  // Ulna shape
  // vector values
  vector<REAL> rR_v = vector<REAL>(n_R);  //(mm)
  vector<Matrix3> JR_v = vector<Matrix3>(n_R);
  vector<Matrix3> BR_v = vector<Matrix3>(n_R - 1);
  vector<Matrix3> SR_v = vector<Matrix3>(n_R);

  // middle shaft, index goes from 2 to 10 for providing morphology
  for (unsigned int i = 2; i < 10; i++) {
    rR_v[i] = 0.7 * r0;
    JR_v[i] = 0.49 * 0.49 * J0;
    if (i < 9) BR_v[i] = 0.49 * 0.49 * B0;
    SR_v[i] = 0.49 * S0;
  }

  // Heads, same iteration as Humerus above
  for (unsigned int i = 0; i < 2; i++) {
    rt = (9 - 3.12 * (i * i) / 4.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * 15.0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rR_v[i] = rt;
    rR_v[11 - i] = rt;
    JR_v[i] = Jt;
    JR_v[11 - i] = Jt;
    BR_v[i] = Bt;
    BR_v[10 - i] = Bt;
    SR_v[i] = St;
    SR_v[11 - i] = St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  Rod *rod1 = RodInitialConfigurations::straightRod_vscale(
      n_R, density / 2.0, rR_v, JR_v, BR_v, SR_v, L0_R, initialTotalTwist,
      radiusorigin, radiusdirection, radiusnormal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod1);
  rod1->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- RADIUS

  // Driving parameters
  const int n_U = 12;
  const REAL L0_U = 180;

  // Define rod
  const Vector3 ulnadirection = radiusdirection;
  const Vector3 ulnanormal = radiusnormal;
  // 10 is the last index of Humerus
  const Vector3 ulnaorigin = rodPtrs[0]->x[10] + Vector3(-1.1 * Rdis, 0.0, 0.0);

  // Radius shape
  // vector values
  vector<REAL> rU_v = vector<REAL>(n_U);  //(mm)
  vector<Matrix3> JU_v = vector<Matrix3>(n_U);
  vector<Matrix3> BU_v = vector<Matrix3>(n_U - 1);
  vector<Matrix3> SU_v = vector<Matrix3>(n_U);

  const REAL rU = 3.0;  // middle radius
  // Second moment of area for disk cross section
  const REAL AU = M_PI * rU * rU;
  const REAL IU_1 = AU * AU / (4.0 * M_PI);
  const REAL IU_2 = IU_1;
  const REAL IU_3 = 2.0 * IU_1;
  const Matrix3 IU = Matrix3(IU_1, 0.0, 0.0, 0.0, IU_2, 0.0, 0.0, 0.0, IU_3);
  // Mass inertia matrix for disk cross section
  const Matrix3 JU = density * 15.0 * IU;

  // Bending matrix
  Matrix3 BU =
      Matrix3(E * IU_1, 0.0, 0.0, 0.0, E * IU_2, 0.0, 0.0, 0.0, G * IU_3);
  // Shear matrix
  Matrix3 SU = Matrix3((4.0 / 3.0) * G * AU, 0.0, 0.0, 0.0,
                       (4.0 / 3.0) * G * AU, 0.0, 0.0, 0.0, E * AU);

  // middle shaft
  for (unsigned int i = 2; i < 10; i++) {
    rU_v[i] = rU;
    JU_v[i] = JU;
    if (i < 9) BU_v[i] = BU;
    SU_v[i] = SU;
  }

  for (unsigned int i = 0; i < 2; i++) {
    rt = (8 - 5 * (i * i) / 4.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * 15.0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rU_v[i] = rt;
    rU_v[11 - i] = rt;
    JU_v[i] = Jt;
    JU_v[11 - i] = Jt;
    BU_v[i] = Bt;
    BU_v[10 - i] = Bt;
    SU_v[i] = St;
    SU_v[11 - i] = St;
  }

  Rod *rod2 = RodInitialConfigurations::straightRod_vscale(
      n_U, density / 2.0, rU_v, JU_v, BU_v, SU_v, L0_U, initialTotalTwist,
      ulnaorigin, ulnadirection, ulnanormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod2);
  rod2->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- Carpometacarpus1

  const int n_C = 10;
  const REAL L0_C = 90;
  // Driving parameters

  // Define rod
  // const Vector3 Carp2direction = Vector3(0.1,-1.0,0.0).unitize();
  const Vector3 Carp2direction =
      (cos(0.7744) * rod2->Q[11][2] - sin(0.7744) * rod2->Q[11][1]).unitize();
  const Vector3 Carp2normal = Vector3(0.0, 0.0, 1.0);
  const Vector3 Carp2origin = rodPtrs[1]->x[11] + 5.0 * radiusdirection;

  // Carp1 shape
  // vector values
  vector<REAL> rC2_v = vector<REAL>(n_C);  //(mm)
  vector<Matrix3> JC2_v = vector<Matrix3>(n_C);
  vector<Matrix3> BC2_v = vector<Matrix3>(n_C - 1);
  vector<Matrix3> SC2_v = vector<Matrix3>(n_C);

  // middle shaft
  for (unsigned int i = 3; i < 7; i++) {
    rC2_v[i] = 0.3 * r0;
    JC2_v[i] = 0.09 * 0.09 * J0;
    if (i < 6) {
      /*
          Here, we harden the bending matrix B. The rationale comes from the
          following:

          If we set the bending stiffness according to the model radius
          above (0.3 * r0), we observe that the bone is prone to bending
          deformations, while in reality it is not (indeed bones are very
          rigid structures).
          This difference may be attributed to us abstracting away key details
         of the internal structure of the bone, which might be of importance. To
         maintain biologically realistic representation of the wing then, we
         harden the bending mode of the bone by a factor of 10 below. We tested
         other factors too, and found that the resulting dynamics remain
         essentially the same (<1% difference)
      */
      BC2_v[i] = 0.09 * B0;
    }
    SC2_v[i] = 0.09 * S0;
  }

  // Heads
  for (unsigned int i = 0; i < 3; i++) {
    rt = (5 - 2.48 * (i * i) / 9.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * 9.0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rC2_v[i] = rt;
    rC2_v[9 - i] = rt;
    JC2_v[i] = Jt;
    JC2_v[9 - i] = Jt;
    BC2_v[i] = Bt;
    BC2_v[8 - i] = Bt;
    SC2_v[i] = St;
    SC2_v[9 - i] = St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  Rod *rod4 = RodInitialConfigurations::straightRod_vscale(
      n_C, density / 2.0, rC2_v, JC2_v, BC2_v, SC2_v, L0_C, initialTotalTwist,
      Carp2origin, Carp2direction, Carp2normal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod4);
  rod4->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- Carpometacarpus2
  // Driving parameters

  // Define rod
  const Vector3 Carp1direction = Carp2direction;
  const Vector3 Carp1normal = Carp2normal;
  const Vector3 Carp1origin = Carp2origin + Vector3(-1.7 * Rdis, 0.0, 0.0);

  // Carp2 shape
  // vector values
  vector<REAL> rC_v = vector<REAL>(n_C);  //(mm)
  vector<Matrix3> JC_v = vector<Matrix3>(n_C);
  vector<Matrix3> BC_v = vector<Matrix3>(n_C - 1);
  vector<Matrix3> SC_v = vector<Matrix3>(n_C);

  // middle shaft
  for (unsigned int i = 3; i < 7; i++) {
    rC_v[i] = 0.7 * r0;
    JC_v[i] = 0.49 * 0.49 * J0;
    if (i < 6) BC_v[i] = 0.49 * 0.49 * B0;
    SC_v[i] = 0.49 * S0;
  }

  // Heads
  for (unsigned int i = 0; i < 3; i++) {
    rt = (10 - 4 * (i * i) / 9.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * 9.0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rC_v[i] = rt;
    rC_v[9 - i] = rt;
    JC_v[i] = Jt;
    JC_v[9 - i] = Jt;
    BC_v[i] = Bt;
    BC_v[8 - i] = Bt;
    SC_v[i] = St;
    SC_v[9 - i] = St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  Rod *rod3 = RodInitialConfigurations::straightRod_vscale(
      n_C, density / 2.0, rC_v, JC_v, BC_v, SC_v, L0_C, initialTotalTwist,
      Carp1origin, Carp1direction, Carp1normal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod3);
  rod3->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- Digits
  // Driving parameters
  const int n_S = 5;
  const REAL L0_S = 75;

  // Define rod
  const Vector3 SDdirection = Carp1direction;
  const Vector3 SDnormal = Carp1normal;
  const Vector3 SDorigin = rodPtrs[3]->x[9] + Vector3(-Rdis, 0.0, 0.0);

  // Digits shape
  // vector values
  vector<REAL> rS_v = vector<REAL>(n_S);  //(mm)
  vector<Matrix3> JS_v = vector<Matrix3>(n_S);
  vector<Matrix3> BS_v = vector<Matrix3>(n_S - 1);
  vector<Matrix3> SS_v = vector<Matrix3>(n_S);

  // Heads
  for (unsigned int i = 0; i < 5; i++) {
    rt = (12 - 2.5 * i);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * 15.0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rS_v[i] = rt;
    JS_v[i] = Jt;
    if (i < 4) BS_v[i] = Bt;
    SS_v[i] = St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  Rod *rod5 = RodInitialConfigurations::straightRod_vscale(
      n_S, density / 2.0, rS_v, JS_v, BS_v, SS_v, L0_S, initialTotalTwist,
      SDorigin, SDdirection, SDnormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod5);
  rod5->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Feathers-Primary Remiges : Rachis

  // Driving parameters
  const int n_F = 8;
  const REAL density_F = 0.0025e-3;
  // Scaled Length - same for all dimensions below
  // We note that applying the above scaling, this length is equivalent to
  // 18 cm as in a real pigeon
  const REAL L_F = 460.0;
  const REAL r_F = 2.5;
  const REAL totalMass_F = density_F * M_PI * r_F * r_F * L_F;

  // E_F and r_F are scaled to match bending stiffness according to rule 2.
  const REAL E_F = 2.5e6;
  const REAL G_F = E_F / (poissonRatio + 1.0);

  vector<Vector3> F1direction;
  vector<Vector3> F1normal;
  vector<Vector3> F1origin;

  // 10 represents number of feathers
  for (unsigned int i = 0; i < 10; i++) {
    const Vector3 directionF1 =
        Vector3((cos((i * 7.7) / 180 * M_PI) * SDdirection[0] -
                 sin((i * 7.7) / 180 * M_PI) * SDdirection[1]),
                (sin((i * 7.7) / 180 * M_PI) * SDdirection[0] +
                 cos((i * 7.7) / 180 * M_PI) * SDdirection[1]),
                0.0);
    F1direction.push_back(directionF1);

    const Vector3 normalF1 = Vector3(0.0, 0.0, 1.0);
    F1normal.push_back(normalF1);

    if (i < 5) {
      const Vector3 originF1 = rodPtrs[5]->x[5 - i];
      F1origin.push_back(originF1);
    } else {
      const Vector3 originF1 = rodPtrs[4]->x[10 - 2 * (i - 4)];
      F1origin.push_back(originF1);
    }
  }

  // Second moment of area for disk cross section for rachis
  const REAL A_F = M_PI * r_F * r_F;
  const REAL IF_1 = A_F * A_F / (4.0 * M_PI);
  const REAL IF_2 = IF_1;
  const REAL IF_3 = 2.0 * IF_1;
  const Matrix3 IF = Matrix3(IF_1, 0.0, 0.0, 0.0, IF_2, 0.0, 0.0, 0.0, IF_3);

  // Mass inertia matrix for disk cross section
  const Matrix3 J_F = density_F * (L_F / n_F) * IF;

  // Bending matrix
  const Matrix3 B_F =
      Matrix3(E_F * IF_1, 0.0, 0.0, 0.0, E_F * IF_2, 0.0, 0.0, 0.0, G_F * IF_3);

  // Shear matrix
  const Matrix3 S_F =
      Matrix3((4.0 / 3.0) * G_F * A_F, 0.0, 0.0, 0.0, (4.0 / 3.0) * G_F * A_F,
              0.0, 0.0, 0.0, E_F * A_F);

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist_F = 0.0;
  const REAL nu_F = 0.005;
  const REAL relaxationNu_F = 0.0;
  const bool useSelfContact_F = false;

  for (unsigned int i = 0; i < 10; i++) {
    Rod *rod6 = RodInitialConfigurations::straightRod(
        n_F, totalMass_F, r_F, J_F, B_F, S_F, L_F, initialTotalTwist_F,
        F1origin[i], F1direction[i], F1normal[i], nu_F, relaxationNu_F,
        useSelfContact_F);
    rodPtrs.push_back(rod6);
    rod6->update(0.0);
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Feathers-Primary Remiges : barbs

  // Driving parameters
  const int n_Fb = 4;
  const REAL density_Fb = 0.0025e-3;
  const REAL L_Fb = 60.0;
  const REAL r_Fb = 1.5;
  const REAL totalMass_Fb = density_Fb * M_PI * r_Fb * r_Fb * L_Fb;

  // E_Fb and r_Fb are scaled to match bending stiffness
  const REAL E_Fb = 2.5e6;
  const REAL G_Fb = E_Fb / (poissonRatio + 1.0);

  // Second moment of area for disk cross section for tendon
  const REAL A_Fb = M_PI * r_Fb * r_Fb;
  const REAL IFb_1 = A_Fb * A_Fb / (4.0 * M_PI);
  const REAL IFb_2 = IFb_1;
  const REAL IFb_3 = 2.0 * IFb_1;
  const Matrix3 IFb =
      Matrix3(IFb_1, 0.0, 0.0, 0.0, IFb_2, 0.0, 0.0, 0.0, IFb_3);

  // Mass inertia matrix for disk cross section
  const Matrix3 J_Fb = density_Fb * (L_Fb / n_Fb) * IFb;

  // Bending matrix
  const Matrix3 B_Fb = Matrix3(E_Fb * IFb_1, 0.0, 0.0, 0.0, E_Fb * IFb_2, 0.0,
                               0.0, 0.0, G_Fb * IFb_3);

  // Shear matrix
  const Matrix3 S_Fb =
      Matrix3((4.0 / 3.0) * G_Fb * A_Fb, 0.0, 0.0, 0.0,
              (4.0 / 3.0) * G_Fb * A_Fb, 0.0, 0.0, 0.0, E_Fb * A_Fb);

  // Initialize straight rod and pack it into a vector of pointers to rod

  const REAL nu_Fb = 0.005;

  for (unsigned int i = 0; i < 10; i++) {  // 10

    const Vector3 Fdirection = F1direction[i];
    const Vector3 Forigin = F1origin[i];

    for (unsigned int j = 0; j < 77; j++) {  // 77
      const Vector3 directionFb1 =
          (Fdirection + Vector3(-cos((55.6798 + i * 7.7) / 180 * M_PI),
                                -1 * sin((55.6798 + i * 7.7) / 180 * M_PI),
                                -0.08))
              .unitize();  // 5.710593137
      // cout<<directionFb1<<"direction"<<endl;
      const Vector3 directionFb1_2 =
          (Fdirection + Vector3(cos((55.6798 + i * 7.7) / 180 * M_PI),
                                1 * sin((55.6798 + i * 7.7) / 180 * M_PI),
                                0.08))
              .unitize();

      const Vector3 normalFb1 = Fdirection * directionFb1;
      // cout<<normalFb1<<"normal"<<endl;

      const Vector3 originFb1 = Forigin + Fdirection * (70.0 + j * 4.0);

      if (j < 10) {
        Rod *rod7 = RodInitialConfigurations::straightRod(
            n_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * totalMass_Fb,
            r_Fb, J_Fb, B_Fb, S_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * L_Fb,
            initialTotalTwist_F, originFb1, directionFb1, normalFb1, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod7);
        rod7->update(0.0);
        Rod *rod8 = RodInitialConfigurations::straightRod(
            n_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * totalMass_Fb,
            r_Fb, J_Fb, B_Fb, S_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * L_Fb,
            initialTotalTwist_F, originFb1, directionFb1_2, normalFb1, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod8);
        rod8->update(0.0);
      } else {
        Rod *rod7 = RodInitialConfigurations::straightRod(
            n_Fb, totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb, L_Fb,
            initialTotalTwist_F, originFb1, directionFb1, normalFb1, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod7);
        rod7->update(0.0);
        Rod *rod8 = RodInitialConfigurations::straightRod(
            n_Fb, totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb, L_Fb,
            initialTotalTwist_F, originFb1, directionFb1_2, normalFb1, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod8);
        rod8->update(0.0);
      }
    }

    for (unsigned int j = 0; j < 19; j++) {
      const Vector3 directionFb1 =
          (Fdirection +
           Vector3(-cos((55.6798 + i * 7.7 + 5 * j) / 180 * M_PI),
                   -1 * sin((55.6798 + i * 7.7 + 5 * j) / 180 * M_PI),
                   -0.08 * (19 - j) / 19))
              .unitize();
      const Vector3 directionFb1_2 =
          (Fdirection +
           Vector3(cos((55.6798 + i * 7.7 - 5 * j) / 180 * M_PI),
                   1 * sin((55.6798 + i * 7.7 - 5 * j) / 180 * M_PI),
                   0.08 * (19 - j) / 19))
              .unitize();

      const Vector3 normalFb1 = Fdirection * directionFb1;

      const Vector3 originFb1 = Forigin + Fdirection * 378.0;

      Rod *rod9 = RodInitialConfigurations::straightRod(
          n_Fb, (1 + j / 20.0) * totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb,
          (1 + j / 20.0) * L_Fb, initialTotalTwist_F, originFb1, directionFb1,
          normalFb1, nu_F, relaxationNu_F, useSelfContact_F);
      rodPtrs.push_back(rod9);
      rod9->update(0.0);

      const Vector3 normalFb2 = -Fdirection * directionFb1_2;

      Rod *rod10 = RodInitialConfigurations::straightRod(
          n_Fb, (1 + j / 20.0) * totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb,
          (1 + j / 20.0) * L_Fb, initialTotalTwist_F, originFb1, directionFb1_2,
          normalFb2, nu_F, relaxationNu_F, useSelfContact_F);
      rodPtrs.push_back(rod10);
      rod10->update(0.0);
    }
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Feathers-Secondary Remiges : Rachis

  // Driving parameters
  vector<Vector3> F2direction;
  vector<Vector3> F2normal;
  vector<Vector3> F2origin;

  for (unsigned int i = 0; i < 9; i++) {
    const Vector3 directionF2 =
        Vector3((cos((122 + i * 6.9) / 180 * M_PI) * ulnadirection[0] -
                 sin((122 + i * 6.9) / 180 * M_PI) * ulnadirection[1]),
                (sin((122 + i * 6.9) / 180 * M_PI) * ulnadirection[0] +
                 cos((122 + i * 6.9) / 180 * M_PI) * ulnadirection[1]),
                0.0);
    F2direction.push_back(directionF2);

    const Vector3 normalF2 = Vector3(0.0, 0.0, 1.0);
    F2normal.push_back(normalF2);

    if (i < 7) {
      const Vector3 originF2 = rodPtrs[2]->x[11 - i];
      F2origin.push_back(originF2);
    } else {
      const Vector3 originF2 = rodPtrs[2]->x[3 - 2 * (i - 7)];
      F2origin.push_back(originF2);
    }
  }

  for (unsigned int i = 0; i < 9; i++) {
    Rod *rod11 = RodInitialConfigurations::straightRod(
        n_F, totalMass_F, r_F, J_F, B_F, S_F, (1.0 - i / 15.0) * L_F,
        initialTotalTwist_F, F2origin[i], F2direction[i], F2normal[i], nu_F,
        relaxationNu_F, useSelfContact_F);
    rodPtrs.push_back(rod11);
    rod11->update(0.0);
  }
  //-----------------------------------------------------------------------------------------------------------------
  // Feathers-Secondary Remiges : Barbs

  for (unsigned int i = 0; i < 9; i++) {
    const Vector3 FFdirection = F2direction[i].unitize();
    const Vector3 FForigin = F2origin[i];
    const Vector3 Binormal = -Vector3(0.0, 0.0, 1.0) * FFdirection;

    Vector3 directionFb2 =
        (FFdirection + Vector3(Binormal[0], Binormal[1], -0.08)).unitize();
    Vector3 directionFb2_2 =
        (FFdirection + Vector3(-Binormal[0], -Binormal[1], 0.08)).unitize();

    Vector3 normalFb2 = (FFdirection * directionFb2).unitize();

    for (unsigned int j = 0; j < (77 - i * 7); j++) {
      const Vector3 originFb2 = FForigin + FFdirection * (70.0 + j * 4.0);

      if (j < 10) {
        Rod *rod12 = RodInitialConfigurations::straightRod(
            n_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * totalMass_Fb,
            r_Fb, J_Fb, B_Fb, S_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * L_Fb,
            initialTotalTwist_F, originFb2, directionFb2, normalFb2, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod12);
        rod12->update(0.0);
        Rod *rod13 = RodInitialConfigurations::straightRod(
            n_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * totalMass_Fb,
            r_Fb, J_Fb, B_Fb, S_Fb,
            (1 - (20 - j * 2 - 2) * (20 - j * 2 - 2) / 400.0) * L_Fb,
            initialTotalTwist_F, originFb2, directionFb2_2, normalFb2, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod13);
        rod13->update(0.0);
      } else {
        Rod *rod12 = RodInitialConfigurations::straightRod(
            n_Fb, totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb, L_Fb,
            initialTotalTwist_F, originFb2, directionFb2, normalFb2, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod12);
        rod12->update(0.0);
        Rod *rod13 = RodInitialConfigurations::straightRod(
            n_Fb, totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb, L_Fb,
            initialTotalTwist_F, originFb2, directionFb2_2, normalFb2, nu_Fb,
            relaxationNu_F, useSelfContact_F);
        rodPtrs.push_back(rod13);
        rod13->update(0.0);
      }
    }

    for (unsigned int j = 0; j < 19; j++) {
      Vector3 NewBi = Vector3((cos((j * 5.0) / 180 * M_PI) * Binormal[0] -
                               sin((j * 5.0) / 180 * M_PI) * Binormal[1]),
                              (sin((j * 5.0) / 180 * M_PI) * Binormal[0] +
                               cos((j * 5.0) / 180 * M_PI) * Binormal[1]),
                              0.0);
      const Vector3 directionFb2 =
          (FFdirection + Vector3(NewBi[0], NewBi[1], -0.08 * (19 - j) / 19))
              .unitize();
      NewBi = Vector3((cos((180 - j * 5.0) / 180 * M_PI) * Binormal[0] -
                       sin((180 - j * 5.0) / 180 * M_PI) * Binormal[1]),
                      (sin((180 - j * 5.0) / 180 * M_PI) * Binormal[0] +
                       cos((180 - j * 5.0) / 180 * M_PI) * Binormal[1]),
                      0.0);
      const Vector3 directionFb2_2 =
          (FFdirection + Vector3(NewBi[0], NewBi[1], 0.08 * (19 - j) / 19))
              .unitize();

      const Vector3 normalFb2 = FFdirection * directionFb2;

      const Vector3 originFb2 = FForigin + FFdirection * (378.0 - 28.0 * i);

      Rod *rod14 = RodInitialConfigurations::straightRod(
          n_Fb, totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb, (1 + j / 20.0) * L_Fb,
          initialTotalTwist_F, originFb2, directionFb2, normalFb2, nu_F,
          relaxationNu_F, useSelfContact_F);
      rodPtrs.push_back(rod14);
      rod14->update(0.0);

      Rod *rod15 = RodInitialConfigurations::straightRod(
          n_Fb, totalMass_Fb, r_Fb, J_Fb, B_Fb, S_Fb, (1 + j / 20.0) * L_Fb,
          initialTotalTwist_F, originFb2, directionFb2_2, normalFb2, nu_F,
          relaxationNu_F, useSelfContact_F);
      rodPtrs.push_back(rod15);
      rod15->update(0.0);
    }
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers --- Tendon and muscle in one rod
  // Driving parameters
  const int nm = 7;
  const REAL densitym = 1.06e-3;
  const REAL Em = 2e4;

  const REAL poissonRatiom = 0.5;  // Incompressible
  const REAL Gm = Em / (poissonRatiom + 1.0);
  const REAL SfactorM = sqrt(2.0);

  Vector3 originmuscle = Vector3(20.0, 20.0, 200.0);
  Vector3 directionmuscle = (rodPtrs[0]->x[8] - originmuscle).unitize();
  Vector3 Projectionmuscle =
      (rodPtrs[0]->x[8] - Vector3(originmuscle[0], originmuscle[1], 0.0))
          .unitize();
  Vector3 normalmuscle = Projectionmuscle * directionmuscle;

  REAL Lm = (rodPtrs[0]->x[8] - originmuscle).length();
  REAL dLm = Lm / (double)nm;  // length of cross-section element
  // Second moment of area for disk cross section for tendon
  const REAL rm_t = 5.0 / SfactorM;
  const REAL Am_t = M_PI * rm_t * rm_t;
  const REAL I0_1m_t = Am_t * Am_t / (4.0 * M_PI);
  const REAL I0_2m_t = I0_1m_t;
  const REAL I0_3m_t = 2.0 * I0_1m_t;
  const Matrix3 I0m_t =
      Matrix3(I0_1m_t, 0.0, 0.0, 0.0, I0_2m_t, 0.0, 0.0, 0.0, I0_3m_t);

  // Mass inertia matrix for disk cross section
  /*
  Here we simplify the wing structure by only including four single muscles. No
  connective tissues are considered, owing to lack of available biological data.
  Nevertheless, these tissues effectively constrain the lateral muscle motion.

  To account for these tissues then, we increase the inertia (J) and bending (B)
  matrices for the tendon and muscle elements below.
   */
  const Matrix3 Jm_t = 4.0 * densitym * dLm * I0m_t;

  // Bending matrix
  const Matrix3 Bm_t = 2.0 * Matrix3(Em * I0_1m_t, 0.0, 0.0, 0.0, Em * I0_2m_t,
                                     0.0, 0.0, 0.0, Gm * I0_3m_t);

  /*
    The tendons in case of the pigeon are short, flat and conforming to the
    joint, and only
    help to connect the muscles with the bones. Hence at leading order, we
    expect them to play a negligible role in the dynamics of the wing.

    Given these observations, we simplify our model by
    lumping the muscle with the tendon, which now occupies the first and last
    two elements of the muscular rod below. To mimic the tendon's ability to
    resist extension but undergo bending, we use a lower Young's  Modulus but
    harden the stiffness of the shear and stretch modes only, similar to the
    case of the MuscularSnake (see that case for more information).

  */
  // Shear matrix
  const Matrix3 Sm_t =
      50000 * Matrix3((4.0 / 3.0) * Gm * Am_t, 0.0, 0.0, 0.0,
                      (4.0 / 3.0) * Gm * Am_t, 0.0, 0.0, 0.0, Em * Am_t);

  // Second moment of area for disk cross section for muscle
  const REAL rm = 15.0 / SfactorM;
  const REAL Am = M_PI * rm * rm;
  const REAL I0_1m = Am * Am / (4.0 * M_PI);
  const REAL I0_2m = I0_1m;
  const REAL I0_3m = 2.0 * I0_1m;
  const Matrix3 I0m =
      Matrix3(I0_1m, 0.0, 0.0, 0.0, I0_2m, 0.0, 0.0, 0.0, I0_3m);
  // Mass inertia matrix for disk cross section
  // Here we use extra increment to the muscle J for the same reason.
  const Matrix3 Jm = 4.0 * 1.75e-3 * 15.0 * I0m;

  // Bending matrix
  const Matrix3 Bm = 2.0 * Matrix3(Em * I0_1m, 0.0, 0.0, 0.0, Em * I0_2m, 0.0,
                                   0.0, 0.0, Gm * I0_3m);
  // Shear matrix
  const Matrix3 Sm = Matrix3((4.0 / 3.0) * Gm * Am, 0.0, 0.0, 0.0,
                             (4.0 / 3.0) * Gm * Am, 0.0, 0.0, 0.0, Em * Am);

  vector<REAL> rm_v = vector<REAL>(nm);  //(mm)
  vector<Matrix3> Jm_v = vector<Matrix3>(nm);
  vector<Matrix3> Bm_v = vector<Matrix3>(nm - 1);
  vector<Matrix3> Sm_v = vector<Matrix3>(nm);

  // tendon
  for (unsigned int i = 0; i < 2; i++) {
    rm_v[i] = rm_t;
    Jm_v[i] = Jm_t;
    Bm_v[i] = Bm_t;
    Sm_v[i] = Sm_t;
  }
  Bm_v[1] = (Bm_t + Bm) / 2;

  for (unsigned int i = 2; i < 5; i++) {
    rm_v[i] = rm;
    Jm_v[i] = Jm;
    Bm_v[i] = Bm;
    Sm_v[i] = Sm;
  }
  Bm_v[5] = (Bm_t + Bm) / 2;
  // muscle SP-SP
  for (unsigned int i = 5; i < 7; i++) {
    rm_v[i] = rm_t;
    Jm_v[i] = Jm_t;
    if (i == 5) Bm_v[i] = Bm_t;
    Sm_v[i] = Sm_t;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist_m = 0.0;
  const REAL nu_m = 2;  // 0.1....0.5......1
  const REAL relaxationNu_m = 0.0;
  const bool useSelfContact_m = false;

  Rod *rod_m1 = RodInitialConfigurations::straightRod_vscale(
      nm, densitym, rm_v, Jm_v, Bm_v, Sm_v, Lm, initialTotalTwist_m,
      originmuscle, directionmuscle, normalmuscle, nu_m, relaxationNu_m,
      useSelfContact_m);
  rodPtrs.push_back(rod_m1);
  rod_m1->update(0.0);
  rod_m1->L0 = Lm;
  //-----------------------------------------------------------------------------------------------------------------
  // Muscle 2

  originmuscle = Vector3(-15.0, -15.0, 100.0);  // 30.0,0.0,100.0
  directionmuscle = (rodPtrs[0]->x[8] - originmuscle).unitize();
  Projectionmuscle =
      (rodPtrs[0]->x[8] - Vector3(originmuscle[0], originmuscle[1], 0.0))
          .unitize();
  normalmuscle = Projectionmuscle * directionmuscle;

  Lm = (rodPtrs[0]->x[8] - originmuscle).length();
  dLm = Lm / (double)nm;  // length of cross-section element

  Rod *rod_m2 = RodInitialConfigurations::straightRod_vscale(
      nm, densitym, rm_v, Jm_v, Bm_v, Sm_v, Lm, initialTotalTwist_m,
      originmuscle, directionmuscle, normalmuscle, nu_m, relaxationNu_m,
      useSelfContact_m);
  rodPtrs.push_back(rod_m2);
  rod_m2->update(0.0);
  rod_m2->L0 = Lm;
  //-----------------------------------------------------------------------------------------------------------------
  // Muscle 3

  originmuscle = Vector3(0.0, 0.0, 150.0);  //-30.0,-10.0,150.0
  directionmuscle = (rodPtrs[1]->x[3] - originmuscle).unitize();
  normalmuscle = Vector3(0.0, 0.0, 1.0);

  Lm = (rodPtrs[1]->x[3] - originmuscle).length();
  dLm = Lm / (double)nm;  // length of cross-section element

  Rod *rod_m3 = RodInitialConfigurations::straightRod_vscale(
      nm, densitym, rm_v, Jm_v, Bm_v, Sm_v, Lm, initialTotalTwist_m,
      originmuscle, directionmuscle, normalmuscle, nu_m, relaxationNu_m,
      useSelfContact_m);
  rodPtrs.push_back(rod_m3);
  rod_m3->update(0.0);
  rod_m3->L0 = Lm;
  //-----------------------------------------------------------------------------------------------------------------
  // Muscle 4

  originmuscle = Vector3(30.0, -10.0, 150.0);  // 30 -10 150
  Vector3 destination = rodPtrs[1]->x[0] - 15.0 * rodPtrs[1]->Q[0][2];
  directionmuscle = (destination - originmuscle).unitize();
  normalmuscle = Vector3(0.0, 0.0, 1.0);

  Lm = (destination - originmuscle).length();
  dLm = Lm / (double)nm;  // length of cross-section element

  Rod *rod_m4 = RodInitialConfigurations::straightRod_vscale(
      nm, densitym, rm_v, Jm_v, Bm_v, Sm_v, Lm, initialTotalTwist_m,
      originmuscle, directionmuscle, normalmuscle, nu_m, relaxationNu_m,
      useSelfContact_m);
  rodPtrs.push_back(rod_m4);
  rod_m4->update(0.0);
  rod_m4->L0 = Lm;

  //-----------------------------------------------------------------------------------------------------------------
  // Pack boundary conditions

  const int size = rodPtrs.size();
  vector<RodBC *> boundaryConditionsPtrs;
  // BC for Humerus
  HingeBC hingeBC0 = HingeBC(rodPtrs[0]);
  boundaryConditionsPtrs.push_back(&hingeBC0);
  // BCs for all other bones and feathers (primary and secondary)
  FreeBC freeBC = FreeBC();
  for (unsigned int i = 1; i < (size - 4); i++) {
    boundaryConditionsPtrs.push_back(&freeBC);
  }
  // Muscles
  HingeBC hingeBCM1 = HingeBC(rodPtrs[size - 4]);
  boundaryConditionsPtrs.push_back(&hingeBCM1);
  HingeBC hingeBCM2 = HingeBC(rodPtrs[size - 3]);
  boundaryConditionsPtrs.push_back(&hingeBCM2);
  HingeBC hingeBCM3 = HingeBC(rodPtrs[size - 2]);
  boundaryConditionsPtrs.push_back(&hingeBCM3);
  HingeBC hingeBCM4 = HingeBC(rodPtrs[size - 1]);
  boundaryConditionsPtrs.push_back(&hingeBCM4);
  // Pack all forces together (no forces applied)
  vector<ExternalForces *> externalForcesPtrs;

  // Gravity -- no gravity in this case, muscle force in ExternalContact
  MultipleForces multipleForces1;
  GravityForce gravity = GravityForce(Vector3(0.0, 0.0, 0.0));
  multipleForces1.add(&gravity);
  MultipleForces *multipleForcesPtr1 = multipleForces1.get();
  for (unsigned int i = 0; i < (size); i++) {
    externalForcesPtrs.push_back(multipleForcesPtr1);
  }

  // Not used, kept for legacy purposes
  vector<Interaction *> substrateInteractionsPtrs;

  // Set up External Contact
  vector<pair<int, int>> attachpoint;
  vector<ExternalContact *> externalcontactPtrs;
  /* The second and third argument are unimportant, but
         are preserved here for legacy purposes. Hence we simply
         set it to 0.0

   Externalcontact takes care of various connections and muscle
   actuations within the wing case.
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
  poly.setWindowStats(0.1, 0.5);

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
    fwdAvgVel.push_back(avgVel[i] % humerusdirection);
  }
  // const REAL latAvgVel = fabs(avgVel % Vector3(0.0,1.0,0.0));
  const vector<REAL> fitness = fwdAvgVel;

  return (fitness);
}

void Wing::run() {
  const vector<REAL> fitness = _wingRun();
  exit(0);
}
