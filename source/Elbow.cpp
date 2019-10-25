#include "Elbow.h"

Elbow::Elbow(const int argc, const char **argv) : amp(0.0), w(0.0), v(0.0) {}

// Units in this case are m/kg/s

vector<REAL> Elbow::_elbowRun() {
  vector<Rod *> rodPtrs;

  // Dumping frequencies (number of frames/dumps per unit time)
  const REAL diagPerUnitTime = 120;
  const REAL povrayPerUnitTime = 50;
  const REAL dt = 0.5e-7;  // artificial muscle + load
  const REAL timeSimulation = (4.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- HUMERUS
  // Driving parameters
  const int n = 20;
  const REAL density = 1750.0;  // 1.75g/cm^3
  const REAL L0 = 0.34;
  const REAL E = 15e9;
  const REAL g = 9.81;

  const REAL dL0 = L0 / (double)n;  // length of cross-section element
  const REAL poissonRatio = 0.5;    // Incompressible
  const REAL G = E / (poissonRatio + 1.0);

  // Define rod
  const Vector3 humerusdirection = Vector3(0.0, 0.0, -1.0);
  const Vector3 humerusnormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 humerusorigin = Vector3(0.0, 0.0, 1.0);

  // Humerus shape
  // vector values
  vector<REAL> r0_v = vector<REAL>(n);  //(mm)
  vector<Matrix3> J0_v = vector<Matrix3>(n);
  vector<Matrix3> B0_v = vector<Matrix3>(n - 1);
  vector<Matrix3> S0_v = vector<Matrix3>(n);

  const REAL r0 = 0.0105;  // middle radius 10.5mm
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

  // middle shaft
  for (unsigned int i = 2; i < 18; i++) {
    r0_v[i] = r0;
    J0_v[i] = J0;
    if (i < 17) B0_v[i] = B0;
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

  for (unsigned int i = 0; i < 2; i++) {
    rt = 0.0231 - 0.0126 * (i * i) / 4.0;
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * dL0 * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    r0_v[i] = rt;
    r0_v[19 - i] = rt;
    J0_v[i] = Jt;
    J0_v[19 - i] = Jt;
    B0_v[i] = Bt;
    B0_v[18 - i] = Bt;
    S0_v[i] = St;
    S0_v[19 - i] = St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist = 0.0;
  const REAL nu = 0.5;
  const REAL relaxationNu = 0.0;
  const bool useSelfContact = false;

  Rod *rod0 = RodInitialConfigurations::straightRod_v(
      n, density, r0_v, J0_v, B0_v, S0_v, L0, initialTotalTwist, humerusorigin,
      humerusdirection, humerusnormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod0);
  rod0->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- RADIUS
  // Driving parameters
  const int n_R = 15;
  const REAL L0_R = 0.255;
  const REAL dL0_R = L0_R / (double)n_R;  // length of cross-section element

  // Define rod
  const Vector3 radiusdirection = Vector3(0.0, 0.0, -1.0);
  const Vector3 radiusnormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 radiusorigin = Vector3(-r0, 0.0, 0.66);

  // Radius shape
  // vector values
  vector<REAL> rR_v = vector<REAL>(n_R);  //(mm)
  vector<Matrix3> JR_v = vector<Matrix3>(n_R);
  vector<Matrix3> BR_v = vector<Matrix3>(n_R - 1);
  vector<Matrix3> SR_v = vector<Matrix3>(n_R);

  // middle shaft
  for (unsigned int i = 2; i < 13; i++) {
    rR_v[i] = r0;
    JR_v[i] = J0;
    if (i < 12) BR_v[i] = B0;
    SR_v[i] = S0;
  }

  // Heads
  for (unsigned int i = 0; i < 2; i++) {
    rt = 0.015 - 0.0045 * (i * i) / 4.0;
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * dL0_R * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rR_v[i] = rt;
    rR_v[14 - i] = rt;
    JR_v[i] = Jt;
    JR_v[14 - i] = Jt;
    BR_v[i] = Bt;
    BR_v[13 - i] = Bt;
    SR_v[i] = St;
    SR_v[14 - i] = St;
  }

  // Initialize straight rod and pack it into a vector of pointers to rod
  Rod *rod1 = RodInitialConfigurations::straightRod_v(
      n_R, density, rR_v, JR_v, BR_v, SR_v, L0_R, initialTotalTwist,
      radiusorigin, radiusdirection, radiusnormal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod1);
  rod1->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Bones --- ULNA

  // Driving parameters
  const int n_U = 15;
  const REAL L0_U = 0.255;
  const REAL dL0_U = L0_U / (double)n_U;

  // Define rod
  const Vector3 ulnadirection = Vector3(0.0, 0.0, -1.0);
  const Vector3 ulnanormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 ulnaorigin = Vector3(1.2 * r0, 0.0, 0.66);

  // Ulna shape
  // vector values
  vector<REAL> rU_v = vector<REAL>(n_U);  //(mm)
  vector<Matrix3> JU_v = vector<Matrix3>(n_U);
  vector<Matrix3> BU_v = vector<Matrix3>(n_U - 1);
  vector<Matrix3> SU_v = vector<Matrix3>(n_U);

  const REAL rU = 0.0055;  // middle radius 10.5mm
  // Second moment of area for disk cross section
  const REAL AU = M_PI * rU * rU;
  const REAL IU_1 = AU * AU / (4.0 * M_PI);
  const REAL IU_2 = IU_1;
  const REAL IU_3 = 2.0 * IU_1;
  const Matrix3 IU = Matrix3(IU_1, 0.0, 0.0, 0.0, IU_2, 0.0, 0.0, 0.0, IU_3);
  // Mass inertia matrix for disk cross section
  const Matrix3 JU = density * dL0 * IU;

  // Bending matrix
  Matrix3 BU =
      Matrix3(E * IU_1, 0.0, 0.0, 0.0, E * IU_2, 0.0, 0.0, 0.0, G * IU_3);
  // Shear matrix
  Matrix3 SU = Matrix3((4.0 / 3.0) * G * AU, 0.0, 0.0, 0.0,
                       (4.0 / 3.0) * G * AU, 0.0, 0.0, 0.0, E * AU);

  // middle shaft
  for (unsigned int i = 2; i < 13; i++) {
    rU_v[i] = rU;
    JU_v[i] = JU;
    if (i < 12) BU_v[i] = BU;
    SU_v[i] = SU;
  }

  for (unsigned int i = 0; i < 2; i++) {
    rt = 0.01 - 0.0045 * (i * i) / 4.0;
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density * dL0_U * It;
    Bt = Matrix3(E * It_1, 0.0, 0.0, 0.0, E * It_2, 0.0, 0.0, 0.0, G * It_3);
    St = Matrix3((4.0 / 3.0) * G * At, 0.0, 0.0, 0.0, (4.0 / 3.0) * G * At, 0.0,
                 0.0, 0.0, E * At);
    rU_v[i] = rt;
    rU_v[14 - i] = rt;
    JU_v[i] = Jt;
    JU_v[14 - i] = Jt;
    BU_v[i] = Bt;
    BU_v[13 - i] = Bt;
    SU_v[i] = St;
    SU_v[14 - i] = St;
  }

  Rod *rod2 = RodInitialConfigurations::straightRod_v(
      n_U, density, rU_v, JU_v, BU_v, SU_v, L0_U, initialTotalTwist, ulnaorigin,
      ulnadirection, ulnanormal, nu, relaxationNu, useSelfContact);
  rodPtrs.push_back(rod2);
  rod2->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Tendons --- Long Head
  // Driving parameters
  const int n_LT = 7;
  const REAL density_T = 1670.0;  // 1.67g/cm^3
  const REAL L0_LT = 0.017 * 7;
  const REAL dL0_LT = L0_LT / (double)n_LT;

  // const REAL totalMass = density*M_PI*r0*r0*L0;
  const REAL E_T = 0.5e9;
  const REAL G_T = E_T / (poissonRatio + 1.0);

  // Define rod
  const Vector3 LongTdirection = Vector3(0.0, 0.0, -1.0);
  const Vector3 LongTnormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 LongTorigin = Vector3(-r0, -3 * r0, 1.0);

  // Long Head Tendon shape
  // vector values
  vector<REAL> rLT_v = vector<REAL>(n_LT);  //(mm)
  vector<Matrix3> JLT_v = vector<Matrix3>(n_LT);
  vector<Matrix3> BLT_v = vector<Matrix3>(n_LT - 1);
  vector<Matrix3> SLT_v = vector<Matrix3>(n_LT);

  for (unsigned int i = 0; i < n_LT; i++) {
    rt = 0.005 + (n_LT - 1) * 0.0029 * (i * i / 36.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density_T * dL0_LT * It;
    Bt = Matrix3(E_T * It_1, 0.0, 0.0, 0.0, E_T * It_2, 0.0, 0.0, 0.0,
                 G_T * It_3);
    St = Matrix3(1e1 * (4.0 / 3.0) * G_T * At, 0.0, 0.0, 0.0,
                 1e1 * (4.0 / 3.0) * G_T * At, 0.0, 0.0, 0.0, E_T * At);
    rLT_v[i] = rt;
    JLT_v[i] = Jt;
    if (i < (n_LT - 1)) BLT_v[i] = Bt;
    SLT_v[i] = St;
  }

  Rod *rod3 = RodInitialConfigurations::straightRod_v(
      n_LT, density_T, rLT_v, JLT_v, BLT_v, SLT_v, L0_LT, initialTotalTwist,
      LongTorigin, LongTdirection, LongTnormal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod3);
  rod3->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Tendons --- Short Head
  // Driving parameters
  const int n_ST = 8;
  const REAL L0_ST = 0.017 * 8;
  const REAL dL0_ST = L0_ST / (double)n_ST;

  // Define rod
  const Vector3 ShortTdirection =
      Vector3(-sin(12.5 * M_PI / 180), 0.0, -cos(12.5 * M_PI / 180));
  const Vector3 ShortTnormal = Vector3(0.0, 1.0, 0.0);
  const Vector3 ShortTorigin = Vector3(6 * r0, -3 * r0, 1.0);

  // Long Head Tendon shape
  // vector values
  vector<REAL> rST_v = vector<REAL>(n_ST);
  vector<Matrix3> JST_v = vector<Matrix3>(n_ST);
  vector<Matrix3> BST_v = vector<Matrix3>(n_ST - 1);
  vector<Matrix3> SST_v = vector<Matrix3>(n_ST);

  for (unsigned int i = 0; i < n_ST; i++) {
    rt = 0.005 + i * 0.0024;
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density_T * dL0_ST * It;
    Bt = Matrix3(E_T * It_1, 0.0, 0.0, 0.0, E_T * It_2, 0.0, 0.0, 0.0,
                 G_T * It_3);
    St = Matrix3(1e1 * (4.0 / 3.0) * G_T * At, 0.0, 0.0, 0.0,
                 1e1 * (4.0 / 3.0) * G_T * At, 0.0, 0.0, 0.0, E_T * At);
    rST_v[i] = rt;
    JST_v[i] = Jt;
    if (i < (n_ST - 1)) BST_v[i] = Bt;
    SST_v[i] = St;
  }

  Rod *rod4 = RodInitialConfigurations::straightRod_v(
      n_ST, density_T, rST_v, JST_v, BST_v, SST_v, L0_ST, initialTotalTwist,
      ShortTorigin, ShortTdirection, ShortTnormal, nu, relaxationNu,
      useSelfContact);
  rodPtrs.push_back(rod4);
  rod4->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Tendons --- Long Head End
  // Driving parameters
  const int n_LET = 7;
  const REAL L0_LET = 0.017 * 7;
  const REAL dL0_LET = L0_LET / (double)n_LET;

  // Define rod
  const Vector3 LongETdirection =
      Vector3(0.0, -sin(17.99 * M_PI / 180), cos(17.99 * M_PI / 180));
  const Vector3 LongETnormal = Vector3(1.0, 0.0, 0.0);
  const Vector3 LongETorigin = Vector3(-r0, 0.0, 0.592);

  // Long Head Tendon shape
  // vector values
  vector<REAL> rLET_v = vector<REAL>(n_LET);  //(mm)
  vector<Matrix3> JLET_v = vector<Matrix3>(n_LET);
  vector<Matrix3> BLET_v = vector<Matrix3>(n_LET - 1);
  vector<Matrix3> SLET_v = vector<Matrix3>(n_LET);

  for (unsigned int i = 0; i < n_LET; i++) {
    rt = 0.005 + (n_LT - 1) * 0.0031 * (i * i / 36.0);
    At = M_PI * rt * rt;
    It_1 = At * At / (4.0 * M_PI);
    It_2 = It_1;
    It_3 = 2.0 * It_1;
    It = Matrix3(It_1, 0.0, 0.0, 0.0, It_2, 0.0, 0.0, 0.0, It_3);

    Jt = density_T * dL0_LET * It;
    Bt = Matrix3(E_T * It_1, 0.0, 0.0, 0.0, E_T * It_2, 0.0, 0.0, 0.0,
                 G_T * It_3);
    St = Matrix3(1e1 * (4.0 / 3.0) * G_T * At, 0.0, 0.0, 0.0,
                 1e1 * (4.0 / 3.0) * G_T * At, 0.0, 0.0, 0.0, E_T * At);
    rLET_v[i] = rt;
    JLET_v[i] = Jt;
    if (i < (n_LET - 1)) BLET_v[i] = Bt;
    SLET_v[i] = St;
  }

  Rod *rod5 = RodInitialConfigurations::straightRod_v(
      n_LET, density_T, rLET_v, JLET_v, BLET_v, SLET_v, L0_LET,
      initialTotalTwist, LongETorigin, LongETdirection, LongETnormal, nu,
      relaxationNu, useSelfContact);
  rodPtrs.push_back(rod5);
  rod5->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Define Muscle Fibers

  // number of muscle fibers
  NOR = 18;

  // Define Muscle
  // Driving parameters
  const int n_m = 14;
  const REAL density_m = 1060.0;
  const REAL L_m = 0.017 * 14;
  const REAL r_m = 0.004;
  const REAL totalMass_m = density_m * M_PI * r_m * r_m * L_m;
  const REAL dL0_m = L_m / (double)n_m;

  const REAL E_m = 1e4;
  const REAL G_m = E_m / (poissonRatio + 1.0);

  // Long head
  const Vector3 muscledirection = Vector3(0.0, 0.0, -1.0);
  const Vector3 musclenormal = Vector3(1.0, 0.0, 0.0);
  vector<Vector3> muscleorigin;

  // Arrangement of muscle fibers
  const Vector3 center = LongTorigin + 5 * dL0 * Vector3(0.0, 0.0, -1.0);
  for (unsigned int i = 0; i < NOR; i++) {
    if (i < 6) {
      const Vector3 origin =
          center + 0.0082 * Vector3(cos((15 + i * 60) * M_PI / 180),
                                    sin((15 + i * 60) * M_PI / 180), 0.0);
      muscleorigin.push_back(origin);
    } else {
      const Vector3 origin =
          center + 0.016 * Vector3(cos(((i - 6) * 30) * M_PI / 180),
                                   sin(((i - 6) * 30) * M_PI / 180), 0.0);
      muscleorigin.push_back(origin);
    }
  }

  // Short head
  const int n_ms = 14;
  const REAL L_ms = 0.017 * 14;
  const REAL totalMass_ms = density_m * M_PI * r_m * r_m * L_ms;

  const Vector3 muscledirection_s = ShortTdirection;
  const Vector3 musclenormal_s = Vector3(0.0, 1.0, 0.0);
  vector<Vector3> muscleorigin_s;

  // Arrangement of muscle fibers
  const Vector3 center_s = rodPtrs[4]->x[6];
  for (unsigned int i = 0; i < NOR; i++) {
    if (i < 6) {
      const Vector3 origin_s =
          center_s + 0.0083 * Vector3(cos((15 + i * 60) * M_PI / 180),
                                      sin((15 + i * 60) * M_PI / 180), 0.0);
      muscleorigin_s.push_back(origin_s);
    } else {
      const Vector3 origin_s =
          center_s + 0.017 * Vector3(cos(((i - 6) * 30) * M_PI / 180),
                                     sin(((i - 6) * 30) * M_PI / 180), 0.0);
      muscleorigin_s.push_back(origin_s);
    }
  }

  // Second moment of area for disk cross section for muscles
  const REAL A_m = M_PI * r_m * r_m;
  const REAL Im_1 = A_m * A_m / (4.0 * M_PI);
  const REAL Im_2 = Im_1;
  const REAL Im_3 = 2.0 * Im_1;
  const Matrix3 Im = Matrix3(Im_1, 0.0, 0.0, 0.0, Im_2, 0.0, 0.0, 0.0, Im_3);

  // Mass inertia matrix for disk cross section
  const Matrix3 J_m = density_m * dL0_m * Im;

  // Bending matrix
  const Matrix3 B_m =
      Matrix3(E_m * Im_1, 0.0, 0.0, 0.0, E_m * Im_2, 0.0, 0.0, 0.0, G_m * Im_3);

  // Shear matrix
  const Matrix3 S_m =
      Matrix3(1e4 * (4.0 / 3.0) * G_m * A_m, 0.0, 0.0, 0.0,
              1e4 * (4.0 / 3.0) * G_m * A_m, 0.0, 0.0, 0.0, E_m * A_m);

  // Initialize straight rod and pack it into a vector of pointers to rod
  const REAL initialTotalTwist_m = 0.0;
  const REAL nu_m = 0.0;
  const REAL relaxationNu_m = 0.0;
  const bool useSelfContact_m = false;

  for (unsigned int i = 0; i < NOR; i++) {
    Rod *rod6 = RodInitialConfigurations::straightRod(
        n_m, totalMass_m, r_m, J_m, B_m, S_m, L_m, initialTotalTwist_m,
        muscleorigin[i], muscledirection, musclenormal, nu_m, relaxationNu_m,
        useSelfContact_m);
    rodPtrs.push_back(rod6);
    rod6->update(0.0);
  }

  for (unsigned int i = 0; i < NOR; i++) {
    Rod *rod7 = RodInitialConfigurations::straightRod(
        n_ms, totalMass_ms, r_m, J_m, B_m, S_m, L_ms, initialTotalTwist_m,
        muscleorigin_s[i], muscledirection_s, musclenormal_s, nu_m,
        relaxationNu_m, useSelfContact_m);
    rodPtrs.push_back(rod7);
    rod7->update(0.0);
  }

  //-----------------------------------------------------------------------------------------------------------------
  // Artificial Muscle --- Helical, twisted filament

  const REAL density_am = 1670.0;  // 1.67g/cm^3
  const REAL r0_am = 0.005;
  const REAL E_am = 3e10;
  const REAL poissonRatio_am = 0.5;
  const REAL G_am = E_am / (poissonRatio_am + 1.0);

  const REAL initialTotalTwist_am = 0.0;
  const REAL nu_am = 5;
  const REAL relaxationNu_am = 0.0;
  const bool useSelfContact_am = false;

  // Total length of rod
  // = 2*pi*n_turns*sqrt(pitch_factor^2 + helix_radius^2)
  const REAL helix_radius = 0.01;
  const REAL n_turns = 9.0;
  // Each turn needs to have at least 10 points say
  const int n_am = (int)(12 * n_turns);
  // centerline z is 1.0
  const REAL pitch_factor = 0.425 / 2.0 / M_PI / n_turns;
  // const REAL pitch_factor =
  //     sqrt(pow(1.0 / 2.0 / M_PI / n_turns, 2) - pow(helix_radius, 2));

  // Direction and origin
  REAL dir_s = 0.5 * (2 * M_PI * n_turns) / n_am;
  Vector3 origin_am = Vector3(-r0, -0.09, 0.575);
  Vector3 direction_am = Vector3(-cos(dir_s), -sin(dir_s), 0.0);
  // REAL d3_prefac =
  //     1. / sqrt(helix_radius * helix_radius + pitch_factor * pitch_factor);
  // Vector3 direction = (-Vector3(-cos(dir_s), -sin(dir_s), 0.0) * d3_prefac *
  //                      Vector3(-helix_radius * sin(dir_s),
  //                              helix_radius * cos(dir_s), pitch_factor))
  //                         .unitize();
  // Vector3 direction = Vector3(0.0, 0.0, 1.0);

  Rod *rod8 = RodInitialConfigurations::helicalRod(
      n_am, density_am, r0_am, pitch_factor, helix_radius, n_turns, E_am, G_am,
      initialTotalTwist_am, origin_am, direction_am, nu_am, relaxationNu_am,
      useSelfContact_am);

  rodPtrs.push_back(rod8);
  rod8->update(0.0);

  //-----------------------------------------------------------------------------------------------------------------
  // Pack boundary conditions
  vector<RodBC *> boundaryConditionsPtrs;
  // Humerus
  /*
    The upper arm is fixed, as usually done in a standard laboratory isokinetic
    test
  */
  FixedBC fixed = FixedBC(rodPtrs[0]);
  boundaryConditionsPtrs.push_back(&fixed);
  // Radius & Ulna
  FreeBC freeBC = FreeBC();
  boundaryConditionsPtrs.push_back(&freeBC);  // Radius
  boundaryConditionsPtrs.push_back(&freeBC);  // Ulna

  // Tendons up
  HingeBC hingeBC3 = HingeBC(rodPtrs[3]);
  boundaryConditionsPtrs.push_back(&hingeBC3);  // Tendon 1
  HingeBC hingeBC4 = HingeBC(rodPtrs[4]);
  boundaryConditionsPtrs.push_back(&hingeBC4);  // Tendon 2
  // Tendon bottom
  boundaryConditionsPtrs.push_back(&freeBC);
  // muscles 6-23
  for (unsigned int i = 0; i < 2 * NOR; i++) {
    boundaryConditionsPtrs.push_back(&freeBC);
  }

  // Boundary condition for the artificial muscle
  SolenoidsBCTeja tejaBC = SolenoidsBCTeja(rodPtrs[42], 0.1, 0.1);
  boundaryConditionsPtrs.push_back(&tejaBC);

  // Pack all forces together (no forces applied)
  vector<ExternalForces *> externalForcesPtrs;

  // Gravity
  MultipleForces multipleForces1;
  GravityForce gravity = GravityForce(Vector3(0.0, 0.0, -g));
  multipleForces1.add(&gravity);
  MultipleForces *multipleForcesPtr1 = multipleForces1.get();
  // Only for three bones and three tendons
  for (unsigned int i = 0; i < 6; i++) {
    externalForcesPtrs.push_back(multipleForcesPtr1);
  }
  // Muscle forces
  MultipleForces multipleForces2;
  GravityForce gravity1 = GravityForce(Vector3(
      0.0, 0.0, 0.0));  // No gravity on muscle fibers, since we do not model
                        // the various tissues wrapping around the muscle to
                        // hold the shape and position of the muscle.

  /*
  50% strength corresponds to a value of a peak force of 21.69N
  If you want 100% strength, please set peak force to 43.38N

 The second argument below is kept for legacy purposes, but is not used
 Hence we set it to 0.
 */
  MuscleContraction muscleC = MuscleContraction(21.69, 0.0);
  multipleForces2.add(&gravity1);
  multipleForces2.add(&muscleC);
  MultipleForces *multipleForcesPtr2 = multipleForces2.get();
  // For all the muscle fibers, undergo same actuation
  for (unsigned int i = 0; i < 2 * NOR; i++) {
    externalForcesPtrs.push_back(multipleForcesPtr2);
  }
  MultipleForces multipleForces3;
  multipleForces3.add(&gravity1);
  MultipleForces *multipleForcesPtr3 = multipleForces3.get();
  externalForcesPtrs.push_back(multipleForcesPtr3);

  // No substrate in this case, so the substrateInteractionsPtrs is
  // empty, for legacy purposes.
  vector<Interaction *> substrateInteractionsPtrs;

  // Set up External Contact
  vector<pair<int, int>> attachpoint;
  vector<ExternalContact *> externalcontactPtrs;
  /* The second and third argument are unimportant, but
         are preserved here for legacy purposes. Hence we simply
         set it to 0.0

    This function includes the assembly of the elbow and the
    contact between muscle fibers
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
  // case we don't collect any data, but it is dumped in the velocity.txt
  // file below.
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
    fwdAvgVel.push_back(avgVel[i] % humerusdirection);
  }
  const vector<REAL> fitness = fwdAvgVel;

  return (fitness);
}

void Elbow::run() {
  const vector<REAL> fitness = _elbowRun();
  exit(0);
}
