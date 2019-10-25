/*
Detect Rod-Rod collision
Rod-Rod connection
Detect Rod-Object collision
*/
#ifndef EXTERNALCONTACT_H_
#define EXTERNALCONTACT_H_

#include "GeometryFunctions.h"
#include "Matrix3.h"
#include "Rod.h"
#include "UsualHeaders.h"
#include "Vector3.h"

using namespace std;

class ExternalContact {
  // Not all variables are used
 protected:
  typedef std::vector<Rod *> Vrodptr;
  Vrodptr rodptrs;
  const REAL isotropicFCK;
  const REAL isotropicFCS;
  vector<pair<int, int>> attachpoint;
  REAL angle_old = 0.0;
  REAL angle_new = 0.0;
  REAL angleVelocity = 0.0;
  REAL ml_old = 0.0;
  REAL ml_new = 0.0;
  REAL mVelocity = 0.0;
  REAL angleVelocity_old = 0.0;
  Vector3 externalforce = Vector3(0.0, 0.0, 0.0);
  REAL Localtime = -0.1;
  int SwitchFlag = 0;
  int timecountL = 0;
  int timecountT = 0;

  vector<int> Point_Rod;
  vector<int> Point_Rod_Direction;
  int nor;

  // Used for receiving spikes
  vector<REAL> Spike_time;

  // for the oscillators theta_i
  REAL theta[6][1000];
  REAL omega[6][1000];
  REAL radius[6];
  REAL Angleold[3];
  REAL Lengthold[3] = {0.4, 0.4, 0.34};
  // snake trajectory
  REAL Lift_ratio = 0.0;

  const REAL _linear(const REAL _vel, const REAL _velThresold) const {
    const REAL vel = fabs(_vel);
    const REAL velThresold = fabs(_velThresold);

    const REAL width = 0.5 * velThresold;
    const REAL velDiff = vel - velThresold;
    REAL f = 1.0;

    if (vel > (velThresold)) f = fabs(1.0 - min(1.0, velDiff / width));

    if (vel > (velThresold + width)) f = 0.0;

    return f;
  }

 public:
  ExternalContact(Vrodptr rodptrs, const REAL _isotropicFCK,
                  const REAL _isotropicFCS, vector<pair<int, int>> _attachpoint)
      : rodptrs(rodptrs),
        isotropicFCK(_isotropicFCK),
        isotropicFCS(_isotropicFCS),
        attachpoint(_attachpoint) {}
  ~ExternalContact(){};

  // This function implements the collision model in walker case. Contact check
  // and repulsive force are performed here between pillars and the muscle
  // tissue.
  void RodRodCollision_Simple(const REAL time) {
    REAL Forcesum = 0.0;

    if (time > 0.1) {
      const REAL Kpush = 1e8;
      const REAL Nupush = 0.001;

      const REAL distance = 2.04;

      // Left and Right
      for (unsigned int k = 0; k < 2; k++) {
        // muscles on each side
        for (unsigned int j = 10; j < 12; j++) {
          Rod *muscle = rodptrs[j + 2 * k];

          Vector3 target =
              (k == 0)
                  ? ((rodptrs[3]->x[3 + (j - 10)] - rodptrs[5]->x[3 + (j - 10)])
                         .unitize())
                  : ((rodptrs[2]->x[3 + (j - 10)] - rodptrs[4]->x[3 + (j - 10)])
                         .unitize());
          Vector3 base = (k == 0) ? (rodptrs[5]->x[3 + (j - 10)])
                                  : (rodptrs[4]->x[3 + (j - 10)]);
          Vector3 Bodyvelocity = (k == 0)
                                     ? (0.5 * (rodptrs[3]->v[3 + (j - 10)] +
                                               rodptrs[5]->v[3 + (j - 10)]))
                                     : (0.5 * (rodptrs[2]->v[3 + (j - 10)] +
                                               rodptrs[4]->v[3 + (j - 10)]));

          const Vector3 forcedirection =
              (rodptrs[3]->x[3 + (j - 10)] - rodptrs[2]->x[3 + (j - 10)])
                  .unitize();
          const int bodyindex = (k == 0) ? (3) : (2);

          // Check edge one -- front and back
          for (unsigned int i = 0; i < 8; i++) {
            const Vector3 position = muscle->x[i] - base;
            const Vector3 direction = position - (position % target) * target;

            if ((direction.length() - distance) < 0) {
              // Bouncing Force
              const REAL pushing =
                  Kpush *
                  (distance - direction.length());  // should be positive

              // Compute damping force
              const Vector3 vInterpenetration = Bodyvelocity - muscle->v[i];
              const REAL vNorm = vInterpenetration % (direction.unitize());
              const REAL damping = Nupush * vNorm;

              // Compute collision force (the 0.5 in fron is because it is
              // distributed onto the two masses at the end of an edge)
              const Vector3 cForce =
                  (pushing + damping) * (direction.unitize());

              rodptrs[bodyindex]->externalForces[3 + (j - 10)] +=
                  ((-0.5 * cForce) % forcedirection) * forcedirection;
              rodptrs[bodyindex + 2]->externalForces[3 + (j - 10)] +=
                  ((-0.5 * cForce) % forcedirection) * forcedirection;
              if (k == 0) {
                Forcesum += -cForce % forcedirection;
              }
              muscle->externalForces[i] += Vector3(cForce[0], cForce[1], 0.0);
            }  // if
          }

          // Check edge two -- side
          target = (rodptrs[5]->x[3 + (j - 10)] - rodptrs[4]->x[3 + (j - 10)])
                       .unitize();
          base = rodptrs[4]->x[3 + (j - 10)];

          // elements on each muscle------Shortside
          for (unsigned int i = 8; i < 12; i++) {
            const Vector3 position = muscle->x[i] - base;
            const Vector3 direction = position - (position % target) * target;

            if ((direction.length() - distance) < 0) {
              // Bouncing Force
              const REAL pushing =
                  Kpush *
                  (distance - direction.length());  // should be positive

              // Compute damping force
              const Vector3 vInterpenetration = -muscle->v[i];
              const REAL vNorm = vInterpenetration % (direction.unitize());
              const REAL damping = Nupush * vNorm;

              // Compute collision force (the 0.5 in fron is because it is
              // distributed onto the two masses at the end of an edge)
              const Vector3 cForce =
                  (pushing + damping) * (direction.unitize());

              muscle->externalForces[i] += Vector3(cForce[0], cForce[1], 0.0);
            }  // if
          }

          // Check edge three -- side
          target = (rodptrs[3]->x[3 + (j - 10)] - rodptrs[2]->x[3 + (j - 10)])
                       .unitize();
          base = rodptrs[2]->x[3 + (j - 10)];

          // elements on each muscle------Shortside
          for (unsigned int i = 18; i < 22; i++) {
            const Vector3 position = muscle->x[i] - base;
            const Vector3 direction = position - (position % target) * target;

            if ((direction.length() - distance) < 0) {
              // Bouncing Force
              const REAL pushing =
                  Kpush *
                  (distance - direction.length());  // should be positive

              // Compute damping force
              const Vector3 vInterpenetration = -muscle->v[i];
              const REAL vNorm = vInterpenetration % (direction.unitize());
              const REAL damping = Nupush * vNorm;

              // Compute collision force (the 0.5 in fron is because it is
              // distributed onto the two masses at the end of an edge)
              const Vector3 cForce =
                  (pushing + damping) * (direction.unitize());

              muscle->externalForces[i] += Vector3(cForce[0], cForce[1], 0.0);
            }  // if
          }
        }
      }
      REAL SumL = 0.0;
      for (unsigned int i = 0; i < 22; i++) {
        SumL += rodptrs[10]->l[i];
      }
      // setup for data output---optional
      rodptrs[0]->Tforce = Forcesum;
      rodptrs[0]->MaxHeight = SumL;
    }
  }

//****************________________________________****************________________________________
//****************———————————————————————————————-****************________________________________
// Elbow-Elbow-Elbow-Elbow-Elbow
#ifdef FLAGELBOW
  // This function takes care of the various connections and muscle fibers
  // collision checks within the elbow case.
  void RodRodAttach(const REAL time) {
    //-----------------------------------------------------------------------------------------------------------------
    // Connections
    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);
    Vector3 check3 = Vector3(0.0, 0.0, 0.0);
    Vector3 check4 = Vector3(0.0, 0.0, 0.0);

    const REAL k = 5e5;

    Rod *Humerus = rodptrs[0];
    Rod *Radius = rodptrs[1];
    Rod *Ulna = rodptrs[2];
    Rod *Long = rodptrs[3];
    Rod *Short = rodptrs[4];
    Rod *LongEnd = rodptrs[5];

    const int nom = 18;

    // Force1 —- Humerus—Radius

    /*
        Here, we utilize a simpler version of the Hinge Joint (in
       HingeJoint.cpp).

       Here the radius moves only in the y–z plane. Instead of using restoring
       torques to correct out-of-plane motions (in the x direction), we  simply
       project all forces due to the connection onto the y–z plane.
    */
    const int size_H = Humerus->n;
    check1 = Humerus->x[size_H] + 0.017 * (Humerus->edge[size_H - 1]).unitize();
    check2 = Radius->x[1];

    Vector3 force_H_R_y = k * ((check2 - check1) % Vector3(0.0, 1.0, 0.0)) *
                          Vector3(0.0, 1.0, 0.0);
    Humerus->externalForces[size_H] += force_H_R_y;
    Radius->externalForces[1] -= force_H_R_y;

    Vector3 force_H_R_z = k * ((check2 - check1) % Vector3(0.0, 0.0, 1.0)) *
                          Vector3(0.0, 0.0, 1.0);
    Humerus->externalForces[size_H] += force_H_R_z;
    Radius->externalForces[1] -= force_H_R_z;

    // Force2 —- Radius—Ulna
    const int size_UR = Radius->n;
    check1 = Radius->x[0];
    check2 = Ulna->x[0];

    Vector3 force_U_R_1 = k * ((check2 - Vector3(check2[0], 0.0, 0.0)) -
                               (check1 - Vector3(check1[0], 0.0, 0.0)));
    Radius->externalForces[0] += force_U_R_1;
    Ulna->externalForces[0] -= force_U_R_1;

    check1 = Radius->x[size_UR];
    check2 = Ulna->x[size_UR];

    Vector3 force_U_R_2 = k * ((check2 - Vector3(check2[0], 0.0, 0.0)) -
                               (check1 - Vector3(check1[0], 0.0, 0.0)));
    Radius->externalForces[size_UR] += force_U_R_2;
    Ulna->externalForces[size_UR] -= force_U_R_2;

    // Force3 —- Radius—Tendon
    check1 = LongEnd->x[0];
    check2 = Radius->x[4];

    Vector3 force_R_T = k * (check2 - check1);
    LongEnd->externalForces[0] += force_R_T;
    Radius->externalForces[4] -= Vector3(0.0, force_R_T[1], force_R_T[2]);

    // Force4 —- Tendon—Muscles
    for (unsigned int i = 0; i < nom; i++) {
      Rod *Muscle = rodptrs[i + 6];
      if (i < 6) {
        check1 =
            Long->x[5] + 0.0082 * Vector3(cos((15 + i * 60) * M_PI / 180),
                                          sin((15 + i * 60) * M_PI / 180), 0.0);
        check3 = LongEnd->x[6] + Vector3(0.0, 0.0, -0.012013) +
                 0.0082 * Vector3(cos((15 + i * 60) * M_PI / 180),
                                  sin((15 + i * 60) * M_PI / 180), 0.0);
      } else {
        check1 =
            Long->x[5] + 0.016 * Vector3(cos(((i - 6) * 30) * M_PI / 180),
                                         sin(((i - 6) * 30) * M_PI / 180), 0.0);
        check3 = LongEnd->x[6] + Vector3(0.0, 0.0, -0.012013) +
                 0.016 * Vector3(cos(((i - 6) * 30) * M_PI / 180),
                                 sin(((i - 6) * 30) * M_PI / 180), 0.0);
      }
      check2 = Muscle->x[0];
      const int size_Muscle = Muscle->n;
      check4 = Muscle->x[size_Muscle];

      Vector3 force_L_M = 0.5 * k * (check2 - check1);
      Long->externalForces[5] += force_L_M;
      Muscle->externalForces[0] -= force_L_M;

      Vector3 force_M_L = 0.5 * k * (check4 - check3);
      LongEnd->externalForces[6] += force_M_L;
      Muscle->externalForces[size_Muscle] -= force_M_L;
    }

    // Force4 —- Tendon—Muscles
    for (unsigned int i = 0; i < nom; i++) {
      Rod *Muscle = rodptrs[i + 24];
      if (i < 6) {
        check1 = Short->x[6] + 0.0083 * Vector3(cos((15 + i * 60) * M_PI / 180),
                                                sin((15 + i * 60) * M_PI / 180),
                                                0.0);
        check3 = LongEnd->x[5] + Vector3(-0.00008947, -0.0052477, -0.004785) +
                 0.0083 * Vector3(cos((15 + i * 60) * M_PI / 180),
                                  sin((15 + i * 60) * M_PI / 180), 0.0);
      } else {
        check1 = Short->x[6] + 0.017 * Vector3(cos(((i - 6) * 30) * M_PI / 180),
                                               sin(((i - 6) * 30) * M_PI / 180),
                                               0.0);
        check3 = LongEnd->x[5] + Vector3(-0.00008947, -0.0052477, -0.004785) +
                 0.017 * Vector3(cos(((i - 6) * 30) * M_PI / 180),
                                 sin(((i - 6) * 30) * M_PI / 180), 0.0);
      }
      check2 = Muscle->x[0];
      const int size_Muscle = Muscle->n;
      check4 = Muscle->x[size_Muscle];

      Vector3 force_L_M = 0.5 * k * (check2 - check1);
      // cout<<check3<<endl;
      // cout<<check4<<endl;
      Short->externalForces[6] += force_L_M;
      Muscle->externalForces[0] -= force_L_M;

      Vector3 force_M_L = 0.5 * k * (check4 - check3);
      LongEnd->externalForces[5] += force_M_L;
      Muscle->externalForces[size_Muscle] -= force_M_L;
    }

    // Force5 —- Radius—Helices
    check1 = rodptrs[42]->x[0];
    check2 = Radius->x[5] + 0.08 * Radius->Q[4][1];

    Vector3 force_R_H = k * (check2 - check1);
    rodptrs[42]->externalForces[0] += force_R_H;

    const Vector3 direction_arm = check2 - Radius->x[1];
    force_R_H = -Vector3(0.0, force_R_H[1], force_R_H[2]);
    Radius->externalTorques[1] += Radius->Q[1] * (direction_arm * force_R_H);

    // Calculating artificial muscle force
    const REAL AmForce =
        sqrt(force_R_H[1] * force_R_H[1] + force_R_H[2] * force_R_H[2]);

    //-----------------------------------------------------------------------------------------------------------------
    // External loads
    // A constant load at the tip
    REAL Addforce = time * 2000.0;
    Addforce = (Addforce > 200.0) ? 200.0 : Addforce;
    Radius->externalForces[15] += Addforce * Vector3(0.0, 0.0, -1.0);

    // A small initial force to activate the joint flexion
    if (time < 0.1) {
      Humerus->v[20] += Vector3(0.0, 0.0000008, 0.0);
    }
    //-----------------------------------------------------------------------------------------------------------------
    // Muscle fiber interactions
    const Vector3 TBone = LongEnd->x[7] - Humerus->x[17];
    if (TBone.length() < 0.04) {
      LongEnd->externalForces[7] +=
          1 * (0.04 - TBone.length()) * Vector3(0.0, -1000, 0.0);
      Humerus->externalForces[17] +=
          1 * (0.04 - TBone.length()) * Vector3(0.0, 1000, 0.0);
    }

    const REAL zeta = 1.0e1;

    for (unsigned int i = 0; i < 6; i++) {
      Rod *current = rodptrs[i + 6];
      for (unsigned int j = 0; j < 3; j++) {
        Rod *adjacent = current;
        if (j == 0) {
          adjacent = rodptrs[i + 7];
          if (i == 5) adjacent = rodptrs[6];
        }
        if (j == 1) adjacent = rodptrs[(i + 6) * 2];
        if (j == 2) adjacent = rodptrs[(i + 6) * 2 + 1];

        const Vector3 x1 = current->x[7];
        const Vector3 x2 = adjacent->x[7];
        const REAL r1 = current->r[7];
        const REAL r2 = adjacent->r[7];
        const REAL sum_r1r2 = r1 + r2;

        // find the shortest line segment between the two centerline segments
        const Vector3 edge1 = current->edge[7];
        const Vector3 edge2 = adjacent->edge[7];
        const vector<Vector3> minVectors =
            findMinDistVectors(x1, edge1, x2, edge2);
        const Vector3 dVector = minVectors[1];
        const Vector3 dVectorDirection = dVector.unitize();

        // Within the range———>reaction forces
        if ((dVector.length() - sum_r1r2) < 0) {
          // Compute forces
          const Vector3 cForce =
              zeta * (sum_r1r2 - dVector.length()) * dVectorDirection;

          for (unsigned int j = 1; j < 7; j++) {
            current->externalForces[j] -= j / 7.0 * cForce;
            current->externalForces[14 - j] -= j / 7.0 * cForce;
            adjacent->externalForces[j] += j / 7.0 * cForce;
            adjacent->externalForces[14 - j] += j / 7.0 * cForce;
          }
          current->externalForces[7] -= cForce;
          adjacent->externalForces[7] += cForce;
        }
      }
    }

    for (unsigned int i = 0; i < 12; i++) {
      Rod *current = rodptrs[i + 12];

      Rod *adjacent = rodptrs[i + 13];
      if (i == 11) adjacent = rodptrs[12];

      const Vector3 x1 = current->x[7];
      const Vector3 x2 = adjacent->x[7];
      const REAL r1 = current->r[7];
      const REAL r2 = adjacent->r[7];
      const REAL sum_r1r2 = r1 + r2;

      // find the shortest line segment between the two centerline segments
      const Vector3 edge1 = current->edge[7];
      const Vector3 edge2 = adjacent->edge[7];
      const vector<Vector3> minVectors =
          findMinDistVectors(x1, edge1, x2, edge2);
      const Vector3 dVector = minVectors[1];
      const Vector3 dVectorDirection = dVector.unitize();

      // Within the range———>reaction forces
      if ((dVector.length() - sum_r1r2) < 0) {
        const Vector3 cForce =
            zeta * (sum_r1r2 - dVector.length()) * dVectorDirection;

        current->externalForces[7] -= cForce;
        adjacent->externalForces[7] += cForce;
      }
    }

    // cout<<rodptrs[12]->shearInternalForces0[7]<<endl;

    //——————————————————

    for (unsigned int i = 0; i < 6; i++) {
      Rod *current = rodptrs[i + 24];
      for (unsigned int j = 0; j < 3; j++) {
        Rod *adjacent = current;
        if (j == 0) {
          adjacent = rodptrs[i + 25];
          if (i == 5) adjacent = rodptrs[24];
        }
        if (j == 1) adjacent = rodptrs[(i + 24) + i + 6];
        if (j == 2) adjacent = rodptrs[(i + 24) + i + 7];

        const Vector3 x1 = current->x[7];
        const Vector3 x2 = adjacent->x[7];
        const REAL r1 = current->r[7];
        const REAL r2 = adjacent->r[7];
        const REAL sum_r1r2 = r1 + r2;

        // find the shortest line segment between the two centerline segments
        const Vector3 edge1 = current->edge[7];
        const Vector3 edge2 = adjacent->edge[7];
        const vector<Vector3> minVectors =
            findMinDistVectors(x1, edge1, x2, edge2);
        const Vector3 dVector = minVectors[1];
        const Vector3 dVectorDirection = dVector.unitize();

        // Within the range———>reaction forces
        if ((dVector.length() - sum_r1r2) < 0) {
          // Compute forces -- only elastic force here, damping force ignored
          const Vector3 cForce =
              zeta * (sum_r1r2 - dVector.length()) * dVectorDirection;

          for (unsigned int j = 1; j < 7; j++) {
            current->externalForces[j] -= j / 7.0 * cForce;
            current->externalForces[14 - j] -= j / 7.0 * cForce;
            adjacent->externalForces[j] += j / 7.0 * cForce;
            adjacent->externalForces[14 - j] += j / 7.0 * cForce;
          }
          current->externalForces[7] -= cForce;
          adjacent->externalForces[7] += cForce;
        }
      }
    }

    for (unsigned int i = 0; i < 12; i++) {
      Rod *current = rodptrs[i + 30];

      Rod *adjacent = rodptrs[i + 31];
      if (i == 11) adjacent = rodptrs[30];

      const Vector3 x1 = current->x[7];
      const Vector3 x2 = adjacent->x[7];
      const REAL r1 = current->r[7];
      const REAL r2 = adjacent->r[7];
      const REAL sum_r1r2 = r1 + r2;

      // find the shortest line segment between the two centerline segments
      const Vector3 edge1 = current->edge[7];
      const Vector3 edge2 = adjacent->edge[7];
      const vector<Vector3> minVectors =
          findMinDistVectors(x1, edge1, x2, edge2);
      const Vector3 dVector = minVectors[1];
      const Vector3 dVectorDirection = dVector.unitize();

      // Within the range———>reaction forces
      if ((dVector.length() - sum_r1r2) < 0) {
        const Vector3 cForce =
            zeta * (sum_r1r2 - dVector.length()) * dVectorDirection;

        current->externalForces[7] -= cForce;
        adjacent->externalForces[7] += cForce;
      }
    }

    // Printing values
    const Vector3 Rdirection = (Radius->x[15] - Radius->x[0]).unitize();
    angle_new = angle(Rdirection, Vector3(0.0, 0.0, -1.0));
    angleVelocity = (angle_new - angle_old) / (0.5e-7);
    angle_old = angle_new;

    ml_new = (rodptrs[10]->x[0] - rodptrs[10]->x[14]).length();
    mVelocity = (ml_new - ml_old) / (0.5e-7);
    ml_old = ml_new;

    // Calculating actuating torque
    const Vector3 TorqArm = 0.017 * 3.0 * (Radius->edge[1]).unitize();
    const Vector3 TorqForce = Vector3(0.0, Radius->externalForces[4][1],
                                      Radius->externalForces[4][2]);
    const Vector3 ActTorq = TorqArm * TorqForce;

    if (Localtime < 0.0) {
      ofstream control;
      control.open("velocity.txt",
                   (time == 2.5e-8) ? std::ofstream::out : std::ofstream::app);
      control << time << " " << angle_new << " " << mVelocity << " "
              << ActTorq[0] << " " << AmForce << endl;

      Localtime = 0.002;
    }
    Localtime -= 0.5e-7;

    angleVelocity_old = angleVelocity;
  }
#endif

//****************________________________________****************________________________________
//****************———————————————————————————————-****************________________________________
#ifdef FLAGFLAGELLA
  // This function takes care of the various connections within the flagella
  // case.

  void RodRodAttach(const REAL time) {
    REAL kposition = 3.8e2;

    nor = (int)rodptrs.size();

    Rod *body = rodptrs[0];
    Rod *muscle = rodptrs[1];
    Vector3 normal = Vector3(0.0, 0.0, 1.0);

    const int size = muscle->n;
    const int span = 1;
    int start = attachpoint[0].first;
    int direction = attachpoint[0].second;
    int end = start + span;

    const REAL armlength = 0.0053;

    // Attach two end points and generate torque to the body
    const Vector3 armdirection1 = normal * body->edge[start];
    const Vector3 armposition1 =
        direction * armlength * (armdirection1.unitize());
    const Vector3 startposition =
        (body->x[start] + body->x[start + 1]) / 2.0 + armposition1;
    Vector3 armdirection2 = normal * body->edge[end];
    const Vector3 armposition2 =
        direction * armlength * (armdirection2.unitize());
    const Vector3 endposition =
        (body->x[end] + body->x[end + 1]) / 2.0 + armposition2;

    Vector3 forcestart = kposition * (muscle->x[0] - startposition);
    Vector3 forceend = kposition * (muscle->x[size] - endposition);

    muscle->externalForces[0] -= forcestart;
    muscle->externalForces[size] -= forceend;

    Vector3 Torque1 =
        armposition1 * (forcestart - Vector3(0.0, 0.0, forcestart[2]));
    Vector3 Torque2 =
        armposition2 * (forceend - Vector3(0.0, 0.0, forceend[2]));

    REAL Torqueaverage2 = Torque1[2] > 0 ? (0.5 * (Torque1[2] - Torque2[2]))
                                         : (-0.5 * (-Torque1[2] + Torque2[2]));
    Vector3 Torqueaverage = Vector3(0.0, 0.0, Torqueaverage2);

    body->externalTorques[start] += body->Q[start] * (Torqueaverage) / 2;
    body->externalTorques[end] += body->Q[end] * (-1 * Torqueaverage) / 2;
  }
#endif

//****************________________________________****************________________________________
//****************———————————————————————————————-****************________________________________
#ifdef FLAGWALKER
  // This function takes care of the various connections, ground setup and
  // muscle actuations within the walker case.

  void RodRodAttach(const REAL time) {
    //-----------------------------------------------------------------------------------------------------------------
    // Connections
    const REAL LegLength =
        rodptrs[2]
            ->MaxHeight;  // This is passed from WAlker.cpp as a proxy variable

    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);
    Vector3 check3 = Vector3(0.0, 0.0, 0.0);
    Vector3 check4 = Vector3(0.0, 0.0, 0.0);

    const REAL k = 5e8;  // 5e8

    const REAL width = 3.5;

    // In the Top Rectangle
    // 1111111111111111
    /* This follows from a simple fixed joint, see FixedJoint.cpp for more
     * details*/
    // Position Connection
    check1 = rodptrs[0]->x[2];
    check2 = rodptrs[6]->x[0];

    Vector3 force_R_R = k * (check2 - check1);

    /* Equal and opposite restoring torques */
    rodptrs[0]->externalForces[2] += force_R_R;
    rodptrs[6]->externalForces[0] -= force_R_R;

    // Perpendicular
    ///////////////////////////////////////////////////////////////////////////
    /*
    In this block of code, we demonstrate how to form a fixed joint connection
    for two rods that are aligned at 90 degrees to one another like so

                                          Rod6
                                           |            ^
                                           |            |
                                           |            | width (in meters)
                                           |            |
                                           |            v
            ============================= (f)
                        Rod0         Fixed joint

    This is a prototypical way of doing such a joint and we explain the details
    in this instance only. Other joints in this function are done similarly.

    Also see FixedJoint.cpp for more details.
    */

    /*
      In this example, we are constraining rod 0 and rod 6 (see ASCII figure
      above).
    */
    /* check3 expresses the desired position of the last element of rod6
    based on connectivity */
    check3 = rodptrs[0]->x[2] + width * rodptrs[0]->Q[2][0];
    /* check4 represents the current position of the last element of rod6 */
    check4 = rodptrs[6]->x[5];

    /* Force is calculated based on a factor multiplied by the difference in
    target position and current position, akin to a linear spring*/
    force_R_R = 0.0001 * k * (check4 - check3);

    /* Equal and opposite restoring torques are now applied on each rod. The
    force is calculated in the lab frame, and the torques need to be calculated
    in the material frame, hence the multiplication by Q matrices below.*/
    rodptrs[0]->externalTorques[2] +=
        rodptrs[0]->Q[2] * ((rodptrs[0]->Q[2][0]) * force_R_R);
    rodptrs[6]->externalTorques[0] +=
        rodptrs[6]->Q[0] * (rodptrs[6]->Q[0][2] * (-force_R_R));
    ///////////////////////////////////////////////////////////////////////////

    // 22222222222222222
    // Position Connection
    check1 = rodptrs[0]->x[14];
    check2 = rodptrs[7]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[0]->externalForces[14] += force_R_R;
    rodptrs[7]->externalForces[0] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[0]->x[14] + width * rodptrs[0]->Q[14][0];
    check4 = rodptrs[7]->x[5];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[0]->externalTorques[14] +=
        rodptrs[0]->Q[14] * ((rodptrs[0]->Q[14][0]) * force_R_R);
    rodptrs[7]->externalTorques[0] +=
        rodptrs[7]->Q[0] * (rodptrs[7]->Q[0][2] * (-force_R_R));

    // 33333333333333333
    // Position Connection
    check1 = rodptrs[1]->x[2];
    check2 = rodptrs[6]->x[5];

    force_R_R = k * (check2 - check1);
    rodptrs[1]->externalForces[2] += force_R_R;
    rodptrs[6]->externalForces[5] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[1]->x[2] - width * rodptrs[1]->Q[2][0];
    check4 = rodptrs[6]->x[0];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[1]->externalTorques[2] +=
        rodptrs[1]->Q[2] * (-rodptrs[1]->Q[2][0] * force_R_R);
    rodptrs[6]->externalTorques[4] +=
        rodptrs[6]->Q[4] * ((-rodptrs[6]->Q[4][2]) * (-force_R_R));

    // 44444444444444444
    // Position Connection
    check1 = rodptrs[1]->x[14];
    check2 = rodptrs[7]->x[5];

    force_R_R = k * (check2 - check1);
    rodptrs[1]->externalForces[14] += force_R_R;
    rodptrs[7]->externalForces[5] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[1]->x[14] - width * rodptrs[1]->Q[14][0];
    check4 = rodptrs[7]->x[0];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[1]->externalTorques[14] +=
        rodptrs[1]->Q[14] * (-rodptrs[1]->Q[14][0] * force_R_R);
    rodptrs[7]->externalTorques[4] +=
        rodptrs[7]->Q[4] * ((-rodptrs[7]->Q[4][2]) * (-force_R_R));

    // Connections between beam and leg
    // LegLegLegLegLeg1111111111111111111
    check1 = rodptrs[0]->x[2];
    check2 = rodptrs[2]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[0]->externalForces[2] += force_R_R;
    rodptrs[2]->externalForces[0] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[0]->x[2] - LegLength * rodptrs[0]->Q[2][1];
    check4 = rodptrs[2]->x[6];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[0]->externalTorques[2] +=
        rodptrs[0]->Q[2] * ((-rodptrs[0]->Q[2][1]) * force_R_R -
                            0.00001 * k * rodptrs[0]->Q[2][0] *
                                (rodptrs[0]->w[2][0] - rodptrs[2]->w[0][0]));
    rodptrs[2]->externalTorques[0] +=
        rodptrs[2]->Q[0] * (rodptrs[2]->Q[0][2] * (-force_R_R) +
                            0.00001 * k * rodptrs[0]->Q[2][0] *
                                (rodptrs[0]->w[2][0] - rodptrs[2]->w[0][0]));

    // LegLegLegLegLeg2222222222222222222
    check1 = rodptrs[0]->x[14];
    check2 = rodptrs[3]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[0]->externalForces[14] += force_R_R;
    rodptrs[3]->externalForces[0] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[0]->x[14] - 3.4 * rodptrs[0]->Q[14][1];
    check4 = rodptrs[3]->x[6];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[0]->externalTorques[14] +=
        rodptrs[0]->Q[14] * ((-rodptrs[0]->Q[14][1]) * force_R_R -
                             0.00001 * k * rodptrs[0]->Q[14][0] *
                                 (rodptrs[0]->w[14][0] - rodptrs[3]->w[0][0]));
    rodptrs[3]->externalTorques[0] +=
        rodptrs[3]->Q[0] * (rodptrs[3]->Q[0][2] * (-force_R_R) +
                            0.00001 * k * rodptrs[0]->Q[14][0] *
                                (rodptrs[0]->w[14][0] - rodptrs[3]->w[0][0]));

    // LegLegLegLegLeg33333333333333333333
    check1 = rodptrs[1]->x[2];
    check2 = rodptrs[4]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[1]->externalForces[2] += force_R_R;
    rodptrs[4]->externalForces[0] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[1]->x[2] - LegLength * rodptrs[1]->Q[2][1];
    check4 = rodptrs[4]->x[6];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[1]->externalTorques[2] +=
        rodptrs[1]->Q[2] * ((-rodptrs[1]->Q[2][1]) * force_R_R -
                            0.00001 * k * rodptrs[1]->Q[2][0] *
                                (rodptrs[1]->w[2][0] - rodptrs[4]->w[0][0]));
    rodptrs[4]->externalTorques[0] +=
        rodptrs[4]->Q[0] * (rodptrs[4]->Q[0][2] * (-force_R_R) +
                            0.00001 * k * rodptrs[1]->Q[2][0] *
                                (rodptrs[1]->w[2][0] - rodptrs[4]->w[0][0]));

    // LegLegLegLegLeg44444444444444444444
    check1 = rodptrs[1]->x[14];
    check2 = rodptrs[5]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[1]->externalForces[14] += force_R_R;
    rodptrs[5]->externalForces[0] -= force_R_R;

    // Perpendicular
    check3 = rodptrs[1]->x[14] - 3.4 * rodptrs[1]->Q[14][1];
    check4 = rodptrs[5]->x[6];

    force_R_R = 0.0001 * k * (check4 - check3);
    rodptrs[1]->externalTorques[14] +=
        rodptrs[1]->Q[14] * ((-rodptrs[1]->Q[14][1]) * force_R_R -
                             0.00001 * k * rodptrs[0]->Q[14][0] *
                                 (rodptrs[1]->w[14][0] - rodptrs[5]->w[0][0]));
    rodptrs[5]->externalTorques[0] +=
        rodptrs[5]->Q[0] * (rodptrs[5]->Q[0][2] * (-force_R_R) +
                            0.00001 * k * rodptrs[0]->Q[14][0] *
                                (rodptrs[1]->w[14][0] - rodptrs[5]->w[0][0]));

    // Connections of the crosses
    check1 = rodptrs[2]->x[4];
    check2 = rodptrs[8]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[2]->externalForces[4] += force_R_R;
    rodptrs[8]->externalForces[0] -= force_R_R;

    check1 = rodptrs[4]->x[4];
    check2 = rodptrs[8]->x[5];

    force_R_R = k * (check2 - check1);
    rodptrs[4]->externalForces[4] += force_R_R;
    rodptrs[8]->externalForces[5] -= force_R_R;

    check1 = rodptrs[3]->x[4];
    check2 = rodptrs[9]->x[0];

    force_R_R = k * (check2 - check1);
    rodptrs[3]->externalForces[4] += force_R_R;
    rodptrs[9]->externalForces[0] -= force_R_R;

    check1 = rodptrs[5]->x[4];
    check2 = rodptrs[9]->x[5];

    force_R_R = k * (check2 - check1);
    rodptrs[5]->externalForces[4] += force_R_R;
    rodptrs[9]->externalForces[5] -= force_R_R;

    //-----------------------------------------------------------------------------------------------------------------
    // Setup the ground and the friction model

    const REAL kPlane = 1e6;
    const REAL etaPlane = 1e4;
    const REAL muKinetic = 0.3;
    const REAL muStaticF = 0.42;
    const REAL vStatic = 1e-4;

    rodptrs[0]->Tforce = 0.0;
    rodptrs[0]->MaxHeight = 0.0;
    rodptrs[1]->Tforce = 0.0;
    rodptrs[1]->MaxHeight = 0.0;

    // Four legs in our simulation -- two pillars in experiment
    for (unsigned int j = 2; j < 6; j++) {
      Rod *Leg = rodptrs[j];

      // Each leg may contact with ground at two corners
      for (unsigned int i = 0; i < 2; i++) {
        const REAL factor = (i == 0) ? -1.0 : 1.0;
        const Vector3 CheckPoint = Leg->x[6] + factor * 1.75 * (Leg->Q[5][1]);

        if (CheckPoint[2] < 1e-4) {
          const Vector3 elementV = Leg->v[6];
          const Vector3 normalPlane = Vector3(0.0, 0.0, 1.0);

          const Vector3 overallRodForces =
              Leg->totalInternalForces[6] + Leg->externalForces[6];
          const REAL currentRodForcesInNormaldirectionSign =
              (overallRodForces % normalPlane);
          const Vector3 currentRodForcesInNormaldirection =
              (currentRodForcesInNormaldirectionSign < 0.0)
                  ? currentRodForcesInNormaldirectionSign * normalPlane
                  : 0.0 * normalPlane;
          const Vector3 elasticPlaneResponse =
              -kPlane * min(CheckPoint[2], 0.0) * normalPlane;
          const Vector3 dampingForcePlane =
              -etaPlane * (elementV % normalPlane) * normalPlane;
          const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
          const Vector3 totalforce =
              (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
          // Ground surface pushing force
          Leg->externalForces[6] += totalforce;

          const REAL normalForce = normalForcePlane.length();

          // Compute axial velocity
          const Vector3 slipVel =
              elementV - (elementV % normalPlane) * normalPlane;

          // Friction in axial direction --- here walker goes in a straight line
          // so we ignore the forces in lateral directions
          {
            const REAL vel = fabs(slipVel.length());
            const REAL width = 0.5 * vStatic;
            const REAL velDiff = vel - vStatic;
            REAL f = 1.0;

            if (vel > vStatic) f = 1.0 - velDiff / width;
            if (vel > (vStatic + width)) f = 0.0;

            const Vector3 kineticFrictionForce =
                -(1.0 - f) * muKinetic * normalForce * slipVel.unitize();
            const Vector3 staticFrictionForceF =
                f * muStaticF * normalForce * Vector3(0.0, 1.0, 0.0);

            Leg->kineticFrictionsForce[6] += kineticFrictionForce;
            if ((j == 2) || (j == 4))
              rodptrs[0]->Tforce +=
                  kineticFrictionForce % Vector3(0.0, 1.0, 0.0);

            if ((j == 3) || (j == 5))
              rodptrs[1]->Tforce +=
                  kineticFrictionForce % Vector3(0.0, 1.0, 0.0);

            // same friction coefficients for forward and backward static
            // frictions.
            Leg->staticFrictionsAxialForceForward[6] += staticFrictionForceF;
            Leg->staticFrictionsAxialForceBackward[6] += staticFrictionForceF;
          }  // Friction in axial direction

        }  // if within range
      }    // for 2 corners
    }      // for 4 legs

    //-----------------------------------------------------------------------------------------------------------------
    // Muscle actuation

    // Muscle Actuation Muscle Actuation Muscle Actuation Muscle Actuation
    // Muscle Actuation
    if (time > 0.1) {
      Vector3 check1 = Vector3(0.0, 0.0, 0.0);
      Vector3 check2 = Vector3(0.0, 0.0, 0.0);

      /*
        Here we initialize the muscle ring to circumscribe the scaffold, with
        some initial gap. To make the muscle conform to the scaffold then, we
        include a small constant, uniform, compressive base force of 200 microN
        (compared to the active muscle force) force on the muscle.

        We ramp this force up to account for initialization effects.
      */
      REAL Base = 1000.0 * (time - 0.1);
      Base = (Base > 200.0) ? 200.0 : Base;

      REAL realtime = time - (0.1);
      realtime = (realtime < 0.0) ? 0.0 : realtime;

      // Muscle Contraction
      const REAL T = 1.0;
      const REAL w = 2.0 * M_PI / T;
      const REAL Kattach = 1.0e8;
      /*
      The muscle contractile force =  Muscle contractile stress sigma_m *
      PI * r^2, where sigma_m = 3.706 kPa and r = 0.227 mm. The derivation
      of sigma_m can be found in Supplementary Note 7.2, which we summarize
      below for the sake of clarity.

// clang-format off
      sigma_m = \frac{1}{\gamma} \left( F_m/A_m + E_m * \epsilon/ (1 -\epsilon)
\right)
// clang-format on

       where ratio of active-total muscle cross section area \gamma = 1,
       muscle output force F_m = 1.4 mN, muscle total area A_m = 0.648 mm^2
       and muscle strain \epsilon = 0.134 are characterized from experiments
       (see ref. 28 in SI). Muscle’s Young’s modulus Em = 10 kPa is found in
       Supplementary Table 1.

       Thus based on this, the contractile force given to each individual
       muscle fiber is 600 mN. We note that the muscle rings (four in number)
       have a total cross-sectional area 1.66 times larger than the cross
       sectional area of the straight fibers (i.e. in the muscle strip).
       On the other hand, being thicker, they contain some necrotic tissue. This
       means that \gamma is some value between 0.6 and 0.8 based on experimental
       estimates (following from ref. 28 in SI). Here for simplicity, we set the
       gamma to be 0.6 so that the contractile stresses in the ring is still 600
       mN.
      */
      const REAL Amp = 600.0;  // This amplitude will make the muscle output the
                               // same force measured in experiment

      for (unsigned int j = 10; j < 18; j++) {
        Rod *muscle = rodptrs[j];
        const int size = (j < 14) ? (muscle->n - 2) : muscle->n;
        REAL Actuation2 = Amp * abs(sin(w * realtime));

        muscle->externalForces[0] =
            (muscle->edge[0].unitize()) * (Actuation2 + Base);
        muscle->externalForces[size] =
            -1 * (muscle->edge[size - 1].unitize()) * (Actuation2 + Base);
        for (unsigned int i = 1; i < (size); i++) {
          muscle->externalForces[i] +=
              (muscle->edge[i].unitize()) * (Actuation2 + Base);
          muscle->externalForces[i] +=
              -1 * (muscle->edge[i - 1].unitize()) * (Actuation2 + Base);
        }  // for elements

        // Muscle End Connections
        if (j < 14) {
          const Vector3 direction = (muscle->x[size] - muscle->x[0]).unitize();
          muscle->externalForces[0] +=
              direction * (Actuation2 + Base) +
              Kattach * (muscle->x[size] - muscle->x[0]);
          muscle->externalForces[size] -=
              direction * (Actuation2 + Base) +
              Kattach * (muscle->x[size] - muscle->x[0]);

          if ((j == 10) || (j == 12)) {
            for (unsigned int i = 12; i < 18; i++) {
              check1 = rodptrs[j]->x[i];
              check2 = rodptrs[j + 1]->x[i];
              Vector3 force_ = Kattach * (check2 - check1);
              rodptrs[j]->externalForces[i] +=
                  Vector3(force_[0], force_[1], 0.0);
              rodptrs[j + 1]->externalForces[i] -=
                  Vector3(force_[0], force_[1], 0.0);
            }
          }

        } else {
          check1 = rodptrs[j]->x[0];
          check2 = rodptrs[12]->x[j - 1];

          Vector3 force_ = Kattach * (check2 - check1);
          rodptrs[j]->externalForces[0] += force_;
          rodptrs[12]->externalForces[j - 1] -= force_;

          check1 = rodptrs[j]->x[size];
          check2 = rodptrs[10]->x[j - 1];

          force_ = Kattach * (check2 - check1);
          rodptrs[j]->externalForces[size] += force_;
          rodptrs[10]->externalForces[j - 1] -= force_;
        }  // connection
      }    // for muscles
    }
  }
#endif

//****************________________________________****************________________________________
//****************———————————————————————————————-****************________________________________
#ifdef FLAGMUSCULARSNAKE
  // This function takes care of the various connections and muscle actuations
  // within the muscularsnake case.

  void RodRodAttach(const REAL time) {
    //-----------------------------------------------------------------------------------------------------------------
    // Connections
    REAL kposition = 5.0 * 1e9;
    REAL kattach = 5 * 1e6;

    nor = (int)rodptrs.size();

    Rod *snake = rodptrs[0];
    // These are torque amplitude values. We convert these to muscle
    // forces later on
    REAL Amp[8] = {22.96, 22.96, 20.95, 20.95, 9.51, 9.51, 13.7, 13.7};
    for (unsigned int k = 1; k < 9; k++) {
      Rod *muscle = rodptrs[k];
      const int size = muscle->n;
      const int start = muscle->mindex;
      const int direction = ((k % 2) == 1) ? 1 : -1;
      const int end = start + 3 * size;

      // Attach two end points and generate torque to the snake body
      const Vector3 startposition =
          snake->x[start] + direction * snake->Q[start][1] * 0.025 / 2.0;
      const Vector3 endposition =
          snake->x[end] + direction * snake->Q[end][1] * 0.025 / 2.0;

      const Vector3 forcestart = kposition * (muscle->x[0] - startposition);
      const Vector3 forceend = kposition * (muscle->x[size] - endposition);

      muscle->externalForces[0] -= forcestart;
      muscle->externalForces[size] -= forceend;

      const Vector3 Torque1 = (direction * snake->Q[start][1] * 0.025 / 2.0) *
                              (forcestart - Vector3(0.0, 0.0, forcestart[2]));
      const Vector3 Torque2 = (direction * snake->Q[end][1] * 0.025 / 2.0) *
                              (forceend - Vector3(0.0, 0.0, forceend[2]));

      /*
        The Torque1 and Torque2 experienced at the end of the muscle groups are
        (as expected) consistently found to be almost the same. Nevertheless,
        to filter out potential numerical noise effects, we average the two
        torques and apply this averaged torque to the snake body, as we would
        expect from a biological snake.
      */
      REAL Torqueaverage2 = Torque1[2] > 0
                                ? (0.5 * (Torque1[2] - Torque2[2]))
                                : (-0.5 * (-Torque1[2] + Torque2[2]));
      Vector3 Torqueaverage = Vector3(0.0, 0.0, Torqueaverage2);

      snake->externalTorques[start] += snake->Q[start] * Torqueaverage;
      snake->externalTorques[end] += snake->Q[end] * (-1 * Torqueaverage);

      // Glue the muscle to the snake body
      // This means that the muscle experiences an attractive force towards
      // the body, but we ignore the effect of the opposite, equivalent force
      // applied on the body
      for (unsigned int j = 1; j < 5; j++) {
        int indexsnake = start + j * 3;

        const Vector3 snakeposition =
            snake->x[indexsnake] +
            direction * 0.025 / 4 *
                (snake->Q[indexsnake][1] + snake->Q[indexsnake - 1][1]);
        const Vector3 forcedirection =
            0.5 * (muscle->Q[j - 1][1] + muscle->Q[j][1]);
        Vector3 force = kattach * ((muscle->x[j] - snakeposition) %
                                   forcedirection * forcedirection);
        muscle->externalForces[j] -= force;
      }
      for (unsigned int j = 9; j < size; j++) {
        int indexsnake = start + j * 3;

        const Vector3 snakeposition =
            snake->x[indexsnake] +
            direction * 0.025 / 4 *
                (snake->Q[indexsnake][1] + snake->Q[indexsnake - 1][1]);
        const Vector3 forcedirection =
            0.5 * (muscle->Q[j - 1][1] + muscle->Q[j][1]);
        Vector3 force = kattach * ((muscle->x[j] - snakeposition) %
                                   forcedirection * forcedirection);
        muscle->externalForces[j] -= force;
      }

      //-----------------------------------------------------------------------------------------------------------------
      // Muscle activations
      const int index = muscle->mindex;
      const REAL timeD = index * 1.0 / 101.76;

      const REAL Realtime = time - timeD;

      REAL ramp = (Realtime <= 1.0) ? (Realtime / 1.0) : 1.0;
      ramp = (Realtime <= 0.0) ? 0.0 : ramp;

      const REAL w = 2.0 * M_PI / 1.0;
      /*
        For legacy purposes, the input from the optimizer is given in terms of
        torque amplitudes (see Gazzola et al. Royal Society Open Source, 2018).

        We then use the same set up and map those values into muscle
        contractile forces by dividing them by the arm length. This is captured
        through the parameter "factor". This also enables a meaningful
        comparison with the continuum snake case of the above reference.
      */

      // Convert from torques (Amp) to muscle forces
      // 80 = 1./0.0125, 0.0125 is the arm length
      // This division is to map the torques optimized by CMA into the
      // contraction forces of our muscles
      const REAL factor = 80.0;

      // left and right muscles
      if ((k % 2) == 1) {
        for (unsigned int i = 4; i < 9; i++) {
          muscle->externalForces[i] += ramp * (muscle->edge[i].unitize()) *
                                       factor * Amp[k - 1] *
                                       (0.5 * sin(w * Realtime) + 0.5);
        }
        for (unsigned int i = 5; i < 10; i++) {
          muscle->externalForces[i] +=
              -1 * ramp * (muscle->edge[i - 1].unitize()) * factor *
              Amp[k - 1] * (0.5 * sin(w * Realtime) + 0.5);
        }
      }
      if ((k % 2) == 0) {
        for (unsigned int i = 4; i < 9; i++) {
          muscle->externalForces[i] += ramp * (muscle->edge[i].unitize()) *
                                       factor * Amp[k - 1] *
                                       (-0.5 * sin(w * Realtime) + 0.5);
        }
        for (unsigned int i = 5; i < 10; i++) {
          muscle->externalForces[i] +=
              -1 * ramp * (muscle->edge[i - 1].unitize()) * factor *
              Amp[k - 1] * (-0.5 * sin(w * Realtime) + 0.5);
        }
      }
    }
  }
#endif

//****************________________________________****************________________________________
//****************———————————————————————————————-****************________________________________
#ifdef FLAGWING
  // This function takes care of the various connections and muscle actuations
  // within the wing case.

  void RodRodAttach(const REAL time) {
    REAL Liftforce = 0.0;
    //-----------------------------------------------------------------------------------------------------------------
    // Connections Muscle and Bones
    const int Rodsize = rodptrs.size();
    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);
    Vector3 check3 = Vector3(0.0, 0.0, 0.0);
    Vector3 check4 = Vector3(0.0, 0.0, 0.0);

    const REAL k = 5e8;

    Rod *Humerus = rodptrs[0];
    Rod *Radius = rodptrs[1];
    Rod *Ulna = rodptrs[2];
    Rod *Carp1 = rodptrs[3];
    Rod *Carp2 = rodptrs[4];
    Rod *Carp3 = rodptrs[5];

    Rod *LiftMu = rodptrs[Rodsize - 4];
    Rod *PullMu = rodptrs[Rodsize - 3];
    Rod *BiMu = rodptrs[Rodsize - 2];
    Rod *TriMu = rodptrs[Rodsize - 1];

    check1 = Humerus->x[8];
    check2 = LiftMu->x[7];

    Vector3 force_M_B = 0.01 * k * (check2 - check1);
    Humerus->externalForces[8] += force_M_B;
    LiftMu->externalForces[7] -= force_M_B;

    check1 = Humerus->x[8];
    check2 = PullMu->x[7];

    force_M_B = 0.01 * k * (check2 - check1);
    Humerus->externalForces[8] += force_M_B;
    PullMu->externalForces[7] -= force_M_B;

    check1 = Radius->x[3];
    check2 = BiMu->x[7];

    force_M_B = 0.01 * k * (check2 - check1);
    Radius->externalTorques[0] +=
        4.0 * 15.0 * Radius->Q[0] * (Radius->Q[0][2] * force_M_B);
    BiMu->externalForces[7] -= force_M_B;

    check1 = Radius->x[0] - 15.0 * Radius->Q[0][2];
    check2 = TriMu->x[7];

    force_M_B = 0.01 * k * (check2 - check1);
    Radius->externalTorques[0] +=
        3.0 * 15.0 * Radius->Q[0] * (-Radius->Q[0][2] * force_M_B);
    TriMu->externalForces[7] -= force_M_B;

    //-----------------------------------------------------------------------------------------------------------------
    // Muscle actuation pattern -- these are hard tuned for initiation phase and
    // one full cycle Peak forces are scaled according to rule 1.

    REAL PhaseForce = 0.0;
    // Lifts
    REAL Amp = 0.0;

    if (time < 0.125) {
      Amp = 15e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * time);
    } else if (time < 0.28) {
      PhaseForce = 0.0;
    } else if (time < (0.28 + 0.125 / 2.0)) {
      Amp = 18.0e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - 0.28));
    } else if (time < (0.405)) {
      Amp = 18e6 * (1.0 - 0.6 * ((time - 0.3425) / 0.0625));
      PhaseForce = 1.0;
    } else if (time < (0.405 + 0.125 / 2.0)) {
      Amp = 0.4 * 18.0e6;
      PhaseForce = 1.0;
    } else {
      PhaseForce = 0.0;
    }
    for (unsigned int i = 2; i < 5; i++) {
      LiftMu->externalForces[i] +=
          (LiftMu->edge[i].unitize()) * Amp * PhaseForce;
    }
    for (unsigned int i = 3; i < 6; i++) {
      LiftMu->externalForces[i] +=
          -1 * (LiftMu->edge[i - 1].unitize()) * Amp * PhaseForce;
    }

    // Pulling Pectoralis
    Amp = 0.0e6;
    if (time < 0.11) {
      PhaseForce = 0.0;
    } else if (time < (0.11 + 0.125 / 2.0)) {
      Amp = 18e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - 0.11));
    } else if (time < (0.235)) {
      Amp = 18e6 * (1.0 - 0.5 * ((time - 0.1725) / 0.0625));
      PhaseForce = 1.0;
    } else if (time < (0.235 + 0.125 / 2.0)) {
      Amp = 0.5 * 18.0e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - (0.11 + 0.125 / 2.0)));
    } else {
      PhaseForce = 0.0;
    }
    for (unsigned int i = 2; i < 5; i++) {
      PullMu->externalForces[i] +=
          (PullMu->edge[i].unitize()) * Amp * PhaseForce;
    }
    for (unsigned int i = 3; i < 6; i++) {
      PullMu->externalForces[i] +=
          -1 * (PullMu->edge[i - 1].unitize()) * Amp * PhaseForce;
    }

    // Triceps
    Amp = 0.0;
    PhaseForce = 0.0;
    if (time < (0.125 / 2.0)) {
      Amp = 27.0e6;
      PhaseForce = sin(4 * M_PI / 0.25 * time);
    } else if (time < 0.27) {
      PhaseForce = 0.0;
    } else if (time < (0.27 + 0.125 / 2.0)) {
      Amp = 8.0e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - 0.27));
    } else if (time < (0.395)) {
      Amp = 8.0e6 * (1.0 - 0.35 * ((time - 0.3325) / 0.0625));
      PhaseForce = 1.0;
    } else if (time < (0.395 + 0.125 / 2.0)) {
      Amp = 0.65 * 8.0e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - 0.3325));
    } else {
      PhaseForce = 0.0;
    }
    for (unsigned int i = 2; i < 5; i++) {
      TriMu->externalForces[i] += (TriMu->edge[i].unitize()) * Amp * PhaseForce;
    }
    for (unsigned int i = 3; i < 6; i++) {
      TriMu->externalForces[i] +=
          -1 * (TriMu->edge[i - 1].unitize()) * Amp * PhaseForce;
    }

    // Biceps
    Amp = 0.0;

    if (time < 0.11) {
      Amp = 0.0e6;
      PhaseForce = time / 0.1;
    } else if (time < (0.11 + 0.125 / 2.0)) {
      Amp = 13.5e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - 0.11));
    } else if (time < (0.235)) {
      Amp = 24e6 * (1.0 - 0.7 * ((time - 0.1725) / 0.0625));
      PhaseForce = 1.0;
    } else if (time < (0.235 + 0.125 / 2.0)) {
      Amp = 0.3 * 24e6;
      PhaseForce = sin(2.0 * M_PI / 0.25 * (time - (0.11 + 0.125 / 2.0)));
    } else {
      PhaseForce = 0.0;
    }
    for (unsigned int i = 2; i < 5; i++) {
      BiMu->externalForces[i] += (BiMu->edge[i].unitize()) * Amp * PhaseForce;
    }
    for (unsigned int i = 3; i < 6; i++) {
      BiMu->externalForces[i] +=
          -1 * (BiMu->edge[i - 1].unitize()) * Amp * PhaseForce;
    }

    const int size_H = Humerus->n;

    //-----------------------------------------------------------------------------------------------------------------
    // Get joint angles and muscle strains
    // Get joint angles
    // Dorsal ventral
    Humerus->MaxHeight = atan((Humerus->x[size_H][2] - Humerus->x[0][2]) /
                              abs(Humerus->x[size_H][1])) *
                         180.0 / M_PI;
    // Antero posterior
    Humerus->Tforce =
        atan(-Humerus->x[size_H][0] / abs(Humerus->x[size_H][1])) * 180.0 /
        M_PI;
    // elbow
    Radius->MaxHeight =
        90.0 - 180.0 / M_PI * angle(Humerus->Q[size_H - 1][2], Radius->Q[0][2]);

    // Get Muscle strain
    REAL musclelength = 0.0;
    for (unsigned int i = 0; i < 7; i++) musclelength += LiftMu->l[i];

    LiftMu->MaxHeight = musclelength / (LiftMu->L0) - 1.0;

    musclelength = 0.0;
    for (unsigned int i = 0; i < 7; i++) musclelength += PullMu->l[i];

    PullMu->MaxHeight = musclelength / (PullMu->L0) - 1.0;

    musclelength = 0.0;
    for (unsigned int i = 0; i < 7; i++) musclelength += BiMu->l[i];

    BiMu->MaxHeight = musclelength / (BiMu->L0) - 1.0;

    musclelength = 0.0;
    for (unsigned int i = 0; i < 7; i++) musclelength += TriMu->l[i];

    TriMu->MaxHeight = musclelength / (TriMu->L0) - 1.0;

    //-----------------------------------------------------------------------------------------------------------------
    // Connections Bones

    // Force1 —- Humerus—Radius
    check1 = Humerus->x[size_H] +
             8.4 * (0.8 * Humerus->Q[size_H - 1][2] - Humerus->Q[size_H - 1][1])
                       .unitize();
    check2 = Radius->x[0];

    Vector3 force_H_R = k * (check2 - check1);
    Humerus->externalForces[size_H] += force_H_R;
    Radius->externalForces[0] -= force_H_R;

    const Vector3 normalH = Humerus->Q[size_H - 1][0];
    check3 = check1 +
             90 * (0.84 * Humerus->Q[size_H - 1][2] + Humerus->Q[size_H - 1][1])
                      .unitize();
    check4 = Radius->x[6];

    force_H_R = 10 * k * ((check4 - check3) % normalH) * normalH;
    Humerus->externalForces[size_H] += Vector3(force_H_R[0], 0.0, force_H_R[2]);
    Radius->externalForces[6] -= Vector3(force_H_R[0], 0.0, force_H_R[2]);

    // Force2 —- Radius—Ulna
    const int size_UR = Radius->n;
    check1 = Radius->x[0] +
             2.1 * 8.4 * (-0.2 * Radius->Q[0][2] + Radius->Q[0][1]).unitize();
    check2 = Ulna->x[0];

    Vector3 force_U_R = k * (check2 - check1);
    Radius->externalForces[0] += force_U_R;
    Ulna->externalForces[0] -= force_U_R;

    check3 = Radius->x[size_UR] +
             2.1 * 8.4 *
                 (-0.2 * Radius->Q[size_UR - 1][2] + Radius->Q[size_UR - 1][1])
                     .unitize();
    check4 = Ulna->x[size_UR];

    force_U_R = k * (check4 - check3);
    Radius->externalForces[size_UR] += force_U_R;
    Ulna->externalForces[size_UR] -= force_U_R;

    const REAL kk = 0.8e8;
    // Force3 —- Radius-Carp1
    const int size_Cp1 = Carp1->n;
    check1 =
        Radius->x[size_UR - 1] + 5.0 * (Radius->Q[size_UR - 2][2]).unitize();
    check2 = Carp1->x[0];

    Vector3 force_U_C = 10 * kk * (check2 - check1);
    Radius->externalForces[size_UR - 1] += force_U_C;
    Carp1->externalForces[0] -= force_U_C;

    check3 = check1 + 90 * (cos(0.7744) * Radius->Q[11][2] -
                            sin(0.7744) * Radius->Q[11][1])
                               .unitize();
    check4 = Carp1->x[size_Cp1];

    force_U_C = kk * (check4 - check3);
    Radius->externalForces[size_UR] += 0.1 * force_U_C;
    Carp1->externalForces[size_Cp1] -= 0.1 * force_U_C;

    // Force4 —- Carp1—Carp2
    check1 = Carp1->x[0] +
             1.7 * 8.4 *
                 (-sin(0.9718) * Carp1->Q[0][2] + cos(0.9718) * Carp1->Q[0][1])
                     .unitize();
    check2 = Carp2->x[0];

    Vector3 force_C_C = kk * (check2 - check1);
    Carp1->externalForces[0] += 0.1 * force_C_C;
    Carp2->externalForces[0] -= 0.1 * force_C_C;

    check3 =
        Carp1->x[size_Cp1] + 1.7 * 8.4 *
                                 (-sin(0.9718) * Carp1->Q[size_Cp1 - 1][2] +
                                  cos(0.9718) * Carp1->Q[size_Cp1 - 1][1])
                                     .unitize();
    check4 = Carp2->x[size_Cp1];

    force_C_C = kk * (check4 - check3);
    Carp1->externalForces[size_Cp1] += 0.1 * force_C_C;
    Carp2->externalForces[size_Cp1] -= 0.1 * force_C_C;

    // Force5 —- Carp1-Carp3
    const int size_Cp3 = Carp3->n;
    check1 = Carp1->x[size_Cp1 - 1] +
             8.4 * (-sin(0.9718) * Carp1->Q[size_Cp1 - 2][2] +
                    cos(0.9718) * Carp1->Q[size_Cp1 - 2][1])
                       .unitize();
    check2 = Carp3->x[0];

    Vector3 force_C_C3 = kk * (check2 - check1);
    Carp1->externalForces[size_Cp1 - 1] += force_C_C3;
    Carp3->externalForces[0] -= force_C_C3;

    check3 = check1 + 75 * (Carp1->Q[size_Cp1 - 2][2]).unitize();
    check4 = Carp3->x[size_Cp3];

    // We constrain the digit with a one way joint similar to SI Note 5.3.
    // The rationale follows from the fact that the digit is tiny and does
    // not play a significant effect on the dynamics
    force_C_C3 = kk * (check4 - check3);
    Carp3->externalForces[size_Cp3] -= 0.1 * force_C_C3;

    //-----------------------------------------------------------------------------------------------------------------
    // Connections Bone and Feathe
    // For simplicity, the orientation of the first
    // element in all rachis and barbs are constrained by boundary conditions,
    // not dynamic corrections. See SI note 5.3 for rationale.
    for (unsigned int i = 0; i < 10; i++) {
      // connection between rachis and bones
      /*
      We ignore the aerodynamic force on the rachis. In a biological wing
      structure, part of the rachis (close to the base) may not be subject to
      aeroforce since they are stick to the skin/flesh and will be covered by
      the top layers of feathers. Meanwhile, the actual radius of the rachis is
      smaller than what we model here, thus the actual aerodynamic force due to
      rachis is very small.
      */
      Rod *FearM = rodptrs[6 + i];
      if (i < 5) {
        check1 = Carp3->x[5 - i];
      } else {
        check1 = Carp2->x[10 - 2 * (i - 4)];
      }
      check2 = FearM->x[0];
      Vector3 force_B_F = k * (check2 - check1);
      if (i < 5) {
        Carp3->externalForces[5 - i] += force_B_F;
        FearM->externalForces[0] -= force_B_F;

        Matrix3 Temp =
            Matrix3(1.0, 0.0, 0.0, 0.0, cos((i * 7.7) / 180 * M_PI),
                    sin((i * 7.7) / 180 * M_PI), 0.0,
                    -sin((i * 7.7) / 180 * M_PI), cos((i * 7.7) / 180 * M_PI));

        FearM->Q[0] = Temp * Carp3->Q[4 - i];
      } else {
        Carp2->externalForces[10 - 2 * (i - 4)] += force_B_F;
        FearM->externalForces[0] -= force_B_F;

        Matrix3 Temp =
            Matrix3(1.0, 0.0, 0.0, 0.0, cos((i * 7.7) / 180 * M_PI),
                    sin((i * 7.7) / 180 * M_PI), 0.0,
                    -sin((i * 7.7) / 180 * M_PI), cos((i * 7.7) / 180 * M_PI));

        FearM->Q[0] = Temp * Carp2->Q[10 - 2 * (i - 4)];
      }

      // Connection between rachis and barbs
      /*
      To avoid numerical stiffness, the reacting connection forces from
      barbs to rachis are not applied directly on the rachis, instead they are
      lumped into one vector and applied on the bones at the node which the base
      of each rachis attaches to. Therefore, the inertial and aerodynamic force
      from the barbs are transfered to the bones and will effect the joint
      dynamics. See SI note 5.3 for detailed rationale. and justification for
      the same.
      */
      for (unsigned int j = 0; j < 77; j++) {
        const REAL xaxis = 70.0 + j * 4.0;
        const int index = (int)(xaxis / 57.5);
        const REAL distance = xaxis - index * 57.5;
        check1 = FearM->x[index] + distance * (FearM->Q[index][2]).unitize();
        Rod *Feather = rodptrs[16 + i * 192 + j * 2];

        check2 = Feather->x[0];
        Vector3 force_F_F = 0.1 * k * (check2 - check1);
        /*
        Due to the elasticity of the rachis, the force transfered from barbs to
        bones will experience a decay. We assume the force decays linearly with
        the axial location of the corresponding barb. We have validated this
        assumption in a separate study, given in SI note 5.3.
        */
        Vector3 featherforce = 1e-1 * (460.0 - xaxis) / 460.0 * force_F_F;

        // Here is where the loads due to the barbs are transferred at the base
        // of the rachis at the corresponding bone which is the carp. See
        // rationale in the SI note 5.3.
        if (i < 5) {
          Carp3->externalForces[5 - i] += featherforce;
        } else {
          Carp2->externalForces[10 - 2 * (i - 4)] += featherforce;
        }

        Feather->externalForces[0] -= 1e-1 * force_F_F;

        Vector3 temp2 = (FearM->Q[index][2] + FearM->Q[index][1] -
                         0.08 * (FearM->Q[index][0]))
                            .unitize();
        Vector3 temp0 = (FearM->Q[index][2] * temp2).unitize();
        Vector3 temp1 = (temp2 * temp0).unitize();

        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);
        for (unsigned int m = 1; m < 4; m++) {
          // The first element of each barb is fixed due to orientation
          // constraint. Thus doesn't contribute to lift. Again, rationale
          // follows from SI note 5.3.
          REAL factor = (time > 0.32) ? (1.2 - 3.0 * (time - 0.32)) : (1.2);
          factor = (factor < 0.9) ? (0.9) : (factor);
          factor = 1.2;
          const Vector3 Direction = -Feather->v[m].unitize();
          const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                    (Feather->v[m] % Feather->v[m]) * factor *
                                    (3.0 * Feather->l[m]) * Direction;
          Feather->externalForces[m] += AeroForce;
          Liftforce += AeroForce[2];
        }

        Feather = rodptrs[16 + i * 192 + j * 2 + 1];
        check2 = Feather->x[0];

        force_F_F = 0.1 * k * (check2 - check1);
        featherforce = 1e-1 * (460.0 - xaxis) / 460.0 * force_F_F;

        // Here is where the loads due to the barbs are transferred at the base
        // of the rachis at the corresponding bone which is the carp. See
        // rationale in the SI note 5.3.

        if (i < 5) {
          Carp3->externalForces[5 - i] += featherforce;
        } else {
          Carp2->externalForces[10 - 2 * (i - 4)] += featherforce;
        }

        Feather->externalForces[0] -= 1e-1 * force_F_F;

        temp2 = (FearM->Q[index][2] - FearM->Q[index][1] +
                 0.08 * (FearM->Q[index][0]))
                    .unitize();
        temp0 = (temp2 * FearM->Q[index][2]).unitize();
        temp1 = (temp2 * temp0).unitize();

        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);

        for (unsigned int m = 1; m < 4; m++) {
          REAL factor = (time > 0.32) ? (1.2 - 3.0 * (time - 0.32)) : (1.2);
          factor = (factor < 0.9) ? (0.9) : (factor);
          factor = 1.2;
          const Vector3 Direction = -Feather->v[m].unitize();
          const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                    (Feather->v[m] % Feather->v[m]) * factor *
                                    (3.0 * Feather->l[m]) * Direction;
          Feather->externalForces[m] += AeroForce;
          Liftforce += AeroForce[2];
        }

      }  // for j

      // Tip
      /*
      Here we ignore the aerodynamic load from the tip. The rationale is
      that in a biological feather, the tip barbs are much
      tinier and fewer, and have negligible contribution to lift force.
      However, they play an important role in creating tip vortices, that
      in the case of a pigeon display turbulent behavior. We ignore such
      effects, as explained in SI Note 5.2.
      */

      for (unsigned int j = 0; j < 19; j++) {
        check1 = FearM->x[6] + 33.0 * (FearM->Q[6][2]).unitize();

        Rod *Feather = rodptrs[16 + i * 192 + 154 + j * 2];
        check2 = Feather->x[0];

        Vector3 force_F_F = 0.1 * k * (check2 - check1);
        Feather->externalForces[0] -= 1e-1 * force_F_F;

        const REAL angle = (90.0 - j * 5.0) / 180 * M_PI;
        Vector3 Unit =
            cos(angle) * FearM->Q[6][2] + sin(angle) * FearM->Q[6][1];

        Vector3 temp2 =
            (FearM->Q[6][2] + Unit - 0.08 * (19 - j) / 19 * (FearM->Q[6][0]))
                .unitize();
        Vector3 temp0 = (FearM->Q[6][2] * temp2).unitize();
        Vector3 temp1 = (temp2 * temp0).unitize();

        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);
        for (unsigned int m = 1; m < 4; m++) {
          const Vector3 Direction = -Feather->v[m].unitize();
          const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                    (Feather->v[m] % Feather->v[m]) * 1.2 *
                                    (3.0 * Feather->l[m]) * Direction;
          Feather->externalForces[m] += AeroForce;
          Liftforce += AeroForce[2];
        }

        Feather = rodptrs[16 + i * 192 + 154 + j * 2 + 1];
        check2 = Feather->x[0];

        force_F_F = 0.1 * k * (check2 - check1);
        Feather->externalForces[0] -= 1e-1 * force_F_F;

        Unit = cos(angle) * FearM->Q[6][2] - sin(angle) * FearM->Q[6][1];

        temp2 =
            (FearM->Q[6][2] + Unit + 0.08 * (19 - j) / 19 * (FearM->Q[6][0]))
                .unitize();
        temp0 = (temp2 * FearM->Q[6][2]).unitize();
        temp1 = (temp2 * temp0).unitize();
        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);
        for (unsigned int m = 1; m < 4; m++) {
          const Vector3 Direction = -Feather->v[m].unitize();
          const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                    (Feather->v[m] % Feather->v[m]) * 1.2 *
                                    (3.0 * Feather->l[m]) * Direction;
          Feather->externalForces += AeroForce;
          Liftforce += AeroForce[2];
        }
      }  // for j
    }    // for i

    //-----------------------------------------------------------------------------------------------------------------
    // Connections Bone and Feather
    // Secondary Remige
    const int Size = rodptrs.size();
    for (unsigned int i = 0; i < 9; i++) {
      // Rachis and bones
      Rod *FearM = rodptrs[6 + 1930 + i];  // Full
      if (i < 7) {
        check1 = Ulna->x[11 - i];
      } else {
        check1 = Ulna->x[3 - 2 * (i - 7)];
      }
      check2 = FearM->x[0];

      Vector3 force_B_F = k * (check2 - check1);
      if (i < 7) {
        Ulna->externalForces[11 - i] += force_B_F;
        FearM->externalForces[0] -= force_B_F;

        Matrix3 Temp =
            Matrix3(1.0, 0.0, 0.0, 0.0, cos((122 + i * 6.9) / 180 * M_PI),
                    sin((122 + i * 6.9) / 180 * M_PI), 0.0,
                    -sin((122 + i * 6.9) / 180 * M_PI),
                    cos((122 + i * 6.9) / 180 * M_PI));

        FearM->Q[0] = Temp * Ulna->Q[11 - i];
      } else {
        Ulna->externalForces[3 - 2 * (i - 7)] += force_B_F;
        FearM->externalForces[0] -= force_B_F;

        Matrix3 Temp =
            Matrix3(1.0, 0.0, 0.0, 0.0, cos((122 + i * 6.9) / 180 * M_PI),
                    sin((122 + i * 6.9) / 180 * M_PI), 0.0,
                    -sin((122 + i * 6.9) / 180 * M_PI),
                    cos((122 + i * 6.9) / 180 * M_PI));

        FearM->Q[0] = Temp * Ulna->Q[3 - 2 * (i - 7)];
      }
      for (unsigned int m = 0; m < 8; m++) {
        const Vector3 Direction = -FearM->v[m].unitize();
        const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                  (FearM->v[m] % FearM->v[m]) * 1.2 *
                                  (5.0 * FearM->l[0]) * Direction;
        Liftforce += AeroForce[2];
      }

      // Rachis and Barbs
      for (unsigned int j = 0; j < (77 - i * 7); j++) {
        const REAL xaxis = 70.0 + j * 4.0;
        const REAL featherlength = (1.0 - i / 15.0) * 460.0;
        const int index = (int)(xaxis / (57.5 * (1.0 - i / 15.0)));
        const REAL distance = xaxis - index * 57.5 * (1.0 - i / 15.0);
        check1 = FearM->x[index] + distance * (FearM->Q[index][2]).unitize();

        Rod *Feather =
            rodptrs[1945 + (i * (192 + (192 - (i - 1) * 14))) / 2 + 2 * j];

        check2 = Feather->x[0];
        Vector3 force_F_F = 0.1 * k * (check2 - check1);
        Vector3 featherforce =
            1e-1 * (featherlength - xaxis) / featherlength * force_F_F;

        // Here is where the loads due to the barbs are transferred at the base
        // of the rachis at the corresponding bone which is the ulna. See
        // rationale in the SI note 5.3

        if (i < 7) {
          Ulna->externalForces[11 - i] += featherforce;
        } else {
          Ulna->externalForces[3 - 2 * (i - 7)] += featherforce;
        }
        Feather->externalForces[0] -= 1e-1 * force_F_F;

        Vector3 temp2 = (FearM->Q[index][2] + FearM->Q[index][1] -
                         0.08 * (FearM->Q[index][0]))
                            .unitize();
        Vector3 temp0 = (FearM->Q[index][2] * temp2).unitize();
        Vector3 temp1 = (temp2 * temp0).unitize();

        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);
        if (j > 0) {
          for (unsigned int m = 1; m < 4; m++) {
            REAL factor = (time > 0.32) ? (1.2 - 3.0 * (time - 0.32)) : (1.2);
            factor = (factor < 0.9) ? (0.9) : (factor);
            factor = 1.2;
            const Vector3 Direction = -Feather->v[m].unitize();
            const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                      (Feather->v[m] % Feather->v[m]) * factor *
                                      (3.0 * Feather->l[m]) * Direction;
            Feather->externalForces[m] += AeroForce;
            Liftforce += AeroForce[2];
          }
        }

        Feather =
            rodptrs[1945 + (i * (192 + (192 - (i - 1) * 14))) / 2 + 2 * j + 1];
        check2 = Feather->x[0];

        force_F_F = 0.1 * k * (check2 - check1);
        featherforce =
            1e-1 * (featherlength - xaxis) / featherlength * force_F_F;
        // Here is where the loads due to the barbs are transferred at the base
        // of the rachis at the corresponding bone which is the ulna. See
        // rationale in the SI note 5.3

        if (i < 7) {
          Ulna->externalForces[11 - i] += featherforce;
        } else {
          Ulna->externalForces[3 - 2 * (i - 7)] += featherforce;
        }

        Feather->externalForces[0] -= 1e-1 * force_F_F;

        temp2 = (FearM->Q[index][2] - FearM->Q[index][1] +
                 0.08 * (FearM->Q[index][0]))
                    .unitize();
        temp0 = (temp2 * FearM->Q[index][2]).unitize();
        temp1 = (temp2 * temp0).unitize();

        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);

        if (j > 0) {
          for (unsigned int m = 1; m < 4; m++) {
            REAL factor = (time > 0.32) ? (1.2 - 3.0 * (time - 0.32)) : (1.2);
            factor = (factor < 0.9) ? (0.9) : (factor);
            factor = 1.2;
            const Vector3 Direction = -Feather->v[m].unitize();
            const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                      (Feather->v[m] % Feather->v[m]) * factor *
                                      (3.0 * Feather->l[m]) * Direction;
            Feather->externalForces[m] += AeroForce;
            Liftforce += AeroForce[2];
          }
        }
      }

      // Tip feathers
      for (unsigned int j = 0; j < 19; j++) {
        const REAL xaxis = 378.0 - 28.0 * i;
        const int index = (int)(xaxis / (57.5 * (1.0 - i / 15.0)));
        const REAL distance = xaxis - index * 57.5 * (1.0 - i / 15.0);
        check1 = FearM->x[index] + distance * (FearM->Q[index][2]).unitize();

        Rod *Feather = rodptrs[1945 + (i * (192 + (192 - (i - 1) * 14))) / 2 +
                               (77 - i * 7) * 2 + 2 * j];
        check2 = Feather->x[0];

        Vector3 force_F_F = 0.1 * k * (check2 - check1);
        Feather->externalForces[0] -= force_F_F;

        const REAL angle = (90.0 - j * 5.0) / 180 * M_PI;
        Vector3 Unit =
            cos(angle) * FearM->Q[index][2] + sin(angle) * FearM->Q[index][1];

        Vector3 temp2 = (FearM->Q[index][2] + Unit -
                         0.08 * (19 - j) / 19 * (FearM->Q[index][0]))
                            .unitize();
        Vector3 temp0 = (FearM->Q[index][2] * temp2).unitize();
        Vector3 temp1 = (temp2 * temp0).unitize();

        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);
        for (unsigned int m = 1; m < 4; m++) {
          const Vector3 Direction = -Feather->v[m].unitize();
          const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                    (Feather->v[m] % Feather->v[m]) * 1.2 *
                                    (3.0 * Feather->l[m]) * Direction;
          Feather->externalForces += AeroForce;
          Liftforce += AeroForce[2];
        }
        Feather = rodptrs[1945 + (i * (192 + (192 - (i - 1) * 14))) / 2 +
                          (77 - i * 7) * 2 + 2 * j + 1];
        check2 = Feather->x[0];

        force_F_F = 0.1 * k * (check2 - check1);
        Feather->externalForces[0] -= force_F_F;

        Unit =
            cos(angle) * FearM->Q[index][2] - sin(angle) * FearM->Q[index][1];

        temp2 = (FearM->Q[index][2] + Unit +
                 0.08 * (19 - j) / 19 * (FearM->Q[index][0]))
                    .unitize();
        temp0 = (temp2 * FearM->Q[index][2]).unitize();
        temp1 = (temp2 * temp0).unitize();
        Feather->Q[0] =
            Matrix3(temp0[0], temp0[1], temp0[2], temp1[0], temp1[1], temp1[2],
                    temp2[0], temp2[1], temp2[2]);
        for (unsigned int m = 1; m < 4; m++) {
          const Vector3 Direction = -Feather->v[m].unitize();
          const Vector3 AeroForce = 0.5 * (1.225e-6) *
                                    (Feather->v[m] % Feather->v[m]) * 1.2 *
                                    (3.0 * Feather->l[m]) * Direction;
          Feather->externalForces += AeroForce;
          Liftforce += AeroForce[2];
        }
      }
    }  // for i
  }
#endif
};
#endif
