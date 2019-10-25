/*
Simple Rod-Rod connection for toy examples
*/
#include "GeometryFunctions.h"
#include "Matrix3.h"
#include "Rod.h"
#include "UsualHeaders.h"
#include "Vector3.h"

using namespace std;

class SimpleConnection {
 protected:
  typedef std::vector<Rod *> Vrodptr;
  Vrodptr rodptrs;

 public:
  SimpleConnection(Vrodptr rodptrs) : rodptrs(rodptrs) {}
  ~SimpleConnection(){};

//-----------------------------------------------------------------------------------------------------------------
/*This function exemplifies how to connect rods with a spherical joint.
In this example, we will connect the last node of link one with the first node
of link two. Link one is completely fixed, while link two is allowed to rotate
in any direction.
*/
#ifdef FLAGSPHERICALJOINT
  void RodRodSimpleConnection(const REAL time) {
    // Connection between link one and link two
    // set up check points
    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);

    // stiffness of the joint connection -- tuned empirically
    const REAL k = 1e8;

    // Select rods
    Rod *One = rodptrs[0];
    Rod *Two = rodptrs[1];

    // Force —- Rod One-Two Connection
    const int size = One->n;

    // The checks points are the nodes sopposed to be connected.
    // Here, they are the last node of rod one and first node of rod two
    check1 = One->x[size];
    check2 = Two->x[0];

    // Compute their position difference, and then connecting force
    Vector3 force = k * (check2 - check1);
    // Force will be applied on each node with opposite direction to bring them
    // together.
    One->externalForces[size] += force;
    Two->externalForces[0] -= force;

    // Some load
    // This is the external load to test out our joint connection.
    const REAL Amp =
        5e3;  // 5 mN --- a very small force since we assume no gravity
    // Forces are applied on the last node of link two.
    // Initially, force is applied in negative z direction
    if (time < 0.2) {
      Two->externalForces[size] += Vector3(0.0, 0.0, -2 * Amp);
    }
    // After a short time period, force will remain in x-z plane,
    // but start changing direction to rotate link two around link one.
    else {
      Two->externalForces[size] +=
          Vector3(Amp * cos(0.5 * M_PI * (time - 0.2)), 0.0,
                  Amp * sin(0.5 * M_PI * (time - 0.2)));
    }
  }
#endif

//-----------------------------------------------------------------------------------------------------------------
/*This function exemplifies how to connect rods with a hinge joint.
A hinge connection is based on a spherical joint with additional constraint to
allow only in plane rotation. In this example, we will still connect the last
node of link one with the first node of link two, and link two is only allowed
to rotation in y-z plane.
*/
#ifdef FLAGHINGEJOINT
  void RodRodSimpleConnection(const REAL time) {
    // Connection between link one and link two
    // Same with the spherical joint connection
    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);

    const REAL k = 1e8;

    Rod *One = rodptrs[0];
    Rod *Two = rodptrs[1];

    // Force —- Rod One-Two Connection
    const int size = One->n;
    check1 = One->x[size];
    check2 = Two->x[0];

    Vector3 force = k * (check2 - check1);
    One->externalForces[size] += force;
    Two->externalForces[0] -= force;

    // Additional in-plane constraint through restoring torque
    // Stiffness of the restoring constraint -- tuned empirically
    const REAL kt = 5e6;
    // Find the normal direction of the contraint plane (y-z plane)
    const Vector3 normaldirection = Vector3(1.0, 0.0, 0.0);
    // Current direction of the first element of link two
    const Vector3 linkdirection = Two->x[1] - Two->x[0];
    // Projection of the linkdirection onto the plane normal
    const Vector3 Forcedirection =
        -(linkdirection % normaldirection) * normaldirection;
    // Compute the restoring torque
    Vector3 Torque = kt * linkdirection * Forcedirection;
    // The opposite torque will be applied on link one (no effect in this case
    // since link one is completely fixed)
    One->externalTorques[size - 1] -= One->Q[size - 1] * Torque;
    Two->externalTorques[0] += Two->Q[0] * Torque;

    // Some load
    // Same load as the previous case.
    const REAL Amp = 5e3;
    if (time < 0.2) {
      Two->externalForces[size] += Vector3(0.0, 0.0, -2 * Amp);
    } else {
      Two->externalForces[size] +=
          Vector3(Amp * cos(0.5 * M_PI * (time - 0.2)), 0.0,
                  Amp * sin(0.5 * M_PI * (time - 0.2)));
    }
  }
#endif

//-----------------------------------------------------------------------------------------------------------------
/*This function exemplifies how to connect rods with a fixed joint.
The concept of fixed connection is same with a hinge joint, but now the
restoring torque will resist any motion that changes the orientation of link
two. In this example, we will still connect the last node of link one with the
first node of link two, and link two is constrained to have the same orientation
with link one.
*/
#ifdef FLAGFIXEDJOINT
  void RodRodSimpleConnection(const REAL time) {
    // Connection between link one and link two
    // Same with the spherical joint connection
    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);
    Vector3 check3 = Vector3(0.0, 0.0, 0.0);
    Vector3 check4 = Vector3(0.0, 0.0, 0.0);

    const REAL k = 1e8;

    Rod *One = rodptrs[0];
    Rod *Two = rodptrs[1];

    // Force —- Rod One-Two Connection
    const int size = One->n;
    check1 = One->x[size];
    check2 = Two->x[0];

    Vector3 force = k * (check2 - check1);
    One->externalForces[size] += force;
    Two->externalForces[0] -= force;

    // Additional orientation constraint through restoring torque
    // Stiffness of the restoring constraint -- tuned empirically
    const REAL kt = 5e6;
    // Segment length of link two
    const REAL dl = 200.0 / 10.0;
    // Current direction of the first element of link two
    const Vector3 linkdirection = Two->x[1] - Two->x[0];
    /*To constrain the orientation of link two, the second node of link two
    should align with the direction of link one. Thus, we compute the desired
    position of the second node of link two as check 3, and the current
    position of the second node of link two as check 4. Check 3 and check 4
    should overlap.*/
    check3 = One->x[size] + dl * One->edge[size - 1].unitize();
    check4 = Two->x[1];
    // Compute the restoring torque
    const Vector3 Forcedirection = -kt * (check4 - check3);
    Vector3 Torque = linkdirection * Forcedirection;
    // The opposite torque will be applied on link one (no effect in this case
    // since link one is completely fixed)
    One->externalTorques[size - 1] -= One->Q[size - 1] * Torque;
    Two->externalTorques[0] += Two->Q[0] * Torque;

    // Some load
    // Same load as the previous case.
    const REAL Amp = 5e3;
    if (time < 0.2) {
      Two->externalForces[size] += Vector3(0.0, 0.0, -2 * Amp);
    } else {
      Two->externalForces[size] +=
          Vector3(Amp * cos(0.5 * M_PI * (time - 0.2)), 0.0,
                  Amp * sin(0.5 * M_PI * (time - 0.2)));
    }
  }
#endif

//-----------------------------------------------------------------------------------------------------------------
// This function exemplifies how to connect muscle with substrate, and how to
// setup muscle force.
#ifdef FLAGPULLINGMUSCLE
  void RodRodSimpleConnection(const REAL time) {
    // Connection
    Vector3 check1 = Vector3(0.0, 0.0, 0.0);
    Vector3 check2 = Vector3(0.0, 0.0, 0.0);
    Vector3 check3 = Vector3(0.0, 0.0, 0.0);
    Vector3 check4 = Vector3(0.0, 0.0, 0.0);

    // stiffness of the joint connection -- tuned empirically
    const REAL k = 1e8;

    // One is the left side substrate and Two is right side substrate that the
    // muscle will be connected to and pulling against
    Rod *One = rodptrs[0];
    Rod *Two = rodptrs[1];
    Rod *Muscle = rodptrs[2];

    // Spherical joint for both ends -- same implementation as before
    // Left side
    const int size = One->n;
    check1 = One->x[6];
    check2 = Muscle->x[0];

    Vector3 force = k * (check2 - check1);
    One->externalForces[6] += force;
    Muscle->externalForces[0] -= force;

    // Right side
    check3 = Two->x[6];
    check4 = Muscle->x[10];

    force = k * (check4 - check3);
    Two->externalForces[6] += force;
    Muscle->externalForces[10] -= force;

    // Muscle Force
    // Amplitude of the muscle force -- 2N
    const REAL Amp = 2e6;
    // Frequency of the muscle actuation -- 1Hz after taking absolute value of
    // sin
    const REAL w = 2.0 * M_PI / 2.0;

    const int size_m = Muscle->n;
    const int Tendon_n = 3;
    const int Muscle_n = size_m - 2 * Tendon_n;
    /*Muscle force will be applied only on the muscle segments
    Here we assume a simple sinusoidal signal for muscle activation, but it can
    be replaced by any experimental/ theoretical muscle functions. Muscle forces
    are applied at the two ends of each muscle segment, with opposite direction.
    Thus the muscle contraction will be uniform through out the entire muscle.
    No force on tendons.
    */
    for (unsigned int i = Tendon_n; i < (Tendon_n + Muscle_n); i++) {
      Muscle->externalForces[i] +=
          Amp * (Muscle->edge[i].unitize()) * abs(sin(w * time));
    }
    for (unsigned int i = (Tendon_n + 1); i < (Tendon_n + Muscle_n + 1); i++) {
      Muscle->externalForces[i] +=
          -1 * Amp * (Muscle->edge[i - 1].unitize()) * abs(sin(w * time));
    }
  }
#endif
};