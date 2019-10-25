#ifndef RODBOUNDARYCONDITIONS_H_
#define RODBOUNDARYCONDITIONS_H_

#include "ArithmeticPrecision.h"
#include "Matrix3.h"
#include "Rod.h"
#include "Tolerance.h"
#include "UsualHeaders.h"
#include "Vector3.h"

using namespace std;

// Boundary condition base class: strict separation between dirichlet and
// neumann boundary conditions since they are enforced at different location
// during symplectic integration
class RodBC {
 protected:
  typedef enum { DIRICHLET, NEUMANN } TypeBC;
  TypeBC status;

 public:
  RodBC() : status(TypeBC::NEUMANN) {}
  virtual ~RodBC() {}
  virtual void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                         const REAL time) = 0;
  virtual void neumann(Rod &rod, const unsigned int step, const REAL dt,
                       const REAL time) = 0;

  virtual void operator()(Rod &rod, const unsigned int step, const REAL dt,
                          const REAL time) {
    switch (status) {
      case TypeBC::NEUMANN:
        neumann(rod, step, dt, time);
        status = TypeBC::DIRICHLET;
        break;

      case TypeBC::DIRICHLET:
        dirichlet(rod, step, dt, time);
        status = TypeBC::NEUMANN;
        break;

      default:
        break;
    };
  }
};

// Pack multiple boundary conditions together
class MultipleBC : public RodBC {
 public:
  vector<RodBC *> bcs;

  MultipleBC(vector<RodBC *> &bcs) : RodBC(), bcs(bcs) {}

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time){};
  void neumann(Rod &rod, const unsigned int step, REAL dt, const REAL time){};

  void operator()(Rod &rod, const unsigned int step, const REAL dt,
                  const REAL time) {
    assert(step >= 0);

    for (unsigned int i = 0; i < bcs.size(); i++)
      (*bcs[i])(rod, step, dt, time);
  }
};

// No external forces, moments, or constraints
class FreeBC : public RodBC {
 public:
  FreeBC() : RodBC() {}

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time){};
  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time){};
};

class LongitudinalWavesBC : public RodBC {
 public:
  const REAL A;
  const REAL omega;
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  LongitudinalWavesBC(vector<Rod *> &rodptrs, const REAL _A, const REAL _omega)
      : RodBC(), A(_A), omega(_omega) {
    assert(rodptrs.size() == 1);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp first end position
    rod.x.front() = startX;

    // Oscillate second end
    rod.x.back() = endX + direction * A * sin(omega * time);

    // Clamp both ends material frame
    rod.Q.front() = startQ;
    rod.Q.back() = endQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set initial point velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set initial point angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Let last point slide in the x-direction
    rod.v.back().x =
        A * omega *
        cos(omega * time);  // Derivative in time of the displacement
    rod.v.back().y = 0.0;
    rod.v.back().z = 0.0;

    // Set last point angular velocity to zero
    rod.w.back().x = 0.0;
    rod.w.back().y = 0.0;
    rod.w.back().z = 0.0;
  }
};

class LongitudinalWavesLoadBC : public RodBC {
 public:
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  LongitudinalWavesLoadBC(vector<Rod *> &rodptrs) : RodBC() {
    assert(rodptrs.size() == 1);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp first end position
    rod.x.front() = startX;

    rod.x.front().y = 0.0;
    rod.x.front().z = 0.0;

    // Clamp both ends material frame
    rod.Q.front() = startQ;
    rod.Q.back() = endQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set initial point velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set initial point angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Let last point slide in the x-direction
    rod.v.back().y = 0.0;
    rod.v.back().z = 0.0;

    // Set last point angular velocity to zero
    rod.w.back().x = 0.0;
    rod.w.back().y = 0.0;
    rod.w.back().z = 0.0;
  }
};

class TorsionallWavesBC : public RodBC {
 public:
  const REAL angularAmplitude;
  const REAL omega;
  Vector3 axisOfRotInMaterialFrame;
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  TorsionallWavesBC(vector<Rod *> &rodptrs, const REAL _angularAmplitude,
                    const REAL _omega)
      : RodBC(), angularAmplitude(_angularAmplitude), omega(_omega) {
    assert(rodptrs.size() == 1);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();

    // Axis of rotation in material frame
    axisOfRotInMaterialFrame = (endQ * direction).unitize();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp both ends positions
    rod.x.front() = startX;
    rod.x.back() = endX;

    // Clamp first end material frame
    rod.Q.front() = startQ;

    // Oscillate second hand material frame
    rod.Q.back() =
        exp(angularAmplitude * sin(omega * time) * axisOfRotInMaterialFrame) *
        endQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set first end velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set first end angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Set second end velocities to zero
    rod.v.back().x = 0.0;
    rod.v.back().y = 0.0;
    rod.v.back().z = 0.0;

    // Set last point angular velocity
    rod.w.back() =
        angularAmplitude * omega * cos(omega * time) *
        axisOfRotInMaterialFrame;  // Derivative in time of the displacement
  }
};

class TorsionallWavesCoupleBC : public RodBC {
 public:
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  TorsionallWavesCoupleBC(vector<Rod *> &rodptrs) : RodBC() {
    assert(rodptrs.size() == 1);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp both ends positions
    rod.x.front() = startX;
    rod.x.back() = endX;

    // Clamp first end material frame
    rod.Q.front() = startQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set first end velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set first end angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Set second end velocities to zero
    rod.v.back().x = 0.0;
    rod.v.back().y = 0.0;
    rod.v.back().z = 0.0;

    // Allow twist only at the second end
    // rod.w.front().x = 0.0;
    // rod.w.front().y = 0.0;
  }
};

class TimoshenkoBeamBC : public RodBC {
 public:
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  TimoshenkoBeamBC(vector<Rod *> &rodptrs) : RodBC() {
    assert(rodptrs.size() == 1);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp first ends position
    rod.x.front() = startX;

    // Clamp first end material frame
    rod.Q.front() = startQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set first end velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set first end angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;
  }
};

class TorsionallWavesCoupleStretchBC : public RodBC {
 public:
  const REAL extension;
  const REAL timeExtension;
  Vector3 extensionVel;
  Vector3 direction;
  Vector3 startX;
  Vector3 endXBeginning;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  TorsionallWavesCoupleStretchBC(vector<Rod *> &rodptrs, const REAL _extension,
                                 const REAL _timeExtension)
      : RodBC(), extension(_extension), timeExtension(_timeExtension) {
    assert(rodptrs.size() == 1);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endXBeginning = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endXBeginning - startX).unitize();

    // Final position
    endX = startX + extension * (endXBeginning - startX);

    extensionVel = (endX - endXBeginning) / timeExtension;

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp first end
    rod.x.front() = startX;

    // Clamp second end after time extension
    if (time > timeExtension) rod.x.back() = endX;

    // Clamp first end material frame
    rod.Q.front() = startQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set first end velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set first end angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Allow twist only at the second end
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;

    if (time > timeExtension) {
      // Set second end velocities to zero
      rod.v.back().x = 0.0;
      rod.v.back().y = 0.0;
      rod.v.back().z = 0.0;
    } else {
      // Set second end velocities to zero
      rod.v.back() = extensionVel;
    }
  }
};

// Clamped on one end, longitudinal displacement only allowed on the other end
class StretchReleaseBC : public RodBC {
 public:
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  StretchReleaseBC(Rod *rodptr) {
    // Initial positions of the two end of the rod
    startX = rodptr->x.front();
    endX = rodptr->x.back();

    // Direction of the rod
    // direction = (endX-startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptr->Q.front();
    endQ = rodptr->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp first end position
    rod.x.front() = startX;

    // Clamp both ends material frame
    rod.Q.front() = startQ;
    // rod.Q.back() = endQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set initial point velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set initial point angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Let last point slide in the x-direction
    // rod.v.back().y = 0.0;
    // rod.v.back().z = 0.0;

    // Set last point angular velocity to zero
    // rod.w.back().x = 0.0;
    // rod.w.back().y = 0.0;
    // rod.w.back().z = 0.0;
  }
};

// Restrains a rod's final point to not move in the x or y directions
// restrains a rod's initial point to not move at all
// for Euler Buckling testing
class RestrictEndBC : public RodBC {
 public:
  RestrictEndBC() : RodBC() {}

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Completely fix one end
    rod.x[0].x = 0;
    rod.x[0].y = 0;
    rod.x[0].z = 0;

    // Let the other end slide in the z-direction
    rod.x[rod.n].x = 0;
    rod.x[rod.n].y = 0;
  };

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Completely fix one end
    rod.v[0].x = 0;
    rod.v[0].y = 0;
    rod.v[0].z = 0;

    // Let the other end slide in the z-direction
    rod.v[rod.n].x = 0;
    rod.v[rod.n].y = 0;
  };
};

// Hinge fix_start
class HingeBC : public RodBC {
 public:
  Vector3 startX;

  HingeBC(Rod *rodptr) : RodBC() { startX = rodptr->x.front(); }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Completely fix one end
    rod.x.front() = startX;
  };

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Completely fix one end
    rod.v[0].x = 0;
    rod.v[0].y = 0;
    rod.v[0].z = 0;
  };
};
// Hinge fix_end
class HingeBC_End : public RodBC {
 public:
  Vector3 endX;

  HingeBC_End(Rod *rodptr) : RodBC() { endX = rodptr->x.back(); }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Completely fix one end
    rod.x.back() = endX;
  };

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Completely fix one end
    const int size = rod.n;
    rod.v[size].x = 0;
    rod.v[size].y = 0;
    rod.v[size].z = 0;
  };
};

// Flapping end
class FlappingBC : public RodBC {
 public:
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Matrix3 startQ;
  Matrix3 endQ;

  FlappingBC(Rod *rodptr) : RodBC() {
    // Initial material frames of the two end of the rod
    startQ = rodptr->Q.front();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Clamp first ends position
    // const Vector3 frontposition =
    // Vector3(0.0,0.0,0.01*sin(120*2.0*M_PI*time)); rod.x.front() =
    // frontposition;
    // Clamp first end material frame
    rod.Q.front() = startQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    REAL frontvelocity = 0.0;

    if (time < 0.2) {
      frontvelocity = 1e-1;
    } else if (time < 0.21) {
      frontvelocity = 1e-1 * (1.0 - 2.0 * (time - 0.2) / 0.01);
    } else if (time < 0.6) {
      frontvelocity = -1e-1;
    } else if (time < 0.61) {
      frontvelocity = -1e-1 * (1.0 - 2.0 * (time - 0.6) / 0.01);
    } else {
      frontvelocity = 1e-1;
    }

    // Set first end velocities to zero
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = frontvelocity;

    // Set first end angular velocity to zero
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;
  }
};

// constant angular velocity
class ConAVBC : public RodBC {
 public:
  int SwitchFlag = 0;
  const REAL setting = 5.236;
  ConAVBC() : RodBC() {}

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time){};

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // const Vector3 direction = (rod.x[15]-rod.x[0]).unitize();
    // const REAL angle_ = angle(direction, Vector3(0.0,0.0,-1.0));

    if (-1.0 * rod.w[2][0] > setting) SwitchFlag = 1;
    if (SwitchFlag == 1) {
      rod.w[2] = Vector3(-1.0 * setting, 0.0, 0.0);
      //  if(angle_ > (40.0/180.0*M_PI))
      //    rod.w[2] = Vector3(-1.0*setting,0.0,0.0) + (setting/0.349)*(angle_
      //    - 40.0/180.0*M_PI)*Vector3(1.0,0.0,0.0);
    }
    if (rod.w[2][0] > 0) rod.w[2] = Vector3(0.0, 0.0, 0.0);
  };
};

// one end completely fixed
class FixOneBC : public RodBC {
 public:
  Vector3 startX;
  Matrix3 startQ;

  FixOneBC(Rod *rodptr) : RodBC() {
    startX = rodptr->x.front();
    startQ = rodptr->Q.front();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Completely fix both ends
    rod.x.front() = startX;
    rod.Q.front() = startQ;
  };

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Completely fix both ends
    rod.v[0].x = 0;
    rod.v[0].y = 0;
    rod.v[0].z = 0;
  };
};

class FixedBC : public RodBC {
 public:
  Vector3 startX;
  Vector3 endX;

  FixedBC(Rod *rodptr) : RodBC() {
    startX = rodptr->x.front();
    endX = rodptr->x.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // Completely fix both ends
    rod.x.front() = startX;
    rod.x.back() = endX;
  };

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Completely fix both ends
    rod.v[0].x = 0;
    rod.v[0].y = 0;
    rod.v[0].z = 0;
    const int size = rod.n;
    rod.v[size].x = 0;
    rod.v[size].y = 0;
    rod.v[size].z = 0;
  };
};

// Restrains a rod's final point to not move in the x or y directions
// restrains a rod's initial point to not move at all
// for Euler Buckling testing
class HelicalBucklingWithForcesBC : public RodBC {
 public:
  HelicalBucklingWithForcesBC() : RodBC() {}

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time){};

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Let the first end slide in the z-direction
    rod.v[0].x = 0;
    rod.v[0].y = 0;

    // Let the other end slide in the z-direction
    rod.v[rod.v.size() - 1].x = 0;
    rod.v[rod.v.size() - 1].y = 0;

    // Let the first end rotate only about the z-direction
    rod.w[0].x = 0;
    rod.w[0].y = 0;

    // Let the other end rotate only about the z-direction
    rod.w[rod.w.size() - 1].x = 0;
    rod.w[rod.w.size() - 1].y = 0;
  };
};

class HelicalBucklingBC : public RodBC {
 public:
  const REAL twistingTime;  // Time during which the twisting is performed
                            // before relaxing the system
  const REAL D;
  const REAL R;
  Vector3 angVel;
  Vector3 shrinkVel;
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Vector3 finalStartX;
  Vector3 finalEndX;
  Matrix3 startQ;
  Matrix3 endQ;
  Matrix3 finalStartQ;
  Matrix3 finalEndQ;

  HelicalBucklingBC(vector<Rod *> &rodptrs, const REAL _twistingTime,
                    const REAL _D, const REAL _R)
      : RodBC(), twistingTime(_twistingTime), D(_D), R(_R) {
    assert(rodptrs.size() == 1);

    // Sanity checks
    assert(twistingTime > 0);
    assert(R > 0);
    assert(D > 0);

    // Compute rotation and shrinkage velocities. Note that I divide by 2
    // because I am going to apply the rotation at both ends of the rod
    const REAL angVelScalar = (R * 2.0 * M_PI / twistingTime) / 2.0;
    const REAL shrinkVelScalar = (D / twistingTime) / 2.0;

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Final position due to displacement boundary conditions
    finalStartX = startX + D / 2.0 * direction;
    finalEndX = endX - D / 2.0 * direction;

    // Note that in this case the axis of rotation and shrinkage coincides with
    // the one of the material frame but it is not always the case!
    angVel = angVelScalar * direction;
    shrinkVel = shrinkVelScalar * direction;

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();

    // Final material frames of the two end of the rod
    finalStartQ = exp(R * M_PI * direction) * startQ;
    finalEndQ = exp(-R * M_PI * direction) * endQ;
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    if (time > twistingTime) {
      rod.x.front() = finalStartX;
      rod.x.back() = finalEndX;

      rod.Q.front() = finalStartQ;
      rod.Q.back() = finalEndQ;
    }
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    if (time > twistingTime) {
      // Set initial point velocities to zero
      rod.v.front().x = 0.0;
      rod.v.front().y = 0.0;
      rod.v.front().z = 0.0;

      // Set initial point angular velocity to zero
      rod.w.front().x = 0.0;
      rod.w.front().y = 0.0;
      rod.w.front().z = 0.0;

      // Set last point velocities to zero
      rod.v.back().x = 0;
      rod.v.back().y = 0;
      rod.v.back().z = 0;

      // Set last point angular velocity to zero
      rod.w.back().x = 0.0;
      rod.w.back().y = 0.0;
      rod.w.back().z = 0.0;
    } else {
      // Set initial point velocities
      rod.v.front() = shrinkVel;

      // Set initial point angular velocity
      rod.w.front() = angVel;

      // Set last point velocities
      rod.v.back() = -shrinkVel;

      // Set last point angular velocity
      rod.w.back() = -angVel;
    }
  }
};

class QuasiStaticEndpointBC : public RodBC {
 public:
  vector<Vector3> xStart;  // position of x[0] in time
  vector<Vector3> xEnd;    // position of x[n] in time
  vector<Matrix3> QStart;  // material frame at x[0] in time
  vector<Matrix3> QEnd;    // material frame at x[n] in time

  QuasiStaticEndpointBC(vector<Vector3> xStart, vector<Vector3> xEnd,
                        vector<Matrix3> QStart, vector<Matrix3> QEnd)
      : RodBC(), xStart(xStart), xEnd(xEnd), QStart(QStart), QEnd(QEnd) {
    // Sanity checks
    assert(xStart.size() > 0);
    assert(xEnd.size() > 0);
    assert(QStart.size() > 0);
    assert(QEnd.size() > 0);

    assert(xStart.size() == xEnd.size());
    assert(xStart.size() == QStart.size());
    assert(xStart.size() == QEnd.size());
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    assert(step >= 0 && step < xStart.size());

    // Set first and last point positions
    rod.x.front() = xStart[step];
    rod.x.back() = xEnd[step];

    // Set first and last material frames
    rod.Q.front() = QStart[step];
    rod.Q.back() = QEnd[step];
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set initial point velocities to zero (quasi-static assumption)
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;
    rod.v.front().z = 0.0;

    // Set initial point angular velocity to zero (quasi-static assumption)
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Set last point velocities to zero (quasi-static assumption)
    rod.v.back().x = 0;
    rod.v.back().y = 0;
    rod.v.back().z = 0;

    // Set last point angular velocity to zero (quasi-static assumption)
    rod.w.back().x = 0.0;
    rod.w.back().y = 0.0;
    rod.w.back().z = 0.0;
  }
};

class SolenoidQuasiStaticEndpointBC : public RodBC {
 public:
  vector<Vector3> xEnd;    // position of x[n] in time
  vector<Matrix3> QStart;  // material frame at x[0] in time
  vector<Matrix3> QEnd;    // material frame at x[n] in time

  SolenoidQuasiStaticEndpointBC(vector<Vector3> _xEnd, vector<Matrix3> _QStart,
                                vector<Matrix3> _QEnd)
      : RodBC(), xEnd(_xEnd), QStart(_QStart), QEnd(_QEnd) {
    // Sanity checks
    assert(xEnd.size() > 0);
    assert(QStart.size() > 0);
    assert(QEnd.size() > 0);

    assert(xEnd.size() == QStart.size());
    assert(xEnd.size() == QEnd.size());
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    assert(step >= 0 && step < xEnd.size());

    // Pin x y start
    rod.x.front().x = 0.0;
    rod.x.front().y = 0.0;

    // Set last point positions
    rod.x.back() = xEnd[step];

    // Set first and last material frames
    rod.Q.front() = QStart[step];
    rod.Q.back() = QEnd[step];
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set initial point velocities to zero (quasi-static assumption)
    rod.v.front().x = 0.0;
    rod.v.front().y = 0.0;

    // Set initial point angular velocity to zero (quasi-static assumption)
    rod.w.front().x = 0.0;
    rod.w.front().y = 0.0;
    rod.w.front().z = 0.0;

    // Set last point velocities to zero (quasi-static assumption)
    rod.v.back().x = 0;
    rod.v.back().y = 0;
    rod.v.back().z = 0;

    // Set last point angular velocity to zero (quasi-static assumption)
    rod.w.back().x = 0.0;
    rod.w.back().y = 0.0;
    rod.w.back().z = 0.0;
  }
};

class SolenoidsBC : public RodBC {
 public:
  const REAL twistingTime;  // Time during which the twisting is performed
                            // before relaxing the system
  const REAL R;
  Vector3 angVel;
  Vector3 shrinkVel;
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Vector3 axisOfRotInMaterialFrame;
  Matrix3 startQ;
  Matrix3 endQ;
  Matrix3 finalEndQ;

  SolenoidsBC(vector<Rod *> &rodptrs, const REAL _twistingTime, const REAL _R)
      : RodBC(), twistingTime(_twistingTime), R(_R) {
    assert(rodptrs.size() == 1);

    // Sanity checks
    assert(twistingTime > 0);
    assert(R > 0);

    // Compute total rotation angle
    const REAL totalRotAngle = R * 2.0 * M_PI;

    // Compute rotation velocities of the hanging end
    const REAL angVelScalar = (totalRotAngle / twistingTime);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();

    // Direction of the rod
    direction = -(endX - startX).unitize();

    // Axis of rotation in material frame
    axisOfRotInMaterialFrame = (endQ * direction).unitize();

    // Note that in this case the axis of rotation and shrinkage coincides with
    // the one of the material frame but it is not always the case!
    angVel = angVelScalar * axisOfRotInMaterialFrame;

    // Final material frames of the two end of the rod
    finalEndQ = exp(-totalRotAngle * axisOfRotInMaterialFrame) * endQ;
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    if (time > twistingTime) {
      rod.x.front() = startX;
      rod.x.back().y = 0.0;
      rod.x.back().z = 0.0;

      rod.Q.front() = startQ;
      rod.Q.back() = finalEndQ;
    }
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    if (time > twistingTime) {
      // Set initial point velocities to zero
      rod.v.front().x = 0.0;
      rod.v.front().y = 0.0;
      rod.v.front().z = 0.0;

      // Set initial point angular velocity to zero
      rod.w.front().x = 0.0;
      rod.w.front().y = 0.0;
      rod.w.front().z = 0.0;

      // Set last point velocities to zero
      rod.v.back().x = 0;
      rod.v.back().y = 0;
      // rod.v.back().z = 0;

      // Set last point angular velocity to zero
      rod.w.back().x = 0.0;
      rod.w.back().y = 0.0;
      rod.w.back().z = 0.0;
    } else {
      // Set initial point velocities
      rod.v.front() = Vector3(0.0, 0.0, 0.0);

      // Set initial point angular velocity
      rod.w.front() = Vector3(0.0, 0.0, 0.0);

      // Set last point velocities
      rod.v.back().y = 0.0;
      rod.v.back().z = 0.0;

      // Set last point angular velocity
      rod.w.back() = -angVel;
    }
  }
};

class SolenoidsBCTeja : public RodBC {
 public:
  const REAL twistingTime;  // Time during which the twisting is performed
                            // before relaxing the system
  const REAL R;
  Vector3 angVel;
  Vector3 shrinkVel;
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Vector3 axisOfRotInMaterialFrame;
  Matrix3 startQ;
  Matrix3 endQ;
  Matrix3 finalEndQ;

  SolenoidsBCTeja(Rod *rodptrs, const REAL _twistingTime, const REAL _R)
      : RodBC(), twistingTime(_twistingTime), R(_R) {
    // Sanity checks
    assert(twistingTime > 0);
    assert(R > 0);

    // Initial positions of the two end of the rod
    endX = rodptrs->x.back();

    // Initial material frames of the two end of the rod
    endQ = rodptrs->Q.back();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    rod.x.back() = endX;
    rod.Q.back() = endQ;
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    // Set initial point velocities to zero
    rod.v.back().x = 0.0;
    rod.v.back().y = 0.0;
    rod.v.back().z = 0.0;

    // Set initial point angular velocity to zero
    rod.w.back().x = 0.0;
    rod.w.back().y = 0.0;
    rod.w.back().z = 0.0;
  }
};

class SolenoidsBC_JCP : public RodBC {
 public:
  const REAL twistingTime;  // Time during which the twisting is performed
                            // before relaxing the system
  const REAL R;
  Vector3 angVel;
  Vector3 shrinkVel;
  Vector3 direction;
  Vector3 startX;
  Vector3 endX;
  Vector3 axisOfRotInMaterialFrame;
  Matrix3 startQ;
  Matrix3 endQ;
  Matrix3 finalEndQ;

  SolenoidsBC_JCP(vector<Rod *> &rodptrs, const REAL _twistingTime,
                  const REAL _R)
      : RodBC(), twistingTime(_twistingTime), R(_R) {
    assert(rodptrs.size() == 1);

    // Sanity checks
    assert(twistingTime > 0);
    assert(R > 0);

    // Compute total rotation angle
    const REAL totalRotAngle = R * 2.0 * M_PI;

    // Compute rotation velocities of the hanging end
    const REAL angVelScalar = (totalRotAngle / twistingTime);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();
    endX = rodptrs[0]->x.back();

    // Direction of the rod
    direction = (endX - startX).unitize();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
    endQ = rodptrs[0]->Q.back();

    // Axis of rotation in material frame
    axisOfRotInMaterialFrame = (endQ * direction).unitize();

    // Note that in this case the axis of rotation and shrinkage coincides with
    // the one of the material frame but it is not always the case!
    angVel = angVelScalar * axisOfRotInMaterialFrame;

    // Final material frames of the two end of the rod
    finalEndQ = exp(-totalRotAngle * axisOfRotInMaterialFrame) * endQ;
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    if (time > twistingTime) {
      rod.x.front() = startX;
      rod.x[100].x = 0.0;
      rod.x[100].y = 0.0;

      rod.Q[0] = startQ;
      rod.Q[99] = finalEndQ;
    }
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    /*
//additional constraints for making multiple plectonemes
for(unsigned int i = 60;i<70;i++)
{
      rod.v[i].x = 0;
      rod.v[i].y = 0;
}
for(unsigned int i = 130;i<140;i++)
{
      rod.v[i].x = 0;
      rod.v[i].y = 0;
}
    */

    if (time < (twistingTime + 25.0)) {
      if (time > twistingTime) {
        // Set initial point velocities to zero
        rod.v.front().x = 0.0;
        rod.v.front().y = 0.0;
        rod.v.front().z = 0.0;

        // Set initial point angular velocity to zero
        rod.w.front().x = 0.0;
        rod.w.front().y = 0.0;
        rod.w.front().z = 0.0;

        // Set last point velocities to zero
        rod.v.back().x = 0;
        rod.v.back().y = 0;
        rod.v.back().z = 0.02;

        // Set last point angular velocity to zero
        rod.w.back().x = 0.0;
        rod.w.back().y = 0.0;
        rod.w.back().z = 0.0;
      } else {
        // Set initial point velocities
        rod.v.front() = Vector3(0.0, 0.0, 0.0);

        // Set initial point angular velocity
        rod.w.front() = Vector3(0.0, 0.0, 0.0);

        // Set last point velocities
        rod.v.back().x = 0.0;
        rod.v.back().y = 0.0;
        rod.v.back().z = 0.0;
        // Set last point angular velocity
        rod.w.back() = -angVel;
      }
    } else {
      // Set initial point velocities to zero
      rod.v.front().x = 0.0;
      rod.v.front().y = 0.0;
      rod.v.front().z = 0.0;

      // Set initial point angular velocity to zero
      rod.w.front().x = 0.0;
      rod.w.front().y = 0.0;
      rod.w.front().z = 0.0;

      // Set last point velocities to zero
      rod.v.back().x = 0;
      rod.v.back().y = 0;

      // Set last point angular velocity to zero
      rod.w.back().x = 0.0;
      rod.w.back().y = 0.0;
      rod.w.back().z = 0.0;
    }
  }
};

class BucklingCageBC : public RodBC {
 public:
  const REAL compressingTime;  // Time during which the twisting is performed
                               // before relaxing the system
  const REAL Dis;
  Vector3 compressVel;
  Vector3 direction;
  REAL Height;
  Vector3 Position1;
  Vector3 Position2;

  BucklingCageBC(Rod *rodptr, const REAL _compressingTime, const REAL _Dis)
      : RodBC(), compressingTime(_compressingTime), Dis(_Dis) {
    // Sanity checks
    assert(compressingTime > 0);
    assert(Dis > 0);

    // Compute compressing velocities
    const REAL VelScalar = (Dis / compressingTime);

    // Direction of the rod
    /*
    direction = (rodptr->x[27]-rodptr->x[11]).unitize();
Height = rodptr->x[27][2];
*/
    direction = (rodptr->x[30] - rodptr->x[8]).unitize();
    Height = rodptr->x[30][2];
    compressVel = VelScalar * direction;
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    // if ( time < twistingTime )
    //{
    if (time < compressingTime) {
      for (unsigned int i = 8; i < 12; i++) {
        rod.x[i].z = Height;
      }
      for (unsigned int i = 27; i < 31; i++) {
        rod.x[i].z = Height;
      }
    }

    if ((time > compressingTime) && (time < 2.0)) {
      rod.x[14] = Position1;
      rod.x[24] = Position2;
    }
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    if (time < compressingTime) {
      for (unsigned int i = 0; i < 12; i++) {
        rod.v[i] = compressVel;
        rod.v[i].z = 0.0;
      }
      for (unsigned int i = 27; i < 31; i++) {
        rod.v[i] = -compressVel;
        rod.v[i].z = 0.0;
      }
      Position1 = rod.x[14];
      Position2 = rod.x[24];
    }

    if ((time > compressingTime) && (time < 2.0)) {
      rod.v[14].x = 0.0;
      rod.v[14].y = 0.0;
      rod.v[14].z = 0.0;

      rod.v[24].x = 0.0;
      rod.v[24].y = 0.0;
      rod.v[24].z = 0.0;
    }
  }
};

class KnotBC : public RodBC {
 public:
  const REAL shakingTime;  // Time during which the shaking is performed before
                           // relaxing
  const REAL A;            // Shake magnitude
  const REAL N;            // Number of shakes
  Vector3 startX;
  Matrix3 startQ;

  KnotBC(vector<Rod *> &rodptrs, const REAL _shakingTime, const REAL _A,
         const REAL _N)
      : RodBC(), shakingTime(_shakingTime), A(_A), N(_N) {
    assert(rodptrs.size() == 1);

    // Sanity checks
    assert(shakingTime > 0);
    assert(A > 0);
    assert(N > 0);

    // Initial positions of the two end of the rod
    startX = rodptrs[0]->x.front();

    // Initial material frames of the two end of the rod
    startQ = rodptrs[0]->Q.front();
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    if (time > shakingTime) {
      rod.x.front() = startX;

      // rod.Q.front() = startQ;
    }
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    if (time > shakingTime) {
      // Set initial point velocities to zero
      rod.v.front().x = 0.0;
      rod.v.front().y = 0.0;
      rod.v.front().z = 0.0;

      // Set initial point angular velocity to zero
      // rod.w.front().x = 0.0;
      // rod.w.front().y = 0.0;
      // rod.w.front().z = 0.0;

    } else {
      // Set initial point velocities
      // rod.v.front() =
      // Vector3(.1*A*cos(.1*2.0*M_PI*N*time/shakingTime),.1*A*cos(.1*2.0*M_PI*N*time/shakingTime),A*cos(2.0*M_PI*N*time/shakingTime));
      rod.v.front() =
          Vector3(0.0, 0.0,
                  A * cos(2.0 * M_PI * N * time / shakingTime) /
                      abs(cos(2.0 * M_PI * N * time / shakingTime)));

      // Set initial point angular velocity
      // rod.w.front() = Vector3(0.0,0.0,0.0);
    }
  }
};

class SelfContactBC : public RodBC {
 public:
  const REAL moveTime;  // Time during which moving is performed
  const REAL D;         // Distance over which endpoints move
  REAL speed;
  Vector3 initialStartX, initialEndX, finalStartX, finalEndX, direction;
  Matrix3 initialQ;

  SelfContactBC(vector<Rod *> &rodptrs, const REAL _moveTime, const REAL _D)
      : RodBC(), moveTime(_moveTime), D(_D) {
    assert(rodptrs.size() == 1);

    // Sanity checks
    assert(_moveTime > 0);
    assert(D > 0);

    // Initial positions of the two end of the rod
    initialStartX = rodptrs[0]->x.front();
    initialEndX = rodptrs[0]->x.back();

    // Compute rod direction
    direction = (initialEndX - initialStartX).unitize();

    // Endpoint speed
    speed = D / moveTime;

    // Initial material frames of the two end of the rod
    initialQ = rodptrs[0]->Q.front();

    // Compute final positions of endpoints
    finalStartX = initialStartX + D * direction;
    finalEndX = initialEndX - D * direction;
  }

  void dirichlet(Rod &rod, const unsigned int step, const REAL dt,
                 const REAL time) {
    if (time > moveTime) {
      rod.x.front() = finalStartX;
      rod.x.back() = finalEndX;

      rod.Q.front() = initialQ;
      rod.Q.back() = initialQ;
    }
  }

  void neumann(Rod &rod, const unsigned int step, const REAL dt,
               const REAL time) {
    if (time > moveTime) {
      // Set initial point velocities to zero
      rod.v.front() = Vector3(0.0, 0.0, 0.0);

      // Set initial point velocities to zero
      rod.v.back() = Vector3(0.0, 0.0, 0.0);

      // Set angular velocities at endpoints to zero
      rod.w.front() = Vector3(0.0, 0.0, 0.0);
      rod.w.back() = Vector3(0.0, 0.0, 0.0);

    } else {
      // Set initial point velocities
      rod.v.front() = speed * direction;
      rod.v.back() = -speed * direction;
      // rod.v.front() =
      // Vector3(D*sin(3.0/2.0*M_PI*time/moveTime),D*cos(3.0/2.0*M_PI*time/moveTime),0.0);
      // rod.v.back() =
      // Vector3(D*sin(3.0/2.0*M_PI*time/moveTime),-D*cos(3.0/2.0*M_PI*time/moveTime),0.0);

      // Set angular velocities at endpoints to zero
      rod.w.front() = Vector3(0.0, 0.0, 0.0);
      rod.w.back() = Vector3(0.0, 0.0, 0.0);
    }
  }
};

/*
class EndpointBC : public RodBC
{
public:
        VV3 xStart; // position of x[0] in time
        VV3 xEnd; // position of x[n] in time
        VM3 QStart; // value of phi[0] in time
        VM3 QEnd; // value of phi[n-1] in time

        bool xStartProscribed;
        bool xEndProscribed;
        bool QStartProscribed;
        bool QEndProscribed;

        EndpointBC(VV3 xStart, VV3 xEnd, VM3 QStart, VM3 QEnd) : RodBC(),
xStart(xStart), xEnd(xEnd), QStart(QStart), QEnd(QEnd)
        {
                xStartProscribed = (xStart.size() > 0) ? true : false;
                xEndProscribed = (xEnd.size() > 0) ? true : false;
                QStartProscribed = (QStart.size() > 0) ? true : false;
                QEndProscribed = (QEnd.size() > 0) ? true : false;
        }

        // Free BCs, ie this function does nothing
        EndpointBC() : RodBC()
        {
                xStartProscribed = false;
                xEndProscribed = false;
                QStartProscribed = false;
                QEndProscribed = false;
        }

        void operator()(Rod& rod, int step, REAL h)
        {
#ifndef ANDREW
                assert(step>=0);

                if (xStartProscribed)
                {
                        assert(step<(int)xStart.size());

                        const int stepPlus = min(step+1, (int)xStart.size()-1);
                        rod.x[0] = xStart[step];
                        rod.v[0] = (xStart[stepPlus] - xStart[step]) / h;
                }

                if (xEndProscribed)
                {
                        assert(step<(int)xEnd.size());

                        const int stepPlus = min(step+1, (int)xEnd.size()-1);
                        rod.x[rod.n] = xEnd[step];
                        rod.v[rod.n] = (xEnd[stepPlus] - xEnd[step]) / h;
                }

                if (QStartProscribed)
                {
                        assert(step<(int)QStart.size());

                        rod.Q[0] = QStart[step];

                        const int stepPlus = min(step+1, (int)QStart.size()-1);
                        const Matrix3 Q = QStart[step];
                        const Matrix3 Q1 = QStart[stepPlus];

#ifndef NDEBUG
                        // Make sure these are nice matrices and that I can
compute R safely if ( !Q.goodNumerics() )
                        {
                                cout << "Problem in EndpointBC -->
QStartProscribed: an element of Q is NaN" << endl; abort();
                        }

                        if (fabs(Q.det())<Tolerance::tol())
                        {
                                cout << "Problem in EndpointBC -->
QStartProscribed: the matrix Q is singular" << endl; abort();
                        }

                        if ( !Q1.goodNumerics() )
                        {
                                cout << "Problem in EndpointBC -->
QStartProscribed: an element of Q1 is NaN" << endl; abort();
                        }
#endif

                        // Compute roatation matrix R: Q --> Q1
                        const Matrix3 R = Q1 * Q.I();

                        // Make sure that the argument of acos lies in the
domain [-1 +1] const REAL argumentACos = (R.tr() - 1.0)/2.0; assert(
argumentACos>=-1.0-std::numeric_limits<REAL>::epsilon() &&
argumentACos<=1.0+std::numeric_limits<REAL>::epsilon() ); const REAL
argumentACosClamped = max(-1.0,min(argumentACos,1.0));

                        // Compute rotation angle associated with R
                        const REAL angle = acos( argumentACosClamped );
                        assert( !(angle!=angle) );

                        // Compute frame rotation velocity vector
                        if (angle <= Tolerance::tol())
                        {
                                rod.w[0] = Vector3();
                        }
                        else
                        {
                                // Compute axis of rotation in material frame of
reference and assign angular velocity
                                // Note that we have to put a minus in front
thanks to that genius of Andrew assert( sin(angle)!=0 ); const Vector3 k =
-1.0/(2*sin(angle))*Vector3( (R.r3c2-R.r2c3), (R.r1c3-R.r3c1), (R.r2c1-R.r1c2)
);

                                rod.w[0] = (angle/h)*k = k.unitize();

                                // Extra care here we really want to make sure
we are not assigning garbage here assert( !(rod.w[0][0]!=rod.w[0][0]) ); assert(
!(rod.w[0][1]!=rod.w[0][1]) ); assert( !(rod.w[0][2]!=rod.w[0][2]) );
                        }

                        rod.w[0] = Vector3();
                }

                if (QEndProscribed)
                {
                        assert(step<(int)QEnd.size());

                        rod.Q[rod.n-1] = QEnd[step];

                        const int stepMinus = max(step-1, 0);
                        const Matrix3 Q = QEnd[stepMinus];
                        const Matrix3 Q1 = QEnd[step];

#ifndef NDEBUG
                        // Make sure these are nice matrices and that I can
compute R safely if ( !Q.goodNumerics() )
                        {
                                cout << "Problem in EndpointBC -->
QEndProscribed: an element of Q is NaN" << endl; abort();
                        }

                        if (fabs(Q.det())<Tolerance::tol())
                        {
                                cout << "Problem in EndpointBC -->
QEndProscribed: the matrix Q is singular" << endl; abort();
                        }

                        if ( !Q1.goodNumerics() )
                        {
                                cout << "Problem in EndpointBC -->
QEndProscribed: an element of Q1 is NaN" << endl; abort();
                        }
#endif

                        // Compute roatation matrix R: Q --> Q1
                        const Matrix3 R = Q1 * Q.I();

                        // Make sure that the argument of acos lies in the
domain [-1 +1] const REAL argumentACos = (R.tr() - 1.0)/2.0; assert(
argumentACos>=-1.0-std::numeric_limits<REAL>::epsilon() &&
argumentACos<=1.0+std::numeric_limits<REAL>::epsilon() ); const REAL
argumentACosClamped = max(-1.0,min(argumentACos,1.0));

                        // Compute rotation angle associated with R
                        const REAL angle = acos( argumentACosClamped );
                        assert( !(angle!=angle) );

                        // Compute frame rotation velocity vector
                        if (angle <= Tolerance::tol())
                        {
                                rod.w[rod.n-1] = Vector3();
                        }
                        else
                        {
                                // Compute axis of rotation in material frame of
reference and assign angular velocity
                                // Note that we have to put a minus in front
thanks to that genius of Andrew assert( sin(angle)!=0 ); const Vector3 k =
-1.0/(2*sin(angle))*Vector3( (R.r3c2-R.r2c3), (R.r1c3-R.r3c1), (R.r2c1-R.r1c2)
);

                                rod.w[rod.n-1] = (angle/h)*k = k.unitize();

                                // Extra care here we really want to make sure
we are not assigning garbage here assert(
!(rod.w[rod.n-1][0]!=rod.w[rod.n-1][0]) ); assert(
!(rod.w[rod.n-1][1]!=rod.w[rod.n-1][1]) ); assert(
!(rod.w[rod.n-1][2]!=rod.w[rod.n-1][2]) );
                        }

                        rod.w[rod.n-1] = Vector3();
                }
#else
                if (xStartProscribed)
                {
                        rod.x[0] = xStart[step];
                        rod.v[0] = (xStart[step] - xStart[step-1]) / h;
                }
                if (xEndProscribed)
                {
                        rod.x[rod.n] = xEnd[step];
                        rod.v[rod.n] = (xEnd[step] - xEnd[step-1]) / h;
                }
                if (QStartProscribed)
                {
                        rod.Q[0] = QStart[step];
                        rod.w[0] = Vector3(); //sorta HAX
                }
                if (QEndProscribed)
                {
                        rod.Q[rod.n-1] = QEnd[step];
                        rod.w[rod.n-1] = Vector3(); //sorta HAX
                }
#endif
        }
};


// restrains a rod's final and initional point to not move in the x or y
directions class ZClampedBC : public RodBC
{
public:
        ZClampedBC(): RodBC() {}

        void operator()(Rod& rod, int step, REAL h)
        {
                switch (status)
                {
                        case TypeBC::NEUMANN:
                                {
                                        // Set initial point velocities to zero
                                        rod.v.front().x = 0.0;
                                        rod.v.front().y = 0.0;

                                        // Set initial point angular velocity to
zero rod.w.front().x = 0.0; rod.w.front().y = 0.0;

                                        // Set last point velocity
                                        rod.v.back().x = 0;
                                        rod.v.back().y = 0;

                                        // Set last point angular velocity
                                        rod.w.back().x = 0.0;
                                        rod.w.back().y = 0.0;

                                        status = TypeBC::DIRICHLET;
                                }
                                break;

                        case TypeBC::DIRICHLET:

                                status = TypeBC::NEUMANN;
                                break;

                        default:
                                break;
                }
        }
};
 */

#endif
