#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <assert.h>
#include <limits>
#include "MathFunctions.h"
#include "Rod.h"

// Class defining interaction forces between rod and rod, substrate, etc
class Interaction {
 protected:
  typedef std::vector<Vector3> VV3;
  typedef std::vector<Matrix3> VM3;
  typedef std::vector<REAL> VREAL;
  typedef std::vector<bool> VBOOL;
  typedef std::vector<int> VINT;
  typedef std::vector<Rod *> Vrodptr;
  typedef std::vector<Interaction *> Vinterptr;

  Vrodptr rodptrs;

 public:
  Interaction(Vrodptr rodptrs) : rodptrs(rodptrs) {}
  virtual ~Interaction() {}

  virtual void applyForces(const REAL time) = 0;
};

// A frictional plane, for the moment set at y=0, rod free for y>r.

// A container with ground, wall, top surface
// Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————
// Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————Nest—————

class FrictionPlaneInteraction : public Interaction {
 protected:
  const REAL kPlane;  // Elastic response to prevent interpenetration rod-plane
  const REAL
      etaPlane;  // Dissipative term appliet to the elastice response rod-plane
  const REAL muKineticPlane;
  const REAL muStaticPlane;
  const REAL vStaticPlane;
  const Vector3 normalPlane;
  const Vector3 originPlane;
  REAL MaxH = 0.0;
  REAL depth = 0.0;
  const REAL timeT = 1.5;
  REAL Strains[8] = {80.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0};
  int index = 0;

  const REAL _linear(const REAL _vel, const REAL _velThresold) const {
    const REAL vel = fabs(_vel);
    const REAL velThresold = fabs(_velThresold);

    const REAL width = 0.5 * velThresold;
    const REAL velDiff = vel - velThresold;
    REAL f = 1.0;

    if (vel > (velThresold)) f = 1.0 - velDiff / width;

    if (vel > (velThresold + width)) f = 0.0;

    return f;
  }

  const REAL _sigmoid(const REAL _vel, const REAL _velThresold) const {
    const REAL vel = fabs(_vel);
    const REAL velThresold = fabs(_velThresold);

    const REAL velDiff = vel - 1.5 * velThresold;
    return 0.5 + 0.5 * erf(-(8.0 / velThresold) * velDiff);
  }

 public:
  FrictionPlaneInteraction(Vrodptr &rodptrs, const Vector3 _normalPlane,
                           const Vector3 _originPlane, const REAL _kPlane,
                           const REAL _etaPlane, const REAL _muKineticPlane,
                           const REAL _muStaticPlane,
                           const REAL _vStatic = 1e-4)
      : Interaction(rodptrs),
        normalPlane(_normalPlane),
        originPlane(_originPlane),
        kPlane(_kPlane),
        etaPlane(_etaPlane),
        muStaticPlane(_muStaticPlane),
        muKineticPlane(_muKineticPlane),
        vStaticPlane(_vStatic) {}

  void applyForces(const REAL time) {
    const int nor = rodptrs.size();
    REAL TotalF = 0.0;

    for (unsigned int j = 0; j < nor; j++) {
      Rod *Rod = rodptrs[j];

      if ((time > (timeT - 0.001)) && (time < timeT))
        MaxH = ((Rod->x[8][2]) > MaxH) ? Rod->x[8][2] : MaxH;

      //************************************************************
      // Surface Above Compressing force ——— assume top surface to be horizontal
      // for now Check only the end point Ground and wall and cap surface
      // pushing force and friction

      for (unsigned int i = 0; i < (Rod->n + 1); i = i + (Rod->n)) {
        // For ground=======

        if (Rod->x[i][2] < 0.5) {
          const Vector3 elementV = Rod->v[i];

          const Vector3 overallRodForces =
              Rod->totalInternalForces[i] + Rod->externalForces[i];
          const REAL currentRodForcesInNormaldirectionSign =
              (overallRodForces % normalPlane);
          const Vector3 currentRodForcesInNormaldirection =
              (currentRodForcesInNormaldirectionSign < 0.0)
                  ? currentRodForcesInNormaldirectionSign * normalPlane
                  : 0.0 * normalPlane;
          const Vector3 elasticPlaneResponse =
              -kPlane * min(Rod->x[i][2], 0.0) * normalPlane;
          const Vector3 dampingForcePlane =
              -etaPlane * (elementV % normalPlane) * normalPlane;
          const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
          const Vector3 totalforce =
              (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
          Rod->externalForces[i] += totalforce;

          const REAL normalForce = normalForcePlane.length();

          // Compute axial velocity
          const Vector3 slipVel =
              elementV - (elementV % normalPlane) * normalPlane;

          // Friction in axial direction
          {
            const REAL vel = fabs(slipVel.length());
            const REAL f = _linear(vel, vStaticPlane);

            const Vector3 kineticFrictionForce =
                -(1.0 - f) * muKineticPlane * normalForce * slipVel.unitize();
            const Vector3 staticFrictionForce =
                f * muStaticPlane * normalForce * Vector3(0.0, 1.0, 0.0);

            Rod->kineticFrictionsForce[i] += kineticFrictionForce;
            Rod->staticFrictionsAxialForceForward[i] += staticFrictionForce;

          }  // Friction in axial direction
        }    // if within range

        // For wall=======

        const REAL distanceFromAxes =
            sqrt(Rod->x[i][0] * Rod->x[i][0] + Rod->x[i][1] * Rod->x[i][1]);
        if (distanceFromAxes > (85.0 - 0.5)) {
          const Vector3 elementV = Rod->v[i];
          const Vector3 normalPlane =
              -Vector3(Rod->x[i][0], Rod->x[i][1], 0.0).unitize();
          const Vector3 overallRodForces =
              Rod->totalInternalForces[i] + Rod->externalForces[i];

          const REAL currentRodForcesInNormaldirectionSign =
              (overallRodForces % normalPlane);
          const Vector3 currentRodForcesInNormaldirection =
              (currentRodForcesInNormaldirectionSign < 0.0)
                  ? currentRodForcesInNormaldirectionSign * normalPlane
                  : 0.0 * normalPlane;
          const Vector3 elasticPlaneResponse =
              -kPlane * min((85.0 - distanceFromAxes), 0.0) * normalPlane;
          const Vector3 dampingForcePlane =
              -etaPlane * (elementV % normalPlane) * normalPlane;
          const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
          const Vector3 totalforce =
              (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
          Rod->externalForces[i] += totalforce;
        }

        // For Cap=======
        if (time > timeT) {
          REAL speed = 15.0 - 8.0 * (time - timeT) / (3 * timeT);  //(mm/s)
          speed = (speed < 8.0) ? 8.0 : speed;
          const REAL startingheight = MaxH + 0.01;  //(m)
          if ((j == 0) && (i == 0)) {
            depth += speed * (5.7e-7);
            if (depth > (2 * Strains[index])) {
              depth = 0;
              index++;
            }
          }

          REAL realdepth = (depth > Strains[index])
                               ? (2 * Strains[index] - depth)
                               : depth;  //(m)

          const REAL currentheight = startingheight - realdepth;
          if (j == 0) rodptrs[0]->MaxHeight = currentheight;

          if (Rod->x[i][2] > currentheight) {
            const Vector3 elementV =
                (Rod->v[i] + speed * Vector3(0.0, 0.0, 1.0));
            const Vector3 normalPlane = Vector3(0.0, 0.0, -1.0);

            const Vector3 overallRodForces =
                Rod->totalInternalForces[i] + Rod->externalForces[i];
            const REAL currentRodForcesInNormaldirectionSign =
                (overallRodForces % normalPlane);
            const Vector3 currentRodForcesInNormaldirection =
                (currentRodForcesInNormaldirectionSign < 0.0)
                    ? currentRodForcesInNormaldirectionSign * normalPlane
                    : 0.0 * normalPlane;
            const Vector3 elasticPlaneResponse =
                -kPlane * min((currentheight - Rod->x[i][2]), 0.0) *
                normalPlane;
            const Vector3 dampingForcePlane =
                -etaPlane * (elementV % normalPlane) * normalPlane;
            const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
            const Vector3 totalforce =
                (normalForcePlane + elasticPlaneResponse + dampingForcePlane);

            Rod->externalForces[i] += totalforce;
            TotalF += (totalforce % Vector3(0.0, 0.0, -1.0));
          }  // if within range

        }  // if time >timeT
      }    // for loop all elements

    }  // for loop all rods

    rodptrs[0]->Tforce = TotalF;

  }  // apply force
};

/*

//Crab———Crab———Crab———Crab———Crab———Crab———Crab———Crab———


//Simple friction plane
class FrictionPlaneInteraction: public Interaction
{
protected:
        const REAL kPlane;		// Elastic response to prevent
interpenetration rod-plane const REAL etaPlane;	// Dissipative term appliet to
the elastice response rod-plane const REAL muKineticPlane; const REAL
muStaticPlane; const REAL vStaticPlane; const Vector3 normalPlane; const Vector3
originPlane;


        const REAL _linear(const REAL _vel, const REAL _velThresold) const
        {
                const REAL vel = fabs(_vel);
                const REAL velThresold = fabs(_velThresold);

                const REAL width = 0.5*velThresold;
                const REAL velDiff = vel - velThresold;
                REAL f = 1.0;

                if ( vel>(velThresold) )
                        f = 1.0-velDiff/width;

                if ( vel>(velThresold+width) )
                        f = 0.0;

                return f;
        }

        const REAL _sigmoid(const REAL _vel, const REAL _velThresold) const
        {
                const REAL vel = fabs(_vel);
                const REAL velThresold = fabs(_velThresold);

                const REAL velDiff = vel - 1.5*velThresold;
                return 0.5 + 0.5*erf( -(8.0/velThresold)*velDiff );
        }

public:
        FrictionPlaneInteraction(Vrodptr& rodptrs, const Vector3 _normalPlane,
const Vector3 _originPlane, const REAL _kPlane, const REAL _etaPlane, const REAL
_muKineticPlane, const REAL _muStaticPlane, const REAL _vStatic = 1e-4) :
                Interaction(rodptrs), normalPlane(_normalPlane),
originPlane(_originPlane), kPlane(_kPlane), etaPlane(_etaPlane),
muStaticPlane(_muStaticPlane), muKineticPlane(_muKineticPlane),
                vStaticPlane(_vStatic){}

        void applyForces(const REAL time)
        {
           const int nor = rodptrs.size();

           for(unsigned int j=0; j<nor; j++)
            {
                Rod* rod = rodptrs[j];


//Ground surface pushing force and friction

                vector<Vector3> overallRodForces =
vector<Vector3>(rod->totalInternalForces.size()-1);
                vFromPointsToElements(rod->totalInternalForces +
rod->externalForces, overallRodForces);

                for (int i=(rod->n-1); i < rod->n; i++)
                {
                        // Fetch initial useful quantities
                        const REAL r = rod->r[i];
                        const Vector3 elementX = 0.5*(rod->x[i] + rod->x[i+1]);
                        const REAL distanceFromPlane = (elementX - originPlane)
% normalPlane;

                        // Numerical fluctuations may cause the distnace from
plane to be larger than r.
                        // We ameng this problem by allowing a tolerance
proportinal (arbitrary constant) to the radius.
                        // For more robustness check the SUGGESTION in
Rod.cpp->applyFrictions if ( (distanceFromPlane-r) < 1e-4 )
                        {
                                // Fetch remaining quantites
                                const Vector3 elementV = 0.5*(rod->v[i] +
rod->v[i+1]); const Vector3 axialDir = (rod->edge[i] /
rod->edge[i].length()).unitize(); const Vector3 w = rod->w[i]; const Matrix3 Q =
rod->Q[i];

                                // Compute plane tangent and bitangent (rolling)
direction const Vector3 rollingDirection = (normalPlane * axialDir).unitize();
                                const Vector3 axialDirection = (rollingDirection
* normalPlane).unitize();

                                const REAL currentRodForcesInNormaldirectionSign
= (overallRodForces[i] % normalPlane); const Vector3
currentRodForcesInNormaldirection = (currentRodForcesInNormaldirectionSign<0.0)?
currentRodForcesInNormaldirectionSign*normalPlane : 0.0*normalPlane; const
Vector3 elasticPlaneResponse = - kPlane * min(distanceFromPlane-r,0.0) *
normalPlane; const Vector3 dampingForcePlane = - etaPlane * (elementV %
normalPlane)*normalPlane; const Vector3 normalForcePlane =
-currentRodForcesInNormaldirection; rod->externalForces[i]   += 0.5 *
(normalForcePlane + elasticPlaneResponse + dampingForcePlane);
                                rod->externalForces[i+1] += 0.5 *
(normalForcePlane + elasticPlaneResponse + dampingForcePlane); const REAL
normalForce = normalForcePlane.length();

                                // Compute slip velocity decomposed in rolling
and tangential directions const Vector3 armTorque = -r * normalPlane; const
Vector3 translationVelInRollingDirection = (elementV %
rollingDirection)*rollingDirection; const Vector3
rotationalVelInRollingDirection = ((Q.T()*(w*(Q*armTorque))) %
rollingDirection)*rollingDirection; const Vector3 slipVelInRollingDirection =
translationVelInRollingDirection + rotationalVelInRollingDirection;

                                // Friction in rolling direction
                                {
                                        const REAL f =
_linear(slipVelInRollingDirection.length(), vStaticPlane);

                                        const Vector3 kineticFrictionForce =
-(1.0-f) * muKineticPlane * normalForce * slipVelInRollingDirection.unitize();
                                        const Vector3 staticFrictionForce  = f *
muStaticPlane  * normalForce * rollingDirection.unitize();

                                        rod->kineticFrictionsForce[i]	 += 0.5
* kineticFrictionForce; rod->kineticFrictionsForce[i+1] += 0.5 *
kineticFrictionForce; rod->kineticFrictionsTorque[i]  += Q*(armTorque *
kineticFrictionForce);

                                        // Note that the static friction
genereta torques as well.
                                        // Since they are annoying to compute
and require more knowledge of the rod's configuration
                                        // their calculation is in the rod's
class itself. But they are accounted for!! rod->staticFrictionsRollingForce[i]
+= 0.5 * staticFrictionForce; rod->staticFrictionsRollingForce[i+1] += 0.5 *
staticFrictionForce; rod->staticFrictionsNormalPlane[i]    += normalPlane;

                                        rod->rollingf[i] = f;
                                }


                                // Compute axial velocity
                                const Vector3 slipVelInAxiallDirection =
(elementV % axialDirection) * axialDirection;

                                // Friction in axial direction
                                {
                                        const REAL f =
_linear(slipVelInAxiallDirection.length(), vStaticPlane);

                                        const Vector3 kineticFrictionForce =
-(1.0-f) * muKineticPlane * normalForce * slipVelInAxiallDirection.unitize();
                                        const Vector3 staticFrictionForce  = f *
muStaticPlane  * normalForce * axialDirection.unitize();

                                        rod->kineticFrictionsForce[i]	 += 0.5
* kineticFrictionForce; rod->kineticFrictionsForce[i+1] += 0.5 *
kineticFrictionForce;
                                        //rod->kineticFrictionsTorque[i]  +=
Q*(armTorque * kineticFrictionForce);

                                        // Note that the static friction
genereta torques as well.
                                        // Since they are annoying to compute
and require more knowledge of the rod's configuration
                                        // their calculation is in the rod's
class itself. But they are accounted for!!
                                        rod->staticFrictionsAxialForceForward[i]
+= 0.5 * staticFrictionForce; rod->staticFrictionsAxialForceForward[i+1] += 0.5
* staticFrictionForce;
                                }// Friction in axial direction

                        }//if within range
                  }//for loop all elements
              }//for loop all rods

        }//apply force
};



*/

// Anisotrpic frictional plane, for the moment set at y=0, rod free for y>r.
class AnisotropicFrictionPlaneInteraction : public Interaction {
 protected:
  const REAL kPlane;  // Elastic response to prevent interpenetration rod-plane
  const REAL
      etaPlane;  // Dissipative term appliet to the elastice response rod-plane
  const REAL muKineticForward;
  const REAL muKineticBackward;
  const REAL muKineticSideways;
  const REAL muStaticForward;
  const REAL muStaticBackward;
  const REAL muStaticSideways;
  const REAL vStaticPlane;
  const Vector3 normalPlane;
  const Vector3 originPlane;

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

  const REAL _sigmoid(const REAL _vel, const REAL _velThresold) const {
    const REAL vel = fabs(_vel);
    const REAL velThresold = fabs(_velThresold);

    const REAL velDiff = vel - 1.5 * velThresold;
    return 0.5 + 0.5 * erf(-(8.0 / velThresold) * velDiff);
  }

 public:
  AnisotropicFrictionPlaneInteraction(
      Vrodptr &rodptrs, const Vector3 _normalPlane, const Vector3 _originPlane,
      const REAL _kPlane, const REAL _etaPlane, const REAL _muKineticForward,
      const REAL _muKineticBackward, const REAL _muKineticSideways,
      const REAL _muStaticForward, const REAL _muStaticBackward,
      const REAL _muStaticSideways, const REAL _vStatic = 1e-4)
      : Interaction(rodptrs),
        normalPlane(_normalPlane),
        originPlane(_originPlane),
        kPlane(_kPlane),
        etaPlane(_etaPlane),
        muKineticForward(_muKineticForward),
        muKineticBackward(_muKineticBackward),
        muKineticSideways(_muKineticSideways),
        muStaticForward(_muStaticForward),
        muStaticBackward(_muStaticBackward),
        muStaticSideways(_muStaticSideways),
        vStaticPlane(_vStatic) {}

  void applyForces(const REAL time) {
    const int nor = rodptrs.size();
    for (unsigned int j = 0; j < nor; j++) {
      Rod *rod = rodptrs[j];

      vector<Vector3> overallRodForces =
          vector<Vector3>(rod->totalInternalForces.size() - 1);
      vFromPointsToElements(rod->totalInternalForces + rod->externalForces,
                            overallRodForces);

      for (unsigned int i = 0; i < (unsigned int)rod->n; i++) {
        // Fetch initial useful quantities
        const REAL r = rod->r[i];
        const Vector3 elementX = 0.5 * (rod->x[i] + rod->x[i + 1]);
        const REAL distanceFromPlane = (elementX - originPlane) % normalPlane;

        // Numerical fluctuations may cause the distnace from plane to be larger
        // than r. We ameng this problem by allowing a tolerance proportinal
        // (arbitrary constant) to the radius. For more robustness check the
        // SUGGESTION in Rod.cpp->applyFrictions
        if ((distanceFromPlane - r) < 1e-4) {
          // Fetch remaining quantites
          const Vector3 elementV = 0.5 * (rod->v[i] + rod->v[i + 1]);
          const Vector3 axialDir =
              (rod->edge[i] / rod->edge[i].length()).unitize();
          const Vector3 w = rod->w[i];
          const Matrix3 Q = rod->Q[i];

          // Compute plane tangent and bitangent (rolling) direction
          const Vector3 rollingDirection = (normalPlane * axialDir).unitize();
          const Vector3 axialDirection =
              (rollingDirection * normalPlane).unitize();

          const REAL currentRodForcesInNormaldirectionSign =
              (overallRodForces[i] % normalPlane);
          const Vector3 currentRodForcesInNormaldirection =
              (currentRodForcesInNormaldirectionSign < 0.0)
                  ? currentRodForcesInNormaldirectionSign * normalPlane
                  : 0.0 * normalPlane;
          const Vector3 elasticPlaneResponse =
              -kPlane * min(distanceFromPlane - r, 0.0) * normalPlane;
          const Vector3 dampingForcePlane =
              -etaPlane * (elementV % normalPlane) * normalPlane;
          const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
          rod->externalForces[i] +=
              0.5 *
              (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
          rod->externalForces[i + 1] +=
              0.5 *
              (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
          const REAL normalForce = normalForcePlane.length();

          // Compute slip velocity decomposed in rolling and tangential
          // directions
          const Vector3 armTorque = -r * normalPlane;
          const Vector3 translationVelInRollingDirection =
              (elementV % rollingDirection) * rollingDirection;
          const Vector3 rotationalVelInRollingDirection =
              ((Q.T() * (w * (Q * armTorque))) % rollingDirection) *
              rollingDirection;
          const Vector3 slipVelInRollingDirection =
              translationVelInRollingDirection +
              rotationalVelInRollingDirection;

          // Friction in rolling direction
          {
            const REAL f =
                _linear(slipVelInRollingDirection.length(), vStaticPlane);

            const Vector3 kineticFrictionForce =
                -(1.0 - f) * muKineticSideways * normalForce *
                slipVelInRollingDirection.unitize();
            const Vector3 staticFrictionForce =
                f * muStaticSideways * normalForce * rollingDirection.unitize();

            rod->kineticFrictionsForce[i] += 0.5 * kineticFrictionForce;
            rod->kineticFrictionsForce[i + 1] += 0.5 * kineticFrictionForce;
            rod->kineticFrictionsTorque[i] +=
                Q * (armTorque * kineticFrictionForce);

            // Note that the static friction genereta torques as well.
            // Since they are annoying to compute and require more knowledge of
            // the rod's configuration their calculation is in the rod's class
            // itself. But they are accounted for!!
            rod->staticFrictionsRollingForce[i] += 0.5 * staticFrictionForce;
            rod->staticFrictionsRollingForce[i + 1] +=
                0.5 * staticFrictionForce;
            rod->staticFrictionsNormalPlane[i] += normalPlane;
          }

          // Compute axial velocity
          const REAL slipVelInAxiallDirectionSign = elementV % axialDirection;
          const Vector3 slipVelInAxiallDirection =
              slipVelInAxiallDirectionSign * axialDirection;

          // Friction in axial direction
          {
            const REAL f =
                _linear(slipVelInAxiallDirection.length(), vStaticPlane);

            // Here we assume that the head of a snake is locate at the first
            // point x[0]!! The axial direction is computed x[i+1]-x[i], that is
            // why there is a (slipVelInAxiallDirectionSign<0.0), and NOT
            // (slipVelInAxiallDirectionSign>0.0)
            const REAL kineticFrictionCoeff =
                (slipVelInAxiallDirectionSign < 0.0) ? muKineticForward
                                                     : muKineticBackward;
            const Vector3 kineticFrictionForce =
                -(1.0 - f) * kineticFrictionCoeff * normalForce *
                slipVelInAxiallDirection.unitize();

            // Again here we assume that the head of a snake is locate at the
            // first point x[0]!! That is why there is a minus, it is supposed
            // to flip the vector "axialDirection" so that it points in the
            // direction tail-to-head and not head-to-tail
            const Vector3 staticFrictionForceForward =
                -f * muStaticForward * normalForce * axialDirection.unitize();
            const Vector3 staticFrictionForceBackward =
                -f * muStaticBackward * normalForce * axialDirection.unitize();

            rod->kineticFrictionsForce[i] += 0.5 * kineticFrictionForce;
            rod->kineticFrictionsForce[i + 1] += 0.5 * kineticFrictionForce;
            // No torque for axial case
            // rod->kineticFrictionsTorque[i]  += Q*(armTorque *
            // kineticFrictionForce);

            // Note that the static friction genereta torques as well.
            // Since they are annoying to compute and require more knowledge of
            // the rod's configuration their calculation is in the rod's class
            // itself. But they are accounted for!!
            rod->staticFrictionsAxialForceForward[i] +=
                0.5 * staticFrictionForceForward;
            rod->staticFrictionsAxialForceForward[i + 1] +=
                0.5 * staticFrictionForceForward;
            rod->staticFrictionsAxialForceBackward[i] +=
                0.5 * staticFrictionForceBackward;
            rod->staticFrictionsAxialForceBackward[i + 1] +=
                0.5 * staticFrictionForceBackward;
          }
        }
        // else
        //{
        //	cout << "detach " << (distanceFromPlane-r) << endl;
        // exit(0);
        //}
      }
    }
  }
};

// Anisotrpic frictional plane, for the moment set at y=0, rod free for y>r.
class SlenderBodyTheoryEnvironment : public Interaction {
 protected:
  const REAL dynamicViscosity;

 public:
  SlenderBodyTheoryEnvironment(Vrodptr &rodptrs, const REAL _dynamicViscosity)
      : Interaction(rodptrs), dynamicViscosity(_dynamicViscosity) {}

  void applyForces(const REAL time) {
    const int rodsize = rodptrs.size();
    for (unsigned int j = 0; j < rodsize; j++) {
      Rod *rod = rodptrs[j];

      // Compute total current length
      REAL totalLength = 0.0;
      vector<REAL>::iterator iterA = rod->l.begin();
      for (iterA = rod->l.begin(); iterA != rod->l.end(); iterA++)
        totalLength += (*iterA);

      // Compute total average current radius
      REAL sumRadii = 0.0;
      vector<REAL>::iterator iterB = rod->r.begin();
      for (iterB = rod->r.begin(); iterB != rod->r.end(); iterB++)
        sumRadii += (*iterB);

      sumRadii /= rod->r.size();

      for (unsigned int i = 0; i < (unsigned int)rod->n; i++) {
        REAL factorSBT =
            -4.0 * M_PI * dynamicViscosity / (log(totalLength / rod->r[i]));
        const Vector3 v = 0.5 * (rod->v[i] + rod->v[i + 1]);
        const Vector3 t = (rod->edge[i] / rod->edge[i].length()).unitize();
        Vector3 stokesForce = Vector3();
        Vector3::projectionSBT(v, t, stokesForce);
        stokesForce *= factorSBT * rod->l[i];
        rod->externalForces[i] += 0.5 * stokesForce;
        rod->externalForces[i + 1] += 0.5 * stokesForce;
      }
    }
  }
};

// Interaction between a rod (snake) and substrate
class SnakePlane : public Interaction {
 public:
  REAL k, nu, vStatic, muForwardKinetic, muSidewaysKinetic, muBackwardKinetic,
      muForwardStatic, muSidewaysStatic, muBackwardStatic;

  SnakePlane(Vrodptr rodptrs, REAL k, REAL nu, REAL vStatic,
             REAL muForwardKinetic, REAL muSidewaysKinetic,
             REAL muBackwardKinetic, REAL muForwardStatic,
             REAL muSidewaysStatic, REAL muBackwardStatic)
      : Interaction(rodptrs),
        k(k),
        nu(nu),
        vStatic(vStatic),
        muForwardKinetic(muForwardKinetic),
        muSidewaysKinetic(muSidewaysKinetic),
        muBackwardKinetic(muBackwardKinetic),
        muForwardStatic(muForwardStatic),
        muSidewaysStatic(muSidewaysStatic),
        muBackwardStatic(muBackwardStatic) {}

  void applyForces(const REAL time) {
    /*
    assert(rodptrs.size()==1);

    Rod* rp = rodptrs[0];

    for (int i=0; i < rp->n; i++)
    {
            REAL r = rp->r[i];
            Vector3 mx = .5*(rp->x[i] + rp->x[i+1]);
            Vector3 mv = .5*(rp->v[i] + rp->v[i+1]);
            Vector3 t = rp->Q[i][2];
            //Vector3 w = rp->w[i];
            Matrix3 Q = rp->Q[i];

            if (mx.y < r)
            {
                    REAL nForce = (r - mx.y) * k;
                    rp->externalForces[i].y +=  (nForce - nu * mv.y)/2;
                    rp->externalForces[i+1].y +=  (nForce - nu * mv.y)/2;
                    rp->w[i] *= (1-nu/1000); // haxxx

                    Vector3 mvp = mv - Vector3(0, 1, 0) * (mv % Vector3(0,1,0));
    // velocity in plane; REAL vTot = mvp.length();

                    if (vTot > vStatic)
                    {
                            // kinetic friction first
                            rp->collided[i] = 1;
                            Vector3 fForce = Vector3();

                            REAL vAxial =  mvp % t;

                            if (vAxial > 0) // front/back friction
                                    fForce -= t * vAxial / vTot * nForce *
    muForwardKinetic; else fForce -= t * vAxial / vTot * nForce *
    muBackwardKinetic;

                            fForce -= (mvp - t * (mvp % t)) / vTot * nForce *
    muSidewaysKinetic; // sideways friction

                            rp->externalForces[i]     += fForce/2;
                            rp->externalForces[i + 1] += fForce/2;
                            rp->externalTorques[i]    -= Q * (Vector3(0, -r, 0)
    * fForce);
                    }
                    else
                    {
                            // rp->v[i] = Vector3();
                            // rp->v[i + 1] = Vector3();

                            REAL sForce = muForwardStatic * nForce;
                            rp->staticFrictions[i]     += sForce / 2;
                            rp->staticFrictions[i + 1] += sForce / 2;

                            // rp->collided[i] = 2;
                    }
            }
    }
    */
  }
};

#endif
