// -*- coding: utf-8 -*-
// @Author: Tejaswin Parthasarathy
// @Date:   2019-09-14 23:49:37
// @Last modified by:   Tejaswin Parthasarathy
// @Last Modified time: 2019-10-07 11:02:23

/*
 * Rod.cpp
 *
 *  Created on: Jun 26, 2014
 *      Author: mgazzola
 */

#include "Rod.h"
#include "Tolerance.h"

Rod::Rod(const int _n, vector<Vector3> _x, vector<Vector3> _v,
         vector<Matrix3> _Q, vector<Vector3> _w, vector<REAL> _l0,
         vector<Vector3> _intrinsic_k0, vector<Vector3> _intrinsicShearStrain0,
         const vector<REAL> _m, vector<REAL> _V, const REAL _density,
         vector<Matrix3> _J0, vector<Matrix3> _B, vector<Matrix3> _S,
         const REAL _nu, const REAL _relaxationNu, const bool _useSelfContact)
    : n(_n),
      x(_x),
      v(_v),
      Q(_Q),
      w(_w),
      l0(_l0),
      intrinsic_k0(_intrinsic_k0),
      intrinsicShearStrain0(_intrinsicShearStrain0),
      m(_m),
      V(_V),
      density(_density),
      J0(_J0),
      B0(_B),
      S0(_S),
      nu(_nu),
      time(0.0),
      dt(0.0),
      relaxationNu(_relaxationNu),
      COUNTER(0),
      selfKineticMu(0.0),
      selfStaticMu(0.0),
      selfKineticThreshold(0.0),
      zetaSoft(100.0),  // Not used, kept for legacy purposes
      zetaHard(0.1),    // Not used, kept for legacy purposes
      useSelfContact(_useSelfContact) {
  // Chack masses consistency, that is the first and last element set to m/2
  assert(m.size() >= 3);
  // assert( m.back()==m.front() );
  // assert( m.front()==m[1]/2.0 );

  // Check sizes consistencies
  assert((int)x.size() == n + 1);
  assert((int)v.size() == n + 1);
  assert((int)Q.size() == n);
  assert((int)w.size() == n);
  assert((int)l0.size() == n);
  assert((int)intrinsic_k0.size() == n - 1);
  assert((int)intrinsicShearStrain0.size() == n);
  assert((int)m.size() == n + 1);
  assert((int)V.size() == n);
  assert((int)J0.size() == n);
  assert((int)B0.size() == n - 1);
  assert((int)S0.size() == n);

  /*
    These variables are used for printing out
    and is not used in the core of the Cosserat rod solver
  */
  MaxHeight = 200.0;
  Tforce = 0.0;
  True_State = vector<REAL>(6);     // Not used, kept for legacy purposes
  Observe_State = vector<REAL>(6);  // Not used, kept for legacy purposes
  LC_State = vector<REAL>(6);       // Not used, kept for legacy purposes
  LC_State_rate = vector<REAL>(6);  // Not used, kept for legacy purposes

  // Accelerations
  a = vector<Vector3>(n + 1);  // --> (n+1)
  wDot = vector<Vector3>(n);   // --> (n)

  // Inverse inertia matrix
  J0inv = vector<Matrix3>(n);  // --> (n)
  vDiagI(J0, J0inv);           // --> (n)

  // Radius, length, voronaoi domain in reference configuration
  d0 = vector<REAL>(n - 1);  // --> (n-1)
  vMidAvgInterior(l0, d0);   // (n) --> (n-1)

  // Current radius, length, voronaoi domain
  r = vector<REAL>(n);      // --> (n)
  l = vector<REAL>(n);      // --> (n)
  d = vector<REAL>(n - 1);  // --> (n-1)
  d = d0;                   // --> (n-1)
  // Dilatations

  e = vector<REAL>(n);       // --> (n)
  e_old = vector<REAL>(n);   // --> (n)
  ed = vector<REAL>(n - 1);  // --> (n-1)

  // Rate of dilatation
  deldt = vector<REAL>(n);  // --> (n)

  // Edges
  edge = vector<Vector3>(n);      // --> (n)
  vDiff(x, edge);                 // edge=vDiff(x)
  edge_old = vector<Vector3>(n);  // --> (n)
  vDiff(x, edge_old);             // edge=vDiff(x)

  // Bending
  k0 = vector<Vector3>(n - 1);  // --> (n-1)
  k0_old = intrinsic_k0;
  kDiff0 = vector<Vector3>(n - 1);  // --> (n-1)
  k0_rate = vector<Vector3>(n - 1);
  bendingInternalTorques0 = vector<Vector3>(n - 1);  // --> (n-1)
  bendingTorques = vector<Vector3>(n);               // --> (n)

  // Shear
  shearStrain0 = vector<Vector3>(n);  // --> (n)
  shearStrain0_old = intrinsicShearStrain0;
  shearStrainDiff0 = vector<Vector3>(n);  // --> (n)
  shearStrain_rate = vector<Vector3>(n);
  shearInternalForces0 = vector<Vector3>(n);  // --> (n)
  shearTorques = vector<Vector3>(n);          // --> (n)
  shearForces = vector<Vector3>(n + 1);       // --> (n+1)

  //—————————————
  // Here we only include axial strain rate damping
  /*
  NOTE : Only the last entry in the shear strain matrix
  is non-zero and is 10^-3. This number will be used later on.
  Also note that this strain-rate damping matrix is going to be
  used only in three cases : ELBOW, WALKER and WING.
  */

  D0 = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3);
  //—————————————
  // Self collisions
  selfCollisionForces = vector<Vector3>(n + 1);  // --> (n+1)
  collided = vector<int>(n);                     // --> (n)

  // Total internal forces and torques
  totalInternalForces = vector<Vector3>(n + 1);  // --> (n+1)
  totalInternalTorques = vector<Vector3>(n);     // --> (n)

  // External forces and torques
  externalForces = vector<Vector3>(n + 1);  // --> (n+1)
  externalTorques = vector<Vector3>(n);     // --> (n)

  // Total forces and torques
  totalForces = vector<Vector3>(n + 1);  // --> (n+1)
  totalTorques = vector<Vector3>(n);     // --> (n)

  // Damping
  dampingForces = vector<Vector3>(n + 1);  // --> (n+1)
  dampingTorques = vector<Vector3>(n);     // --> (n)

  // Transport of w
  transportW = vector<Vector3>(n);  // --> (n)

  // Friction
  staticFrictionsAxialForceForward = vector<Vector3>(n + 1);   // --> (n+1)
  staticFrictionsAxialForceBackward = vector<Vector3>(n + 1);  // --> (n+1)
  staticFrictionsRollingForce = vector<Vector3>(n + 1);        // --> (n+1)
  staticFrictionsNormalPlane = vector<Vector3>(n);             // --> (n)

  // All these are not used here and are kept for legacy reasons.
  staticFrictionsAxialForceForward_R = vector<Vector3>(n + 1);   // --> (n+1)
  staticFrictionsAxialForceBackward_R = vector<Vector3>(n + 1);  // --> (n+1)
  staticFrictionsRollingForce_R = vector<Vector3>(n + 1);        // --> (n+1)
  staticFrictionsNormalPlane_R = vector<Vector3>(n);             // --> (n)

  kineticFrictionsForce = vector<Vector3>(n + 1);  // --> (n+1)
  kineticFrictionsTorque = vector<Vector3>(n);     // --> (n)
  rollingf = vector<REAL>(n);                      // --> (n)

  // external contact check
  // This is not used and is kept for legacy reasons.
  globalcheck = vector<int>(n, 0);

  // Allocate temporary vectors
  tempVV3_n = vector<Vector3>(n);
  tempVV3_n2 = vector<Vector3>(n);
  tempVV3_nplus = vector<Vector3>(n + 1);
  tempVV3_nminus = vector<Vector3>(n - 1);
  tempVV3_nminus2 = vector<Vector3>(n - 1);
  tempVREAL_nminus = vector<REAL>(n - 1);
  tempVREAL_n = vector<REAL>(n);
  tempVREAL_nplus = vector<REAL>(n + 1);
  tempVM3_nminus = vector<Matrix3>(n - 1);
  tempVM3_n = vector<Matrix3>(n);

  // Get ready
  update(0.0);  // Update takes simulation time, which in the constructor is 0.0
  computeEnergies();
}

Vector3 Rod::computeVelocityCenterOfMass() {
  Vector3 avgVel = Vector3();

  for (unsigned int i = 0; i < v.size(); i++) avgVel += v[i];

  avgVel /= (REAL)v.size();

  return avgVel;
}

Vector3 Rod::computeAngularVelocityCenterOfMass() {
  Vector3 avgAngVel = Vector3();

  for (unsigned int i = 0; i < w.size(); i++) avgAngVel += Q[i].T() * w[i];

  avgAngVel /= (REAL)w.size();

  return avgAngVel;
}

void Rod::computeShearStrain0() {
  /*
  In biological systems, many components (especially muscles) exhibit
  viscoelastic behavior. This requires the computed force to depend on
  strain-rates in addition to strains (corresponding to the pure elastic case).
  Our framework was not initially coded for computations of shear strain rates.
  In order to skip this legacy issue, here we hard-code dt (the time-step of
  the simulation), depending on the case, for only computing strain rates.
  This dt value is the same as the one set in the corresponding case source
  files. This legacy version will soon be replaced by a more elegant,
  user-friendly implementation.

  In the examples we show, this strain-rate computation is necessary only
  for the elbow, walker and the wing, and the flags below correspond to these
  cases.
  */

  /*
    This dt is the same for all test cases except walker and flag wing and is
    the same one for the ELBOW. Note that except for the cases of the elbow,
    walker and the wing, no other cases use this hard-coded dt and the
    corresponding strain rate.
  */
  REAL dt = 0.5e-7;

#ifdef FLAGWALKER
  dt = 0.75e-7;  // walker
#endif
#ifdef FLAGWING
  /*
     wing -- a remark here: the actual dt in wing is 0.5e-7. We
     use dt=0.75e-7 here, which was set when testing the
     numerical fidelity of the wing simulation. Nevertheless, this
     mismatch in the dt value is factored in the variable
    "factor" as detailed in the
    computeInternalShearForces0 function
  */
  dt = 0.75e-7;
#endif

  // Compute the shear strain in the material frame of reference and with
  // respect to the reference configuration
  const Vector3 onez = Vector3(0, 0, 1);
  v_a_times_b_divide_c_minus_d_equal_c(
      Q, edge, l0, onez,
      shearStrain0);  // shearStrain0 = Q * d/dS{r} - onez = Q * edge / l0 -
                      // onez; (note: edge must be normalized by l0 and not l if
                      // you want to capture extension! tangent=e/l0)

  // Compute the shear strain rate
  // Note that the sign is consistent as we use the definition of
  // shear strain above which is signed
  v_a_minus_b_equal_c(shearStrain0, shearStrain0_old, tempVV3_n);
  v_a_times_b_equal_c(tempVV3_n, 1.0 / dt, shearStrain_rate);
  shearStrain0_old = shearStrain0;
}

void Rod::computeCurvature0() {
  // dt value not used in our code, kept for legacy reasons. See curvature
  // rate calculation below
  const REAL dt = 0.75e-7;

  // Compute the curvature expressed in the material frame k = kx*d1 + ky*d2 +
  // kz*d3, and with respect to the reference configuration k = d/dS{Q}, that is
  // that the derivative is taken by dividing for d0!!
  //
  // Discrtetization:
  //		k[i] = log( Q[i]*Q^T[i-1] ) / d0
  //
  // The next three lines are equivalent to k0 = vLog(vRotDiff(Q)) / d0 --> but
  // no memory allocation (in-place)
  vRotDiff(Q, tempVM3_nminus);
  vLog(tempVM3_nminus, tempVV3_nminus);
  v_a_divide_b_equal_c(tempVV3_nminus, d0, k0);

  /*
    Here, we do not include any viscous (damping) effects while considering the
    bending and twist modes. We however calculate the curvature/twist rate here
    as k0_rate, but do not use it (we multiply it by a 0.0 in
    computeInternalBendingTorques0 function). This code allows
    for computation of curvature strain-rate if needed, depending upon the
    material properties. We once again stress that in all cases shown in the
    paper, the viscous effects due to bending/twist are unimportant and hence we
    do not include it.
  */
  v_a_minus_b_equal_c(k0, k0_old, tempVV3_nminus);
  v_a_times_b_equal_c(tempVV3_nminus, 1.0 / dt, k0_rate);  // the curvature rate
  k0_old = k0;
}

void Rod::computeInternalShearForces0(const REAL time) {
  // Compute the internal shear forces n in the material frame and with respect
  // to the reference configuration! It is important to notice that here we use
  // a linear contitutive relation based on the definition of engineering
  // stress-strain curve, opposite to the true strain-stress curve! Engineering
  // stress-strain --> stress measured experimentally as stress = F/(E*A0) where
  // A0 is the area PRIOR applying the load True stress-strain --> stress
  // measured experimentally as stress = F/(E*A) where A is the area AFTER
  // applying the load Therefore since we use the engineering stress-strain, we
  // do NOT use the current S = S0/e, but we use the S0 PRIOR the load

  // EVERYTHING INCLUDED BETWEEN (1) AND (2) IS KIND OF A LENGTHY TECHNICALITY
  // THAT CAN BE SKIPPED. WHAT IS IMPORTANT ARE THE 10 LINES AFTER (2) (1)
  // ----------------------------------------------------------------------------------

  // Strain rate damping factor for different cases -- not all cases are needed

  /*
    Here the factor controls the viscous behavior of the material in
    stretch modes. This is important as in biological systems, many
    components (especially muscles) exhibit viscoelastic behavior. For more
    details on the computation of this factor and how it is used, see comments
    that follow.

    The use of strain rate damping effects has been originally introduced for
    the elbow characterization (Figure 1(d) of Nat. Comm. 2019 paper) and for
    the wing model (in which due to lack of available data we employed the elbow
    muscle model as reported in the paper.) In this legacy version of the code,
    the strain-rate damping is also used for the walker although its effect is
    negligible as detailed in the comments appropriate to WALKER section below.
    */
  REAL factor =
      0.0;  // The default setting is that the strain rate damping is not used!

///////////////////////////////////////////////////////////////////////////////
/*
This block of code represents the ramp-up in the viscous factor from the initial
time of simulation, varying based on the case. The ramp-up is necessary to avoid
numerical artifacts due to simulation-start up.

Please note that we limit the factor by a maximal amount in the next block of
code, which will be reached in less than 0.1 seconds after the start of the
simulation (for reference, elbow takes 4 seconds for the entire simulation).
*/

/*
THESE NUMBERS ARE SEEMINGLY LARGE, BUT AS DETAILED BELOW, THIS IS
BECAUSE OF THE UNITS (kg,m,s,N) or (g, mm, s, microN) USED IN THE CORRESPONDING
TESTCASES.
*/
#ifdef FLAGELBOW
  factor = (time - 0.01) * 330000;  // Elbow
#endif
#ifdef FLAGWALKER
  factor = (time - 0.01) * 200000.0;  // walker
#endif
#ifdef FLAGWING
  factor = (time - 0.01) * 2e10;  // wing
#endif
  factor = (factor < 0.0) ? 0.0 : factor;

  ///////////////////////////////////////////////////////////////////////////////
/*
This block of code limits the maximum viscous factor allowed, on a case by case
basis. We explain the rationale of this limit for the different cases below.
Here, we explain how the viscous behavior is calculated in the first place.

                                            Only stretch modes
                                                  |
                                                  |
    Max strain rate damping (\zeta) = factor * D0[2,2]

. For definition of D0 above, see constructor call below:

    D0 = Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 1.0e-3);

  It is important to note that damping only manifests for stretching modes, and
not the shearing modes. Note that all the factors above are multiplied by 1e-3!

    Max strain rate damping has a unit of kg m/s (or) g mm/s depending on the
case. Upon multiplication of \zeta by the strain-rate, we recover the viscous
forces (with appropriate units)

*/
#ifdef FLAGELBOW
  /*
    Elbow:

    In this case, the factor is only applied to muscles, tendons and bones. We
    note that the tendons and bones are very rigid and viscous damping plays a
    negligible role (this was well tested in multiple scenarios).

    In elbow case, the units used are in kg, m and s. Hence in this case \zeta
    (defined above) is calculated as:

               factor    D
                 |       |
       \zeta = 33000 * 1e-3 = 33 kg m/s.

    Note that this is for a single contractile rod in the muscle. In the elbow
    simulations, we employed 36 rods for constructing the biceps, which means
    the aggregate viscous damping factor \zeta is 33 * 36 = 1188 kg m/s. We note
    that this value is consistent with those reported in Figure 1 (d) of our
    Nat. Comm., 2019 paper for biological muscle tissue. This motivates the
    maximum value set for factor below.

    We note that the strain rate damping amounts to ~20-25% of the maximum
    forces exerted by the biceps as shown below. Based on the data from Figure 1
    (d) of our Nat. Comm., 2019 paper. This can be estimated as:

    (To understand the following explanation, please have the figure 1 (d)
    open.)

    Strain_max is 0.15 (At 60 degrees, \eta = 0.85, so that Strain_max = 1 -
    \eta  = 0.15).
    Time of flexion : 60 degrees / (Angular Velocity) = 60/100 = 0.6 s
    (Angular velocity = 100 corresponds to strain rate damping of 1188 kg m/s)

    Strain rate = Strain_max / Time_min = 0.15/0.6 = 0.25.

    So the corresponding damping force = (1188 * 0.25) = 297 N, which is
    ~20—25% of the static bicep force.

  */
  factor = (factor > 33000) ? 33000 : factor;
#endif
#ifdef FLAGWALKER
  /*
  Walker:

  In this case, the factor is applied to muscles and scaffold. The scaffold
  is made of PEGDA material, which has visco-elastic behavior as well. Hence the
  rationale to include viscous effects in our model.

  Given the lack of characterization of the material properties, we set
  the damping forces to be small ( ~1% of the contractile muscle forces, see
  calculations below). Still the code allows one to account for these effects
  in case accurate data become available.

  In the walker case, the units used are in g, mm and s. Hence in this case
  \zeta is calculated as:

             factor    D
               |       |
     \zeta = 10000 * 1e-3 = 10 g mm/s.

  Here we set the maximum factor based on engineering estimates from the
  walker design, as explained below:

  The maximum strain in the case of the walker is Strain_max = 0.13,
  see last point on SI note 7.2 in Nat. Comm. 2019 paper.

  As the walker is actuated at 2Hz (or T = 0.5s), the characteristic relaxation
  time-scale is Time_min = T/2 = 0.25s (i.e. the time for it to reach peak
  contraction). Thus the upper bound on the estimated strain rate is then

    StrainRate_max = Strain_max/Time_min = 0.13/0.25 = 0.52 s^{-1}.

  Upon multiplication by \zeta = 10 g mm/s calculated above, this gives a
  viscous damping force of

  \zeta * StrainRate_max = 10 * 0.52 = 5.2 micro N

  which is small ( ~1% ) when compared to the muscle rod contractile force (600
  micro N).
*/
  factor = (factor > 10000.0) ? 10000.0 : factor;
#endif
#ifdef FLAGWING
  /*
Wing:

In this case, the factor is only applied to muscles and it has been set so as to
approximately match the elbow muscle behavior due to lack of better
characterization.

In the wing case, the units used are in g, mm and s. Hence in this case
\zeta is calculated as:

           factor    D
             |       |
   \zeta = 1e9 * 1e-3 = 1e6 g mm/s

The estimated maximum strain in the case of the wing is Strain_max ~ 0.4-0.5,
(as the wing flaps). As the time-period of actuation of wing is 0.38s,
the characteristic relaxation time-scale is Time_min = T/4 ~ 0.095s (i.e. the
time for it to reach peak contraction). Thus the upper bound on the estimated
strain rate is then

  StrainRate_max ~ Strain_max/Time_min = 0.5/0.095 = 5.26 s^{-1}.

The strain-rate calculation employed in the case of the wing is based on
dt=0.75e-7, but our simulation runs with dt=0.5e-7 (i.e. 1.5 times slower).
To account for this scaling, the overall damping force is now scaled as follows:

scalingfactor * \zeta * StrainRate_max = (1/1.5) * 1e6 * 5.26 = 3.5 N

The total muscle force in this case is ~20 N (see Wing.cpp). This means that the
ratio of viscous damping force to the total force is ~ 20%, similar to the
elbow.

NOTE : We note that the value of the forces calculated here are before
re-scaling the wing. Nevertheless, the above calculated percentage remains the
same.
*/
  factor = (factor > 1e9) ? 1e9 : factor;
#endif
  /* In the case of the elbow, we do not apply viscous strain-rate damping
  for the artificial muscle*/
  if (n == 108) factor = 0.0;
#ifdef FLAGWING
  /* In the case of the wing, we only apply viscous strain-rate damping
  for the muscle which all have n=7*/
  if (n != 7) factor = 0.0;
#endif

  // (2)
  // ----------------------------------------------------------------------------------

  const vector<Matrix3> DS_0 = vector<Matrix3>(n, factor * D0);

  v_a_minus_b_equal_c(shearStrain0, intrinsicShearStrain0,
                      shearStrainDiff0);  // shearStrainDiff0 = shearStrain0 -
                                          // intrinsicShearStrain0;
  v_a_times_b_equal_c(
      S0, shearStrainDiff0,
      tempVV3_n);  // shearInternalForces0 = S0 * shearStrainDiff0;
  v_a_times_b_equal_c(DS_0, shearStrain_rate, tempVV3_n2); /* viscous effects*/
  v_a_plus_b_equal_c(tempVV3_n, tempVV3_n2, shearInternalForces0);
}

void Rod::computeInternalBendingTorques0(const REAL time) {
  const vector<Matrix3> DB_0 = vector<Matrix3>(n - 1, 0.0 * D0);

  // Compute the bending internal torques m in the material frame of reference
  // and with respect to the reference configuration
  v_a_minus_b_equal_c(k0, intrinsic_k0, kDiff0);  // kDiff0 = k0 - intrinsic_k0;
  v_a_times_b_equal_c(
      B0, kDiff0,
      tempVV3_nminus);  // bendingInternalTorques0 = B0 * kDiff0;
  v_a_times_b_equal_c(DB_0, k0_rate,
                      tempVV3_nminus2); /* viscous effects, not used in our
                                           code, kept for legacy reasons*/
  v_a_plus_b_equal_c(tempVV3_nminus, tempVV3_nminus2, bendingInternalTorques0);
}

void Rod::computeInternalForces() {
  // Compute the right hand side of the linear momentum equation (excluding
  // external forces) in the lab frame of reference and with respect to the
  // actual configuration
  //
  // From actual configuration (function of s) to reference configuration
  // (function of S)
  // ---------------------
  //
  // Arc-lengths:
  // 		S = reference configuration arc-lenght
  // 		s = actual configuration arc-lenght
  //
  // Divergence of internal stresses gives the net force (per unit length!)
  // acting on the infinitesimal element of length ds:
  //		F = d/ds{n}*ds
  //		n = n0 -->	this is because we use the engineering
  // stress-strain
  // realtion (see above) 					if we were to be
  // using the true stress-strain relationsin we should apply the rescaling n =
  // n0/e but this introduce an asymmetry in the load-strain relation which
  // makes the system unstable (this is why we stick to the engineering
  // definition)
  //
  // Equalities:
  // 		ds = e*dS
  //		d/ds{} = (1/e)*d/dS{}
  // 		n0 = S0*shearStrain0 --> we use the engineering definition of
  // the stress-strain constitutive law, hence the use of S0 (check also the
  // function computeInternalShearForces() )
  //
  // Linear momentum balance:
  //		d/dt{rho*A*ds*v} = d/ds{n}*ds
  //		d/dt{rho*(A0/e)*(e*dS)*v} = (1/e)d/dS{n0}*(e*dS)
  //		rho*A0*dS*d/dt{v} = d/dS{n0}*dS
  //		rho*A0*dS*d/dt{v} = d{n0} --> note that d/dS{n}*dS = d{n} is
  // simply the difference between adjacent internal shear forces n[i]-n[i-1]
  //
  // --> THEREFORE EVERYTHING IS DONE WITH RESPECT TO REFERENCE CONFIGURATION
  // AND NOTHING CHANGES!!
  //
  // For convenience in updating the rod configuarion we want to express the
  // linear momentum balance in the lab frame of reference. But at the moment
  // the internal shear forces n are expressed in the material frame of
  // reference. How to compute F=d{n} in the lab frame of reference.
  // ---------------------
  //
  // n is computed in the material frame of reference, but we want the force F
  // in the lab frame:
  //		x_{mat} = Q*x_{lab}
  //		x_{lab} = Q^{-1}*x_{mat}
  //		Q*Q^T = Q^T*Q = I --> I is the identity matrix here
  //		Q^{-1} = Q^T
  //		n_{mat} = S*shearStrain_{mat}
  //		n_{lab} = Q^T*n_{mat}
  //		F_{lab} = d{Q^T*n}

  // Next lines are equvalent (but in-place) to: shearForces_{lab} =
  // vDelta(vT(Q)*shearInternalForces0_{mat});

  vT(Q, tempVM3_n);

  v_a_times_b_equal_c(tempVM3_n, shearInternalForces0 / e, tempVV3_n);

  vDelta(tempVV3_n, shearForces);
}

void Rod::computeInternalTorques() {
  // Compute the right hand side of the angular momentum equation (excluding
  // external torques) in the material frame of reference and wihrespect to the
  // reference configuration
  // ---------------------
  //
  // Arc-lengths:
  // 		S = reference configuration arc-lenght
  // 		s = actual configuration arc-lenght
  //
  // Divergence of internal torques plus the tangent cross internal shear forces
  // gives the torques (per unit length!) acting on the infinitesimal element of
  // length ds, and expressed in the lab frame of reference:
  //		d/dt{J*w} = d/ds{m}*ds + (d/ds{r} x n)*ds
  //		d/dt{rho*(A*ds)*(A/4pi)*w} = d/ds{m}*ds + (d/ds{r} x n)*ds
  //
  // Fact:
  //		C = c*d1 = cx*d1 + cy*d2 + cz*d3
  //		d/dt{di} = w x di --> by definition
  //		d/dt{C} = d/dt{c*d} = d/dt{c}*di + c*d/dt{di}
  //		d/dt{C} = d/dt{C} + w x C
  //
  //		C = c*d1 = cx*d1 + cy*d2 + cz*d3
  //		d/ds{di} = k x di --> by definition
  //		d/ds{C} = d/ds{c*d} = d/ds{c}*di + c*d/ds{di}
  //		d/ds{C} = d/ds{C} + k x C
  //
  // From lab frame to material frame of reference:
  //		d/dt{J*w} + w x J*w = d/ds{m}*ds + k x m + (d/ds{r} x n)*ds
  //		d/dt{J*w} = d/ds{m}*ds + (k x m)*ds + (Q*d/ds{r} x n)*ds - w x
  // J*w 		d/dt{J*w} = d/ds{ B*d/ds{Q} }*ds + (d/ds{Q} x B*k) +
  // (Q*d/ds{r} x n)*ds
  //- w x J*w
  //
  // From actual configuration (function of s) to reference configuration
  // (function of S):
  //
  // 		ds = e*dS
  //		d/ds{} = (1/e)*d/dS{}
  // 		A = A0/e
  // 		J ~ rho*(A*ds)*(A/4pi) = J0/e
  //		k = d/ds{Q} = (1/e)d/dS{Q} = k0/e
  //		B ~ G*I = G*A*A = G*A0/e2 = B0/e2
  //		n = n0 -->	this is because we use the engineering
  // stress-strain
  // realtion (see above) 					if we were to be
  // using the true stress-strain relationsin we should apply the rescaling n =
  // n0/e but this introduce an asymmetry in the load-strain relation which
  // makes the system unstable (this is why we stick to the engineering
  // definition) n0 = S0*shearStrain0 --> we use the engineering definition of
  // the stress-strain constitutive law, hence the use of S0 (check also the
  // function computeInternalShearForces() )
  //
  //		d/dt{J*w} = d/ds{ B*d/ds{Q} }*ds + (d/ds{Q} x B*d/ds{Q})*ds +
  //(Q*d/ds{r} x n)*ds - w x J*w --> 	note that the mass inertia terms are NOT
  // multiplied by ds, because it
  // is already contained in its definition J = rho*I*ds (where I is the area
  // moment here) 		d/dt{J0/e*w} =	(1/e)*d/dS{
  //(B0/e2)*(1/e)*d/dS{Q} }*e*dS + ( (1/e)*d/dS{Q} x B0/e2*(1/e)*d/dS{Q} )*e*dS
  //+ 						(Q*(1/e)*d/dS{r} x n0)*e*ds -
  // w x J0/e*w 		J0*d/dt{w/e} =	d{B0*k0/e3} + (k0 x B0*k0/e3)*dS
  // + Q*d{r} x n0 - w x (J0/e)*w
  //

  // Compute the terms d{B0*k0/e3} + (k0 x B0*k0/e3)*dS
  // Next lines are equivalent (but in-place) to bendingTorques =
  // vDelta(bendingInternalTorques0/ed^3) +
  // vMidAvg(d0*k0*bendingInternalTorques0/ed^3);
  bendingInternalTorques0 = bendingInternalTorques0 / (ed * ed * ed);
  vDelta(bendingInternalTorques0, tempVV3_n);
  v_a_times_b_cross_c_equal_d(d0, k0, bendingInternalTorques0,
                              tempVV3_nminus);  // d = a * b * c
  vMidAvg(tempVV3_nminus, bendingTorques);
  bendingTorques += tempVV3_n;  // in-place operator

  // Compute the term +Q*d{r} x n0 = -n0 x Q*d{r}
  // Next lines are equvalent (but in-place) to shearTorques =
  // -shearInternalForces0*(Q*edge); // here there is a minus because in my
  // derivation is +Q*r' x S*shearStrain = - S*shearStrain x Q*r'
  v_a_times_b_equal_c(Q, edge, tempVV3_n);  // edge = et*dS
  v_a_times_b_cross_c_equal_d(-1.0, shearInternalForces0 / e, tempVV3_n,
                              shearTorques);

  // Compute the term - w x J0/e*w = J0*w/e x w
  // Term relative to the transport of angular velocity in the angular momentum
  // balance due to the use of material frame of reference (J0*w/e) x w
  transportW = (J0 * w / e) * w;
}

// Numerical dumping mimicking viscous effects of the material
void Rod::computeDampingForces() {
  // This damping model is in the spirit of Coulomb dissipation model (f=-nu*v),
  // but accounts only for internal dissipation by using relative velocities
  // between elements. This way rigid body motions are not artificially slowed
  // down.

  assert(nu >= 0.0);
  assert(relaxationNu >= 0.0);
  assert(dt > 0.0);

  // if the decaying time is zero than nu is constant, otherwise it decays
  // exponentially.
  /*
      Note that we only use relaxationNu in the case of the walker. We
     initialize the walker slightly above ground and then let it fall, to reach
     its initial rest position before it starts moving. During this transient we
     use this exponentially decaying nu, to damp out vibrations due to the fall.

      In all the other cases, relaxationNu is 0.0 and hence a constant,
     time-invariant nu is used, as set by the user.
  */
  REAL nuNow = (relaxationNu == 0.0 || nu == 0.0)
                   ? nu
                   : (0.13 + nu * exp(-time / relaxationNu));
  nuNow = (nuNow < Tolerance::tol()) ? 0.0 : nuNow;

  // Rescale since nu is per unit length
  const REAL ds0 = vSum(l0) / dampingTorques.size();
  nuNow *= ds0;

  // Think of something better!
  for (unsigned int i = 0; i < dampingForces.size(); i++) {
    const REAL factor = (i == 0 || i == (dampingForces.size() - 1)) ? 0.5 : 1.0;
    dampingForces[i] -= (nuNow * factor) * v[i];
  }

  // Think of something better!
  for (unsigned int i = 0; i < dampingTorques.size(); i++)
    dampingTorques[i] -= nuNow * w[i];
}

void Rod::computeAllInternalResultingForcesAndTorques() {
  computeInternalForces();   // in-place
  computeInternalTorques();  // in-place

  if (useSelfContact) computeSelfCollisionForces();  // empty for the moment

  computeDampingForces();  // in-place

  v_a_plus_b_plus_c_equal_d(
      shearForces, selfCollisionForces, dampingForces,
      totalInternalForces);  // (in-place) totalInternalForces = shearForces +
                             // selfCollisionForces + dampingForces;
  v_a_plus_b_plus_c_equal_d(
      bendingTorques, shearTorques, dampingTorques,
      totalInternalTorques);  // (in-place) totalInternalTorques =
                              // bendingTorques + shearTorques + dampingTorques
  totalInternalTorques += transportW;
}

/* This function is never used in our cases and retained for legacy purposes
 */
void Rod::computeSelfCollisionForces() {
  const REAL zeta = 1e4;
  const REAL zeta_nu = 10.0;

  for (int i = 0; i < (int)n; i++) {
    assert(r[i] > 0.0);
    assert(l[i] > 0.0);

    // don't calculate any global forces arising from other segments within
    // distance ~= PI*R of this one b/c the other energy terms already account
    // for the local curvature bounded regime

    // int skip = ceil(.5*PI*r[i]/l[i]);
    const int skip = ceil(.8 * M_PI * r[i] / l[i]) + 1;
    assert(skip >= 1);
    for (int j = i - skip; j > -1; j--) {
      // if outside bounding box don't bother
      const Vector3 x1 = x[i];
      const Vector3 x2 = x[j];
      const REAL r1 = r[i];
      const REAL r2 = r[j];
      const REAL sum_r1r2 = r1 + r2;
      if ((x1 - x2).length() <= (sum_r1r2 + l[i] + l[j])) {
        // find the shortest line segment between the two centerline segments
        const Vector3 edge1 = edge[i];
        const Vector3 edge2 = edge[j];
        const vector<Vector3> minVectors =
            findMinDistVectors(x1, edge1, x2, edge2);
        const Vector3 dVector = minVectors[1];
        const Vector3 dVectorDirection = dVector.unitize();

        // gamma tells you whether ther is overlap
        // const REAL gamma = (sum_r1r2 - dVector.length())/sum_r1r2;
        const REAL gamma = (sum_r1r2 - dVector.length());
        const bool yesCollision = (gamma > 0.0);

        // Set collision flags
        collided[i] = collided[j] = yesCollision;

        // Compute contact force
        const REAL cForceContact = zeta * gamma;

        // Compute damping force
        const Vector3 vInterpenetration =
            ((v[i] + v[i + 1]) - (v[j] + v[j + 1])) / 2.0;
        const REAL vNorm = vInterpenetration % dVectorDirection;
        const REAL cForceDamping = zeta_nu * vNorm;

        // Compute collision force (the 0.5 in fron is because it is distributed
        // onto the two masses at the end of an edge)
        const Vector3 cForce = yesCollision *
                               (0.5 * cForceContact + cForceDamping) *
                               dVectorDirection;

        // Remenber that the first and last points have half the mass! Hence the
        // strange coefficient in fron of cForce
        selfCollisionForces[i] -= ((i == 0) ? 0.5 : 1.0) * cForce;
        selfCollisionForces[i + 1] -= ((i + 1 == n) ? 0.5 : 1.0) * cForce;
        selfCollisionForces[j] += ((j == 0) ? 0.5 : 1.0) * cForce;
        selfCollisionForces[j + 1] += ((j + 1 == n) ? 0.5 : 1.0) * cForce;
      }
    }
  }
}

void Rod::applyStaticFrictions() {
  // AXIAL STATIC
  // -------------------------------

  // Apply static friction forces in the axial direction
  for (unsigned int i = 0; i < totalInternalForces.size(); i++) {
    // Rod surface friction
    // Get forward axial direction in the tail-to-head sense
    Vector3 maxFrictionDirection =
        staticFrictionsAxialForceForward[i].unitize();

    // Compute projection of total forces in the axial (or tangential to the
    // plane) direction
    REAL projectionSign =
        ((totalInternalForces[i] + externalForces[i]) % maxFrictionDirection);
    Vector3 projection = projectionSign * maxFrictionDirection;

    // If total forces are pushing forward, so (projectionSign>0.0) then the max
    // friction is given by staticFrictionsAxialForceForward, otherwise it is
    // given by staticFrictionsAxialForceBackward
    Vector3 maxFriction = (projectionSign > 0.0)
                              ? staticFrictionsAxialForceForward[i]
                              : staticFrictionsAxialForceBackward[i];

    // Compute static friction force, which opposes to the pushing forces
    Vector3 forceStatic =
        -projection.unitize() * min(projection.length(), maxFriction.length());

    // Remove friction component from total force
    externalForces[i] += forceStatic;
  }

  /*
          // This part compute the torque induced by the friction at the bottom
          // In principle this is counter balancd by the rection from the
  substrate and thereofre it is commented out
          // Also it couses instabilities to the code

          // Store force static temporarily to be used to compute torques
          tempVV3_nplus[i] = forceStatic;
  }

  // Apply static friction torques in the axial direction
  vFromPointsToElements(tempVV3_nplus, tempVV3_n);
  for (unsigned int i=0; i<totalInternalTorques.size(); i++)
  {
          // Correct torques
          const Vector3 normalPlaneDirection =
  staticFrictionsNormalPlane[i].unitize(); const Vector3 arm =
  (-r[i]*normalPlaneDirection); externalTorques[i] += Q[i]*(arm * tempVV3_n[i]);
  }
   */

  // ROLLING STATIC
  // -------------------------------

  // From point-wise forces to element forces and masses
  vFromPointsToElements(totalInternalForces + externalForces, tempVV3_n);

  // The static forces have been distributed over the mass points that
  // constitute the rod in the Interaction.h class. Here we translate those
  // forces to the forces acting on the cylindrical elements. SUGGESTION: To
  // ensure robustness, whenever one mass point for numerical reasons one might
  // consider to take the max of the mass point static elements
  vector<Vector3> overallStaticRollingForces =
      vector<Vector3>(staticFrictionsRollingForce.size() - 1);
  vFromPointsToElements(staticFrictionsRollingForce,
                        overallStaticRollingForces);

  // Not used in our cases, kept for legacy reasons
  vector<Vector3> overallStaticRollingForces_R =
      vector<Vector3>(staticFrictionsRollingForce_R.size() - 1);
  vFromPointsToElements(staticFrictionsRollingForce_R,
                        overallStaticRollingForces_R);

  // The point mass velocities are translated into the velocities of cylindrical
  // elements
  vector<Vector3> overallVel = vector<Vector3>(v.size() - 1);
  fill(overallVel.begin(), overallVel.end(), Vector3(0.0, 0.0, 0.0));
  vMidAvgInterior(v, overallVel);

  // Compute force generated by current torque on the surface and store them in
  // the temporary container tempVV3_nplus
  for (unsigned int i = 0; i < totalInternalTorques.size(); i++) {
    // Rod surface static
    // Extract useful quantities
    const Vector3 maxRollFriction = overallStaticRollingForces[i];
    const Vector3 maxRollFrictionDirection = maxRollFriction.unitize();
    const Vector3 normalPlaneDirection =
        staticFrictionsNormalPlane[i].unitize();
    assert(fabs(normalPlaneDirection % maxRollFrictionDirection) <
           Tolerance::tol());
    const Vector3 tangentPlaneDirection =
        (normalPlaneDirection * maxRollFrictionDirection).unitize();
    const REAL R = r[i];
    const Vector3 arm = (-r[i] * normalPlaneDirection);

    // Compute torque in the plane tangetial direction
    const Matrix3 materialFrame = Q[i];
    const Vector3 torqueInMaterialFrame =
        totalInternalTorques[i] + externalTorques[i];
    const Vector3 torqueInLabFrame = materialFrame.T() * torqueInMaterialFrame;
    const Vector3 torqueInTangentDirection =
        (torqueInLabFrame % tangentPlaneDirection) * tangentPlaneDirection;

    // Project total force acting on the element in the rolling direction
    const Vector3 forceElementRolling =
        (tempVV3_n[i] % maxRollFrictionDirection) * maxRollFrictionDirection;

    // Compute static force
    const REAL F = (forceElementRolling % maxRollFrictionDirection);
    const REAL M = (torqueInTangentDirection % tangentPlaneDirection);

    const Vector3 forceNoSlip =
        -(R * F - 2.0 * M) / (3.0 * R) * maxRollFrictionDirection;

    // Compute static friction
    const Vector3 forceStatic =
        forceNoSlip.unitize() *
        min(forceNoSlip.length(), maxRollFriction.length());

    // Correct forces
    externalForces[i] += 0.5 * forceStatic;
    externalForces[i + 1] += 0.5 * forceStatic;

    // Correct torques
    externalTorques[i] += materialFrame * (arm * forceStatic);

    // This part enforces a correction on velocities to exactly meet no slip
    // when the right conditions are satisfied. It is not necessary and it is
    // actually potentially harful, since it abruptly modifies the dynamics of
    // an element inducing strong reaction forces. Nevertheless I leave this
    // piec of code here, just in case, but I do not advice to use it
    // if ( (forceNoSlip.length() - maxRollFriction.length() ) <=
    // Tolerance::tol() )
    //{
    //	const Vector3 translationVelInRollingDirection = (overallVel[i] %
    // maxRollFrictionDirection)*maxRollFrictionDirection; 	const Vector3
    // rotationalVelInRollingDirection =
    //((materialFrame.T()*(w[i]*(materialFrame*arm))) %
    // maxRollFrictionDirection)*maxRollFrictionDirection; 	const Vector3
    // slipVelInRollingDirection = translationVelInRollingDirection +
    // rotationalVelInRollingDirection;
    //	//const REAL factor = min(1.0,max(rollingf[i], 0.0));
    //	//const Vector3 slipVelCorrection =
    //-factor*slipVelInRollingDirection/2.0; 	const Vector3 slipVelCorrection
    //= -slipVelInRollingDirection/2.0;

    //	const Vector3 correction_v_i = ((i>0)?0.5:1.0)*slipVelCorrection;
    //	const Vector3 correction_v_i_plus =
    //(((int)i<(n-1))?0.5:1.0)*slipVelCorrection; 	const Vector3
    // correction_torque = materialFrame*(slipVelCorrection *
    // normalPlaneDirection)/R;

    //	v[i] += correction_v_i;
    //	v[i+1] += correction_v_i_plus;
    //	w[i] += correction_torque;
    //}
  }

  // KINETIC FRICTION (BOTH AXIAL AND ROLLING)
  // -------------------------------

  externalForces += kineticFrictionsForce;
  externalTorques += kineticFrictionsTorque;
}

void Rod::update(const REAL time) {
  // Update edges
  const vector<Vector3> onezVec = vector<Vector3>(n, Vector3(0, 0, 1));
  vDiff(x, edge);  // edge=vDiff(x)

  // Update current lengths
  vLength(edge, l);  // l=vLength(edge)

  // Update dilatation edges with respect to reference configuration
  e_old = e;
  e = l / l0;

  // Update radius accounting for stretching under the assumption of
  // incompressibility and cylindrical cross sections
  v_sqrt_a_dividepar_b_times_c_equal_d(V, l, M_PI,
                                       r);  // r = vSqrt( V/(l * M_PI) )

  // Update current voronoi domains
  vMidAvgInterior(l, d);  // size = n-1

  // Update voronoi dilatation
  ed = d / d0;

  // Compute curvature and shear strain with respect to the reference
  // configuration
  computeCurvature0();
  computeShearStrain0();

  // Compute internal shear forces and bending torques with respect to the
  // reference configuration
  computeInternalBendingTorques0(time);
  computeInternalShearForces0(time);
}

// This is not used in our cases and retained for legacy reasons
void Rod::computeBendingEnergy() {
  bendingEnergy = .5 * vSum((kDiff0 % bendingInternalTorques0) * d0);
}

// This is not used in our cases and retained for legacy reasons
void Rod::computeShearEnergy() {
  shearEnergy = 0.5 * vSum(shearStrainDiff0 % (S0 * shearStrainDiff0) * l0);
}

// This is not used in our cases and retained for legacy reasons
void Rod::computeTranslationalEnergy() {
  translationalEnergy = .5 * vSum(m * (v % v));
}

// This is not used in our cases and retained for legacy reasons
void Rod::computeRotationalEnergy() {
  rotationalEnergy = .5 * vSum(w % (J0 * w / e));
}

// This is not used in our cases and retained for legacy reasons
void Rod::computeEnergies() {
  computeBendingEnergy();
  computeShearEnergy();
  computeTranslationalEnergy();
  computeRotationalEnergy();
  totalInternalEnergy =
      bendingEnergy + rotationalEnergy + translationalEnergy + shearEnergy;
}

bool Rod::nanCheck() {
  for (unsigned int i = 0; i < (unsigned int)n; i++)
    if (x[i].x != x[i].x) {
      cout << i << endl;
      return true;
    }
  return false;
}

void Rod::reset() {
  // some quantities naturally reset; we need only reset those defined via +=
  const Vector3 zero = Vector3(0.0, 0.0, 0.0);
  fill(collided.begin(), collided.end(), 0);
  fill(rollingf.begin(), rollingf.end(), 0);
  fill(a.begin(), a.end(), zero);
  fill(wDot.begin(), wDot.end(), zero);
  fill(dampingForces.begin(), dampingForces.end(), zero);
  fill(dampingTorques.begin(), dampingTorques.end(), zero);
  fill(externalForces.begin(), externalForces.end(), zero);
  fill(externalTorques.begin(), externalTorques.end(), zero);
  fill(selfCollisionForces.begin(), selfCollisionForces.end(), zero);
  fill(staticFrictionsAxialForceForward.begin(),
       staticFrictionsAxialForceForward.end(), zero);
  fill(staticFrictionsAxialForceBackward.begin(),
       staticFrictionsAxialForceBackward.end(), zero);
  fill(staticFrictionsRollingForce.begin(), staticFrictionsRollingForce.end(),
       zero);
  fill(staticFrictionsNormalPlane.begin(), staticFrictionsNormalPlane.end(),
       zero);

  // This is not used in our cases and retained for legacy reasons
  fill(staticFrictionsAxialForceForward_R.begin(),
       staticFrictionsAxialForceForward_R.end(), zero);
  fill(staticFrictionsAxialForceBackward_R.begin(),
       staticFrictionsAxialForceBackward_R.end(), zero);
  fill(staticFrictionsRollingForce_R.begin(),
       staticFrictionsRollingForce_R.end(), zero);
  fill(staticFrictionsNormalPlane_R.begin(), staticFrictionsNormalPlane_R.end(),
       zero);

  fill(kineticFrictionsForce.begin(), kineticFrictionsForce.end(), zero);
  fill(kineticFrictionsTorque.begin(), kineticFrictionsTorque.end(), zero);
}

/*
 NOTE : From now on, we just use functions for I/O and Povray visualization

 LONG and COMPLICATED, but not relevant to the mechanics solver.
 One can modify this as they like.
*/
void Rod::dumpPovray(string filePovray, string fileData, int colorindex,
                     REAL time) {
// This is for walker
#ifdef FLAGWALKER
  if ((colorindex == 0) || (colorindex == 2) || (colorindex == 3)) {
    fstream f;
    f.open(fileData.c_str(), fstream::app);

    f << "#declare Ribbon_Spline" << colorindex << " =\n";
    f << "spline { natural_spline\n";
    for (unsigned int i = 0; i < n; i++) {
      f << (-0.25 + 0.25 * i) << ", "
        << "<" << x[i][0] << "," << x[i][1] << "," << x[i][2] << ">,\n";
    }
    f << (-0.25 + 0.25 * n) << ", "
      << "<" << x[n][0] << "," << x[n][1] << "," << x[n][2] << ">\n";
    f << "}\n";
    f << "union{\n";
    f << "#local Nr = -0.25;\n";
    f << "#local EndNr = " << (-0.25 + 0.25 * n - 0.05) << ";\n";
    f << "#while (Nr <= EndNr)\n";
    if (colorindex == 0) {
      f << "box{ <-1.75,0,0>,\n";
      f << "<5.25,-0.3,0.68>\n";
      f << "pigment{ color rgbt<0.45,0.39,1,0>}\n";
      f << "rotate x*0\n";
    } else {
      const Vector3 direction = (x[6] - x[0]).unitize();
      const REAL Angle =
          angle(direction, Vector3(0.0, 0.0, -1.0)) * 180.0 / M_PI;
      f << "box{ <-1.75,1.75,0>,\n";
      f << "<5.25,-1.75,0.3>\n";
      f << "pigment{ color rgbt<0.45,0.39,1,0>}\n";
      if (x[6][1] < x[0][1]) {
        f << "rotate -x*" << Angle << "\n";
      } else {
        f << "rotate x*" << Angle << "\n";
      }
    }
    f << "translate Ribbon_Spline" << colorindex << "(Nr)\n";
    f << "scale<0.25,0.25,0.25>\n";
    f << "rotate<0,90,90>\n";
    f << "translate<2,0,4>\n";
    f << "}\n";
    f << "#local Nr = Nr + 0.04;\n";
    f << "#end\n";
    f << "}\n";
    f.close();
  }
#endif

  // The start of rendering rods
  int startindex = -1;

#ifdef FLAGWALKER
  startindex = 9;
#endif

  // The start of changing color of rods
  int startindexcolor = 3;
  int startindexcolor_2 = 6;

#ifdef FLAGFLAGELLA
  startindexcolor = 1;
  startindexcolor_2 = 2;
#endif
#ifdef FLAGWALKER
  startindexcolor = 0;
  startindexcolor_2 = 20;
#endif
#ifdef FLAGMUSCULARSNAKE
  startindexcolor = 1;
  startindexcolor_2 = 0;
#endif
#ifdef FLAGWING
  startindexcolor = 6;
  startindexcolor_2 = 5000;
#endif

  if (colorindex > startindex) {
    // header
    fstream f;
    f.open(fileData.c_str(), fstream::app);

    f << "sphere_sweep\n";
    f << "{\n";
    f << "b_spline\n";
    unsigned int numVertex = n + 1;
    if ((n == 19) || (n == 25) || (n == 21) || (n == 9)) {
      numVertex = 2 * n + 1;
    }
    if ((n == 4) || (n == 2) || (n == 1000) || (n == 8)) {
      numVertex = 3 * n + 1;
    }
    if (n == 14) {
      numVertex = n;
    }
    f << numVertex << ",\n";
    f << "<" << x[0][0] << "," << x[0][1] << "," << x[0][2] << ">, " << r[0]
      << ",\n";

    // The main body of the rod
    for (unsigned int i = 1; i < n; i++) {
      if ((n == 19) || (n == 25) || (n == 21)) {
        f << "<" << (x[i - 1][0] + x[i][0]) / 2 << ","
          << (x[i][1] + x[i - 1][1]) / 2 << "," << (x[i][2] + x[i - 1][2]) / 2
          << ">, " << 0.5 * (r[i - 1] + r[i]) << ",\n";
      }
      if (n == 9) {
        REAL r1 = (i == 1) ? r[0] : 0.5 * (r[i - 2] + r[i - 1]);
        REAL r2 = 0.5 * (r[i - 1] + r[i]);
        f << "<" << (x[i - 1][0] + x[i][0]) / 2 << ","
          << (x[i][1] + x[i - 1][1]) / 2 << "," << (x[i][2] + x[i - 1][2]) / 2
          << ">, " << 0.5 * (r1 + r2) << ",\n";
      }
      if ((n == 4) || (n == 2) || (n == 1000) || (n == 8)) {
        f << "<" << (2 * x[i - 1][0] + x[i][0]) / 3 << ","
          << (x[i][1] + 2 * x[i - 1][1]) / 3 << ","
          << (x[i][2] + 2 * x[i - 1][2]) / 3 << ">, "
          << (2 * r[i - 1] + r[i]) / 3 << ",\n";
        f << "<" << (x[i - 1][0] + 2 * x[i][0]) / 3 << ","
          << (2 * x[i][1] + x[i - 1][1]) / 3 << ","
          << (2 * x[i][2] + x[i - 1][2]) / 3 << ">, "
          << (r[i - 1] + 2 * r[i]) / 3 << ",\n";
      }
      if (n == 14) {
        if (i < (n - 1)) {
          f << "<" << x[i][0] << "," << x[i][1] << "," << x[i][2] << ">, "
            << 0.5 * (r[i - 1] + r[i]) << ",\n";
        } else {
          f << "<" << x[i][0] << "," << x[i][1] << "," << x[i][2] << ">, "
            << 0.5 * (r[i - 1] + r[i]) << "\n";
        }
      } else {
        f << "<" << x[i][0] << "," << x[i][1] << "," << x[i][2] << ">, "
          << 0.5 * (r[i - 1] + r[i]) << ",\n";
      }
    }

    // The last point of the rod
    if ((n == 4) || (n == 2) || (n == 1000) || (n == 8)) {
      f << "<" << (x[n][0] + 2 * x[n - 1][0]) / 3 << ","
        << (x[n][1] + 2 * x[n - 1][1]) / 3 << ","
        << (x[n][2] + 2 * x[n - 1][2]) / 3 << ">, "
        << (r[n - 2] + 2 * r[n - 1]) / 3 << ",\n";
      f << "<" << (2 * x[n][0] + x[n - 1][0]) / 3 << ","
        << (2 * x[n][1] + x[n - 1][1]) / 3 << ","
        << (2 * x[n][2] + x[n - 1][2]) / 3 << ">, "
        << (2 * r[n - 2] + r[n - 1]) / 3 << ",\n";
    }
    if ((n == 19) || (n == 25) || (n == 21)) {
      f << "<" << (x[n][0] + x[n - 1][0]) / 2 << ","
        << (x[n][1] + x[n - 1][1]) / 2 << "," << (x[n][2] + x[n - 1][2]) / 2
        << ">, " << 0.5 * (r[n - 2] + r[n - 1]) << ",\n";
    }
    if ((n == 9)) {
      f << "<" << (x[n][0] + x[n - 1][0]) / 2 << ","
        << (x[n][1] + x[n - 1][1]) / 2 << "," << (x[n][2] + x[n - 1][2]) / 2
        << ">, " << 0.5 * (0.5 * (r[n - 2] + r[n - 1]) + r[n - 1]) << ",\n";
    }
    if (n == 14) {
    } else {
      f << "<" << x[n][0] << "," << x[n][1] << "," << x[n][2] << ">, "
        << r[n - 1] << "\n";
    }

    // Tail
    f << "texture\n";
    f << "{\n";
    if (colorindex < startindexcolor) {
      f << "pigment{ color rgbt<0.45,0.39,1,0> }\n";
    } else if (colorindex < startindexcolor_2) {
#if defined FLAGELBOW || FLAGWALKER || FLAGFLAGELLA
      f << "pigment{ color rgb<1,0.65,0> }\n";
#endif

#ifdef FLAGWING
      if ((n == 4) || (n == 8)) {
        f << "pigment{ color MediumSeaGreen }\n";
      } else if (n == 7) {
        REAL color = 0.65 - 0.1 * (r[4] - 10.6066);
        f << "pigment{ color rgb<1," << color << ",0> }\n";
      } else {
        f << "pigment{ color rgb<1,0.65,0> }\n";
      }
#endif
    } else if (colorindex == 42) {
      REAL color = 0.65 - 80 * (r[n / 2] - 0.006);
      f << "pigment{ color rgb<1," << color << ",0> }\n";
    } else {
      REAL color = 0.65 - 700 * (r[n / 2] - 0.004);
      f << "pigment{ color rgb<1," << color << ",0> }\n";
    }
    f << "finish { phong 1 }\n";
    f << "}\n";

#if defined FLAGELBOW || FLAGMUSCULARSNAKE
    f << "scale<4,4,4>\n";  // Elbow,Snake
#endif
#ifdef FLAGFLAGELLA
    f << "scale<8,8,8>\n";  // Flagella
#endif
#ifdef FLAGWALKER
    f << "scale<0.25,0.25,0.25>\n";  // Walker
#endif
#if defined FLAGWING || FLAGSPHERICALJOINT || FLAGHINGEJOINT || \
    FLAGFIXEDJOINT || FLAGPULLINGMUSCLE
    f << "scale<0.07,0.07,0.07>\n";  // wing
#endif
    f << "rotate<0,90,90>\n";
    f << "translate<2,0,4>\n";
    f << "}\n";

    f.close();
  }  // if condition
}  // void Rod::dumpPovray

#ifdef SNAKE_VIZ
void Rod::paint() {
  assert(x.size() > 5);

  const float SPAN = 2.0;
  const float offsetAxis = 0.025;
  const float scalingAxis = 0.75;
  // const float scalingOrigin = 0.25;

  // Draw body top view
  {
    glPushMatrix();

    double color[4] = {1.0, 0.5, 0.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINE_STRIP);
    {
      for (unsigned int i = 0; i < x.size(); i++)
        glVertex2f((0.5 / SPAN) * x[i].x + 0.25, (0.5 / SPAN) * x[i].y + 0.75);
    }

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  // Draw head top view
  {
    glPushMatrix();

    double color[4] = {0.0, 1.0, 0.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINE_STRIP);
    {
      for (unsigned int i = 0; i < (unsigned int)min(5, (int)x.size()); i++)
        glVertex2f((0.5 / SPAN) * x[i].x + 0.25, (0.5 / SPAN) * x[i].y + 0.75);
    }
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  // Draw axis top view
  draw_arrow(offsetAxis, 0.5 + offsetAxis, "x", "y", scalingAxis);
  // draw_origin(0.25, 0.75, scalingOrigin);

  // Draw body lateral view
  {
    glPushMatrix();

    double color[4] = {1.0, 0.5, 0.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINE_STRIP);
    {
      for (unsigned int i = 0; i < x.size(); i++)
        glVertex2f((0.5 / SPAN) * x[i].x + 0.25, (0.5 / SPAN) * x[i].z + 0.25);
    }

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  // Draw head lateral view
  {
    glPushMatrix();

    double color[4] = {0.0, 1.0, 0.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINE_STRIP);
    {
      for (unsigned int i = 0; i < (unsigned int)min(5, (int)x.size()); i++)
        glVertex2f((0.5 / SPAN) * x[i].x + 0.25, (0.5 / SPAN) * x[i].z + 0.25);
    }
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  // Draw axis lateral view
  draw_arrow(offsetAxis, offsetAxis, "x", "z", scalingAxis);
  // draw_origin(0.25, 0.25, scalingOrigin);

  // Draw body front view
  {
    glPushMatrix();

    double color[4] = {1.0, 0.5, 0.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINE_STRIP);
    {
      for (unsigned int i = 0; i < x.size(); i++)
        glVertex2f((0.5 / SPAN) * x[i].y + 0.75, (0.5 / SPAN) * x[i].z + 0.25);
    }

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  // Draw head front view
  {
    glPushMatrix();

    double color[4] = {0.0, 1.0, 0.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_LINE_STRIP);
    {
      for (unsigned int i = 0; i < (unsigned int)min(5, (int)x.size()); i++)
        glVertex2f((0.5 / SPAN) * x[i].y + 0.75, (0.5 / SPAN) * x[i].z + 0.25);
    }
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  // Draw axis front view
  draw_arrow(0.5 + offsetAxis, offsetAxis, "y", "z", scalingAxis);
  // draw_origin(0.75, 0.25, scalingOrigin);

  // Draw subdivisions
  {
    glPushMatrix();

    double color[4] = {1.0, 1.0, 1.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINES);
    {
      glVertex2f(0.5, 0.0);
      glVertex2f(0.5, 1.0);

      glVertex2f(0.0, 0.5);
      glVertex2f(1.0, 0.5);
    }
    glEnd();

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  /*
          // Draw self-contact
          if (useSelfContact)
          {
                  glPushMatrix();

                  double color[4] = {1.0,0.0,0.0,1.0};
                  glColor4f(color[0],color[1],color[2],color[3]);
                  glLineWidth(2.0);
                  glEnable(GL_LINE_SMOOTH);
                  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
                  glEnable(GL_BLEND);
                  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                  glBegin(GL_LINES);
                  {
                          for (unsigned int i=0; i<x.size(); i++)
                                  if((abs(selfCollisionForces[i].x) >
     Tolerance::tol() || abs(selfCollisionForces[i].y) > Tolerance::tol() ||
     abs(selfCollisionForces[i].z) > Tolerance::tol()) &&
     (abs(selfCollisionForces[i+1].x) > Tolerance::tol() ||
     abs(selfCollisionForces[i+1].y) > Tolerance::tol() ||
     abs(selfCollisionForces[i+1].z) > Tolerance::tol()))
                                  {
                                          glVertex2f(x[i][0], x[i][view]);
                                          glVertex2f(x[i+1][0], x[i+1][view]);
                                  }

                  }
                  glEnd();
                  glDisable(GL_LINE_SMOOTH);
                  glDisable(GL_BLEND);

                  glPopMatrix();
          }
  */
}

void Rod::draw_strokestring(void *font, float const size, char const *string) {
  glPushMatrix();
  float const scale = size * 0.01; /* GLUT character base size is 100 units */
  glScalef(scale, scale, scale);

  char const *c = string;
  for (; c && *c; c++) {
    glutStrokeCharacter(font, *c);
  }
  glPopMatrix();
}

void Rod::draw_arrow(const float posx, const float posy,
                     char const *const annotation_hor,
                     char const *const annotation_vert, float scale) {
  const float tx = posx;
  const float ty = posy;
  const float tz = 0.0;
  const float length = scale * 0.1;
  const float fontSize = scale * 0.025;

  {
    // Draw horizontal axis
    glPushMatrix();

    double color[4] = {1.0, 1.0, 1.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINES);
    {
      glVertex2f(tx + 0.0, ty + 0.0);
      glVertex2f(tx + length, ty + 0.0);
    }
    glEnd();

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  {
    // Draw horizontal axis label
    glPushMatrix();
    glTranslatef(tx + 0.75 * length, ty - fontSize, tz);
    draw_strokestring(GLUT_STROKE_MONO_ROMAN, fontSize, annotation_hor);
    glPopMatrix();
  }

  {
    // Draw vertical axis
    glPushMatrix();

    double color[4] = {1.0, 1.0, 1.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINES);
    {
      glVertex2f(tx + 0.0, ty + 0.0);
      glVertex2f(tx + 0.0, ty + length);
    }
    glEnd();

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  {
    // Draw vertical axis label
    glPushMatrix();
    glTranslatef(tx - fontSize, ty + 0.75 * length, tz);
    draw_strokestring(GLUT_STROKE_MONO_ROMAN, fontSize, annotation_vert);
    glPopMatrix();
  }
}

void Rod::draw_origin(const float posx, const float posy, float scale) {
  const float tx = posx;
  const float ty = posy;
  const float length = scale * 0.1;

  {
    // Draw horizontal axis
    glPushMatrix();

    double color[4] = {1.0, 1.0, 1.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINES);
    {
      glVertex2f(tx + 0.0, ty - length / 2.0);
      glVertex2f(tx + 0.0, ty + length / 2.0);
    }
    glEnd();

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }

  {
    // Draw vertical axis
    glPushMatrix();

    double color[4] = {1.0, 1.0, 1.0, 1.0};
    glColor4f(color[0], color[1], color[2], color[3]);
    glLineWidth(2.0);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_LINES);
    {
      glVertex2f(tx - length / 2.0, ty + 0.0);
      glVertex2f(tx + length / 2.0, ty + 0.0);
    }
    glEnd();

    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glPopMatrix();
  }
}
#endif
