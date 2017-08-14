#ifndef ROD_H_
#define ROD_H_

#ifdef SNAKE_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

#include "UsualHeaders.h"
#include "Vector3.h"
#include "Matrix3.h"
#include "VectorFunctions.h"
#include "SpeedFunctions.h"
#include "GeometryFunctions.h"

class Rod
{
protected:
	static const unsigned int PRINTPRECISION = 6;

public:

	// Flags
	const bool useSelfContact;

	// Some utility variable (just for prototyping purposes)
	long unsigned int COUNTER;
	REAL time;
	REAL dt;

	// Number of edges
	const int n; // number of edges, so there are n+1 vertices

	// Material properties
	const REAL density; // rod density --> [kg/m3]

	// Self-contact
	const REAL zetaSoft; // self-contact stiffness, soft core
	const REAL zetaHard; // self-contact stiffness, hard core

	// Masses
	vector<REAL> m; // --> [kg], (n+1)

	// Volumes
	const vector<REAL> V; // edge volume, so for round rod V = PI*r*r*l --> [m3], (n)

	// Position, velocity, acceleration (in lab coordinates)
	vector<Vector3> x;		// centerline position in lab coords --> [m], (n+1)
	vector<Vector3> v;		// centerline velocity in lab coords --> [m/s], (n+1)
	vector<Vector3> a;		// centerline acceleration in lab coords -->[m/s2] (n+1)

	// Frame, angular velocity, angular acceleration (in material frame coordinates)
	vector<Matrix3> Q;		// material frame attached to the edge mapping local coords to lab coords: x_material = Q * x_lab --> [-], (n)
	vector<Vector3> w;		// frame angular velocity in material frame --> [1/s],  (n)
	vector<Vector3> wDot;	// frame angular acceleration in material frame --> [1/s2], (n)

	// Inertia, bending, shear matrices (in reference configuration, they embedd rotatory inertia and constitutive laws)
	const vector<Matrix3> J0;		// mass moment of inertia matrix in reference configuration --> [kg*m2], (n)
	vector<Matrix3> J0inv;	// inverse moment of inertia matrix in reference configuration --> [1/m4], (n)
	const vector<Matrix3> B0;		// bending matrix in reference configuration --> [Nm2], (n-1)
	const vector<Matrix3> S0;		// shear matrix in reference configuration --> [N], (n)

	// Radius, length, voronaoi domain in reference configuration
	const vector<REAL> l0; // lengths in reference configuration--> [m], (n)
	vector<REAL> d0; // voronoi domain in reference configuration --> [m], (n-1)

	// Current radius, length, voronaoi domain
	vector<REAL> r;	// current radii --> [m], (n)
	vector<REAL> l;	// current lengths --> [m], (n)
	vector<REAL> d;	// current voronoi domain --> [m], (n-1)

	// Dilatations
	vector<REAL> e; // dilatation -> l/l0 --> [-], (n)
	vector<REAL> e_old; // dilatation -> l/l0 --> [-], (n)
	vector<REAL> ed; // voronoi dilatation -> d/d0 --> [-], (n-1)

	// Rate of dilatation
	vector<REAL> deldt; // dilatation -> d/dt{l/l0} --> [-], (n)

	// Edges
	vector<Vector3> edge; // edges, ie e_i = x_(i+1) - x_i, expressed in lab frame --> [m], (n-1)
	vector<Vector3> edge_old; // edges, ie e_i = x_(i+1) - x_i, expressed in lab frame --> [m], (n-1)

	// Bending
	vector<Vector3> k0;						// curvature vector in the reference configuration, derivative taken dividing by [d0]!!! Defined at interior vertices --> [1/m], (n-1)
	const vector<Vector3> intrinsic_k0;			// intrinsic curvature vector --> [1/m], (n-1)
	vector<Vector3> kDiff0; 				// difference between k and intrinsic_k0 (k0-intrinsic_k0) --> [1/m], (n-1)
	vector<Vector3> bendingInternalTorques0;// bending internal torques in reference configuration, bendingInternalTorques0=B0*kDiff0 --> [Nm], (n-1)
	vector<Vector3> bendingTorques; 		// resulting bending torques at edges in the current cofiguration (obtain by rescaling appropriately bendingInternalTorques0!) --> [Nm], (n)

	// Shear
	vector<Vector3> shearStrain0;			// shear strain shearStrain0 = Q*d/dS{r} - d3 (expressed in the body frame and with respect to reference configuration) --> [-], (n)
	const vector<Vector3> intrinsicShearStrain0;	// intrinsic shear vector with respect to reference configuration--> [-], (n)
	vector<Vector3> shearStrainDiff0;		// difference between shearStrain0 and intrinsicShearStrain0 --> [-], (n)
	vector<Vector3> shearInternalForces0;	// internat shear forces with respect to the reference configuration S0*shearStrainDiff0 --> [N], (n)
	vector<Vector3> shearForces;			// resulting shear forces with respect to actual configuration at the vertices (obtained by appropriately rescaling)! --> [N], (n+1)
	vector<Vector3> shearTorques;			// resulting shear torques with respect to actual configuration at edges (obtained by appropriately rescaling)! --> [Nm], (n)

	// Damping
	const REAL nu; // linear drag term, ie fdrag = - nu * v
	const REAL relaxationNu;
	vector<Vector3> dampingForces; // --> [N], (n+1)
	vector<Vector3> dampingTorques; // --> [Nm], (n)

	// Total forces and torques
	vector<Vector3> totalForces; // --> [N], (n+1)
	vector<Vector3> totalTorques; // --> [Nm], (n)

	// Total internal forces and torques
	vector<Vector3> totalInternalForces; // --> [N], (n+1)
	vector<Vector3> totalInternalTorques; // --> [Nm], (n)

	// External forces and torques
	vector<Vector3> externalForces; // --> [N], (n+1)
	vector<Vector3> externalTorques; // --> [Nm], (n)

	// Self collision forces
	vector<Vector3> selfCollisionForces; // --> [N], (n+1)

	// Term relative to the transport of angular velocity in the angular momentum balance due to the use of material frame of reference (J0*w/e) x w
	vector<Vector3> transportW;

	// Energies
	REAL translationalEnergy;	// we mean the energy associated with velocity of centerline only here --> [J]
	REAL rotationalEnergy;		// --> [J]
	REAL bendingEnergy;			// --> [J]
	REAL shearEnergy;			// --> [J]
	REAL totalInternalEnergy;	// --> [J]


	const REAL selfKineticMu;
	const REAL selfStaticMu;
	const REAL selfKineticThreshold;
	vector<int> collided; //whether the segment was in contact with substrate the last time step

	// Friction with substrate
	vector<Vector3> staticFrictionsAxialForceForward;
	vector<Vector3> staticFrictionsAxialForceBackward;
	vector<Vector3> staticFrictionsRollingForce;
	vector<Vector3> staticFrictionsNormalPlane;
	vector<Vector3> kineticFrictionsForce;
	vector<Vector3> kineticFrictionsTorque;
	vector<REAL> rollingf;

	// Temporary vectors
	vector<Vector3> tempVV3_n;			// temporary vector of vector3 for speed (size equal to number of edges)
	vector<Vector3> tempVV3_nplus;		// temporary vector of vector3 for speed (size equal to number of nodes)
	vector<Vector3> tempVV3_nminus;		// temporary vector of vector3 for speed (size equal to number of edges minus one)
	vector<REAL> tempVREAL_nminus;		// temporary vector of reals for speed (size equal to number of edges minus one)
	vector<REAL> tempVREAL_n;			// temporary vector of reals for speed (size equal to number of edges)
	vector<REAL> tempVREAL_nplus;		// temporary vector of reals for speed (size equal to number of nodes)
	vector<Matrix3> tempVM3_n;			// temporary vector of vector3 for speed (size equal to number of edges)
	vector<Matrix3> tempVM3_nminus;		// temporary vector of vector3 for speed (size equal to number of edges minus one)
  
	Rod(const int _n, vector<Vector3> _x, vector<Vector3> _v, vector<Matrix3> _Q, vector<Vector3> _w, vector<REAL> _l0, vector<Vector3> _intrinsic_k0, vector<Vector3> _intrinsicShearStrain0,
		const vector<REAL> _m, vector<REAL> _V, const REAL _density, vector<Matrix3> _J0, vector<Matrix3> _B, vector<Matrix3> _S, const REAL _nu, const REAL _relaxationNu, const bool _useSelfContact = false);

	void computeShearStrain0();
	void computeCurvature0();

	void computeInternalShearForces0();
	void computeInternalBendingTorques0();

	void computeInternalForces();
	void computeInternalTorques();

	void computeSelfCollisionForces();
	void computeDampingForces();

	void computeAllInternalResultingForcesAndTorques();

	Vector3 computeVelocityCenterOfMass();
	Vector3 computeAngularVelocityCenterOfMass();

	void applyStaticFrictions();

	void update();

	void computeBendingEnergy();
	void computeShearEnergy();
	void computeTranslationalEnergy();
	void computeRotationalEnergy();
	void computeEnergies();

	void setDt(const REAL _dt){ dt = _dt; }
    void setTime(const REAL _time){ time = _time; }
    bool nanCheck();
    void reset();
    void dumpPovray(string filePovray, string fileData);

#ifdef SNAKE_VIZ
    void paint();
    void draw_strokestring(void *font, float const size, char const *string);
    void draw_arrow(const float posx, const float posy, char const * const annotation_vert, char const * const annotation_hor, float annot_size);
    void draw_origin(const float posx, const float posy, float scale);
#endif

    friend ostream& operator<<(ostream&, Rod&);
};

#endif
