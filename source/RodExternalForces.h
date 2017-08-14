#ifndef RODEXTERNALFORCES_H_
#define RODEXTERNALFORCES_H_

#include <limits>
#include "Rod.h"
#include "MathFunctions.h"
#include "Interaction.h"
#include "RodBoundaryConditions.h"
#include "SplineProfileZeroEnds.h"

using namespace std;

// Abstract ExternalForce class representing some external force acting upon a single Rod
class ExternalForces
{
protected:
	typedef std::vector<Vector3> VV3;
	typedef std::vector<Matrix3> VM3;
	typedef std::vector<REAL> VREAL;
	typedef std::vector<bool> VBOOL;
	typedef std::vector<int> VINT;
	typedef std::vector<Rod*> Vrodptr;
	typedef std::vector<ExternalForces*> Vefptr;
	typedef std::vector<RodBC*> Vbcptr;

public:
	virtual ~ExternalForces(){}

	virtual void applyForces(Rod &R, const REAL time) = 0;
};

class MultipleForces : public ExternalForces
{
public:
	Vefptr forces;

	MultipleForces() : ExternalForces() {}
	MultipleForces(Vefptr forces) : ExternalForces(), forces(forces) {}

	void add(ExternalForces* singleForce) { forces.push_back(singleForce); }

	Vefptr get()
	{
		Vefptr container;
		container.push_back(this);
		return container;
	}

	void applyForces(Rod &R, const REAL time)
	{
		for(unsigned int i=0; i<forces.size(); i++)
			forces[i]->applyForces(R, time);
	}
};

class NoForces : public ExternalForces
{
public:
  void applyForces(Rod &R, const REAL time){}
};

// Periodic load is applied at the second hand
class PeriodicLoad : public ExternalForces
{
public:
	const Vector3 F;
	const REAL omega;

	PeriodicLoad(const Vector3 _F, const REAL _omega) : ExternalForces(), F(_F), omega(_omega) {}

	void applyForces(Rod &R, const REAL time)
	{
		R.externalForces.back() += F*sin(omega*time);
	}
};


// Periodic couple is applied at the second hand
class PeriodicCouple : public ExternalForces
{
public:
	typedef enum {D1, D2, D3} Direction;

	const Direction dir;
	Vector3 maskDirection;
	const REAL M;
	const REAL omega;

	PeriodicCouple(const REAL _M, const REAL _omega, const Direction _dir) : ExternalForces(), M(_M), omega(_omega), dir(_dir)
	{
		// Active torque expressed directly in the material frame of reference, hence the naming D1, D2, D3!
		switch (dir)
		{
			case Direction::D1:
				maskDirection = Vector3(1.0, 0.0, 0.0);
				break;
			case Direction::D2:
				maskDirection = Vector3(0.0, 1.0, 0.0);
				break;
			case Direction::D3:
				maskDirection = Vector3(0.0, 0.0, 1.0);
				break;
			default:
				assert(false);
				break;
		}
	}

	void applyForces(Rod &R, const REAL time)
	{
		R.externalTorques.back() += M*maskDirection*sin(omega*time);
	}
};

// Exerts a force proportional to vertex mass on each vertex in a prescribed direction
class GravityForce : public ExternalForces
{
public:
	const Vector3 F;

	GravityForce(const Vector3 F) : ExternalForces(), F(F) {}

	void applyForces(Rod &R, const REAL time)
	{
		R.externalForces += R.m * F;
	}
};


class ConstantMuscleTorques : public ExternalForces
{
public:
	const REAL A, rumpUpTime;
	const Vector3 director;

	ConstantMuscleTorques(const REAL _A, const Vector3 _director, const REAL _rumpUpTime) :
		ExternalForces(), A(_A), director(_director), rumpUpTime(_rumpUpTime)
	{
	}

	void applyForces(Rod &R, const REAL time)
	{
		const VREAL s = vCumSum(R.l0);

		assert(s.size()==(unsigned int)R.n);
		assert(s.size()==(unsigned int)R.externalTorques.size());

		for(unsigned int i = 1; i < (unsigned int)R.n; i++)
		{
			const REAL factorTime = (time<=rumpUpTime)?(time/rumpUpTime):1.0;
			const REAL torque = factorTime*A;

			R.externalTorques[i - 1] -= torque*director;
			R.externalTorques[i]     += torque*director;
		}
	}
};


// Muscles actually exert torques, not forces, we just name it this way following convention
// T = A sin(wx - vt), traveling wave with Vector3 A, magnitude |A|
// A muscle acts at a vertex, imposing equal and opposite torques on each adjacent edge
class MuscleForce : public ExternalForces
{
public:
	const Vector3 A;
	const REAL w, v;

	MuscleForce(const Vector3 A, const REAL w, const REAL v) : ExternalForces(), A(A), w(w), v(v) {}

	void applyForces(Rod &R, const REAL time)
	{
		const VREAL s = vCumSum(R.l0);
		for(unsigned int i = 1; i < (unsigned int)R.n; i++)
		{
			const Vector3 torque = A * sin( w*time - v*s[i] );
			R.externalTorques[i - 1] -= torque;
			R.externalTorques[i]     += torque;
		}
	}
};

// Spline muscles actually exert torques, not forces, we just name it this way following convention
// T = A * sin(wx - vt), traveling wave with Vector3 A, magnitude |A|. A is a spline profile, so A(s).
// A muscle acts at a vertex, imposing equal and opposite torques on each adjacent edge
class SplineMuscleTorques : public ExternalForces
{
public:
	SplineProfileZeroEnds spline;
	const REAL w, v, rumpUpTime;
	const Vector3 direction;

	SplineMuscleTorques(vector<double> A, const REAL w, const REAL v, const Vector3 _direction, const REAL _rumpUpTime) : ExternalForces(), spline(A), w(w), v(v), direction(_direction), rumpUpTime(_rumpUpTime)
	{
	}

	void applyForces(Rod &R, const REAL time)
	{
		const VREAL s = vCumSum(R.l0);
		for(unsigned int i = 1; i < (unsigned int)R.n; i++)
		{
			const REAL factor = (time<=rumpUpTime)?(time/rumpUpTime):1.0;
			const REAL torque = factor * spline.getProfile(s[i]) * sin( w*time - v*s[i] );
			//R.externalTorques[i] += R.Q[i] * direction * torque * R.l0[i];
			R.externalTorques[i - 1] -= R.Q[i-1] * direction * torque;
			R.externalTorques[i]     += R.Q[i]   * direction * torque;
		}
	}
};

// Instead of a sine wave, this muscular contraction is highly localized, namely a gaussian
// with width w. It travels with speed v down the snake. The torque is determined by the variable direction
// that is expressed in the laboratory frame of reference.
class LocalizedMuscleForce : public ExternalForces
{
public:
	const REAL A, sigma, v, position, rumpUpTime;
	const Vector3 direction;

	LocalizedMuscleForce(const REAL _A, const REAL _sigma, const REAL _v, const REAL _position, const Vector3 _direction, const REAL _rumpUpTime) : ExternalForces(), A(_A), sigma(_sigma), v(_v), position(_position), direction(_direction), rumpUpTime(_rumpUpTime)
	{
	}

	void applyForces(Rod &R, const REAL time)
	{
		const VREAL s = vCumSum(R.l0);

		for(unsigned int i = 1; i < (unsigned int)R.n; i++)
		{
			const REAL factorTime = (time<=rumpUpTime)?(time/rumpUpTime):1.0;
			const REAL torque = factorTime * A * exp( -pow(s[i] - position - v*time,2) / (2.0 * sigma * sigma));

			R.externalTorques[i - 1] -= torque*(R.Q[i-1] * direction);
			R.externalTorques[i]     += torque*(R.Q[i]   * direction);
		}
	}
};

// Instead of a sine wave, this muscular contraction is highly localized, namely a gaussian
// with width w. It travels with speed v down the snake. The torque here is directly associated to the directors
// of the material frame
class LocalizedMuscleForceLagrangian_ForcedToBeInPlane : public ExternalForces
{
public:
	const REAL A, sigma, v, position, rumpUpTime;
	const Vector3 director, normalToDesiredPlane;

	LocalizedMuscleForceLagrangian_ForcedToBeInPlane(const REAL _A, const REAL _sigma, const REAL _v, const REAL _position, const Vector3 _director, const Vector3 _normalToDesiredPlane, const REAL _rumpUpTime) :
		ExternalForces(), A(_A), sigma(_sigma), v(_v), position(_position), director(_director), normalToDesiredPlane(_normalToDesiredPlane), rumpUpTime(_rumpUpTime)
	{
	}

	void applyForces(Rod &R, const REAL time)
	{
		const VREAL s = vCumSum(R.l0);
		//const REAL totalLength = s.back();

		//const REAL widthSigmoid = 20.0;
		//const REAL xSigmoid = 2.0/3.0;

		for(unsigned int i = 1; i < (unsigned int)R.n; i++)
		{
			//const REAL factorSigmoid = 1.0 - 1.0/(1.0+exp(-widthSigmoid*(s[i]/totalLength-xSigmoid))); // normalized between 0 and 1
			const REAL factorTime = (time<=rumpUpTime)?(time/rumpUpTime):1.0;
			const REAL torque = factorTime * A * exp( -pow(s[i] - position - v*time,2) / (2.0 * sigma * sigma));

			R.externalTorques[i - 1] -= torque*director;
			R.externalTorques[i]     += torque*director;

			/*
			const Vector3 normalPlaneInMaterialFrame_minus = R.Q[i-1]*normalToDesiredPlane;
			const Vector3 normalPlaneInMaterialFrame = R.Q[i]*normalToDesiredPlane;

			const Vector3 directorInPlane_minus = (director - (director%normalPlaneInMaterialFrame_minus)*normalPlaneInMaterialFrame_minus );
			const Vector3 directorInPlane = (director - (director%normalPlaneInMaterialFrame)*normalPlaneInMaterialFrame );

			const Vector3 nomralizedDirectorInPlane_minus = directorInPlane_minus.unitize();
			const Vector3 nomralizedDirectorInPlane = directorInPlane.unitize();

			R.externalTorques[i - 1] += torque*nomralizedDirectorInPlane_minus;
			R.externalTorques[i]     -= torque*nomralizedDirectorInPlane;
			*/
		}
	}
};


// Represents constant forces at the endpoints. Note that a fixed endpoint
// will be represented by the bc function, so its tension force should
// be zero -- only free ends make sense with tension forces
class HelicalBucklingForcesAndTorques : public ExternalForces
{
public:
	Vector3 startForce, endForce;
	Vector3 startTorque, endTorque;

	HelicalBucklingForcesAndTorques(Vector3 _startForce, Vector3 _endForce, Vector3 _startTorque, Vector3 _endTorque) :
		ExternalForces(), startForce(_startForce), endForce(_endForce), startTorque(_startTorque), endTorque(_endTorque) {}

	void applyForces(Rod& R, const REAL time)
	{
		R.externalForces.front() = startForce;
		R.externalForces.back() = endForce;
		R.externalTorques.front() = startTorque;
		R.externalTorques.back() = endTorque;
	}
};

// Periodic loads
class PeriodicLoadLongitudinalWaves : public ExternalForces
{
public:
	const REAL frequency;
	const Vector3 load;

	PeriodicLoadLongitudinalWaves(const REAL _frequency, const Vector3 _load) : ExternalForces(), frequency(_frequency), load(_load) {}

	void applyForces(Rod& R, const REAL time)
	{
		R.externalForces.back() = load * sin(frequency*time);
	}
};

// Represents constant forces at the endpoints. Note that a fixed endpoint
// will be represented by the bc function, so its tension force should
// be zero -- only free ends make sense with tension forces
class EndpointForces : public ExternalForces
{
public:
	Vector3 startForce, endForce;

	EndpointForces(Vector3 startForce, Vector3 endForce) : ExternalForces(), startForce(startForce), endForce(endForce) {}

	void applyForces(Rod& R, const REAL time)
	{
		R.externalForces[0] += startForce;
		R.externalForces[R.n] += endForce;
	}
};

// Represents constant forces at the endpoints. Note that a fixed endpoint
// will be represented by the bc function, so its tension force should
// be zero -- only free ends make sense with tension forces
class GradualEndpointForces : public ExternalForces
{
public:
	const REAL rampupTime;
	Vector3 startForce, endForce;

	GradualEndpointForces(const Vector3 _startForce, const Vector3 _endForce, const REAL _rampupTime) : ExternalForces(), startForce(_startForce), endForce(_endForce), rampupTime(_rampupTime) {}

	void applyForces(Rod& R, const REAL time)
	{
		assert(rampupTime>=0.0);
		const REAL fraction = (rampupTime==0.0)?1.0:min(1.0,time/rampupTime);
		R.externalForces[0] += fraction*startForce;
		R.externalForces[R.n] += fraction*endForce;
	}
};

class GradualEndpointTorques : public ExternalForces
{
public:
	const REAL rampupTime;
	Vector3 startTorque, endTorque;

	GradualEndpointTorques(const Vector3 _startTorque, const Vector3 _endTorque, const REAL _rampupTime) : ExternalForces(), startTorque(_startTorque), endTorque(_endTorque), rampupTime(_rampupTime) {}

	void applyForces(Rod& R, const REAL time)
	{
		assert(rampupTime>=0.0);
		const REAL fraction = (rampupTime==0.0)?1.0:min(1.0,time/rampupTime);
		R.externalTorques.front() += fraction*R.Q.front()*startTorque;
		R.externalTorques.back()  += fraction*R.Q.back()*endTorque;
	}
};

// Apply pointwise force on one end for a given amount of time and then release
class OneEndTemporaryForce : public ExternalForces
{
public:
	const REAL timeForceIsApplied;
	const Vector3 force;


	OneEndTemporaryForce(const Vector3 _force, const REAL _timeForceIsApplied) : ExternalForces(), force(_force), timeForceIsApplied(_timeForceIsApplied) {}

	void applyForces(Rod& R, const REAL time)
	{
		if (time<timeForceIsApplied)
			R.externalForces.back() += force;
	}
};

// Represents constant forces at the endpoints. Note that a fixed endpoint
// will be represented by the bc function, so its tension force should
// be zero -- only free ends make sense with tension forces
class EndpointTorques : public ExternalForces
{
public:
	Vector3 startTorque, endTorque;

	EndpointTorques(Vector3 _startTorque, Vector3 _endTorque) : ExternalForces(), startTorque(_startTorque), endTorque(_endTorque) {}

	void applyForces(Rod& R, const REAL time)
	{
		R.externalTorques.front() += (R.Q.front() * startTorque);
		R.externalTorques.back()  += (R.Q.back() * endTorque);
	}
};

class UniformTorques : public ExternalForces
{
public:
	const Vector3 torque;

	UniformTorques(const Vector3 _torque) : ExternalForces(), torque(_torque) {}

	void applyForces(Rod& R, const REAL time)
	{
		R.externalTorques += R.Q*torque;
	}
};

class UniformForces : public ExternalForces
{
public:
	const Vector3 force;

	UniformForces(const Vector3 _force) : ExternalForces(), force(_force) {}

	void applyForces(Rod& R, const REAL time)
	{
		R.externalForces += force;

		// This is bacause, the first and last masses are half the others!
		R.externalForces.front() -= force/2.0;
		R.externalForces.back() -= force/2.0;
	}
};

// A random force with strength A and width sigma
class RandomForces : public ExternalForces
{
public:
	REAL A;

	RandomForces(REAL A) : ExternalForces(), A(A) {}

	void applyForces(Rod& R, const REAL time)
	{
		// A is a force density per unit length
		//const int nSegments = R.n;
		//const REAL L = vSum(R.l);
		//const REAL scaledA = A*L/nSegments;

		const REAL scaledA = A;

		R.externalForces += scaledA*vRandVector3(R.n+1);
	}
};

#endif
