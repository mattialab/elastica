#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <limits>
#include <assert.h>
#include "Rod.h"
#include "MathFunctions.h"

// Class defining interaction forces between rod and rod, substrate, etc
class Interaction
{
protected:
	typedef std::vector<Vector3> VV3;
	typedef std::vector<Matrix3> VM3;
	typedef std::vector<REAL> VREAL;
	typedef std::vector<bool> VBOOL;
	typedef std::vector<int> VINT;
	typedef std::vector<Rod*> Vrodptr;
	typedef std::vector<Interaction*> Vinterptr;

	Vrodptr rodptrs;

public:
	Interaction(Vrodptr rodptrs) : rodptrs(rodptrs) {}
	virtual ~Interaction(){}

	virtual void applyForces(const REAL time) = 0;
};

/*
// A frictional plane, for the moment set at y=0, rod free for y>r.
class FrictionPlaneInteraction: public Interaction
{
protected:
	const REAL kPlane;		// Elastic response to prevent interpenetration rod-plane
	const REAL etaPlane;	// Dissipative term appliet to the elastice response rod-plane
	const REAL muKineticPlane;
	const REAL muStaticPlane;
	const REAL vStaticPlane;
	const Vector3 normalPlane;
	const Vector3 originPlane;


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
	FrictionPlaneInteraction(Vrodptr& rodptrs, const Vector3 _normalPlane, const Vector3 _originPlane, const REAL _kPlane, const REAL _etaPlane, const REAL _muKineticPlane, const REAL _muStaticPlane, const REAL _vStatic = 1e-4) :
		Interaction(rodptrs), normalPlane(_normalPlane), originPlane(_originPlane), kPlane(_kPlane), etaPlane(_etaPlane), muStaticPlane(_muStaticPlane), muKineticPlane(_muKineticPlane),
		vStaticPlane(_vStatic){}

	void applyForces(const REAL time)
	{
		Rod* rod = rodptrs[0];

		vector<Vector3> overallRodForces = vector<Vector3>(rod->totalInternalForces.size()-1);
		vFromPointsToElements(rod->totalInternalForces + rod->externalForces, overallRodForces);

		for (int i=0; i < rod->n; i++)
		{
			// Fetch initial useful quantities
			const REAL r = rod->r[i];
			const Vector3 elementX = 0.5*(rod->x[i] + rod->x[i+1]);
			const REAL distanceFromPlane = (elementX - originPlane) % normalPlane;

			// Numerical fluctuations may cause the distnace from plane to be larger than r.
			// We ameng this problem by allowing a tolerance proportinal (arbitrary constant) to the radius.
			// For more robustness check the SUGGESTION in Rod.cpp->applyFrictions
			if ( (distanceFromPlane-r) < 1e-4 )
			{
				// Fetch remaining quantites
				const Vector3 elementV = 0.5*(rod->v[i] + rod->v[i+1]);
				const Vector3 axialDir = (rod->edge[i] / rod->edge[i].length()).unitize();
				const Vector3 w = rod->w[i];
				const Matrix3 Q = rod->Q[i];

				// Compute plane tangent and bitangent (rolling) direction
				const Vector3 rollingDirection = (normalPlane * axialDir).unitize();
				const Vector3 axialDirection = (rollingDirection * normalPlane).unitize();

				const REAL currentRodForcesInNormaldirectionSign = (overallRodForces[i] % normalPlane);
				const Vector3 currentRodForcesInNormaldirection = (currentRodForcesInNormaldirectionSign<0.0)? currentRodForcesInNormaldirectionSign*normalPlane : 0.0*normalPlane;
				const Vector3 elasticPlaneResponse = - kPlane * min(distanceFromPlane-r,0.0) * normalPlane;
				const Vector3 dampingForcePlane = - etaPlane * (elementV % normalPlane)*normalPlane;
				const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
				rod->externalForces[i]   += 0.5 * (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
				rod->externalForces[i+1] += 0.5 * (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
				const REAL normalForce = normalForcePlane.length();

				// Compute slip velocity decomposed in rolling and tangential directions
				const Vector3 armTorque = -r * normalPlane;
				const Vector3 translationVelInRollingDirection = (elementV % rollingDirection)*rollingDirection;
				const Vector3 rotationalVelInRollingDirection = ((Q.T()*(w*(Q*armTorque))) % rollingDirection)*rollingDirection;
				const Vector3 slipVelInRollingDirection = translationVelInRollingDirection + rotationalVelInRollingDirection;

				// Friction in rolling direction
				{
					const REAL f = _linear(slipVelInRollingDirection.length(), vStaticPlane);

					const Vector3 kineticFrictionForce = -(1.0-f) * muKineticPlane * normalForce * slipVelInRollingDirection.unitize();
					const Vector3 staticFrictionForce  = f * muStaticPlane  * normalForce * rollingDirection.unitize();

					rod->kineticFrictionsForce[i]	 += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsForce[i+1] += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsTorque[i]  += Q*(armTorque * kineticFrictionForce);

					// Note that the static friction genereta torques as well.
					// Since they are annoying to compute and require more knowledge of the rod's configuration
					// their calculation is in the rod's class itself. But they are accounted for!!
					rod->staticFrictionsRollingForce[i]   += 0.5 * staticFrictionForce;
					rod->staticFrictionsRollingForce[i+1] += 0.5 * staticFrictionForce;
					rod->staticFrictionsNormalPlane[i]    += normalPlane;

					rod->rollingf[i] = f;
				}

				// Compute axial velocity
				const Vector3 slipVelInAxiallDirection = (elementV % axialDirection) * axialDirection;

				// Friction in axial direction
				{
					const REAL f = _linear(slipVelInAxiallDirection.length(), vStaticPlane);

					const Vector3 kineticFrictionForce = -(1.0-f) * muKineticPlane * normalForce * slipVelInAxiallDirection.unitize();
					const Vector3 staticFrictionForce  = f * muStaticPlane  * normalForce * axialDirection.unitize();

					rod->kineticFrictionsForce[i]	 += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsForce[i+1] += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsTorque[i]  += Q*(armTorque * kineticFrictionForce);

					// Note that the static friction genereta torques as well.
					// Since they are annoying to compute and require more knowledge of the rod's configuration
					// their calculation is in the rod's class itself. But they are accounted for!!
					rod->staticFrictionsAxialForceForward[i]   += 0.5 * staticFrictionForce;
					rod->staticFrictionsAxialForceForward[i+1] += 0.5 * staticFrictionForce;
				}
			}
			else
			{
				cout << "detach " << (distanceFromPlane-r) << endl;
				exit(0);
			}
		}
	}
};
*/

// Anisotrpic frictional plane, for the moment set at y=0, rod free for y>r.
class AnisotropicFrictionPlaneInteraction: public Interaction
{
protected:
	const REAL kPlane;		// Elastic response to prevent interpenetration rod-plane
	const REAL etaPlane;	// Dissipative term appliet to the elastice response rod-plane
	const REAL muKineticForward;
	const REAL muKineticBackward;
	const REAL muKineticSideways;
	const REAL muStaticForward;
	const REAL muStaticBackward;
	const REAL muStaticSideways;
	const REAL vStaticPlane;
	const Vector3 normalPlane;
	const Vector3 originPlane;


	const REAL _linear(const REAL _vel, const REAL _velThresold) const
	{
		const REAL vel = fabs(_vel);
		const REAL velThresold = fabs(_velThresold);

		const REAL width = 0.5*velThresold;
		const REAL velDiff = vel - velThresold;
		REAL f = 1.0;

		if ( vel>(velThresold) )
			f = fabs(1.0-min(1.0,velDiff/width));

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
	AnisotropicFrictionPlaneInteraction(Vrodptr& rodptrs, const Vector3 _normalPlane, const Vector3 _originPlane, const REAL _kPlane, const REAL _etaPlane,
										const REAL _muKineticForward, const REAL _muKineticBackward, const REAL _muKineticSideways,
										const REAL _muStaticForward, const REAL _muStaticBackward, const REAL _muStaticSideways, const REAL _vStatic = 1e-4) :
		Interaction(rodptrs), normalPlane(_normalPlane), originPlane(_originPlane), kPlane(_kPlane), etaPlane(_etaPlane),
		muKineticForward(_muKineticForward), muKineticBackward(_muKineticBackward), muKineticSideways(_muKineticSideways),
		muStaticForward(_muStaticForward), muStaticBackward(_muStaticBackward), muStaticSideways(_muStaticSideways), vStaticPlane(_vStatic){}

	void applyForces(const REAL time)
	{
		Rod* rod = rodptrs[0];

		vector<Vector3> overallRodForces = vector<Vector3>(rod->totalInternalForces.size()-1);
		vFromPointsToElements(rod->totalInternalForces + rod->externalForces, overallRodForces);

		for (unsigned int i=0; i < (unsigned int)rod->n; i++)
		{
			// Fetch initial useful quantities
			const REAL r = rod->r[i];
			const Vector3 elementX = 0.5*(rod->x[i] + rod->x[i+1]);
			const REAL distanceFromPlane = (elementX - originPlane) % normalPlane;

			// Numerical fluctuations may cause the distnace from plane to be larger than r.
			// We ameng this problem by allowing a tolerance proportinal (arbitrary constant) to the radius.
			// For more robustness check the SUGGESTION in Rod.cpp->applyFrictions
			if ( (distanceFromPlane-r) < 1e-4 )
			{
				// Fetch remaining quantites
				const Vector3 elementV = 0.5*(rod->v[i] + rod->v[i+1]);
				const Vector3 axialDir = (rod->edge[i] / rod->edge[i].length()).unitize();
				const Vector3 w = rod->w[i];
				const Matrix3 Q = rod->Q[i];

				// Compute plane tangent and bitangent (rolling) direction
				const Vector3 rollingDirection = (normalPlane * axialDir).unitize();
				const Vector3 axialDirection = (rollingDirection * normalPlane).unitize();

				const REAL currentRodForcesInNormaldirectionSign = (overallRodForces[i] % normalPlane);
				const Vector3 currentRodForcesInNormaldirection = (currentRodForcesInNormaldirectionSign<0.0)? currentRodForcesInNormaldirectionSign*normalPlane : 0.0*normalPlane;
				const Vector3 elasticPlaneResponse = - kPlane * min(distanceFromPlane-r,0.0) * normalPlane;
				const Vector3 dampingForcePlane = - etaPlane * (elementV % normalPlane)*normalPlane;
				const Vector3 normalForcePlane = -currentRodForcesInNormaldirection;
				rod->externalForces[i]   += 0.5 * (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
				rod->externalForces[i+1] += 0.5 * (normalForcePlane + elasticPlaneResponse + dampingForcePlane);
				const REAL normalForce = normalForcePlane.length();

				// Compute slip velocity decomposed in rolling and tangential directions
				const Vector3 armTorque = -r * normalPlane;
				const Vector3 translationVelInRollingDirection = (elementV % rollingDirection)*rollingDirection;
				const Vector3 rotationalVelInRollingDirection = ((Q.T()*(w*(Q*armTorque))) % rollingDirection)*rollingDirection;
				const Vector3 slipVelInRollingDirection = translationVelInRollingDirection + rotationalVelInRollingDirection;

				// Friction in rolling direction
				{
					const REAL f = _linear(slipVelInRollingDirection.length(), vStaticPlane);

					const Vector3 kineticFrictionForce = -(1.0-f) * muKineticSideways * normalForce * slipVelInRollingDirection.unitize();
					const Vector3 staticFrictionForce  = f * muStaticSideways  * normalForce * rollingDirection.unitize();

					rod->kineticFrictionsForce[i]	 += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsForce[i+1] += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsTorque[i]  += Q*(armTorque * kineticFrictionForce);

					// Note that the static friction genereta torques as well.
					// Since they are annoying to compute and require more knowledge of the rod's configuration
					// their calculation is in the rod's class itself. But they are accounted for!!
					rod->staticFrictionsRollingForce[i]   += 0.5 * staticFrictionForce;
					rod->staticFrictionsRollingForce[i+1] += 0.5 * staticFrictionForce;
					rod->staticFrictionsNormalPlane[i]    += normalPlane;
				}

				// Compute axial velocity
				const REAL slipVelInAxiallDirectionSign = elementV % axialDirection;
				const Vector3 slipVelInAxiallDirection = slipVelInAxiallDirectionSign * axialDirection;

				// Friction in axial direction
				{
					const REAL f = _linear(slipVelInAxiallDirection.length(), vStaticPlane);

					// Here we assume that the head of a snake is locate at the first point x[0]!!
					// The axial direction is computed x[i+1]-x[i], that is why there is a (slipVelInAxiallDirectionSign<0.0), and NOT (slipVelInAxiallDirectionSign>0.0)
					const REAL kineticFrictionCoeff = (slipVelInAxiallDirectionSign<0.0)?muKineticForward:muKineticBackward;
					const Vector3 kineticFrictionForce = -(1.0-f) * kineticFrictionCoeff * normalForce * slipVelInAxiallDirection.unitize();

					// Again here we assume that the head of a snake is locate at the first point x[0]!!
					// That is why there is a minus, it is supposed to flip the vector "axialDirection" so that it points in the direction tail-to-head and not head-to-tail
					const Vector3 staticFrictionForceForward  = - f * muStaticForward  * normalForce * axialDirection.unitize();
					const Vector3 staticFrictionForceBackward  = - f * muStaticBackward  * normalForce * axialDirection.unitize();

					rod->kineticFrictionsForce[i]	 += 0.5 * kineticFrictionForce;
					rod->kineticFrictionsForce[i+1] += 0.5 * kineticFrictionForce;
					// No torque for axial case
					//rod->kineticFrictionsTorque[i]  += Q*(armTorque * kineticFrictionForce);

					// Note that the static friction genereta torques as well.
					// Since they are annoying to compute and require more knowledge of the rod's configuration
					// their calculation is in the rod's class itself. But they are accounted for!!
					rod->staticFrictionsAxialForceForward[i]   += 0.5 * staticFrictionForceForward;
					rod->staticFrictionsAxialForceForward[i+1] += 0.5 * staticFrictionForceForward;
					rod->staticFrictionsAxialForceBackward[i]   += 0.5 * staticFrictionForceBackward;
					rod->staticFrictionsAxialForceBackward[i+1] += 0.5 * staticFrictionForceBackward;
				}
			}
			//else
			//{
			//	cout << "detach " << (distanceFromPlane-r) << endl;
				//exit(0);
			//}
		}
	}
};


// Anisotrpic frictional plane, for the moment set at y=0, rod free for y>r.
class SlenderBodyTheoryEnvironment: public Interaction
{
protected:
	const REAL dynamicViscosity;

public:
		SlenderBodyTheoryEnvironment(Vrodptr& rodptrs, const REAL _dynamicViscosity) : Interaction(rodptrs), dynamicViscosity(_dynamicViscosity) {}

	void applyForces(const REAL time)
	{
		Rod* rod = rodptrs[0];

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

		const REAL factorSBT = -4.0*M_PI*dynamicViscosity/( log(totalLength/sumRadii) );

		for (unsigned int i=0; i < (unsigned int)rod->n; i++)
		{
			const Vector3 v = 0.5*(rod->v[i] + rod->v[i+1]);
			const Vector3 t = (rod->edge[i] / rod->edge[i].length()).unitize();
			Vector3 stokesForce = Vector3();
			Vector3::projectionSBT(v,t,stokesForce);
			stokesForce *= factorSBT * rod->l[i];
			rod->externalForces[i]   += 0.5 * stokesForce;
			rod->externalForces[i+1] += 0.5 * stokesForce;
		}
	}
};



// Interaction between a rod (snake) and substrate
class SnakePlane: public Interaction
{
public:
	REAL k, nu, vStatic, muForwardKinetic, muSidewaysKinetic, muBackwardKinetic, muForwardStatic, muSidewaysStatic, muBackwardStatic;

	SnakePlane(	Vrodptr rodptrs, REAL k, REAL nu, REAL vStatic, REAL muForwardKinetic, REAL muSidewaysKinetic,
			REAL muBackwardKinetic, REAL muForwardStatic, REAL muSidewaysStatic, REAL muBackwardStatic):
				Interaction(rodptrs), k(k), nu(nu), vStatic(vStatic), muForwardKinetic(muForwardKinetic),
				muSidewaysKinetic(muSidewaysKinetic), muBackwardKinetic(muBackwardKinetic),
				muForwardStatic(muForwardStatic), muSidewaysStatic(muSidewaysStatic), muBackwardStatic(muBackwardStatic) {}

	void applyForces(const REAL time)
	{
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

				Vector3 mvp = mv - Vector3(0, 1, 0) * (mv % Vector3(0,1,0)); // velocity in plane;
				REAL vTot = mvp.length();

				if (vTot > vStatic)
				{
					// kinetic friction first
					rp->collided[i] = 1;
					Vector3 fForce = Vector3();

					REAL vAxial =  mvp % t;

					if (vAxial > 0) // front/back friction
						fForce -= t * vAxial / vTot * nForce * muForwardKinetic;
					else
						fForce -= t * vAxial / vTot * nForce * muBackwardKinetic;

					fForce -= (mvp - t * (mvp % t)) / vTot * nForce * muSidewaysKinetic; // sideways friction

					rp->externalForces[i]     += fForce/2;
					rp->externalForces[i + 1] += fForce/2;
					rp->externalTorques[i]    -= Q * (Vector3(0, -r, 0) * fForce);
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

