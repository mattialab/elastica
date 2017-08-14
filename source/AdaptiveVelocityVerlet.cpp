/*
 * AdaptiveVelocityVerlet.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#include "AdaptiveVelocityVerlet.h"

REAL AdaptiveVelocityVerlet::_smoothMod(Vector3 & vec)
{
	return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z + 1e-4);
}


REAL AdaptiveVelocityVerlet::_computeMinLengthOverAcc()
{
	REAL minDt = numeric_limits<REAL>::max();

	for(unsigned int i=0; i<rodptrs.size(); i++)
		for(unsigned int j=0; j<rodptrs[i]->l.size(); j++)
			minDt = min( minDt, rodptrs[i]->l[j]/_smoothMod(rodptrs[i]->a[j]) );

	return minDt;
}

REAL AdaptiveVelocityVerlet::_computeMinOneOverAngAcc()
{
	const REAL oneDegree = (M_PI/180.0);
	REAL minDt = numeric_limits<REAL>::max();

	for(unsigned int i=0; i<rodptrs.size(); i++)
		for(unsigned int j=0; j<rodptrs[i]->wDot.size(); j++)
			minDt = min( minDt, oneDegree/_smoothMod(rodptrs[i]->wDot[j]) );

	return minDt;
}

REAL AdaptiveVelocityVerlet::_computeMinLengthOverVel()
{
	REAL minDt = numeric_limits<REAL>::max();

	for(unsigned int i=0; i<rodptrs.size(); i++)
		for(unsigned int j=0; j<rodptrs[i]->l.size(); j++)
			minDt = min( minDt, rodptrs[i]->l[j]/_smoothMod(rodptrs[i]->v[j]) );

	return minDt;
}

REAL AdaptiveVelocityVerlet::_computeMinOneOverAngVel()
{
	const REAL oneDegree = (M_PI/180.0);
	REAL minDt = numeric_limits<REAL>::max();

	for(unsigned int i=0; i<rodptrs.size(); i++)
		for(unsigned int j=0; j<rodptrs[i]->w.size(); j++)
			minDt = min( minDt, oneDegree/_smoothMod(rodptrs[i]->w[j]) );

	return minDt;
}

REAL AdaptiveVelocityVerlet::_initial_dt(const REAL t, const REAL dt_min, const REAL dt_max)
{
	if(t > 10.0)
		return dt_max;

	const REAL dt_current = dt_min*exp(t);

	return (dt_current<=dt_max)?dt_current:dt_max;
}


REAL AdaptiveVelocityVerlet::_computeTimeReparameterization(const REAL dt0, const REAL maxDt, const REAL minDt, const REAL time)
{
	const REAL dtAccCFL = sqrt( CFL*_computeMinLengthOverAcc() );
	const REAL dtAngAccCFL = sqrt( CFL*_computeMinOneOverAngAcc() );
	const REAL dtVelCFL = CFL*_computeMinLengthOverVel();
	const REAL dtAngVelCFL = CFL*_computeMinOneOverAngVel();

	const REAL current_dt = max(min(min(min(min(dtAccCFL,dtAngAccCFL),dtVelCFL),dtAngVelCFL),maxDt),minDt);

	const REAL RL = 0.01;
	dtRunning = (time==0.0)?current_dt:( (1.0-RL)*dtRunning + RL*current_dt );

	return dt0/dtRunning;
}


REAL AdaptiveVelocityVerlet::integrate(const REAL time, const REAL dt0, const int step)
{
	// Adaptive Velocity Verlet (Note that dt is substituted by n, since the time is reparameterized!)

	// Start:	a(n)
	//			rho(n) = U(n)

	//	1.	x(n+1/2) = x(n) + dt0/(2*rho(n))*v(n)
	//	2.	a(n+1/2)
	//	3.	v(n+1/2) = v(n) + dt0/(2*rho(n))*a(n+1/2)
	//	4.	U(n+1/2)
	//	5.	rho(n+1) = 2*U(n+1/2) - rho(n) ---> explicit update (subject to oscillations, better implicit)
	//	6.	v(n+1) = v(n+1/2) + dt0/(2*rho(n+1))*a(n+1/2)
	//	7.	x(n+1) = x(n+1/2) + dt0/(2*rho(n+1))*v(n+1)
	//	8.	t(n+1) = t(n) + 0.5*dt0/rho(n) + 0.5*dt0/rho(n+1)

	const REAL dtMin = 1e-7;
	const REAL dtMax = _initial_dt(time,1e-5,1e-2);

	const unsigned int numRods = rodptrs.size();

	if( time==0.0 )
	{
		// Enforce boundary conditions (position+velocity)
		for(unsigned int j=0; j<numRods; j++)
		{
			(*bcptrs[j]).neumann(*rodptrs[j], step, dt0, time);
			(*bcptrs[j]).dirichlet(*rodptrs[j], step, dt0, time);
		}

		// Compute acceleration: a(n)
		{
			for(unsigned int j=0; j<numRods; j++)
				rodptrs[j]->update();

			_computeAccelerations(time);
		}

		// Compute time scaling factor: rho(n) = U(n)
		const REAL U = _computeTimeReparameterization(dt0, dtMax, dtMin, time);
		rho = U;
	}

	// Compute half step position: x(n+1/2) = x(n) + dt0/(2*rho(n))*v(n)
	for(unsigned int j=0; j<numRods; j++)
	{
		rodptrs[j]->x += 0.5*dt0/rho * rodptrs[j]->v;
		rodptrs[j]->Q =  vExp(0.5*dt0/rho * rodptrs[j]->w) * rodptrs[j]->Q;

		// Enforce boundary conditions for position
		(*bcptrs[j]).dirichlet(*rodptrs[j], step, dt0, time);
	}

	// Compute half step acceleration: a(n+1/2)
	{
		for(unsigned int j=0; j<numRods; j++)
			rodptrs[j]->update();

		_computeAccelerations(time);
	}

	// Compute half step velocity: v(n+1/2) = v(n) + dt0/(2*rho(n))*a(n+1/2)
	for(unsigned int j=0; j<numRods; j++)
	{
		rodptrs[j]->v += 0.5*dt0/rho * rodptrs[j]->a;
		rodptrs[j]->w += 0.5*dt0/rho * rodptrs[j]->wDot;

		// Enforce boundary conditions for velocity
		(*bcptrs[j]).neumann(*rodptrs[j], step, 0.5*dt0/rho, time);
	}

	// Compute time scaling factor: U(n+1/2) and rho(n+1) = 2*U(n+1/2) - rho(n)
	const REAL U = _computeTimeReparameterization(dt0, dtMax, dtMin, time);
	const REAL rhoOld = rho;
	rho = 2.0*U-rhoOld;
	assert(rho>0.0);

	// Compute final velocity: v(n+1) = v(n+1/2) + dt0/(2*rho(n+1))*a(n+1/2)
	for(unsigned int j=0; j<numRods; j++)
	{
		rodptrs[j]->v += 0.5*dt0/rho * rodptrs[j]->a;
		rodptrs[j]->w += 0.5*dt0/rho * rodptrs[j]->wDot;

		// Enforce boundary conditions for velocity
		(*bcptrs[j]).neumann(*rodptrs[j], step, 0.5*dt0/rho, time);
	}

	// Compute final positions: x(n+1) = x(n+1/2) + dt0/(2*rho(n+1))*v(n+1)
	for(unsigned int j=0; j<numRods; j++)
	{
		rodptrs[j]->x += 0.5*dt0/rho * rodptrs[j]->v;
		rodptrs[j]->Q =  vExp(0.5*dt0/rho * rodptrs[j]->w) * rodptrs[j]->Q;

		// Enforce boundary conditions for position
		(*bcptrs[j]).dirichlet(*rodptrs[j], step, 0.5*dt0/rho, time );
	}

	// Update rod
	for(unsigned int j=0; j<numRods; j++)
		rodptrs[j]->update();

	// Compute actual dt
	const REAL dt = 0.5*dt0/rhoOld + 0.5*dt0/rho;

	return dt;
}
