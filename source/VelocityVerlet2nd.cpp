/*
 * VelocityVerlet2nd.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#include "VelocityVerlet2nd.h"

REAL VelocityVerlet2nd::integrate(const REAL time, const REAL dt, const int step)
{
	// Velocity Verlet
	//	1. v(t+dt/2) = v(t) + 0.5*dt*a(t)
	//	2. x(t+dt) = x(t) + dt*v(t+dt/2)
	//	3. v(t+dt) = v(t+dt/2) + 0.5*dt*a(t+dt)

	const unsigned int numRods = rodptrs.size();

	if( time==0.0 )
	{
		// Enforce boundary conditions (position+velocity)
		for(unsigned int j=0; j<numRods; j++)
		{
			(*bcptrs[j]).neumann(*rodptrs[j], step, 0, time);
			(*bcptrs[j]).dirichlet(*rodptrs[j], step, 0, time);
		}

		// Compute acceleration: a(t)
		{
			for(unsigned int j=0; j<numRods; j++)
				rodptrs[j]->update();

			_computeAccelerations(time);
		}
	}

	// Compute half step velocity: v(t+dt/2) = v(t) + 0.5*dt*a(t)
	for(unsigned int j=0; j<numRods; j++)
	{
		_v_update_v(0.5*dt, rodptrs[j]); // (in-place) rodptrs[j]->v += 0.5 * dt * rodptrs[j]->a;
		_v_update_w(0.5*dt, rodptrs[j]); // (in-place) rodptrs[j]->w += 0.5 * dt * rodptrs[j]->wDot;

		// Enforce boundary conditions for velocity
		(*bcptrs[j]).neumann(*rodptrs[j], step, 0, time);
	}

	// Compute final positions: x(t+dt) = x(t) + dt*v(t+dt/2)
	for(unsigned int j=0; j<numRods; j++)
	{
		_v_update_x(dt, rodptrs[j]); // (in-place) rodptrs[j]->x += dt * rodptrs[j]->v;
		_v_update_Q(dt, rodptrs[j]); // (in-place) rodptrs[j]->Q =  vExp(dt * rodptrs[j]->w) * rodptrs[j]->Q;

		// Enforce boundary conditions for position
		(*bcptrs[j]).dirichlet(*rodptrs[j], step, 0, time);
	}

	// Compute acceleration: a(t+dt)
	{
		for(unsigned int j=0; j<numRods; j++)
			rodptrs[j]->update();

		_computeAccelerations(time);
	}

	// Compute final velocity: v(t+dt) = v(t+dt/2) + 0.5*dt*a(t+dt)
	for(unsigned int j=0; j<numRods; j++)
	{
		_v_update_v(0.5*dt, rodptrs[j]); // (in-place) rodptrs[j]->v += 0.5 * dt * rodptrs[j]->a;
		_v_update_w(0.5*dt, rodptrs[j]); // (in-place) rodptrs[j]->w += 0.5 * dt * rodptrs[j]->wDot;

		// Enforce boundary conditions for velocity
		(*bcptrs[j]).neumann(*rodptrs[j], step, 0, time);
	}

	// Update rod
	for(unsigned int j=0; j<numRods; j++)
		rodptrs[j]->update();

	return dt;
}

