/*
 * PositionVerletGeneral.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#include "PositionVerletGeneral.h"

PositionVerletGeneral::PositionVerletGeneral(Vrodptr& _rodptrs, Vefptr& _efptrs, Vbcptr& _bcptrs, Vinterptr& _interptrs, const unsigned int order):
PolymerIntegrator(_rodptrs, _efptrs, _bcptrs, _interptrs), order(order)
{
	cCoefs = VREAL();
	dCoefs = VREAL();
	REAL c1, c2, d1, d2;
	switch (order)
	{
	case 1:
		cCoefs.push_back(1.0);
		dCoefs.push_back(1.0);
		break;
	case 2:
		cCoefs.push_back(0.5);
		cCoefs.push_back(0.5);
		dCoefs.push_back(0.0);
		dCoefs.push_back(1.0);
		break;
	case 4:
		c1 = 1.0/(2.0*(2.0-pow(2.0,1.0/3.0)));
		c2 = (1.0-pow(2.0,1.0/3.0))*c1;
		d1 = 2.0*c1;
		d2 = -pow(2.0,1.0/3.0)*d1;
		cCoefs.push_back(c1);
		cCoefs.push_back(c2);
		cCoefs.push_back(c2);
		cCoefs.push_back(c1);
		dCoefs.push_back(0.0);
		dCoefs.push_back(d1);
		dCoefs.push_back(d2);
		dCoefs.push_back(d1);
		break;
	default:
		printf("Order %d symplectic integrator not implemented", order);
		throw "SI for this order not implemented!";
	}
}

// At second order it is exactly equivalent to the positionVerlet algorithm (but less efficient)
// I keep it because it offers a compact general way to for position Verlet schemes of different orders
REAL PositionVerletGeneral::integrate(const REAL time, const REAL dt, const int step)
{
	const unsigned int numRods = rodptrs.size();

	for(unsigned int i=0; i<order; i++)
	{
		_computeForces(time);

		// Compute second order dynamics (velocities)
		for(unsigned int j=0; j<numRods; j++)
		{
			rodptrs[j]->v += dCoefs[i] * dt * rodptrs[j]->totalForces / rodptrs[j]->m;
			rodptrs[j]->w += dCoefs[i] * dt * (rodptrs[j]->J0inv * rodptrs[j]->totalTorques);
		}

		// Enforce Neumann boundary conditions
		// Note that an internal flag in (*bcptrs[j]) switches automatically
		// so that Neumann or Dirichlet conditions are applied properly.
		for(unsigned int j=0; j<numRods; j++)
			(*bcptrs[j])(*rodptrs[j], step, dt, time);

		// Compute first order dynamics (positions and frames)
		for(unsigned int j=0; j<numRods; j++)
		{

			rodptrs[j]->x += cCoefs[i] * dt * rodptrs[j]->v;
			rodptrs[j]->Q =  vExp(rodptrs[j]->w * cCoefs[i] * dt) * rodptrs[j]->Q;
		}

		// Enforce Dirichlet boundary conditions
		for(unsigned int j=0; j<numRods; j++)
			(*bcptrs[j])(*rodptrs[j], step, dt, time);

		// Update rod
		for(unsigned int j=0; j<numRods; j++)
			rodptrs[j]->update();
	}

	return dt;
}


