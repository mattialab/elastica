/*
 * AdaptiveVelocityVerlet.h
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#ifndef ADAPTIVEVELOCITYVERLET_H_
#define ADAPTIVEVELOCITYVERLET_H_

#include "PolymerIntegrator.h"

class AdaptiveVelocityVerlet : public PolymerIntegrator
{
protected:
	const REAL CFL;
	REAL dtRunning;
	REAL rho;

	REAL _smoothMod(Vector3 & vec);
	REAL _computeMinLengthOverAcc();
	REAL _computeMinOneOverAngAcc();
	REAL _computeMinLengthOverVel();
	REAL _computeMinOneOverAngVel();
	REAL _initial_dt(const REAL t, const REAL dt_min, const REAL dt_max);
	REAL _computeTimeReparameterization(const REAL dt0, const REAL maxDt, const REAL minDt, const REAL time);

public:
	AdaptiveVelocityVerlet(Vrodptr& _rodptrs, Vefptr& _efptrs, Vbcptr& _bcptrs, Vinterptr& _interptrs, const REAL _CFL) :
		PolymerIntegrator(_rodptrs, _efptrs, _bcptrs, _interptrs), CFL(_CFL), dtRunning(0.0), rho(0.0)
	{
	};

	virtual ~AdaptiveVelocityVerlet(){};

	virtual REAL integrate(const REAL time, const REAL dt, const int step);
};

#endif /* ADAPTIVEVELOCITYVERLET_H_ */
