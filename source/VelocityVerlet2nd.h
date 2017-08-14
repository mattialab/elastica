/*
 * VelocityVerlet2nd.h
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#ifndef VELOCITYVERLET2ND_H_
#define VELOCITYVERLET2ND_H_

#include "PolymerIntegrator.h"

class VelocityVerlet2nd : public PolymerIntegrator
{
public:
	VelocityVerlet2nd(Vrodptr& _rodptrs, Vefptr& _efptrs, Vbcptr& _bcptrs, Vinterptr& _interptrs) : PolymerIntegrator(_rodptrs, _efptrs, _bcptrs, _interptrs){};

	virtual ~VelocityVerlet2nd(){};

	virtual REAL integrate(const REAL time, const REAL dt, const int step);
};

#endif /* VELOCITYVERLET2ND_H_ */

