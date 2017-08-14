/*
 * PositionVerlet2nd.h
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#ifndef POSITIONVERLET2ND_H_
#define POSITIONVERLET2ND_H_

#include "PolymerIntegrator.h"

class PositionVerlet2nd : public PolymerIntegrator
{
public:
	PositionVerlet2nd(Vrodptr& _rodptrs, Vefptr& _efptrs, Vbcptr& _bcptrs, Vinterptr& _interptrs) : PolymerIntegrator(_rodptrs, _efptrs, _bcptrs, _interptrs){};

	virtual ~PositionVerlet2nd(){};

	virtual REAL integrate(const REAL time, const REAL dt, const int step);
};

#endif /* POSITIONVERLET2ND_H_ */
