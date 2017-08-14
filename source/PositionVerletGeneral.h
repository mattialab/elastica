/*
 * PositionVerletGeneral.h
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#ifndef POSITIONVERLETGENERAL_H_
#define POSITIONVERLETGENERAL_H_

#include "PolymerIntegrator.h"

class PositionVerletGeneral : public PolymerIntegrator
{
protected:
	VREAL cCoefs;
	VREAL dCoefs;
	const unsigned int order;

public:
	PositionVerletGeneral(Vrodptr& _rodptrs, Vefptr& _efptrs, Vbcptr& _bcptrs, Vinterptr& _interptrs, const unsigned int order);

	virtual ~PositionVerletGeneral(){};

	virtual REAL integrate(const REAL time, const REAL dt, const int step);
};

#endif /* POSITIONVERLETGENERAL_H_ */
