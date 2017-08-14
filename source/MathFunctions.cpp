/*
 * MathFunctions.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: mgazzola
 */

#include "MathFunctions.h"

// Returns fmod shifted into [0,r]
REAL posMod(const REAL a, const REAL r)
{
	REAL v = fmod(a, r);

	if (v < 0)
		v += r;

	return v;
}

// Calls the normal acos func, but if x>1 returns 0, and if x<-1 returns PI,
// so as to get around retarded precision errors making acos(cos(0)) = nan
REAL arcCos(const REAL x)
{
	const REAL argumentACosClamped = std::max((REAL)-1.0,(REAL)std::min(x,(REAL)1.0));
	return acos(argumentACosClamped);
}

// Gaussian random normal generator without trigonometric calls using Box-Muller transform
double randn_notrig(const double mu, const double sigma)
{
	double var1 = 0.0;
	double var2 = 0.0;
	double rsquared = 0.0;

	// If no deviate has been stored, the polar Box-Muller transformation is
	// performed, producing two independent normally-distributed random
	// deviates.  One is stored for the next round, and one is returned.

	// Choose pairs of uniformly distributed deviates, discarding those
	// that don't fall within the unit circle
	do
	{
		var1 = 2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
		var2 = 2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
		rsquared=var1*var1+var2*var2;
	} while ( rsquared>=1.0 || rsquared == 0.0);

	// Calculate polar tranformation for each deviate
	const double polar = sqrt(-2.0*log(rsquared)/rsquared);

	// Return second deviate
	return var2*polar*sigma + mu;
}

/*
// PREVIOUS VERSION FUCKIN ANDREW
double randn_notrig(const double mu, const double sigma)
{
	bool deviateAvailable = false;	//flag
	float storedDeviate;		//deviate from previous calculation
	double polar, rsquared, var1, var2;

	// If no deviate has been stored, the polar Box-Muller transformation is
	// performed, producing two independent normally-distributed random
	// deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable)
	{
		// Choose pairs of uniformly distributed deviates, discarding those
		// that don't fall within the unit circle
		do
		{
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);

		// Calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);

		// Store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;

		// Return second deviate
		return var2*polar*sigma + mu;
	}
	else
	{
		// If a deviate is available from a previous call to this function, it is
		// returned, and the flag is set to false.
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}
*/



