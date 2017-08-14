/*
 * Vector3.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: mgazzola
 */

#include "Vector3.h"

// constructors
Vector3::Vector3() : x(0), y(0), z(0)
{
}

Vector3::Vector3(REAL x, REAL y, REAL z) : x(x), y(y), z(z)
{
}

Vector3::Vector3(const Vector3& V) : x(V.x), y(V.y), z(V.z)
{
}

// io
std::ostream& operator<< (std::ostream& out, const Vector3& A)
{
	char buffer[1000];
	sprintf(buffer,"%10.10e %10.10e %10.10e", A.x, A.y, A.z);
	out << buffer;
	return out;
}

///// Additional vector functions not defined as member functions //////

Vector3 operator*(REAL a, const Vector3& V)
{
	return V * a;
}

// gets the angle between two vector in range(0, PI)
REAL angle(const Vector3& A, const Vector3& B)
{
	REAL cosTheta = (A % B) / (A.length() * B.length());
	return arcCos(cosTheta);
}

// generates a uniform gaussian random variable in 3 dimensions
Vector3 randVector3()
{
	return Vector3(randn_notrig(0.,1.), randn_notrig(0.,1.), randn_notrig(0.,1.));
}

Vector3 randVector3(REAL mu, REAL sigma)
{
	return Vector3(randn_notrig(mu, sigma), randn_notrig(mu, sigma), randn_notrig(mu, sigma));
}

// random only in x direction
Vector3 randXOnly()
{
	return Vector3(randn_notrig(0.,1.), 0, 0);
}


