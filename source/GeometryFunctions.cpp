/*
 * GeometryFunctions.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: mgazzola
 */

#include "GeometryFunctions.h"

// Given a set of points x[3], it constructs a set of frames
// that are aligned with the tangent of the curve expressed by x[3]
std::vector<Matrix3> alignFrames(const std::vector<Vector3>&x, const Vector3& u0)
{
	// Check that u0 is indeed a unit vector
	assert( fabs(u0%u0-1.0) < Tolerance::tol() );

	// Compute tangent vectors
	std::vector<Vector3> t = vUnitize(vDiff(x));

	// Check that u0 is orthogonal to tangent
	assert( fabs(t[0] % u0) < Tolerance::tol() );

	// Create first body frame
	std::vector<Matrix3> Q;
	Q.push_back(Matrix3(u0, t[0] * u0, t[0]));

	// Index i starts from 1, because the first frame has been filled already
	for(unsigned int i=1; i<t.size(); i++)
	{
		// Compute angle of rotation: in this case, this is a counterclockwise rotation
		const REAL angleRotation = angle(t[i-1], t[i]);

		// Use cross product to compute axis of rotation in the lab frame of reference (LabFOR)
		const Vector3 axisRotationInLabFOR = (t[i-1] * t[i]).unitize();

		// Express the axis of rotation in LabFOR in the material frame of reference:
		// Formula: k_material = Q * (k_lab - x_material), where k_lab is the vector to be transformed,
		// x_material is its position in the material frame of reference
		const Vector3 axisRotationInMaterialFOR = (Q[i-1] * axisRotationInLabFOR).unitize();

		// Apply rotation using Rodrigues' rotation formula: Q_rotated = exp(angle * k) * Q.
		// This formula is based on the definition of skew map S.
		// In general it applies a CLOCKWISE rotation, but here Andrew for no reason decided to
		// implement S_andrew = -S_rest_of_the_world and so in this code the formula applies a
		// COUNTERCLOCKWISE rotation by default, hence no minus in front of 'angleRotation' below
		const Matrix3 rotatedQ = exp(angleRotation * axisRotationInMaterialFOR) * Q[i-1];

		// Store frame
		Q.push_back( rotatedQ );
	}

	return Q;
}

// Apply clockwise twist about the director d3 (the one paired to the tangent vector of the rod)
void applyTwists(std::vector<Matrix3>& Q, const std::vector<REAL>& twists)
{
	// Check consistency in size
	assert(Q.size() == twists.size());

	// Because Andrew had the brilliant idea to set S_andrew=-S_rest_of_the_world,
	// there is a minus in front of the axis of rotation
	for (unsigned int i=0; i<Q.size(); i++)
	{
		const Vector3 axisRotationInMaterialFOR = Vector3(0,0,twists[i]);
		Q[i] = exp(-axisRotationInMaterialFOR) * Q[i];
	}
}

// finds the minimum distance vector between two line segments
// with positions x_i and extent e_i, returning the point on the
// first segment plus the distance vector
std::vector< Vector3 > findMinDistVectors(const Vector3& x1, const Vector3& e1, const Vector3& x2, const Vector3& e2, const REAL tol)
{
	// we assume all e1 and e2 are non-0, problems occur if not
	// if they are within tol of parallel, we use another formula
	// note of course <a,b> = <b, a> for vectors so
	// we only compute one of them and use it for both
	const REAL e1e1 = e1%e1;
	const REAL e1e2 = e1%e2;
	const REAL e1x1 = e1%x1;
	const REAL e1x2 = e1%x2;
	const REAL e2e2 = e2%e2;
	const REAL e2x1 = e2%x1;
	const REAL e2x2 = e2%x2;
	REAL s = 0.0;
	REAL t = 0.0;

	// check if lines are nearly parallel; if so use diff formula
	// than the skew case.
	if (fabs((1 - e1e2*e1e2/(e1e1*e2e2)))<tol)
	{
		t = (e1x2 - e1x1)/e1e1;
		if (t<0) t = 0;
		if (t>1) t = 1;
		s = (e2%(x1+t*e1-x2))/e2e2;
		if (s<0) s = 0;
		if (s>1) s = 1;
	}
	else
	{
		// now for the skew case
		s = (e1e1*e2x1 + e1e2*e1x2 - e1e2*e1x1 - e1e1*e2x2) / (e1e1*e2e2 - e1e2*e1e2);
		t = (e1e2*s + e1x2 -e1x1) / e1e1;

		// if CP given by s and t above is in [0,1]x[0,1]
		// then return it unchanged
		// as it is the minimum. Otherwise the minimum will occur on the border
		// of the unit square.
		if (s<0 || s>1 || t<0 || t>1)
		{
			REAL sTemp, tTemp, dTemp;
			REAL curMinDist=1e20;

			// now check each edge of the square
			// s = 0
			tTemp = (e1x2 - e1x1)/e1e1;

			if (tTemp<0)
				tTemp = 0;

			if (tTemp>1)
				tTemp = 1;

			s = 0;
			t = tTemp;
			curMinDist = (x1 + e1*tTemp - x2).length();

			// s = 1
			tTemp = (e1x2 +e1e2 - e1x1)/e1e1;

			if (tTemp<0)
				tTemp = 0;

			if (tTemp>1)
				tTemp = 1;

			dTemp = (x1 + e1*tTemp - x2 - e1*1).length();

			if ( dTemp < curMinDist)
			{
				curMinDist = dTemp;
				s = 1;
				t = tTemp;
			}

			// t = 0
			sTemp = (e2x1 - e2x2)/e2e2;

			if (sTemp<0)
				sTemp = 0;

			if (sTemp>1)
				sTemp = 1;

			dTemp = (x2 + e2*sTemp - x1).length();

			if ( dTemp < curMinDist)
			{
				curMinDist = dTemp;
				s = sTemp;
				t = 0;
			}

			// t = 1
			sTemp = (e2x1 + e1e2 - e2x2)/e2e2;

			if (sTemp<0)
				sTemp = 0;

			if (sTemp>1)
				sTemp = 1;

			dTemp = (x2 + e2*sTemp - x1 - e1).length();

			if ( dTemp < curMinDist)
			{
				curMinDist = dTemp;
				s = sTemp;
				t = 1;
			}
		}
	}
	// now construct the minimal distance vector using s and t
	std::vector<Vector3> result;
	result.push_back(x1 + t*e1);
	result.push_back(-x1 - t*e1 + x2 + s*e2);
	return result;
}

