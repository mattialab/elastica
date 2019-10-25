/*
 * SplineProfileZeroEnds.h
 *
 *  Created on: Jul 2, 2014
 *      Author: mgazzola
 */

#ifndef SPLINEPROFILEZEROENDS_H_
#define SPLINEPROFILEZEROENDS_H_

#include "UsualHeaders.h"

using namespace std;

class SplineProfileZeroEnds {
  vector<double> ctrlPoints;
  vector<double> xProfile;
  vector<double> yProfile;
  double _coxDeBoorRecursion(const int j, const int d, const double *t,
                             const double u) const;
  void _computeAndStoreProfile();

 public:
  // Number of discretization points for spline.
  // It should be at any time larger the the number of
  // discretization points of a rod - check assert in spline force class)
  const unsigned int N;

  // Total number of control points, two of which are anchored to zero at the
  // extrema therefore the number of actual control parameters is
  // Nw=ctrlPoints.size()+2
  const unsigned int Nw;

  // Spline degree
  const unsigned int ndegree;

  SplineProfileZeroEnds(vector<double> ctrlPoints);
  virtual ~SplineProfileZeroEnds();

  double getProfile(const double ss) const;
};

#endif /* SPLINEPROFILEZEROENDS_H_ */
