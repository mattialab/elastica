/*
 * SplineProfileZeroEnds.cpp
 *
 *  Created on: Jul 2, 2014
 *      Author: mgazzola
 */

#include "SplineProfileZeroEnds.h"

SplineProfileZeroEnds::SplineProfileZeroEnds(vector<double> ctrlPoints)
    : N(1000), Nw(ctrlPoints.size() + 2), ndegree(3) {
  // Clear internal profile
  xProfile.clear();
  yProfile.clear();

  // Store control points
  this->ctrlPoints = ctrlPoints;

  // Check control parameter vector size
  if (this->ctrlPoints.size() < 3) {
    printf("WIDTH vector too small!\n");
    abort();
  }

  _computeAndStoreProfile();
}

SplineProfileZeroEnds::~SplineProfileZeroEnds() {}

double SplineProfileZeroEnds::_coxDeBoorRecursion(const int j, const int d,
                                                  const double *t,
                                                  const double u) const {
  const double thres = 1e-6;
  double fac1 = 0.0;
  double fac2 = 0.0;

  if (fabs((double)d - 1.0) < thres) return (t[j] <= u && u < t[j + 1]) ? 1 : 0;

  if (fabs(t[j + d - 1] - t[j]) < thres)
    fac1 = 0.0;
  else
    fac1 = (u - t[j]) / (t[j + d - 1] - t[j]);

  if (fabs(t[j + d] - t[j + 1]) < thres)
    fac2 = 0.0;
  else
    fac2 = (t[j + d] - u) / (t[j + d] - t[j + 1]);

  const double bOld = _coxDeBoorRecursion(j, d - 1, t, u);
  const double bOldRight = _coxDeBoorRecursion(j + 1, d - 1, t, u);

  return fac1 * bOld + fac2 * bOldRight;
}

void SplineProfileZeroEnds::_computeAndStoreProfile() {
  // Very inefficient implementation, but the filling of the W vector
  // is done once at the beginning and that s it, so no big deal!

  // Constants
  const unsigned int Nx = N;
  const unsigned int n = Nw - 1;
  const unsigned int d = ndegree + 1;
  const unsigned int m = n + d + 1;
  double t[m];

  // Allocate control points
  double Pwx[Nw];
  double Pwy[Nw];

  // Allocate x and y vectors and set them to zero
  double xr[Nx];
  double yr[Nx];
  for (unsigned int i = 0; i < Nx; i++) {
    xr[i] = 0.0;
    yr[i] = 0.0;
  }

  // Define control points, first the ends
  Pwx[0] = 0.0;
  Pwy[0] = 0.0;
  Pwx[Nw - 1] = 1.0;
  Pwy[Nw - 1] = 0.0;

  // Now the interior
  unsigned int counter = 0;
  const double dxw = 1.0 / ((double)Nw - 3.0);
  for (unsigned int i = 1; i < Nw - 1; i++) {
    Pwx[i] = (i - 1) * dxw;
    Pwy[i] = ctrlPoints[counter];
    counter++;
  }

  // Build t vector
  for (unsigned int i = 0; i < m; i++) {
    if (i < d)
      t[i] = 0.0;
    else if (i > n)
      t[i] = n - d + 2.0;
    else
      t[i] = i - d + 1.0;
  }

  for (unsigned int i = 0; i < m; i++) t[i] /= t[m - 1];

  // Build the spline
  const double du = t[m - 1] / ((double)Nx - 1.0);
  for (unsigned int i = 0; i < Nx; i++) {
    const double u = (double)i * du;
    for (unsigned int j = 0; j < n + 1; j++) {
      xr[i] += Pwx[j] * _coxDeBoorRecursion(j, d, t, u);
      yr[i] += Pwy[j] * _coxDeBoorRecursion(j, d, t, u);
    }
  }

  xr[Nx - 1] = Pwx[Nw - 1];
  yr[Nx - 1] = Pwy[Nw - 1];

  for (unsigned int i = 0; i < Nx; i++) {
    xProfile.push_back(xr[i]);
    yProfile.push_back(yr[i]);
  }
}

double SplineProfileZeroEnds::getProfile(const double ss) const {
  assert((ss >= 0.0) && (ss <= 1.0));

  // Generate profile via linear interpolation
  const unsigned int Nx = N;
  const double s = ss;
  double width = 0.0;

  unsigned int idx = 0;
  for (unsigned int i = 0; i < Nx; i++) {
    if (xProfile[i] >= s) {
      idx = std::max(0, (int)i);
      break;
    }
  }

  const unsigned int idxPlus = std::min(int(idx + 1), int(Nx - 1));
  const double denominator =
      (idx == idxPlus) ? 1.0 : (xProfile[idxPlus] - xProfile[idx]);
  width = yProfile[idx] + (yProfile[idxPlus] - yProfile[idx]) *
                              (s - xProfile[idx]) / denominator;

  return width;
}
