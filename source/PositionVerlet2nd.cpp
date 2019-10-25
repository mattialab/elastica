/*
 * PositionVerlet2nd.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: mgazzola
 */

#include "PositionVerlet2nd.h"

REAL PositionVerlet2nd::integrate(const REAL time, const REAL dt,
                                  const int step) {
  // Position Verlet
  //	1. x(t+dt/2) = x(t) + 0.5*dt*v(t)
  //	2. v(t+dt) = v(t) + dt*a(t+dt/2)
  //	3. x(t+dt) = x(t+dt/2) + 0.5*dt*v(t+dt)

  // Get the number of rods in the system
  const unsigned int numRods = rodptrs.size();

  // Set dt
  for (unsigned int j = 0; j < numRods; j++) (*rodptrs[j]).setDt(dt);

  // Enforce boundary conditions (position+velocity)
  if (time == 0.0) {
    for (unsigned int j = 0; j < numRods; j++) {
      (*rodptrs[j]).setTime(time);
      (*bcptrs[j]).neumann(*rodptrs[j], step, dt, time);
      (*bcptrs[j]).dirichlet(*rodptrs[j], step, dt, time);
    }
  }

  // Compute half step position: x(t+dt/2) = x(t) + 0.5*dt*v(t)
  for (unsigned int j = 0; j < numRods; j++) {
    _v_update_x(
        0.5 * dt,
        rodptrs[j]);  // (in-place) rodptrs[j]->x += 0.5 * dt * rodptrs[j]->v;
    _v_update_Q(0.5 * dt,
                rodptrs[j]);  // (in-place) rodptrs[j]->Q =  vExp(0.5 * dt *
                              // rodptrs[j]->w) * rodptrs[j]->Q;

    // Enforce boundary conditions for position
    (*rodptrs[j]).setTime(time + 0.5 * dt);
    (*bcptrs[j]).dirichlet(*rodptrs[j], step, dt, time + 0.5 * dt);
  }

  // Compute half step acceleration: a(t+dt/2)
  {
    for (unsigned int j = 0; j < numRods; j++) rodptrs[j]->update(time);

    _computeAccelerations(time + 0.5 * dt, step);
  }

  // Compute final velocity: v(t+dt) = v(t) + dt*a(t+dt/2)
  for (unsigned int j = 0; j < numRods; j++) {
    _v_update_v(dt,
                rodptrs[j]);  // (in-place) rodptrs[j]->v += dt * rodptrs[j]->a;
    _v_update_w(
        dt, rodptrs[j]);  // (in-place) rodptrs[j]->w += dt * rodptrs[j]->wDot;

    // Enforce boundary conditions for velocity
    (*rodptrs[j]).setTime(time + dt);
    (*bcptrs[j]).neumann(*rodptrs[j], step, dt, time + dt);
  }
  // cout << "Reach4"<< endl;
  // Compute final positions: x(t+dt) = x(t+dt/2) + 0.5*dt*v(t+dt)
  for (unsigned int j = 0; j < numRods; j++) {
    _v_update_x(
        0.5 * dt,
        rodptrs[j]);  // (in-place) rodptrs[j]->x += 0.5 * dt * rodptrs[j]->v;
    _v_update_Q(0.5 * dt,
                rodptrs[j]);  // (in-place) rodptrs[j]->Q =  vExp(0.5 * dt *
                              // rodptrs[j]->w) * rodptrs[j]->Q;

    // Enforce boundary conditions for position
    (*bcptrs[j]).dirichlet(*rodptrs[j], step, dt, time + dt);
  }

  return dt;
}
