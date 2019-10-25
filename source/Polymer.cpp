/*
 * Polymer.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: mgazzola
 */

#include "Polymer.h"

void Polymer::computeEnergies() {
  bendingEnergy = 0.0;
  shearEnergy = 0.0;
  translationalEnergy = 0.0;
  rotationalEnergy = 0.0;
  totalInternalEnergy = 0;

  for (int i = 0; i < numRods; i++) {
    rodptrs[i]->computeEnergies();
    bendingEnergy += rodptrs[i]->bendingEnergy;
    shearEnergy += rodptrs[i]->shearEnergy;
    translationalEnergy += rodptrs[i]->translationalEnergy;
    rotationalEnergy += rodptrs[i]->rotationalEnergy;
    totalInternalEnergy += rodptrs[i]->totalInternalEnergy;
  }
}

bool Polymer::simulate(const REAL simulationTime, const REAL dt0,
                       const REAL diagPerUnitTime,
                       const REAL povrayFramesPerUnitTime,
                       const string diagnostics, const string integrationType,
                       const REAL CFL) {
#ifdef SNAKE_POV
  long unsigned int povrayStep = 0;
#endif

  unsigned int N = 1;
  long unsigned int step = 0;
  double time = 0.0;
  double diagFotoTimer = 0.0;
  double povrayFotoTimer = 0.0;

#ifdef FLAGELBOW
  // Defining parameters change for artificial muscles
  auto temp_intrinsic_curvature = rodptrs[42]->intrinsic_k0;
  REAL initial_density = rodptrs[42]->density;
  REAL initial_radius = rodptrs[42]->r[0];
  REAL final_radius = 2. * initial_radius;
#endif

  while (time <= simulationTime) {
    // Dump Povray files
#ifdef SNAKE_POV
    const double fotoDT =
        1.0 / ((povrayFramesPerUnitTime > 0) ? povrayFramesPerUnitTime : 1e-6);
    if ((povrayFramesPerUnitTime != 0) &&
        (povrayFotoTimer > fotoDT || step == 0)) {
      povrayFotoTimer = 0.0;
      char bufPov[1000];
      sprintf(bufPov, "snake_%07d.pov", int(povrayStep));
      char bufData[1000];
      sprintf(bufData, "data_%07d.inc", int(povrayStep));

#if defined FLAGELBOW || FLAGFLAGELLA || FLAGWALKER || FLAGMUSCULARSNAKE || \
    FLAGWING
      // Dump Povray file
      fstream f;
      f.open(bufPov, fstream::out);
      f << "#include \"scenepovray.inc\""
        << "\n";
      f << "#include \"" << bufData << "\""
        << "\n";
      f << "\n";
      f.close();

      // Dump data to be rendered in the pov file
      for (unsigned int i = 0; i < numRods; i++) {
        rodptrs[i]->dumpPovray(bufPov, bufData, i, time);
      }
#endif
      // Data outputs for each case
      ofstream matlab;
      matlab.open("Flagella.txt",
                  (step == 0) ? std::ofstream::out : std::ofstream::app);
      //-------------Flagella-----------
#ifdef FLAGFLAGELLA
      matlab << time << " " << (rodptrs[0]->x[9][0]) << " "
             << (rodptrs[0]->x[9][1]) << endl;
#endif
      //-------------Walker-----------
#ifdef FLAGWALKER
      const Vector3 direction1 =
          (rodptrs[1]->x[2] - rodptrs[1]->x[0]).unitize();
      const Vector3 direction2 =
          (rodptrs[1]->x[16] - rodptrs[1]->x[14]).unitize();

      const REAL Angle1 = abs(angle(direction1, direction2)) * 180.0 / M_PI;
      const REAL Angle2 =
          abs(angle(direction1, rodptrs[4]->edge[0])) * 180.0 / M_PI;
      const REAL Angle3 =
          abs(angle(direction2, rodptrs[5]->edge[0])) * 180.0 / M_PI;
      const Vector3 position1 = rodptrs[2]->x[6] + 1.75 * rodptrs[2]->Q[5][1];
      const Vector3 position2 = rodptrs[3]->x[6] - 1.75 * rodptrs[3]->Q[5][1];

      rodptrs[0]->MaxHeight = rodptrs[2]->MaxHeight + rodptrs[4]->MaxHeight;
      rodptrs[1]->MaxHeight = rodptrs[3]->MaxHeight + rodptrs[5]->MaxHeight;

      matlab << time << " " << rodptrs[0]->x[8][1] << " "
             << (position1 - position2).length() << " " << rodptrs[0]->Tforce
             << " " << rodptrs[0]->MaxHeight << " " << rodptrs[1]->Tforce << " "
             << rodptrs[1]->MaxHeight << endl;
#endif
      //-------------Snake-----------
      // Velocites
#if defined FLAGMUSCULARSNAKE || FLAGSNAKE
      REAL forwardV = 0.0;
      REAL sideV = 0.0;
      for (unsigned int ii = 0; ii < 101; ii++) {
        forwardV += (rodptrs[0]->v[ii]) % Vector3(-1.0, 0.0, 0.0);
        sideV += (rodptrs[0]->v[ii]) % Vector3(0.0, -1.0, 0.0);
      }
      matlab << time << " " << rodptrs[0]->x[50][0] << " "
             << rodptrs[0]->x[50][1] << " " << forwardV / 101 << " "
             << sideV / 101 << endl;
#endif
#ifdef FLAGWING
      //-------------Wing__Force-----------
      const int Rodsize = rodptrs.size();
      matlab << time << " " << rodptrs[0]->MaxHeight << " "
             << rodptrs[0]->Tforce << " " << rodptrs[1]->MaxHeight << " "
             << rodptrs[Rodsize - 4]->MaxHeight << " "
             << rodptrs[Rodsize - 3]->MaxHeight << " "
             << rodptrs[Rodsize - 2]->MaxHeight << " "
             << rodptrs[Rodsize - 1]->MaxHeight << endl;
#endif
      // matlab<<" "<<rodptrs[1]->x[10][0]<<" "<<rodptrs[1]->x[10][2]<<endl;
      matlab.close();

      // Outputting for visualization of the simple test cases
#if defined FLAGSPHERICALJOINT || FLAGHINGEJOINT || FLAGFIXEDJOINT || \
    FLAGPULLINGMUSCLE || FLAGSNAKE || FLAGHELICALBUCKLING
      char visData[1000];
      for (unsigned int i = 0; i < numRods; i++) {
        sprintf(visData, "rod%01d_%07d.txt", int(i + 1), int(povrayStep));
        ofstream pythonvis;
        pythonvis.open(visData, std::ofstream::app);
        const int rodsize = rodptrs[0]->n;
        for (unsigned int j = 0; j < (rodsize + 1); j++) {
          REAL radius = 0.0;
          if (j == 0)
            radius = rodptrs[i]->r[0];
          else if (j == rodsize)
            radius = rodptrs[i]->r[rodsize - 1];
          else
            radius = 0.5 * (rodptrs[i]->r[j - 1] + rodptrs[i]->r[j]);
          pythonvis << rodptrs[i]->x[j][0] << " " << rodptrs[i]->x[j][1] << " "
                    << rodptrs[i]->x[j][2] << " " << radius << " " << time
                    << endl;
        }
        pythonvis.close();
      }
#endif

      povrayStep++;
    }
#endif

#ifdef SNAKE_VIZ
    if (step % 5000 == 0) {
      _paint(snake);
      cout << "time = " << time << endl;
    }
#endif

    // Integrate
    REAL dt = 0.0;

    dt = pint->integrate(time, dt0, step);
    // cout << rodptrs[0]->x[0] << endl;
    // Update cumulative quantities
    time += dt;
    diagFotoTimer += dt;
    povrayFotoTimer += dt;
    step += 1;

    if (step % 25000 == 0) {
      cout << "Simulated time : " << time << endl;
    }

#ifdef FLAGELBOW
    // Parameters change for artificial muscle.  i=42.
    // REAL current_radius =
    //     initial_radius + 1.0 * (final_radius - initial_radius) * time /
    //                          (simulationTime);  // 1.0=max rate
    REAL current_radius = initial_radius + 1.0 *
                                               (final_radius - initial_radius) *
                                               time / (1.0);  // 1.0=max rate

    for (auto it = std::make_pair(rodptrs[42]->intrinsic_k0.begin(),
                                  temp_intrinsic_curvature.cbegin());
         it.first != rodptrs[42]->intrinsic_k0.end(); ++it.first, ++it.second) {
      it.first->z = it.second->z * initial_radius / current_radius;
    }

    std::fill(rodptrs[42]->r.begin(), rodptrs[42]->r.end(), current_radius);

    rodptrs[42]->density = initial_density * (initial_radius / current_radius) *
                           (initial_radius / current_radius);
#endif

    // Window stats
    if (time >= startTimeStats && time <= endTimeStats) {
      for (unsigned int i = 0; i < numRods; i++) {
        avgVel[i] = (avgVel[i] * N + vMean(rodptrs[i]->v)) / (N + 1);
      }
      // cout << avgVel[0] << endl;
      N += 1;
    }

#ifndef NDEBUG

    if (nanCheck()) {
      cout << "Found a NaN. Ending simulation." << endl;
      return false;
    }

#endif
  }

  // Compute energies
  computeEnergies();

  return true;
}

bool Polymer::nanCheck() {
  bool foundNan = false;
  for (int i = 0; i < numRods; i++)
    if (rodptrs[i]->nanCheck()) foundNan = true;
  return foundNan;
}

void Polymer::printEnergies(const int step, const REAL time) {
  computeEnergies();

  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "rodEnergies_%05d", i);

    FILE *outfile = fopen(buffer, (step == 0) ? "w" : "a");
    fprintf(outfile, "%1.10e %1.10e %1.10e %1.10e %1.10e %1.10e\n", time,
            totalInternalEnergy, translationalEnergy, rotationalEnergy,
            bendingEnergy, shearEnergy);
    fclose(outfile);
  }
}

void Polymer::printX(const int step, const REAL time,
                     const string outfilename) {
  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_rodX_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, (step == 0) ? "w" : "a");
    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].x);
    fprintf(outfile, "\n");
    fclose(outfile);
  }

  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    const Matrix3 Q0 = rodptrs[i]->Q.front();
    const Vector3 axisOfRotation = Vector3(0, 0, 1);

    char buffer[1000];
    sprintf(buffer, "%s_rodX_angle_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, (step == 0) ? "w" : "a");
    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->Q.size(); j++) {
      const Vector3 angvector = (Q0 * rodptrs[i]->Q[j].T()).log();
      const REAL length = angvector.length();
      const int sign =
          floor((angvector % (length * axisOfRotation)) / (length * length));
      const REAL angle = -sign * length;
      fprintf(outfile, "%1.10e ", angle);
    }
    fprintf(outfile, "\n");
    fclose(outfile);
  }

  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    const Matrix3 Q0 = rodptrs[i]->Q.front();
    const Vector3 axisOfRotation = Vector3(1, 0, 0);

    char buffer[1000];
    sprintf(buffer, "%s_rodX_anglebend_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, (step == 0) ? "w" : "a");
    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->Q.size(); j++) {
      const Vector3 angvector = (Q0 * rodptrs[i]->Q[j].T()).log();
      const REAL length = angvector.length();
      const int sign =
          floor((angvector % (length * axisOfRotation)) / (length * length));
      const REAL angle = -sign * length;
      fprintf(outfile, "%1.10e ", angle);
    }
    fprintf(outfile, "\n");
    fclose(outfile);
  }

  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_rodX_sigma_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, (step == 0) ? "w" : "a");
    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->Q.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->shearStrain0[j].length());

    fprintf(outfile, "\n");
    fclose(outfile);
  }
}

void Polymer::printXV(const int step, const REAL time,
                      const string outfilename) {
  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_rodXV_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, (step == 0) ? "w" : "a");

    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].x);
    fprintf(outfile, "\n");

    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].y);
    fprintf(outfile, "\n");

    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].z);
    fprintf(outfile, "\n");

    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->v[j].x);
    fprintf(outfile, "\n");

    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->v[j].y);
    fprintf(outfile, "\n");

    fprintf(outfile, "%1.10e ", time);
    for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
      fprintf(outfile, "%1.10e ", rodptrs[i]->v[j].z);
    fprintf(outfile, "\n");

    fclose(outfile);
  }
}

void Polymer::print_s_internalTorques(const string outfilename) {
  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, "w");

    const vector<REAL> s = vCumSum(rodptrs[i]->l);

    assert(s.size() == rodptrs[i]->bendingInternalTorques0.size() + 1);

    for (unsigned int j = 0; j < rodptrs[i]->bendingInternalTorques0.size();
         j++)
      fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j],
              rodptrs[i]->bendingInternalTorques0[j].x,
              rodptrs[i]->bendingInternalTorques0[j].y,
              rodptrs[i]->bendingInternalTorques0[j].z);
    fprintf(outfile, "\n");

    fclose(outfile);
  }
}

void Polymer::print_s_coordinates(const string outfilename) {
  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, "w");

    const vector<REAL> s = vCumSum(rodptrs[i]->l);

    assert(s.size() + 1 == rodptrs[i]->x.size());
    assert(s.size() > 1);

    const REAL zero = 0.0;
    fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", zero,
            rodptrs[i]->x[0].x, rodptrs[i]->x[0].y, rodptrs[i]->x[0].z);
    for (unsigned int j = 1; j < rodptrs[i]->x.size(); j++) {
      assert((j - 1) < s.size());
      fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j - 1],
              rodptrs[i]->x[j].x, rodptrs[i]->x[j].y, rodptrs[i]->x[j].z);
    }
    fprintf(outfile, "\n");

    fclose(outfile);
  }
}

void Polymer::print_s_internalShears(const string outfilename) {
  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, "w");

    vector<REAL> s = vector<REAL>(rodptrs[i]->shearInternalForces0.size());
    assert(s.size() == rodptrs[i]->shearInternalForces0.size());
    s[0] = rodptrs[i]->l[0] / 2.0;
    for (unsigned int j = 1; j < s.size(); j++)
      s[j] = s[j - 1] + (rodptrs[i]->l[j - 1] / 2.0 + rodptrs[i]->l[j] / 2.0);

    for (unsigned int j = 0; j < rodptrs[i]->shearInternalForces0.size(); j++) {
      const Vector3 shearInLabFrame =
          rodptrs[i]->Q[j].T() * rodptrs[i]->shearInternalForces0[j];
      fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j],
              shearInLabFrame.x, shearInLabFrame.y, shearInLabFrame.z);
    }
    fprintf(outfile, "\n");

    fclose(outfile);
  }
}

void Polymer::print_s_curvatures(const string outfilename) {
  for (unsigned int i = 0; i < rodptrs.size(); i++) {
    char buffer[1000];
    sprintf(buffer, "%s_%05d", outfilename.c_str(), i);

    FILE *outfile = fopen(buffer, "w");

    const vector<REAL> s = vCumSum(rodptrs[i]->l);

    assert(s.size() == rodptrs[i]->k0.size() + 1);

    for (unsigned int j = 0; j < rodptrs[i]->k0.size(); j++) {
      const REAL ed = rodptrs[i]->ed[j];
      fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j],
              rodptrs[i]->k0[j].x / ed, rodptrs[i]->k0[j].y / ed,
              rodptrs[i]->k0[j].z / ed);
    }
    fprintf(outfile, "\n");

    fclose(outfile);
  }
}

#ifdef SNAKE_VIZ
void Polymer::_paint(Rod *snake) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushAttrib(GL_ENABLE_BIT);

  snake->paint();

  glPopAttrib();
  glutSwapBuffers();
}
#endif
