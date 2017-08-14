/*
 * Polymer.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: mgazzola
 */

#include "Polymer.h"

void Polymer::computeEnergies()
{
	bendingEnergy = 0.0;
	shearEnergy = 0.0;
	translationalEnergy = 0.0;
	rotationalEnergy = 0.0;
	totalInternalEnergy = 0;

	for (int i=0; i<numRods; i++)
	{
		rodptrs[i]->computeEnergies();
		bendingEnergy += rodptrs[i]->bendingEnergy;
		shearEnergy += rodptrs[i]->shearEnergy;
		translationalEnergy += rodptrs[i]->translationalEnergy;
		rotationalEnergy += rodptrs[i]->rotationalEnergy;
		totalInternalEnergy += rodptrs[i]->totalInternalEnergy;
	}
}

bool Polymer::simulate(const REAL simulationTime, const REAL dt0, const unsigned int diagPerUnitTime, const unsigned int povrayFramesPerUnitTime, const string diagnostics, const string integrationType, const REAL CFL)
{
#ifdef SNAKE_VIZ
	long unsigned int povrayStep = 0;
#endif

	unsigned int N = 1;
	long unsigned int step = 0;
	double time = 0.0;
	double diagFotoTimer = 0.0;
	double povrayFotoTimer = 0.0;

	Rod *snake = rodptrs[0];

	while(time <= simulationTime)
	{
		// Dump diagnostics
#ifdef SNAKE_VIZ
		{
			const double fotoDT = 1.0/( (diagPerUnitTime>0)?diagPerUnitTime:1e-6 ) ;
			if( (diagPerUnitTime!=0) && (diagFotoTimer > fotoDT || step==0) )
			{
				diagFotoTimer = 0.0;
				printEnergies(step, time);
				printX(step, time, diagnostics);
				printXV(step, time, diagnostics);
			}
		}
#endif

		// Dump Povray files
#ifdef SNAKE_VIZ
		const double fotoDT = 1.0/( (povrayFramesPerUnitTime>0)?povrayFramesPerUnitTime:1e-6 ) ;
		if( (povrayFramesPerUnitTime!=0) && (povrayFotoTimer > fotoDT || step==0) )
		{
			povrayFotoTimer = 0.0;
			char bufPov[1000];
			sprintf(bufPov, "snake_%07d.pov", int(povrayStep));
			char bufData[1000];
			sprintf(bufData, "data_%07d.inc", int(povrayStep));
			snake->dumpPovray(bufPov, bufData);
			povrayStep++;
		}
#endif

#ifdef SNAKE_VIZ
		if (step%5000 == 0)
		{
			_paint(snake);
			cout << "time = " << time << endl;
		}
#endif

		// Integrate
		REAL dt = 0.0;
		dt = pint->integrate(time, dt0, step);

		// Update cumulative quantities
		time += dt;
		diagFotoTimer += dt;
		povrayFotoTimer += dt;
		step += 1;

		// Window stats
		if (time >= startTimeStats && time <= endTimeStats)
		{
			avgVel = (avgVel*N + vMean(snake->v))/(N+1);
			N += 1;
		}

#ifndef NDEBUG
		if (nanCheck())
		{
			cout << "Found a NaN. Ending simulation." << endl;
			return false;
		}
#endif
	}

	// Compute energies
	computeEnergies();

	return true;
}

bool Polymer::nanCheck()
{
	bool foundNan = false;
	for (int i=0; i<numRods; i++)
		if (rodptrs[i]->nanCheck())

	foundNan = true;
	return foundNan;
}

void Polymer::printEnergies(const int step, const REAL time)
{
	computeEnergies();

	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"rodEnergies_%05d",i);

		FILE * outfile = fopen( buffer, (step==0)?"w":"a" );
		fprintf( outfile, "%1.10e %1.10e %1.10e %1.10e %1.10e %1.10e\n", time, totalInternalEnergy, translationalEnergy, rotationalEnergy, bendingEnergy, shearEnergy );
		fclose( outfile );
	}
}

void Polymer::printX(const int step, const REAL time, const string outfilename)
{
	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_rodX_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, (step==0)?"w":"a" );
		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->x[j].x );
		fprintf( outfile, "\n");
		fclose( outfile );
	}

	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		const Matrix3 Q0 = rodptrs[i]->Q.front();
		const Vector3 axisOfRotation = Vector3(0,0,1);

		char buffer[1000];
		sprintf(buffer,"%s_rodX_angle_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, (step==0)?"w":"a" );
		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->Q.size(); j++)
		{
			const Vector3 angvector = (Q0*rodptrs[i]->Q[j].T()).log();
			const REAL length = angvector.length();
			const int sign = floor( (angvector % (length*axisOfRotation) ) / ( length*length ) );
			const REAL angle = -sign*length;
			fprintf( outfile, "%1.10e ", angle );
		}
		fprintf( outfile, "\n");
		fclose( outfile );
	}

	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		const Matrix3 Q0 = rodptrs[i]->Q.front();
		const Vector3 axisOfRotation = Vector3(1,0,0);

		char buffer[1000];
		sprintf(buffer,"%s_rodX_anglebend_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, (step==0)?"w":"a" );
		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->Q.size(); j++)
		{
			const Vector3 angvector = (Q0*rodptrs[i]->Q[j].T()).log();
			const REAL length = angvector.length();
			const int sign = floor( (angvector % (length*axisOfRotation) ) / ( length*length ) );
			const REAL angle = -sign*length;
			fprintf( outfile, "%1.10e ", angle );
		}
		fprintf( outfile, "\n");
		fclose( outfile );
	}

	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_rodX_sigma_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, (step==0)?"w":"a" );
		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->Q.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->shearStrain0[j].length() );

		fprintf( outfile, "\n");
		fclose( outfile );
	}
}

void Polymer::printXV(const int step, const REAL time, const string outfilename)
{
	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_rodXV_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, (step==0)?"w":"a" );

		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->x[j].x );
		fprintf( outfile, "\n");

		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->x[j].y );
		fprintf( outfile, "\n");

		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->x[j].z );
		fprintf( outfile, "\n");

		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->v[j].x );
		fprintf( outfile, "\n");

		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->v[j].y );
		fprintf( outfile, "\n");

		fprintf( outfile, "%1.10e ", time );
		for (unsigned int j=0; j<rodptrs[i]->x.size(); j++)
			fprintf( outfile, "%1.10e ", rodptrs[i]->v[j].z );
		fprintf( outfile, "\n");

		fclose( outfile );
	}
}

void Polymer::print_s_internalTorques(const string outfilename)
{
	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, "w" );

		const vector<REAL> s = vCumSum(rodptrs[i]->l);

		assert(s.size()==rodptrs[i]->bendingInternalTorques0.size()+1);

		for (unsigned int j=0; j<rodptrs[i]->bendingInternalTorques0.size(); j++)
			fprintf( outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j], rodptrs[i]->bendingInternalTorques0[j].x, rodptrs[i]->bendingInternalTorques0[j].y, rodptrs[i]->bendingInternalTorques0[j].z );
		fprintf( outfile, "\n");

		fclose( outfile );
	}
}

void Polymer::print_s_coordinates(const string outfilename)
{
	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, "w" );

		const vector<REAL> s = vCumSum(rodptrs[i]->l);

		assert(s.size()+1==rodptrs[i]->x.size());
		assert(s.size()>1);

		const REAL zero = 0.0;
		fprintf( outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", zero, rodptrs[i]->x[0].x, rodptrs[i]->x[0].y, rodptrs[i]->x[0].z );
		for (unsigned int j=1; j<rodptrs[i]->x.size(); j++)
		{
			assert((j-1)<s.size());
			fprintf( outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j-1], rodptrs[i]->x[j].x, rodptrs[i]->x[j].y, rodptrs[i]->x[j].z );
		}
		fprintf( outfile, "\n");

		fclose( outfile );
	}
}

void Polymer::print_s_internalShears(const string outfilename)
{
	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, "w" );

		vector<REAL> s = vector<REAL>(rodptrs[i]->shearInternalForces0.size());
		assert(s.size()==rodptrs[i]->shearInternalForces0.size());
		s[0] = rodptrs[i]->l[0]/2.0;
		for(unsigned int j=1; j<s.size(); j++)
			s[j] = s[j-1] + (rodptrs[i]->l[j-1]/2.0 + rodptrs[i]->l[j]/2.0);

		for (unsigned int j=0; j<rodptrs[i]->shearInternalForces0.size(); j++)
		{
			const Vector3 shearInLabFrame = rodptrs[i]->Q[j].T()*rodptrs[i]->shearInternalForces0[j];
			fprintf( outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j], shearInLabFrame.x, shearInLabFrame.y, shearInLabFrame.z );
		}
		fprintf( outfile, "\n");

		fclose( outfile );
	}
}

void Polymer::print_s_curvatures(const string outfilename)
{
	for (unsigned int i=0; i<rodptrs.size(); i++)
	{
		char buffer[1000];
		sprintf(buffer,"%s_%05d", outfilename.c_str(), i);

		FILE * outfile = fopen( buffer, "w" );

		const vector<REAL> s = vCumSum(rodptrs[i]->l);

		assert(s.size()==rodptrs[i]->k0.size()+1);

		for (unsigned int j=0; j<rodptrs[i]->k0.size(); j++)
		{
			const REAL ed =  rodptrs[i]->ed[j];
			fprintf( outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j], rodptrs[i]->k0[j].x/ed, rodptrs[i]->k0[j].y/ed, rodptrs[i]->k0[j].z/ed );
		}
		fprintf( outfile, "\n");

		fclose( outfile );
	}
}

#ifdef SNAKE_VIZ
void Polymer::_paint(Rod * snake)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushAttrib(GL_ENABLE_BIT);

	snake->paint();

	glPopAttrib();
	glutSwapBuffers();
}
#endif

