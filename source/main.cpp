/*
 * main.cpp
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#include "UsualHeaders.h"
#include "ArgumentParser.h"
#include "Test.h"
#include "MRAGEnvironment.h"
#include "EulerBuckling.h"
#include "EulerBeam.h"
#include "InstabilityHelical.h"
#include "MassSpringSystem.h"
#include "MitchellBuckling.h"
#include "TorsionalWavesCouple.h"
#include "TorsionalWavesCoupleStretch.h"
#include "RollingFrictionInitialVelocity.h"
#include "RollingFrictionInclinedPlane.h"
#include "RollingFrictionTorque.h"
#include "AxialFriction.h"
#include "QuasistaticTimoshenkoBeam.h"
#include "SlenderBodyStokes.h"
#include "Snake.h"
#include "SolenoidsJCP.h"
#include "Solenoids.h"
#include "StigmaticStart.h"

#ifdef SNAKE_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

using namespace std;

Test * test = NULL;

#ifdef SNAKE_VIZ
struct VisualSupport
{
	static void display()
	{
	}

	static void idle(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		test->run();
		glutSwapBuffers();
	}

	static void run(int argc, const char ** argv)
	{
		static bool bSetup = false;

		if (!bSetup)
		{
			setup(argc, argv);
			bSetup = true;
		}

		glutDisplayFunc(display);
		glutIdleFunc(idle);
		glutMainLoop();
	}

	static void setup(int argc,  const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		//glutInitWindowSize(1024,1024);
		glutInitWindowSize(700,700);
		glutCreateWindow("School");
		glutDisplayFunc(display);
		//glClearColor(1,1,1,1);
		glClearColor(0,0,0,1);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, 1, 0, 1, 0, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
};
#endif

int main(const int argc, const char** argv)
{
	MRAG::ArgumentParser parser(argc, argv);

	MRAG::Environment::setup(max(1, parser("-nthreads").asInt()));

	const string studycase = parser("-study").asString();

	if( studycase == "EULER_BUCKLING" )
		test = new EulerBuckling(argc, argv);
	else if( studycase == "EULER_BEAM" )
		test = new EulerBeam(argc, argv);
	else if( studycase == "SNAKE" )
		test = new Snake(argc, argv);
	else if( studycase == "STIGMATIC_START" )
		test = new StigmaticStart(argc, argv);
	else if( studycase == "MITCHELL_BUCKLING" )
		test = new MitchellBuckling(argc, argv);
	else if( studycase == "LOCALIZED_HELICAL_INSTABILITY" )
		test = new InstabilityHelical(argc, argv);
	else if( studycase == "MASS_SPRING_SYSTEM" )
		test = new MassSpringSystem(argc, argv);
	else if( studycase == "TORSIONAL_WAVES_COUPLE" )
		test = new TorsionalWavesCouple(argc, argv);
	else if( studycase == "TORSIONAL_WAVES_COUPLE_STRETCH" )
		test = new TorsionalWavesCoupleStretch(argc, argv);
	else if( studycase == "QUASISTATIC_TIMOSHENKO_BEAM" )
		test = new QuasistaticTimoshenkoBeam(argc, argv);
	else if( studycase == "ROLLING_FRICTION_INITIAL_VELOCITY" )
		test = new RollingFrictionInitialVelocity(argc, argv);
	else if( studycase == "ROLLING_FRICTION_INCLINED_PLANE" )
		test = new RollingFrictionInclinedPlane(argc, argv);
	else if( studycase == "ROLLING_FRICTION_TORQUE" )
		test = new RollingFrictionTorque(argc, argv);
	else if( studycase == "AXIAL_FRICTION" )
		test = new AxialFriction(argc, argv);
	else if( studycase == "SLENDER_BODY_STOKES" )
		test = new SlenderBodyStokes(argc, argv);
	else if( studycase == "SOLENOIDS_JCP" )
		test = new SolenoidsJCP(argc, argv);
	else if( studycase == "SOLENOIDS" )
		test = new Solenoids(argc, argv);
	else
	{
		printf("Study case not defined!\n");
		abort();
	}


	try
	{
#ifdef SNAKE_VIZ
		VisualSupport::run(argc, argv);
#else
		test->run();
#endif

		return 0;
	}
	catch(string &exc)
	{
		cout << exc << endl;
	}
	catch(char const* exc)
	{
		cout << exc << endl;
	}
}



