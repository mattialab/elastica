/*
 * main.cpp
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#include "ArgumentParser.h"
#include "Elbow.h"
#include "FixedJoint.h"
#include "Flagella.h"
#include "HingeJoint.h"
#include "InstabilityHelical.h"
#include "MRAGEnvironment.h"
#include "MuscularSnake.h"
#include "PullingMuscle.h"
#include "QuasistaticTimoshenkoBeam.h"
#include "Snake.h"
#include "SphericalJoint.h"
#include "Test.h"
#include "UsualHeaders.h"
#include "Walker.h"
#include "Wing.h"

#ifdef SNAKE_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

using namespace std;

Test *test = NULL;

#ifdef SNAKE_VIZ
struct VisualSupport {
  static void display() {}

  static void idle(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    test->run();
    glutSwapBuffers();
  }

  static void run(int argc, const char **argv) {
    static bool bSetup = false;

    if (!bSetup) {
      setup(argc, argv);
      bSetup = true;
    }

    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMainLoop();
  }

  static void setup(int argc, const char **argv) {
    glutInit(&argc, const_cast<char **>(argv));
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    // glutInitWindowSize(1024,1024);
    glutInitWindowSize(700, 700);
    glutCreateWindow("School");
    glutDisplayFunc(display);
    // glClearColor(1,1,1,1);
    glClearColor(0, 0, 0, 1);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  }
};
#endif

int main(const int argc, const char **argv) {
  MRAG::ArgumentParser parser(argc, argv);

  MRAG::Environment::setup(max(1, parser("-nthreads").asInt()));

  const string studycase = parser("-study").asString();

#ifdef FLAGELBOW
  test = new Elbow(argc, argv);
#endif
#ifdef FLAGFLAGELLA
  test = new Flagella(argc, argv);
#endif
#ifdef FLAGWALKER
  test = new Walker(argc, argv);
#endif
#ifdef FLAGMUSCULARSNAKE
  test = new MuscularSnake(argc, argv);
#endif
#ifdef FLAGWING
  test = new Wing(argc, argv);
#endif
#ifdef FLAGSPHERICALJOINT
  test = new SphericalJoint(argc, argv);
#endif
#ifdef FLAGHINGEJOINT
  test = new HingeJoint(argc, argv);
#endif
#ifdef FLAGFIXEDJOINT
  test = new FixedJoint(argc, argv);
#endif
#ifdef FLAGPULLINGMUSCLE
  test = new PullingMuscle(argc, argv);
#endif
#ifdef FLAGQUASISTATICTIMOSHENKOBEAM
  test = new QuasistaticTimoshenkoBeam(argc, argv);
#endif
#ifdef FLAGHELICALBUCKLING
  test = new InstabilityHelical(argc, argv);
#endif
#ifdef FLAGSNAKE
  test = new Snake(argc, argv);
#endif

  try {
#ifdef SNAKE_VIZ
    VisualSupport::run(argc, argv);
#else
    test->run();
#endif

    return 0;
  } catch (string &exc) {
    cout << exc << endl;
  } catch (char const *exc) {
    cout << exc << endl;
  }
}
