/*
 * ValidationMitchellBuckling.cpp
 *
 *  Created on: Aug 4, 2014
 *      Author: mgazzola
 */

#include "MitchellBuckling.h"

MitchellBuckling::MitchellBuckling(const int argc, const char ** argv)
{
}

MitchellBuckling::~MitchellBuckling()
{
}

bool MitchellBuckling::_testMitchellBuckling(const REAL alpha, const REAL beta, const REAL twist, const REAL L, const REAL timeSimulation, const REAL translationalLimit, REAL &translationalEnergy)
{
	// Discretization parameters
	const int nPerUnit = 50;
	const REAL T = timeSimulation;
	const REAL dt = 1e-5;
	const int numSteps = ceil(T/dt);

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 0;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const int n = L*nPerUnit;
	const REAL totalMass = 1.0;
	const REAL deltaL = L/(REAL)n;
	const REAL r = 0.01*L;
	const REAL nu = 0.0;
	const REAL totTwist = 0.0;
	const REAL A = M_PI*r*r;
	const REAL density = totalMass/(A*L);

	// Second moment of area for disk cross section
	const REAL I1 = A*A/(4.0*M_PI);
	const REAL I2 = I1;
	const REAL I3 = 2.0*I1;
	const Matrix3 I = Matrix3(	 I1,	0.0,	0.0,
								0.0,	 I2,	0.0,
								0.0,	0.0,	 I3);

	// Mass inertia matrix for disk cross section
	const Matrix3 J = density*deltaL*I;

	// Bending matrix
    Matrix3 B = Matrix3(	alpha, 0, 0,
							0, alpha, 0,
							0, 0, beta);

    // Shear matrix
    Matrix3 S = 100000 * Matrix3(	1, 0, 0,
    								0, 1, 0,
				 	 	 	 	 	0, 0, 1);

    // Initialize rod and pack it into a vector of pointers to rod.
    // This circle rod lies in the x-y plane
    const bool useSelfContact = false;
    Rod rod = RodInitialConfigurations::circleRod(n, totalMass, r, J, B, S, L, totTwist, nu, useSelfContact);
    Vrodptr rodPtrs;
    rodPtrs.push_back(&rod);

    // Initialize boundary conditions (cut the rod, twist, relax)
    const REAL relaxFraction = 0.5;
    const unsigned int twistSteps = ceil(numSteps*(1.0 - relaxFraction));
    const unsigned int relaxSteps = numSteps - twistSteps;
    VM3 QStart;
    VM3 QEnd;

    // ...cut
    QStart.push_back(rod.Q[0]);
    QEnd.push_back(rod.Q[n-1]);

    // ...twist
    for (unsigned int i=1; i<twistSteps; i++)
    {
		// Compute angle of rotation: in this case, this is a counterclockwise rotation
		const REAL angleRotation = twist/(double)twistSteps;

		// Use cross product to compute axis of rotation in the lab frame of reference (LabFOR)
		const Vector3 axisRotationInLabFOR = (rod.x[1] - rod.x[0]).unitize();

		// Express the axis of rotation in LabFOR in the material frame of reference:
		// Formula: k_material = Q * (k_lab - x_material), where k_lab is the vector to be transformed,
		// x_material is its position in the material frame of reference
		const Vector3 axisRotationInMaterialFOR = (QStart[QStart.size()-1] * axisRotationInLabFOR).unitize();

		// Apply rotation using Rodrigues' rotation formula: Q_rotated = exp(angle * k) * Q.
		// This formula is based on the definition of skew map S.
		// In general it applies a CLOCKWISE rotation, but here Andrew for no reason decided to
		// implement S_andrew = -S_rest_of_the_world and so in this code the formula applies a
		// COUNTERCLOCKWISE rotation by default, hence no minus in front of 'angleRotation' below
		QStart.push_back( exp(angleRotation * axisRotationInMaterialFOR) * QStart[QStart.size()-1] );

    	// The second end is not twisted
    	QEnd.push_back(QEnd[QEnd.size()-1]);
    }

    // ...relax
    for (unsigned int i=0; i<relaxSteps; i++)
    {
    	QStart.push_back(QStart[relaxSteps-1]);
    	QEnd.push_back(QEnd[relaxSteps-1]);
    }

    // Pack boundary consitions
    QuasiStaticEndpointBC twistedEndsBC = QuasiStaticEndpointBC(VV3(numSteps, rod.x[0]), VV3(numSteps, rod.x[n]), QStart, QEnd);
    Vbcptr boundaryConditionsPtrs;
    boundaryConditionsPtrs.push_back(&twistedEndsBC);

    // Radom forces (in this case it doesnt really matter, since the transition is so abrupt)
    RandomForces randomForces = RandomForces(1e-3);

   	// Pack all forces together
   	MultipleForces multipleForces;
   	multipleForces.add(&randomForces);
   	Vefptr externalForcesPtrs = multipleForces.get();

   	// Empty interaction forces (no substrate in this case)
   	Vinterptr substrateInteractionsPtrs;

   	// Set up integrator (define integration order)
   	PolymerIntegrator * integrator = new PositionVerlet2nd(rodPtrs, externalForcesPtrs, boundaryConditionsPtrs, substrateInteractionsPtrs);

   	// Allocate polymer simulator
   	Polymer poly = Polymer(integrator);

   	// Simulate
   	bool goodRun = poly.simulate(T, dt, diagPerUnitTime, povrayPerUnitTime, "diagnostics");

   	// Throw exception if something went wrong
   	if (!goodRun)
   		throw "not good run in Mitchell buckling, what is going on?";

   	// Evaluate whether solution is stable or not
   	translationalEnergy = poly.getTotalTranslationalEnergy();
   	const bool stable = (translationalEnergy>translationalLimit)?false:true;

   	return stable;
}

void MitchellBuckling::_sweepMitchellBuckling(const REAL alpha, const REAL L, string outfileName)
{
	const int numChecks = 12;
	const double minTwist = 0.0;
	const double maxTwist = 25.0;
	const double minBetaAlphaRatio = 0.5;
	const double maxBetaAlphaRatio = 2.0;
	const int nCoarseTwist = 10;
	const int nCoarseBetaAlphaRatio = 10;
	const REAL timeSimulation = 2.0;
	const REAL translationalLimit = 1e-4; // In this case it doesnt really matter, since the transition is so abrupt

	// Build coarse + adaptive grid
	vector< pair<double, double> > adaptiveGrid;

	// Build coarse grid
	{
		REAL betaAlphaRatio = minBetaAlphaRatio;
		for(unsigned int i=0; i<nCoarseBetaAlphaRatio; i++)
		{
			REAL twist = minTwist;

			for(unsigned int j=0; j<nCoarseTwist; j++)
			{
				pair<double, double> coord;
				coord.first = betaAlphaRatio;
				coord.second = twist;
				adaptiveGrid.push_back(coord);

				twist += (maxTwist-minTwist)/double(nCoarseTwist-1);
			}

			betaAlphaRatio += (maxBetaAlphaRatio-minBetaAlphaRatio)/double(nCoarseBetaAlphaRatio-1);
		}
	}

	// Compute probabilites
	for(unsigned int i=0; i<adaptiveGrid.size(); i++)
	{
		ofstream forceLenghtFile;
		forceLenghtFile.open (outfileName.c_str(), (i==0)?std::ofstream::out:std::ofstream::app);

		const REAL betaAlphaRatio = adaptiveGrid[i].first;
		const REAL twist = adaptiveGrid[i].second;

		const double dy = 0.5;
		const double analTwist = 2.0*M_PI*sqrt(3.0)/betaAlphaRatio;
		const double ym = analTwist-dy;
		const double yp = analTwist+dy;

		if( twist>=ym && twist<=yp )
		{
			vector<bool> stableCases;
			for (unsigned int check=0; check<numChecks; check++)
				stableCases.push_back(false);

#ifndef SNAKE_VIZ
#pragma omp parallel for
#endif
			for (unsigned int check=0; check<numChecks; check++)
			{
				const REAL beta = betaAlphaRatio * alpha;
				REAL translationalEnergy = 0.0;
				const bool stable = _testMitchellBuckling(alpha, beta, twist, L, timeSimulation, translationalLimit, translationalEnergy);
				string stableStr = stable?"yes":"no";
				printf("twist=%f, betaAlphaRatio=%f --> translational energy=%e, stable=%s\n", twist, betaAlphaRatio, translationalEnergy, stableStr.c_str());
				stableCases[check] = stable;
			}

			int nBuckles = 0;
			for (unsigned int j=0; j<stableCases.size(); j++)
			{
				nBuckles += stableCases[j]?0:1;
				printf("flags: %d ",stableCases[j]?0:1);
			}
			printf("\n");

			const double probability = (double)nBuckles/(double)numChecks;
			forceLenghtFile << twist << "  " << betaAlphaRatio << " " << probability << endl;
			printf("twist=%f, betaAlphaRatio=%f --> P=%e\n", twist, betaAlphaRatio, probability);
		}
		else if (twist > yp)
		{
			const double probability = 1.0;
			forceLenghtFile << twist << "  " << betaAlphaRatio << " " << probability << endl;
			printf("twist=%f, betaAlphaRatio=%f --> P=%e\n", twist, betaAlphaRatio, probability);
		}
		else if (twist < ym)
		{
			const double probability = 0.0;
			forceLenghtFile << twist << "  " << betaAlphaRatio << " " << probability << endl;
			printf("twist=%f, betaAlphaRatio=%f --> P=%e\n", twist, betaAlphaRatio, probability);
		}
		else
		{
			printf("wtf something wrong!");
			exit(0);
		}

		forceLenghtFile.close();
	}
}

void MitchellBuckling::run()
{
	cout << "Mitchell buckilng validation test" << endl << endl;

	/*
	{
		// One instance example for numerical checks and povray images
		const REAL alpha = 1.0;
		const REAL L = 1.0;
		const REAL beta = 2.0;
		const REAL twist = 10.0;
		const REAL timeSimulation = 0.5;
		const REAL translationalLimit = 0.1;
		REAL translationalEnergy = 0.0;
		_testMitchellBuckling(alpha, beta, twist, L, timeSimulation, translationalLimit, translationalEnergy);
		printf("energy=%e\n\n", translationalEnergy);
	}
	*/

	{
		// Phase space
		const REAL alpha = 1.0;
		const REAL L = 1.0;
		_sweepMitchellBuckling(alpha, L, "mitchell_Alpha1_L1.txt");
	}

	cout << "Mitchell buckling study completed!" << endl << endl;

	exit(0);
}

