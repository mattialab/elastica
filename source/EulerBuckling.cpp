/*
 * ValidationEulerBuckling.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: mgazzola
 */

#include "EulerBuckling.h"

#include <omp.h>

EulerBuckling::EulerBuckling(const int argc, const char ** argv)
{
}

EulerBuckling::~EulerBuckling()
{
}

bool EulerBuckling::_testEulerBuckling(const REAL force, const REAL lengthRod, const REAL alpha, const REAL timeSimulation, const REAL bendLimit, REAL &bendingEnergy)
{
	typedef std::vector<RodBC*> Vbcptr;
	typedef std::vector<ExternalForces*> Vefptr;
	typedef std::vector<Rod*> Vrodptr;
	typedef std::vector<Interaction*> Vinterptr;

	// Discretization paramters
	const int nPerUnit = 100;
	const REAL T = timeSimulation;
	const REAL dt = 1e-5;

	// Dumping frequencies (number of frames/dumps per unit time)
	const unsigned int diagPerUnitTime = 0;
	const unsigned int povrayPerUnitTime = 0;

	// Physical parameters
	const REAL L = lengthRod;
	const int n = L*nPerUnit;
	const REAL totalMass = 1.0;
	const REAL m = totalMass/(double)n;
	const REAL r = 0.025*lengthRod;
	const REAL nu = 0.0;
	const REAL totTwist = 0.0;
	const REAL F = force;
	const Vector3 originRod = Vector3(0.0,0.0,0.0);
	const Vector3 directionRod = Vector3(0.0,0.0,1.0);
	const Vector3 normalRod = Vector3(1.0,0.0,0.0);

	// Inertia matrix (for the discretization cylinders)
	const REAL orthoI = m*r*r/4.0;
	const REAL axialI = m*r*r/2.0;
	Matrix3 I = Matrix3(	orthoI, 0, 0,
							0, orthoI, 0,
							0, 0, axialI);

	// Bending matrix
	Matrix3 B = alpha * Matrix3(	1, 0, 0,
									0, 1, 0,
									0, 0, 2.0/3.0);
	// Shear matrix
	Matrix3 S = 100000 * Matrix3(	1, 0, 0,
									0, 1, 0,
									0, 0, 1);

	// Initialize rod and pack it into a vector of pointers to rods
	const bool useSelfContact = false;
	Rod rod = RodInitialConfigurations::compressedRod(n, totalMass, r, I, B, S, L, F, totTwist, nu, originRod, directionRod, normalRod, useSelfContact);
	Vrodptr rodPtrs;
	rodPtrs.push_back(&rod);

	// Initialize rod boundary conditions and pack it into a vector of pointers to boundary conditions
	RestrictEndBC constrainedEndsBC = RestrictEndBC();
	Vbcptr boundaryConditionsPtrs;
	boundaryConditionsPtrs.push_back(&constrainedEndsBC);

	// Force applied at the extrema of the rod
	EndpointForces endpointsForce = EndpointForces(Vector3(), Vector3(0,0,-F));

	// Random forces
	RandomForces randomForces = RandomForces(0.01);

	// Pack all forces together
	MultipleForces multipleForces;
	multipleForces.add(&endpointsForce);
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
		throw "not good run in Euler buckling, what is going on?";

	// Evaluate whether solution is stable or not
	bendingEnergy = poly.getTotalBendingEnergy();
	const bool stable = (bendingEnergy>bendLimit)?false:true;

	return stable;
}

void EulerBuckling::_sweepEulerBucklingForceLength()
{
	const REAL alpha = 1.0;
	const int numChecks = 12;
	const double minL = 0.5;
	const double maxL = 2.0;
	const double minF = 0.0;
	const double maxF = 20.0;
	const int nCoarseL = 10;
	const int nCoarseF = 10;
	const REAL T = 10.0;
	const REAL bendLimit = 1e-3;

	// Build coarse + adaptive grid
	vector< pair<double, double> > adaptiveGrid;

	// Build coarse grid
	{
		REAL L = minL;
		for(unsigned int i=0; i<nCoarseL; i++)
		{
			REAL F = minF;

			for(unsigned int j=0; j<nCoarseF; j++)
			{
				pair<double, double> coord;
				coord.first = L;
				coord.second = F;
				adaptiveGrid.push_back(coord);

				F += (maxF-minF)/double(nCoarseF-1);
			}

			L += (maxL-minL)/double(nCoarseL-1);
		}
	}

	for(unsigned int i=0; i<adaptiveGrid.size(); i++)
	{
		// Compute probabilites
		ofstream forceLenghtFile;
		forceLenghtFile.open ("forceLength.txt", (i==0)?std::ofstream::out:std::ofstream::app);

		const REAL L = adaptiveGrid[i].first;
		const REAL F = adaptiveGrid[i].second;

		const double dy = 0.5;
		const double analF = M_PI*M_PI*alpha/(L*L);
		const double ym = analF-dy;
		const double yp = analF+dy;

		if( F>=ym && F<=yp )
		{
			vector<bool> stableCases;
			for (unsigned int check=0; check<numChecks; check++)
				stableCases.push_back(false);

#pragma omp parallel for
			for (unsigned int check=0; check<numChecks; check++)
			{
				REAL bendingEnergy = 0.0;
				const bool stable = _testEulerBuckling(F, L, alpha, T, bendLimit, bendingEnergy);
				string stableStr = stable?"yes":"no";
				printf("F=%f, L=%f --> bendingEnergy=%e, stable=%s\n", F, L, bendingEnergy, stableStr.c_str());
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
			forceLenghtFile << F << "  " << L << " " << probability << endl;
			printf("F=%f, L=%f --> P=%e\n", F, L, probability);
		}
		else if (F > yp)
		{
			const double probability = 1.0;
			forceLenghtFile << F << "  " << L << " " << probability << endl;
			printf("F=%f, L=%f --> P=%e\n", F, L, probability);
		}
		else if (F < ym)
		{
			const double probability = 0.0;
			forceLenghtFile << F << "  " << L << " " << probability << endl;
			printf("F=%f, L=%f --> P=%e\n", F, L, probability);
		}

		forceLenghtFile.close();
	}
}

void EulerBuckling::_sweepEulerBucklingForceAlpha()
{
	const REAL L = 1.0;
	const int numChecks = 12;
	const double minAlpha = 0.5;
	const double maxAlpha = 2.0;
	const double minF = 0.0;
	const double maxF = 20.0;
	const int nCoarseAlpha = 10;
	const int nCoarseF = 10;
	const REAL T = 10.0;
	const REAL bendLimit = 1e-3;

	// Build coarse + adaptive grid
	vector< pair<double, double> > adaptiveGrid;

	// Build coarse grid
	{
		REAL alpha = minAlpha;
		for(unsigned int i=0; i<nCoarseAlpha; i++)
		{
			REAL F = minF;

			for(unsigned int j=0; j<nCoarseF; j++)
			{
				pair<double, double> coord;
				coord.first = alpha;
				coord.second = F;
				adaptiveGrid.push_back(coord);

				F += (maxF-minF)/double(nCoarseF-1);
			}

			alpha += (maxAlpha-minAlpha)/double(nCoarseAlpha-1);
		}
	}

	for(unsigned int i=0; i<adaptiveGrid.size(); i++)
	{
		// Compute probabilites
		ofstream forceAlphaFile;
		forceAlphaFile.open ("forceAlpha.txt", (i==0)?std::ofstream::out:std::ofstream::app);

		const REAL alpha = adaptiveGrid[i].first;
		const REAL F = adaptiveGrid[i].second;

		const double dy = 0.5;
		const double analF = M_PI*M_PI*alpha/(L*L);
		const double ym = analF-dy;
		const double yp = analF+dy;

		if( F>=ym && F<=yp )
		{
			vector<bool> stableCases;
			for (unsigned int check=0; check<numChecks; check++)
				stableCases.push_back(false);

#pragma omp parallel for
			for (unsigned int check=0; check<numChecks; check++)
			{
				REAL bendingEnergy = 0.0;
				const bool stable = _testEulerBuckling(F, L, alpha, T, bendLimit, bendingEnergy);
				string stableStr = stable?"yes":"no";
				printf("F=%f, alpha=%f --> bendingEnergy=%e, stable=%s\n", F, alpha, bendingEnergy, stableStr.c_str());
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
			forceAlphaFile << F << "  " << alpha << " " << probability << endl;
			printf("F=%f, alpha=%f --> P=%e\n", F, alpha, probability);
		}
		else if (F > yp)
		{
			const double probability = 1.0;
			forceAlphaFile << F << "  " << alpha << " " << probability << endl;
			printf("F=%f, alpha=%f --> P=%e\n", F, alpha, probability);
		}
		else if (F < ym)
		{
			const double probability = 0.0;
			forceAlphaFile << F << "  " << alpha << " " << probability << endl;
			printf("F=%f, alpha=%f --> P=%e\n", F, alpha, probability);
		}

		forceAlphaFile.close();
	}
}

void EulerBuckling::run()
{
	cout << "Euler buckilng validation test" << endl << endl;

	{
		// Single simulation
		REAL bendingEnergy = 0.0;
		const REAL L = 1.0;
		const REAL F = 10.0;
		const REAL alpha = 1.0;
		const REAL T = 10.0;
		const REAL bendLimit = 1e-3;
		_testEulerBuckling(F, L, alpha, T, bendLimit, bendingEnergy);
	}

	// Swipe loop for phase-space validation
	//_sweepEulerBucklingForceAlpha();
	//_sweepEulerBucklingForceLength();

	cout << "Euler buckling study completed!" << endl << endl;

	exit(0);
}

