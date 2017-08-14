/*
 * InterfaceCma.h
 *
 *  Created on: Jun 18, 2014
 *      Author: mgazzola
 */

#ifndef INTERFACECMA_H
#define INTERFACECMA_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <assert.h>

using namespace std;

class InterfaceCma
{
public:
	InterfaceCma();
	virtual ~InterfaceCma();

	vector<double> parse(string input);
	void printParameters();
	void dumpRandomFitness();
	void dumpFitness(double fitness);

protected:
	std::vector<double> parameters;
	std::mt19937 rng;
	std::normal_distribution<double> normal_dist;
};

#endif
