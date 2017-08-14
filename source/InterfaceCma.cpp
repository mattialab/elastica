/*
 * InterfaceCma.cpp
 *
 *  Created on: Jun 18, 2014
 *      Author: mgazzola
 */

#include "InterfaceCma.h"

InterfaceCma::InterfaceCma()
{
	parameters.clear();

	rng.seed(time(NULL));
}

InterfaceCma::~InterfaceCma()
{
}

vector<double> InterfaceCma::parse(string input)
{
	parameters.clear();

	FILE * ppFile;
	ppFile = fopen(input.c_str(),"r");
	if(ppFile==NULL){ std::cout << "could not open ctrl file " << input << "!" << std::endl;}
	char workDir[10000];
	char jobName[10000];
	unsigned int dim = 0;
	unsigned int ID = 0;
	double dummy = 0.0;
	fscanf(ppFile,"%s",workDir);
	fscanf(ppFile,"%s",jobName);
	fscanf(ppFile,"%d",&dim);
	for(unsigned int i=0; i<dim-1; i++)
	{
		fscanf(ppFile,"%lf",&dummy);
		parameters.push_back(dummy);
	}
	fscanf(ppFile,"%d",&ID);
	fclose(ppFile);

	return parameters;
}

void InterfaceCma::printParameters()
{
	cout << "Cmaes parameters: " << endl;
	for(unsigned int i=0; i<parameters.size(); i++)
		cout << parameters[i] << endl;
}

void InterfaceCma::dumpRandomFitness()
{
	double fitness = normal_dist(rng);

	FILE * fitnessFile = fopen("fitness","w");
	assert(fitnessFile!=NULL);
	fprintf(fitnessFile,"%e", fitness);
	fclose(fitnessFile);
}

void InterfaceCma::dumpFitness(double fitness)
{
	FILE * fitnessFile = fopen("fitness","w");
	assert(fitnessFile!=NULL);
	fprintf(fitnessFile,"%e", fitness);
	fclose(fitnessFile);
}


