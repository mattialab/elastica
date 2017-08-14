/*
 *  MRAGProfiler.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#include "MRAGEnvironment.h"
#include "MRAGProfiler.h"

#ifdef _MRAG_TBB

#include "tbb/tick_count.h"
using namespace tbb;

void MRAG::ProfileAgent::_getTime(tick_count& time)
{
	time = tick_count::now();
}

float MRAG::ProfileAgent::_getElapsedTime(const tick_count& tS, const tick_count& tE)
{
	return (tE - tS).seconds();
}
	
#else
#include <time.h>
void MRAG::ProfileAgent::_getTime(clock_t& time)
{
	time = clock();
}

float MRAG::ProfileAgent::_getElapsedTime(const clock_t& tS, const clock_t& tE)
{
	return (tE - tS)/(double)CLOCKS_PER_SEC;
}

#endif
	
