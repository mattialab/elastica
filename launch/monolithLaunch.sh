#!/bin/bash

NTHREADS=1

SETTINGS=		
SETTINGS+=" -study SOLENOIDS"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -cmaes ctrl.413"
SETTINGS+=" -ncycles 5"
SETTINGS+=" -framesPerUnitTime 0"

OPTIONS=${SETTINGS}

export OMP_NUM_THREADS=${NTHREADS}
export LD_LIBRARY_PATH=/usr/local/Cellar/tbb/4.3-20140724/lib/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/usr/local/Cellar/tbb/4.3-20140724/lib/:$DYLD_LIBRARY_PATH
echo "Numeber of OPENMP threads: "$OMP_NUM_THREADS

rm -fr ../run
mkdir ../run
cp ../makefiles/snake ../run/executable  
cp ctrl* monolithLaunch.sh ../run/
cd ../run
./executable ${OPTIONS}














