#!/bin/bash

NTHREADS=1

SETTINGS=
SETTINGS+=" -study KNOT"
SETTINGS+=" -nthreads ${NTHREADS}"
SETTINGS+=" -cmaes ctrl.413"
SETTINGS+=" -framesPerUnitTime 1"

OPTIONS=${SETTINGS}

rm -fr ../run
mkdir ../run
cp ../makefiles/snake ../run/executable  
cp ctrl* leviLaunch.sh ../run/
cd ../run

./executable ${OPTIONS}














