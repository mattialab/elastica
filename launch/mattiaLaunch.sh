#!/bin/bash

SETTINGS=
SETTINGS+=" -cmaes ctrl.413"
SETTINGS+=" -ncycles 5"

OPTIONS=${SETTINGS}

rm -fr ../run
mkdir ../run
cp ../makefiles/snake ../run/executable  
cp ctrl* mattiaLaunch.sh ../run/
cd ../run
./executable ${OPTIONS}














