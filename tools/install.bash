#!/bin/bash
#

./mkmf -t Makefile.template -p ../lib/libensdam.a ../src/E* ../src/InterpTools ../src/ObsError ../src/StochTools ../src/TranSpHO/

make all

#make mcmc_ensemble_update.x

#make examples

make install

