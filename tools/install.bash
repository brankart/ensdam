#!/bin/bash
#

./mkmf -t Makefile.template -p ../lib/libensdam.a ../src/*

make all

#make mcmc_ensemble_update.x

#make examples

#make install

