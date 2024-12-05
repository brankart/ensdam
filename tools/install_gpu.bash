#!/bin/bash
#

module purge
module load nvidia-compilers/24.3
module load cuda/12.2.0  # not needed before, worked with cuda/12.4.1
module load openmpi/4.0.5-cuda

ln -sf macro/make.jean-zay_gpu Makefile.macro

./mkmf -t Makefile.template -p ../lib/libensdam_gpu.a ../src/E* ../src/InterpTools ../src/ObsError ../src/StochTools ../src/TranSpHO/

make all

#make mcmc_ensemble_update.x

#make examples

make install

