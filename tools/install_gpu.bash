#!/bin/bash
#

cluster='jean-zay'
cluster='adastra'

if [ $cluster = 'jean-zay' ] ; then

  module purge
  module load nvidia-compilers/24.3
  #module load netcdf-fortran/4.5.3-mpi-cuda
  module load cuda/12.2.0  # not needed before, worked with cuda/12.4.1
  module load openmpi/4.0.5-cuda

  ln -sf macro/make.jean-zay_gpu Makefile.macro

elif [ $cluster = 'adastra' ] ; then

  module purge
  module load cpe/24.07
  module load craype-x86-trento
  module load craype-accel-amd-gfx90a
  module load PrgEnv-cray
  module load rocm
  module load cray-mpich

  ln -sf macro/make.adastra_gpu Makefile.macro

else

  echo "Bad cluster name"
  exit

fi

./mkmf -t Makefile.template -p ../lib/libensdam_gpu.a ../src/E* ../src/InterpTools ../src/ObsError ../src/StochTools ../src/TranSpHO/

make all

#make mcmc_ensemble_update.x

#make examples

make install

