#!/bin/bash

# Make sure to have access to: cmake, fortran compiler, python
# and then build the python module as follows.
# To use the MPI version, you also need to install: mpi4py

# To know the default compilation settings: cmake -LAH
# Edit the flags to modify the default

# General options to pass to cmake
flags=""

# Compilation options to pass to the fortran compiler
flags_fortran=""

# Generate the module with or without mpi: pyensdam or pyensdam_mpi
flag_mpi="ON"   # ON or OFF
flag_mpi="OFF"   # ON or OFF

# Uncomment to explicitly specify the python3 exectuable
#flags="$flags -DPython3_EXECUTABLE=/gpfslocalsup/pub/anaconda-py3/2021.05/envs/python-3.10.4/bin/python3"

# Uncomment to explicitly specify the Fortran compiler
flags="$flags -DCMAKE_Fortran_COMPILER=ifort"
if [ ${flag_mpi} = "ON" ]  ;then
  flags="$flags -DCMAKE_Fortran_COMPILER=mpiifort"
fi

# Pass the MPI flag to cmake
flags="$flags -DUSE_MPI=${flag_mpi}"

# Modify project name if mpi
project_name="pyensdam"
if [ ${flag_mpi} = "ON" ]  ;then
  project_tail="_mpi"
  project_name="pyensdam_mpi"
  flags_fortran="-DMPI -DMPI_MODULE $flags_fortran"
fi

flags="$flags -DCMAKE_Fortran_FLAGS=$flags_fortran"

#If needed: remove all files and directories generated by cmake
#rm -fr build build_mpi
#find . -type d -name "CMakeFiles" -exec rm -rf {} +
#find . -name "cmake_install.cmake" -exec rm -f {} +

cmake $flags -B build${project_tail} -S .

cmake --build build${project_tail} --target wheel && pip install --user build${project_tail}/python/dist/${project_name}-0.1.1-py3-none-any.whl --force-reinstall
