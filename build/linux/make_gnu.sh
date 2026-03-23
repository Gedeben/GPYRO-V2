#!/bin/bash


# Gpyro uses several environment variables for compilation. Since the default values
# specified on lines 15 and 20 below probably do not work on your system, add lines
# similar to the following to your ~/.bashrc file:
#

# export GPYRO_FCOMPL_GNU=gfortran
# export GPYRO_FCOMPL_MPI_GNU=/usr/local/openmpi-4.0.3_gnu/bin/mpifort
# export GPYRO_BINARY_DIRECTORY=/usr/local/bin

# Check if environment variables are set and export if not:
if [ -z "$GPYRO_FCOMPL_GNU" ]; then
   export GPYRO_FCOMPL_GNU=gfortran
fi


if [ -z "$GPYRO_FCOMPL_MPI_GNU" ]; then
   export GPYRO_FCOMPL_MPI_GNU=mpifort
fi


mkdir gpyro 2> /dev/null
cd gpyro

rm -f *.o *.mod
make -j6 -f ../../Makefile_gpyro  gnu_linux
rm -f *.o *.mod

rm -f *.o *.mod
make -j6 -f ../../Makefile_gpyro  gnu_linux_debug
rm -f *.o *.mod

cd ..
mkdir gpyro_propest 2> /dev/null
cd gpyro_propest

rm -f *.o *.mod
make -j6 -f ../../Makefile_gpyro_propest  gnu_linux_mpi
rm -f *.o *.mod

rm -f *.o *.mod
make -j6 -f ../../Makefile_gpyro_propest  gnu_linux_mpi_debug
rm -f *.o *.mod
exit 0
