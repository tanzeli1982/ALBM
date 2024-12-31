#! /bin/bash

##-----------------------------------------------------------------------
## Default ALBM Build
##------------
##

module purge 

export NETCDF_HOME=/qfs/people/tanz151/packages/Pnetcdf-1.14.0_gnu
export LD_LIBRARY_PATH=$NETCDF_HOME/lib:$LD_LIBRARY_PATH
module load gcc/4.8.5 mvapich2/2.3.1

export FC=gfortran

arg=$( echo $1 | tr '[:upper:]' '[:lower:]' )
if [ -z "$arg" ]; then
   gmake -j8
elif [ $arg = 'clean' ]; then
   gmake $arg -j8
elif [ $arg = 'debug' ]; then
   gmake DEBUG=1 -j8
else
   echo "Wrong Argument: $1!!!"
fi
