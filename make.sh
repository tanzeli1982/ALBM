#! /bin/bash

##-----------------------------------------------------------------------
## Default ALBM Build
##------------
##

module purge

export NETCDF_HOME=/people/tanz151/packages/Pnetcdf_mvapich_2.3.7
export LD_LIBRARY_PATH=$NETCDF_HOME/lib:$LD_LIBRARY_PATH
module load intel/2022.1.0 mvapich2/2.3.7

export FC=ifort

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
