#! /bin/bash

##-----------------------------------------------------------------------
## Default ALBM Build
##------------
##

module unload netcdf

mach=$( echo $HOSTNAME | cut -f 1 -d . -d 0 )

if [ $mach = 'constance' ]; then

   export NETCDF_HOME=/people/tanz151/packages/Pnetcdf_openmpi_1.8.3
   export LD_LIBRARY_PATH=$NETCDF_HOME/lib:$LD_LIBRARY_PATH
   module load intel/15.0.1 openmpi/1.8.3

elif [ $mach = 'compy' ]; then

   export NETCDF_HOME=/people/tanz151/packages/Pnetcdf_mvapich_2.3.1
   export LD_LIBRARY_PATH=$NETCDF_HOME/lib:$LD_LIBRARY_PATH
   module load intel/19.0.3 mvapich2/2.3.1

fi

arg=$( echo $1 | tr '[:upper:]' '[:lower:]' )
if [ -z "$arg" ]; then
   if [ $mach = 'constance' ]; then
      gmake OPT=0 -j8
   elif [ $mach = 'compy' ]; then
      gmake OPT=1 -j8
   fi
elif [ $arg = 'clean' ]; then
   gmake $arg -j8
elif [ $arg = 'debug' ]; then
   gmake DEBUG=1 TRACE=1 -j8
else
   echo "Wrong Argument: $1!!!"
fi
