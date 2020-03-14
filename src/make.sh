#! /bin/bash

##-----------------------------------------------------------------------
## Default ALBM Build
##------------
##

export NETCDF_HOME=/usr/local
export LD_LIBRARY_PATH=$NETCDF_HOME/lib:$LD_LIBRARY_PATH

arg=$( echo $1 | tr '[:upper:]' '[:lower:]' )
if [ -z "$arg" ]; then
   gmake OPT=1 -j8
elif [ $arg = 'trace' ]; then
   gmake TRACE=1 -j8
elif [ $arg = 'clean' ]; then
   gmake $arg -j8
elif [ $arg = 'debug' ]; then
   gmake DEBUG=1 TRACE=1 -j8
else
   echo "Wrong Argument: $1!!!"
fi
