#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 7: comparing regridding via -boot_file and -regrid_file."
# The list of files to delete when done:
files="foo.nc bar.nc baz.nc"

OPTS="-y 0 -surface constant"

set -e -x

# Create the file to regrid and bootstrap from:
$PISM_PATH/pismv -test G -Lx 4000 -Ly 4000 -Lz 4000 -Mx 41 -My 51 -Mz 31 -y 0 -o foo.nc 

# Bootstrap from this file:
$PISM_PATH/pismr -boot_file foo.nc -Lx 2000 -Ly 2000 -Lz 4000 -Mx 31 -My 41 -Mz 51 $OPTS -o bar.nc 

# Overwrite topg using -regrig_file and save the result to baz.nc:
$PISM_PATH/pismr -i bar.nc -regrid_file foo.nc -regrid_vars topg $OPTS -o baz.nc 

set +e

# Compare:
$PISM_PATH/nccmp.py -v topg bar.nc baz.nc
if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0
