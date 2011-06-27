#!/bin/bash

# minimal spinup of model using SeaRISE-Antarctica plus PIK physics
# run preprocess.sh before this

# compare older at page
#    http://www.pism-docs.org/wiki/doku.php?id=worked-searise-antarctica-uaf

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then
  NN="$1"
fi

#inputs generated by preprocess.sh:   pism_Antarctica_5km.nc pism_dT.nc pism_dSL.nc

# grids
VDIMS="-Lz 5000 -Lbz 2000"
THIRTYKMGRID="${VDIMS} -Mx 200 -My 200 -Mz 71 -Mbz 21"
FIFTEENKMGRID="${VDIMS} -Mx 400 -My 400 -Mz 101 -Mbz 21"
TENKMGRID="${VDIMS} -Mx 600 -My 600 -Mz 151 -Mbz 31"

GRID=$THIRTYKMGRID
SKIP="-skip 10"
if [ $# -gt 1 ] ; then
  if [ $2 -eq "1" ] ; then  # if user says "antspin.sh N 1" then
    echo " 15km grid"
    GRID=$FIFTEENKMGRID
    SKIP="-skip 10"
  fi
  if [ $2 -eq "2" ] ; then  # if user says "antspin.sh N 2" then
    echo " 10km grid: LARGE COMPUTATIONAL TIME AND LARGE MEMORY"
    GRID=$TENKMGRID
    SKIP="-skip 20"
  fi
else
  echo " 30km grid"
fi
echo "  grid = '$GRID $SKIP'"

PIKOPTIONS="-pik -eigen_calving 2.0e18 -calving_at_thickness 150.0"  # parameters preliminary

SSA="-ssa_sliding -ssa_method fd -e_ssa 1.0"
SLIDING="-thk_eff -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -topg_to_phi 5.0,20.0,-300.0,700.0,10.0" # preliminary

COUPLER_SIMPLE="-atmosphere constant -surface simple"
COUPLER_FORCING="-atmosphere constant,dTforcing -surface simple -dTforcing pism_dT.nc -ocean constant,dSLforcing -dSLforcing pism_dSL.nc"

#  bootstrapping plus short SIA run for 100 a
mpiexec -n $NN pismr -e 3 -boot_file pism_Antarctica_5km.nc \
  $GRID $SKIP $COUPLER_SIMPLE -ocean_kill \
  -y 100 -o ant_pre100.nc

# FIXME: very reasonable to add in period of -no_mass; see for example
#   pism-dev/examples/searise-greenland/spinup.sh

#  paleo-climate forcing run from -5000 a BPE to 0a BPE
mpiexec -n $NN pismr -e 3 -i ant_pre100.nc \
  $SKIP $COUPLER_FORCING $PIKOPTIONS $SSA $SLIDING \
  -ts_file ts_ant_m50ka.nc -ts_times -5000:1:0 \
  -extra_file ex_ant_m50ka.nc -extra_times -4975:25:-25 -extra_vars bmelt,tauc,tempicethk_basal,Href,csurf,cbase,mask,IcebergMask,diffusivity,thk \
  -ys -5000 -ye 0 -o ant_m50ka.nc -o_size big
