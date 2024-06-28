#!/bin/bash

SOLVENT_LC="__SOLVENT_LC__"

source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash

if [ "$SOLVENT_LC" == "__SOLVENT_LC__" ]; then
    echo Replace all placeholders!
    exit
fi

mkdir outputs
gmx_mpi editconf -f ../../../parameterize/xff-dash/*/rdk_model.pdb -o outputs/mol-original.gro
cp ../../../parameterize/xff-dash/*/rdk_model.top outputs/mol-original.top
cp ../../../../solvent-models/$SOLVENT_LC/xff-dash/model-correct-density.gro outputs/$SOLVENT_LC.gro
cp ../../../../solvent-models/$SOLVENT_LC/xff-dash/renamed.top outputs/$SOLVENT_LC-original.top
