#!/bin/bash

mkdir outputs/

FORCEFIELD="REPLACEME"
SOLUTE="REPLACEME"
SOLVENT_LC="REPLACEME"
NSTEPS="25000000"

# 4.5e-5 for water, 1.0e-4 for chloroform
# 5.5e-5 for DMSO, see https://doi.org/10.1023/A:1022985208144
if [ "$SOLVENT_LC" == "water" ]; then
  COMPRESSIBILITY="4.5e-5"
  SOLVENT="SOL"
elif [ "$SOLVENT_LC" == "chcl3" ]; then
  COMPRESSIBILITY="1.0e-4"
  SOLVENT="CL3"
elif [ "$SOLVENT_LC" == "dmso" ]; then
  COMPRESSIBILITY="5.5e-5"
  SOLVENT="DSO"
else
  echo Unknown solvent \(SOLVENT_LC\)\; Exiting...
  exit 1
fi

if [[ "$FORCEFIELD" == "REPLACEME" ]] || [[ "$SOLUTE" == "REPLACEME" ]] || [[ "$SOLVENT" == "REPLACEME" ]] || [[ "$NSTEPS" == "REPLACEME" ]] || [[ "$COMPRESSIBILITY" == "REPLACEME" ]]; then
    echo Please replace all placeholders!
    exit 1
fi

# Get input files
sed "s/__SOLUTE__/$SOLUTE/;s/__SOLVENT__/$SOLVENT/;s/__NSTEPS__/$NSTEPS/;s/__COMPRESSIBILITY__/$COMPRESSIBILITY/" \
    hremd-md-template.mdp > outputs/sampling.mdp

cp ../../equilibrate/$FORCEFIELD/outputs/npt.gro outputs/
cp ../../equilibrate/$FORCEFIELD/outputs/topol.top outputs/

# for Chloroform; add extra forcefield files accordingly
cp ../../equilibrate/$FORCEFIELD/outputs/$SOLVENT_LC.top outputs/
