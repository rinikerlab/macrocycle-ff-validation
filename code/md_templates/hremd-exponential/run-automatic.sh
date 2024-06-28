#!/bin/bash

SOLVENT_LC="REPLACEME"
SOLUTE="REPLACEME"
FORCEFIELD="REPLACEME"
NSTEPS="25000000"
NREP="12"
TMIN="300"
TMAX="2400"

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

SEDSTRING="
    s/SOLUTE=\"REPLACEME\"/SOLUTE=\"$SOLUTE\"/;
    s/SOLVENT=\"REPLACEME\"/SOLVENT=\"$SOLVENT\"/;
    s/SOLVENT_LC=\"REPLACEME\"/SOLVENT_LC=\"$SOLVENT_LC\"/;
    s/FORCEFIELD=\"REPLACEME\"/FORCEFIELD=\"$FORCEFIELD\"/;
    s/NSTEPS=\"REPLACEME\"/NSTEPS=\"$NSTEPS\"/;
    s/COMPRESSIBILITY=\"REPLACEME\"/COMPRESSIBILITY=\"$COMPRESSIBILITY\"/;
    s/NREP=\"REPLACEME\"/NREP=\"$NREP\"/;
    s/TMIN=\"REPLACEME\"/TMIN=\"$TMIN\"/;
    s/TMAX=\"REPLACEME\"/TMAX=\"$TMAX\"/;
"
sed -i "$SEDSTRING" 1-get-inputs.sh
sed -i "$SEDSTRING" 2-make-single-topology.sh
sed -i "$SEDSTRING" 3-plumed-prepare-hremd.sh

sed -i "s/__NJOBS__/$NREP/" 4-run-euler.sh
sed -i "s/__NJOBS__/$NREP/" 4-run-local.sh

bash 1-get-inputs.sh
bash 2-make-single-topology.sh
bash 3-plumed-prepare-hremd.sh
