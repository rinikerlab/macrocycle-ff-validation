#!/bin/bash

SOLVENT_LC="REPLACEME"

if [ "$SOLVENT_LC" == "REPLACEME" ]; then
   echo Replace all placeholders!
   exit 1
fi

mkdir outputs
cp ../../../parameterize/amber/outputs/converted.amb2gmx/converted_GMX.gro outputs/mol-original.gro
cp ../../../parameterize/amber/outputs/converted.amb2gmx/converted_GMX.top outputs/mol-original.top
cp ../../../../solvent-models/$SOLVENT_LC/amber/outputs/$SOLVENT_LC-correct-density.gro outputs/$SOLVENT_LC-original.gro
cp ../../../../solvent-models/$SOLVENT_LC/amber/outputs/converted.amb2gmx/converted_GMX.top outputs/$SOLVENT_LC-original.top
