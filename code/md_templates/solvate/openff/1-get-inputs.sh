#!/bin/bash

SOLVENT_LC="REPLACEME"

if [ "$SOLVENT_LC" == "REPLACEME" ]; then
    echo Replace all placeholders!
    exit
fi

mkdir outputs
cp ../../../parameterize/openff/outputs/model.gro outputs/mol-original.gro
cp ../../../parameterize/openff/outputs/model.top outputs/mol-original.top
cp ../../../../solvent-models/$SOLVENT_LC/openff/outputs/model-correct-density.gro outputs/$SOLVENT_LC-original.gro
cp ../../../../solvent-models/$SOLVENT_LC/openff/outputs/renamed.top outputs/$SOLVENT_LC-original.top
