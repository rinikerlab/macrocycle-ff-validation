#!/bin/bash

SOLVENT_LC="REPLACEME"

if [[ "$SOLVENT_LC" == "REPLACEME" ]]; then
    echo "Please replace all placeholders"
    exit 1
fi

mkdir outputs
cp ../../../parameterize/opls/mol.gro outputs/mol-original.gro
cp ../../../parameterize/opls/mol.itp outputs/mol-original.itp
cp ../../../../solvent-models/$SOLVENT_LC/opls/correct-density.gro outputs/$SOLVENT_LC-original.gro
cp ../../../../solvent-models/$SOLVENT_LC/opls/renamed.itp outputs/$SOLVENT_LC-original.itp
