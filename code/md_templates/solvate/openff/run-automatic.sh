#!/bin/bash

SOLVENT_LC="REPLACEME"
SOLUTE="REPLACEME"

sed -i 's/SOLUTE="REPLACEME"/SOLUTE="'$SOLUTE'"/;s/SOLVENT_LC="REPLACEME"/SOLVENT_LC="'$SOLVENT_LC'"/' 1-get-inputs.sh
sed -i 's/SOLUTE="REPLACEME"/SOLUTE="'$SOLUTE'"/;s/SOLVENT_LC="REPLACEME"/SOLVENT_LC="'$SOLVENT_LC'"/' 2-run.sh

bash 1-get-inputs.sh
bash 2-run.sh
