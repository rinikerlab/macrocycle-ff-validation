#!/bin/bash

SOLVENT_LC="$1"
SOLUTE="$2"

source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash

if [ -z "$SOLUTE" ]; then
    echo Usage: $0 SOLVENT_LC SOLUTE
fi

sed -i "s/=\"__SOLUTE__/=\"$SOLUTE/;s/=\"__SOLVENT_LC__/=\"$SOLVENT_LC/" 1-get-inputs.sh
sed -i "s/=\"__SOLUTE__/=\"$SOLUTE/;s/=\"__SOLVENT_LC__/=\"$SOLVENT_LC/" 2-run.sh

bash 1-get-inputs.sh
bash 2-run.sh
