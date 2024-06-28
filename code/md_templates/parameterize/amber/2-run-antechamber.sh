#!/bin/bash

CHARGE="REPLACEME"
NAME="REPLACEME"

if [[ "$CHARGE" == "REPLACEME" || "$NAME" == "REPLACEME" ]]; then
    echo Replace all placeholders!
    exit 1
fi

mkdir outputs
cd outputs || exit
obabel -imol ../rdk_model.mol -xn -opdb -O model.pdb
antechamber -i model.pdb -fi pdb -o model.prep -fo prepi -c bcc -at gaff2 -nc "$CHARGE" -rn "$NAME" || exit
parmchk2 -i model.prep -f prepi -o $NAME.frcmod -s gaff2

#antechamber -i model.pdb -fi pdb -o model.mol2 -fo mol2 -c bcc -at gaff2 -nc "$CHARGE" -rn "$NAME" || exit
#parmchk2 -i model.mol2 -f mol2 -o $NAME.frcmod -s gaff2
