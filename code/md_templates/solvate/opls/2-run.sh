#!/bin/bash

SOLUTE="REPLACEME"
SOLVENT_LC="REPLACEME"

if [[ "$SOLUTE" == "REPLACEME" || "$SOLVENT_LC" == "REPLACEME" ]]; then
    echo "Please replace all placeholders"
    exit 1
fi

GMX="/localhome/fwaibl/programs/gromacs/bin/gmx_mpi"
cd outputs

rm \#*

cp $SOLVENT_LC-original.gro $SOLVENT_LC.gro
cat $SOLVENT_LC-original.itp \
    | bash ../process_gmx_top.sh remove defaults \
    | bash ../process_gmx_top.sh remove atomtypes \
    | bash ../process_gmx_top.sh remove system \
    | bash ../process_gmx_top.sh remove molecules \
    > $SOLVENT_LC.itp
solvent_atomtypes=$( bash ../process_gmx_top.sh extract atomtypes $SOLVENT_LC-original.itp ) 
if [ -z "$solvent_atomtypes" ]; then
    echo "No chloroform atom types found; exiting!"
    exit 1
fi
sed "s/UNK/$SOLUTE/" mol-original.itp \
    | bash ../process_gmx_top.sh append atomtypes "$solvent_atomtypes" \
    > mol.itp

# defaults copied from
# /localhome/fwaibl/programs/gromacs/share/gromacs/top/oplsaa.ff/forcefield.itp
cat << EOF > solvated.top
[ defaults ]
; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ
1		3		yes		0.5	0.5

#include "mol.itp"

#include "$SOLVENT_LC.itp"

[ system ]
mol

[ molecules ]
$SOLUTE   1
EOF

sed "s/UNK/$SOLUTE/" mol-original.gro > mol.gro

$GMX solvate -cp mol.gro -cs $SOLVENT_LC.gro -p solvated.top -o solvated.gro -box 7 7 7

echo
echo
echo Don\'t forget to manually check solvated.top and $SOLVENT_LC.top
echo Example replacement mask to replace atom names in VIM: ':%s/\(C1\|Cl1\|Cl2\|Cl3\|H1\)/\1CL3/'
echo
