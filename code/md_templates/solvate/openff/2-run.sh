#!/bin/bash

SOLUTE="REPLACEME"
SOLVENT_LC="REPLACEME"

if [[ "$SOLUTE" == "REPLACEME" || "$SOLVENT_LC" == "REPLACEME" ]]; then
    echo Replace all placeholders!
    exit 1
fi

GMX="/localhome/fwaibl/programs/gromacs/bin/gmx_mpi"
cd outputs

rm \#*

sed 's/UNK/CL3/;s/x0/CL3/' $SOLVENT_LC-original.gro > $SOLVENT_LC.gro

sed 's/UNK/CL3/;s/x0/CL3/' $SOLVENT_LC-original.top \
    | bash ../process_gmx_top.sh remove defaults \
    | bash ../process_gmx_top.sh remove atomtypes \
    | bash ../process_gmx_top.sh remove system \
    | bash ../process_gmx_top.sh remove molecules \
    > $SOLVENT_LC.top
SOLVENT_atomtypes=$( bash ../process_gmx_top.sh extract atomtypes $SOLVENT_LC-original.top ) 
if [ -z "$SOLVENT_atomtypes" ]; then
    echo "No chloroform atom types found; exiting!"
    exit 1
fi
sed "s/UNK/$SOLUTE/;s/x0/$SOLUTE/" mol-original.top \
    | bash ../process_gmx_top.sh append atomtypes "$SOLVENT_atomtypes" \
    | bash ../process_gmx_top.sh prepend system "#include \"$SOLVENT_LC.top\"" \
    > solvated.top
sed "s/UNK/$SOLUTE/;s/x0/$SOLUTE/" mol-original.gro > mol.gro

$GMX solvate -cp mol.gro -cs $SOLVENT_LC.gro -p solvated.top -o solvated.gro -box 7 7 7

echo
echo
echo Don\'t forget to manually edit solvated.top by adding the atom types from $SOLVENT_LC.top
echo Example replacement mask for VIM: ':%s/\(C1\|Cl1\|Cl2\|Cl3\|H1\)/\1CL3/'
echo Also, remove [defaults] and [atomtypes] from $SOLVENT_LC.top
echo
