#!/bin/bash

SOLUTE="__SOLUTE__"
SOLVENT_LC="__SOLVENT_LC__"

if [[ "$SOLUTE" == "__SOLUTE__" || "$SOLVENT_LC" == "__SOLVENT_LC__" ]]; then
    echo Replace all placeholders!
    exit 1
fi

GMX="/localhome/fwaibl/programs/gromacs/bin/gmx_mpi"
cd outputs

rm \#*

cat $SOLVENT_LC-original.top \
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

sed "s/lig/$SOLUTE/" mol-original.top \
    | bash ../process_gmx_top.sh append atomtypes "$SOLVENT_atomtypes" \
    | bash ../process_gmx_top.sh prepend system "#include \"$SOLVENT_LC.top\"" \
    > solvated_empty_line.top
sed "s/lig/$SOLUTE/" mol-original.gro > mol.gro

$GMX solvate -cp mol.gro -cs $SOLVENT_LC.gro -p solvated_empty_line.top -o solvated.gro -box 7 7 7
bash ../process_gmx_top.sh remove_newlines molecules solvated_empty_line.top > solvated.top

echo
echo
echo Don\'t forget to manually edit solvated.top by adding the atom types from $SOLVENT_LC.top
echo Example replacement mask for VIM: ':%s/\(C1\|Cl1\|Cl2\|Cl3\|H1\)/\1CL3/'
echo Also, remove [defaults] and [atomtypes] from $SOLVENT_LC.top
echo
