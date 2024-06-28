#!/bin/bash
# Create pdb and pqr files without solvent, with CONECT info.
# These files can be used e.g. as "topology" for mdtraj or VMD.

# Important note:
# it seems that cpptraj does not parse the topologies correctly if there are no
# empty lines between blocks. This sometimes leads to errors, which can be
# fixed by adding empty lines between the blocks.

GMX="/localhome/fwaibl/programs/gromacs/bin/gmx"
SOLVENT="REPLACEME"

$GMX editconf -f outputs/npt.gro -o outputs/npt.pdb
sed "/$SOLVENT/d" outputs/npt.pdb > outputs/npt-dry.pdb

cpptraj << EOF 
parm outputs/processed.top
trajin outputs/npt.gro
strip :$SOLVENT
trajout outputs/npt-bonds.pdb conect
trajout outputs/npt-bonds.pqr pdb conect dumpq
EOF
