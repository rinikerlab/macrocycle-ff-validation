#!/bin/bash

# Replace this with your own path
# source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash
( cd outputs; gmx_mpi grompp -f sampling.mdp -c npt.gro -pp )
