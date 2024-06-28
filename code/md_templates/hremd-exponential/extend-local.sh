#!/bin/bash

source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash
cd outputs/

steps=$((25000000*4))

for i in {0..7}; do
	gmx_mpi convert-tpr -nsteps $steps -s topol$i.tpr -o topol$i-2.tpr
	ln -fs ../topol$i-2.tpr remd$i/topol.tpr
done

# FROM ANNA'S SUBMIT SCRIPT
mpirun -np 8 gmx_mpi mdrun -v -multidir remd{0..7} -replex 100 -g logfile.log -plumed ../plumed.dat -hrex -dlb no -cpi state
# mpirun -np 8 gmx_mpi mdrun -v -multi 8 -replex 100 -nsteps 25000000 -g logfile.log -plumed plumed.dat -hrex -nb gpu -dlb no
