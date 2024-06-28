#!/bin/bash

# Replace this with your own path
# source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash
cd outputs/

for i in {0..11}; do
	mkdir remd$i
	ln -fs ../topol$i.tpr remd$i/topol.tpr
done

touch plumed.dat
# FROM ANNA'S SUBMIT SCRIPT
mpirun -oversubscribe -np 12 gmx_mpi mdrun -v -multidir remd{0..11} -replex 100 -g logfile.log -plumed ../plumed.dat -hrex -dlb no
