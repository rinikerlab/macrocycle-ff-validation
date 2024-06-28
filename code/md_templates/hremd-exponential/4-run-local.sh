#!/bin/bash

# Replace this with your own path
# Also set up plumed
# source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash
cd outputs/

N_JOBS=__NJOBS__
let max_job=N_JOBS-1

for i in $( seq 0 $max_job ); do
	mkdir remd$i
	ln -fs ../topol$i.tpr remd$i/topol.tpr
done

touch plumed.dat
# FROM ANNA'S SUBMIT SCRIPT
jobdirs=$( for i in $( seq 0 $max_job ); do echo remd$i; done )
mpirun -oversubscribe -np __NJOBS__ gmx_mpi mdrun -v -multidir $jobdirs -replex 100 -g logfile.log -plumed ../plumed.dat -hrex -dlb no
