#!/bin/bash
#SBATCH --ntasks __NJOBS__
#SBATCH --gpus 1
#SBATCH --cpus-per-task=2
#SBATCH --time=120:00:00

# Replace this with your own path
# Also set up plumed
# source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash

# this needs to be synced also with the SLURM block!
N_JOBS=__NJOBS__
export OMP_NUM_THREADS=2

let max_job=N_JOBS-1

cd outputs/

for i in $( seq 0 $max_job ); do
	mkdir remd$i
	# ln -fs ../topol$i.tpr remd$i/topol.tpr
	cp topol$i.tpr remd$i/topol.tpr
done

touch plumed.dat
# FROM ANNA'S SUBMIT SCRIPT
jobdirs=$( for i in $( seq 0 $max_job ); do echo remd$i; done )
mpirun -np $N_JOBS gmx_mpi mdrun -v -multidir $jobdirs -replex 100 -g logfile.log -plumed ../plumed.dat -hrex -dlb yes -pin on -notunepme
