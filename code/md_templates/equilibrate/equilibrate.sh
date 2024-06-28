#!/bin/bash

FORCEFIELD="REPLACEME"
SOLUTE="REPLACEME"
SOLVENT_LC="REPLACEME"

TOP_FILE_EXT="top"

# 4.5e-5 for water, 1.0e-4 for chloroform
# 5.5e-5 for DMSO, see https://doi.org/10.1023/A:1022985208144
if [ "$SOLVENT_LC" == "water" ]; then
  COMPRESSIBILITY="4.5e-5"
  SOLVENT="SOL"
elif [ "$SOLVENT_LC" == "chcl3" ]; then
  COMPRESSIBILITY="1.0e-4"
  SOLVENT="CL3"
elif [ "$SOLVENT_LC" == "dmso" ]; then
  COMPRESSIBILITY="5.5e-5"
  SOLVENT="DSO"
else
  echo Unknown solvent \(SOLVENT_LC\)\; Exiting...
  exit 1
fi

if [[ "$FORCEFIELD" == "REPLACEME" || "$SOLUTE" == "REPLACEME" || "$SOLVENT_LC" == "REPLACEME" ]]; then
    echo Please replace all placeholders!
    exit 1
fi

# Get input files. You might have to update this part!
mkdir outputs

cp ../../solvate/$FORCEFIELD/outputs/solvated.top outputs/topol.top
cp ../../solvate/$FORCEFIELD/outputs/$SOLVENT_LC.$TOP_FILE_EXT outputs/$SOLVENT_LC.$TOP_FILE_EXT
cp ../../solvate/$FORCEFIELD/outputs/solvated.gro outputs/solvated.gro

GMX='gmx_mpi'

MDP_TEMPLATE_MIN="
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
; nstlist         = 1         ; Frequency to update the neighbor list and long range forces
; cutoff-scheme   = Verlet    ; Buffered neighbor searching
; ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
"

MDP_TEMPLATE_NVT="
title                   = NVT equilibration
; define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
; lincs_iter              = 1         ; accuracy of LINCS
; lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
; cutoff-scheme           = Verlet    ; Buffered neighbor searching
; ns_type                 = grid      ; search neighboring grid cells
; nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
; pme_order               = 4         ; cubic interpolation
; fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = $SOLUTE $SOLVENT   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
"

MDP_TEMPLATE_NPT="
title                   = NPT equilibration 
; define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 500000    ; 2 * 500000 = 1000 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
; lincs_iter              = 1         ; accuracy of LINCS
; lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
; cutoff-scheme           = Verlet    ; Buffered neighbor searching
; ns_type                 = grid      ; search neighboring grid cells
; nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
; pme_order               = 4         ; cubic interpolation
; fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = $SOLUTE $SOLVENT   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = C-rescale     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = $COMPRESSIBILITY                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 
"

cd outputs || exit

function run_gmx() {
    local mdp_string=$1
    local basename=$2
    local topology=$3
    local start_coord_base=$4
    local restart=$5

    local cpt_opt=""
    if [[ "$restart" == "restart" ]]; then
        cpt_opt="-t $start_coord_base.cpt"
    fi
    cat <<< "$mdp_string" > $basename.mdp
    $GMX grompp -f $basename.mdp -c $start_coord_base.gro -r $start_coord_base.gro $cpt_opt -p $topology -o $basename.tpr || exit 1
    $GMX mdrun -deffnm $basename -nb gpu

}

run_gmx "$MDP_TEMPLATE_MIN" min topol.top solvated
run_gmx "$MDP_TEMPLATE_NVT" nvt topol.top min
run_gmx "$MDP_TEMPLATE_NPT" npt topol.top nvt restart
