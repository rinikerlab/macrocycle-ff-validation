This repository contains code and data to reproduce the figures in the paper "Validating Small-Molecule Force Fields for Macrocyclic Compounds Using NMR Data in Different Solvents" (J. Chem. Inf. Model. 2024, 64, 20, 7938–7948, https://doi.org/10.1021/acs.jcim.4c01120).

# Data in this repository
In the paper, we ran REST2 simulations of 11 compounds, 4 force fields, partially in several solvents, and partially with 2 different settings. The nomenclature is as follows:
* 11 compounds: BC1 (`begnini-compound-1`, residue name `BC1`), BC2 (`begnini-compound-2`, `BC2`), G16 (`poongavanam-g16`, `G16`), E2-enant (`poongavanam-e2-enant`, `PE2`), rifampicin (`danelius-rifampicin`, `RIF`), roxithromycin (`danelius-roxithromycin`, `ROX`), telithromycin (`danelius-telithromycin`, `TEL`), spiramycin (`danelius-spiramycin`, `SPI`), lorlatinib (`peng-lorlatinib`, `LOR`), NLeu5R (`comeau-nleu5r`, `N5R`), NLeu5S (`comeau-nleu5s`, `N5S`). Some compounds are protonated in water: rifampicin (`danelius-rifampicin-charged`), roxithromycin (`danelius-roxithromycin-charged`), telithromycin (`danelius-telithromycin-charged`), spiramycin (`danelius-spiramycin-charged`)
* forcefield: OpenFF 2.0 (`openff`), GaFF 2 (`amber`), OPLS/AA (`opls`), XFF with DASH charges (`xff-dash`)
* method: REST2 with quadratic lambda placement (`hremd-quadratic`), REST2 with exponential lambda placement (`hremd-exponential`)
* solvent: chloroform (`chcl3`), water (`water`), DMSO (`dmso`)

# Reproducing the Figures
All figures can be reproduced by following these steps:
* clone this repository: `git clone git@github.com:rinikerlab/macrocycle-ff-validation.git`
* create a new environment:
```
conda env create -f environment.yml
conda activate macrocycle-ff-benchmark
```
* Run the notebook `code/Create-Figures.ipynb`

# Re-running simulations
The data folder contains solvated topologies and .gro files for each combination of compound, solvent, and force-field. To re-run a simulation, follow the following steps:
* create a folder `data/COMPOUND/SOLVENT/equilibrate/FORCEFIELD`, and copy the content of `code/md_templates/equilibrate` there.
* set all placeholders in `equilibration.sh`, and run it. `FORCEFIELD` and `SOLVENT_LC` should be replaced to a name matching the folder structure ("openff" / "amber" / "opls" / "xff-dash" and "water" / "chcl3" / "dmso"), and `SOLUTE` should match the residue name in the topology.
* create a folder `data/COMPOUND/SOLVENT/hremd-quadratic/FORCEFIELD`, and copy the content of `code/md_templates/hremd-quadratic` there.
* set all placeholders in `1-get-inputs.sh`, and execute the scripts in the order 1-4. Use `4-run-local.sh` to run on the current PC. `4-run-euler.sh` is the submit script that was used on the ETH Euler cluster, and might be used as a template to run on other cluster systems.

As an example, you can run the following to start one equilibration + REST2 simulation (after cloning this repository, starting from the base folder)
```
CMP=poongavanam-e2-enant
CMP_NAME=PE2
FF=amber
SOLV=dmso

# Equilibration
eq_dir=data/$CMP/$SOLV/equilibrate/$FF
mkdir -p $eq_dir
cp code/md_templates/equilibrate/equilibrate.sh $eq_dir/ || exit 1
( cd $eq_dir
	sed -i "s/FORCEFIELD=.*/FORCEFIELD=$FF/;s/SOLVENT_LC=.*/SOLVENT_LC=$SOLV/;s/SOLUTE=.*/SOLUTE=$CMP_NAME/;" equilibrate.sh
	bash equilibrate.sh || exit 1
)

# H-REMD / REST2
# Note: if you don't have at least 12 CPU cores (1 per replica), this will oversubscribe and might be inefficient.
hremd_dir=data/$CMP/$SOLV/hremd-quadratic/$FF
mkdir -p $hremd_dir
cp code/md_templates/hremd-quadratic/* $hremd_dir/ || exit 1
( cd $hremd_dir
	sed -i "s/FORCEFIELD=.*/FORCEFIELD=$FF/;s/SOLVENT_LC=.*/SOLVENT_LC=$SOLV/;s/SOLUTE=.*/SOLUTE=$CMP_NAME/;" 1-get-inputs.sh
	bash 1-get-inputs.sh || exit 1
	bash 2-make-single-topology.sh || exit 1
	bash 3-plumed-prepare-hremd.sh || exit 1
	bash 4-run-local.sh || exit 1
)
```

After the simulation is done, run `make-dry-pdb.sh` and then python `run-analysis.py --compound COMPOUND --forcefield FORCEFIELD --method METHOD --solvent SOLVENT`, with parameters matching the folder naming as explained before.

# Re-running the parameterizations
To re-run the parameterization, you can start from SMILES or .mol files.
* the script `code/md_templates/create-initial-structure/smiles-to-structure.py` converts a SMILES code to a conformer.
* the scripts in `code/md_templates/parameterize` can be used to assign force-field parameters for each force field. Note that for OPLS, the molecule must be manually uploaded to LigParGen (https://zarbi.chem.yale.edu/ligpargen/), and for XFF, the molecule must be uploaded to the XFF web server (https://xff.xtalpi.com/)
