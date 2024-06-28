#!/bin/bash

# uncomment this to run clusterin on all compounds:

# COMPOUND_METHOD="begnini-compound-1 chcl3 hremd-quadratic
# begnini-compound-2 chcl3 hremd-quadratic
# comeau-nleu5r chcl3 hremd-quadratic
# comeau-nleu5s chcl3 hremd-quadratic
# danelius-rifampicin chcl3 hremd-quadratic
# danelius-roxithromycin chcl3 hremd-quadratic
# danelius-spiramycin chcl3 hremd-quadratic
# danelius-telithromycin chcl3 hremd-quadratic
# peng-lorlatinib chcl3 hremd-with-angles
# poongavanam-e2-enant chcl3 hremd-quadratic
# poongavanam-g16 chcl3 hremd-quadratic"

# run only on charged roxithromycin:
# COMPOUND_METHOD="danelius-roxithromycin-charged water hremd-quadratic"

# run only on roxithromycin in chloroform:
# COMPOUND_METHOD="danelius-roxithromycin chcl3 hremd-quadratic"

# run on telithromycin in water and chloroform:
COMPOUND_METHOD="danelius-telithromycin chcl3 hremd-quadratic
danelius-telithromycin-charged water hremd-quadratic"

while read compound solvent method; do
	SEL="!@H="
	TRAJ="../data/$compound/$solvent/$method/openff/outputs/remd0/traj_comp.xtc"
	TOP="../data/$compound/$solvent/$method/openff/outputs/npt-bonds.pdb"
	OUT_DIR="outputs/clustering"

	mkdir $OUT_DIR
	cp $TOP $OUT_DIR/$compound-equilibrated.pdb

	center_circular.py -t $TOP $TRAJ $OUT_DIR/$compound-$solvent-centered.nc

	cpptraj <<- EOF | tee $OUT_DIR/cpptraj-$compound-$solvent.log
	parm $TOP
	trajin $OUT_DIR/$compound-$solvent-centered.nc 10000 50000 1
	cluster hieragglo clusters 5 averagelinkage rms $SEL repout $OUT_DIR/$compound-$solvent-cluster info $OUT_DIR/$compound-$solvent-cluster-info.dat out $OUT_DIR/$compound-$solvent-cluster.dat summary $OUT_DIR/$compound-$solvent-cluster-summary.dat
	EOF

	for i in {0..4}; do
		ambpdb -c $OUT_DIR/$compound-$solvent-cluster.c$i.crd -p $TOP > $OUT_DIR/$compound-$solvent-cluster.c$i.pdb
	done
done <<< "$COMPOUND_METHOD"
