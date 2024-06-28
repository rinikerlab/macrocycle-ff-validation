#!/bin/bash

NAME="REPLACEME"

mkdir outputs
cd outputs || exit
# cp ../../../../parameterize/amber/outputs/model.prep .
# cp ../../../../parameterize/amber/outputs/sqm.pdb .
# cp ../../../../parameterize/amber/outputs/$NAME.frcmod .
cat << EOF > leap.in
source leaprc.gaff2
# source leaprc.water.tip3p
loadamberparams $NAME.frcmod
loadamberprep model.prep
# loadoff solvents.lib
pdb = loadpdb sqm.pdb
# translate WAT {0 0 -10}
# pdb_with_wat = sequence { pdb WAT }
# solvateBox pdb CHCL3BOX 22 iso
check pdb
saveamberparm pdb prepared.parm7 prepared.rst7
savepdb pdb prepared.pdb
quit
EOF
tleap -f leap.in
#cat << EOF > leap.in
#source leaprc.gaff2
#loadamberparams $NAME.frcmod
#BC3 = loadmol2 model.mol2
#check BC3
#saveoff BC3 BC3.lib
#quit
#EOF
#tleap -f leap.in
#
#cat << EOF > leap2.in
#source leaprc.gaff2
#loadamberparams $NAME.frcmod
#loadoff BC3.lib
#pdb = loadpdb sqm.pdb
#check pdb
#saveamberparm pdb model.parm7 model-optimized.rst7
#savepdb pdb model-leap.pdb
#quit
#EOF
#tleap -f leap2.in
