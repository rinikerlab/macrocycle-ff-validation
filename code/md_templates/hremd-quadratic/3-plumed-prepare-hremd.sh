#!/bin/bash

# Replace this with your own path
# Also set up plumed
# source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash

AWK_CHECK_UNDERSCORE='
BEGIN{
    has_underscore=0;
    current_flag="none";
}
(NF == 3) && ($1 == "[") && ($3 == "]") {
    current_flag=$2;
    next;
}
current_flag=="atoms" {
    last_letter = substr($2, length($2), 1);
    if (last_letter == "_") {
        has_underscore=1;
    }
}
END{
    if (!has_underscore) {
        print "no underscore!";
        exit 1;
    }
}
'

# PREPARATION WITH GROMPP

# # produce a processed topology
# grompp -pp
# # choose the "hot" atoms
# vi processed.top
# # generate the actual topology
# plumed partial_tempering $scale < processed.top > topol$i.top

# HOW TO CHOOSE "HOT" ATOMS
# * Edit the [atoms] section in the top file
# * add an _ to the atom types (second column)

# CHECKING THE TOPOLOGIES

# grompp -o topol-unscaled.tpr
# grompp -pp
# vi processed.top # choose the "hot" atoms appending "_". You can choose whatever.
# plumed partial_tempering 1.0 < processed.top > topol-scaled.top # scale with factor 1
# grompp -p topol-scaled.top -o topol-scaled.tpr
# # Then do a rerun on a trajectory
# mdrun -s topol-unscaled.tpr -rerun rerun.trr
# mdrun -s topol-scaled.tpr -rerun rerun.trr
# # and compare the resuling energy files. they should be identical


# ASSERT THAT ATOMS FOR TEMPERING WERE CHOSEN!
awk "$AWK_CHECK_UNDERSCORE" outputs/processed.top || exit

lambdas=(1.0000 0.8859 0.7787 0.6785 0.5851 0.4987 0.4191 0.3465 0.2807 0.2219 0.1700 0.1250)

# GENERATE SCALED TOPOLOGIES

# clean directory
rm -fr \#*
rm -fr topol*

cd outputs

for i in ${!lambdas[@]}
do
  lambda=${lambdas[$i]}
  echo $lambda
  # process topology
  # (if you are curious, try "diff topol0.top topol1.top" to see the changes)
  plumed partial_tempering $lambda < processed.top > topol$i.top
  # prepare tpr file
  # -maxwarn is often needed because box could be charged
  # grompp_mpi_d  -maxwarn 0 -o topol$i.tpr -f grompp$i.mdp -p topol$i.top
  gmx_mpi grompp -maxwarn 0 -o topol$i.tpr -f sampling.mdp -p topol$i.top -c npt.gro
done
