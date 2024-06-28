#!/bin/bash

# Replace this with your own path
# source /localhome/fwaibl/programs/gromacs-plumed/bin/GMXRC.bash

AWK_APPEND_UNDERSCORE_TO_FIRST_MOL='
BEGIN{
    current_flag = "none";
    done = 0;
    currently_appending = 0;
}

/^;/ {
    print;
    next;
}

(NF == 0) {
    print;
    next;
}

(NF == 3) && ($1 == "[") && ($3 == "]") {
    if (currently_appending) {
        done = 1;
    }
    current_flag=$2;
    print;
    next;
}

(current_flag == "atoms") {
    if (done) {
        print;
    } else {
        currently_appending = 1;

        # append underscore to second element
        printf "%s %s_", $1, $2;
        for (i=3; i<=NF; i++) {
            printf " %s", $i;
        }
        printf "\n";
    }
    next;
}

{print}
'

(
    cd outputs
    gmx_mpi grompp -f sampling.mdp -c npt.gro -pp || exit
    mv processed.top grompp.top || exit
    awk "$AWK_APPEND_UNDERSCORE_TO_FIRST_MOL" grompp.top > processed.top || exit
)


# AWK_CHECK_UNDERSCORE='
# BEGIN{
#     has_underscore=0;
#     current_flag="none";
# }
# (NF == 3) && ($1 == "[") && ($3 == "]") {
#     current_flag=$2;
#     next;
# }
# current_flag=="atoms" {
#     last_letter = substr($2, length($2), 1);
#     if (last_letter == "_") {
#         has_underscore=1;
#     }
# }
# END{
#     if (!has_underscore) {
#         print "no underscore!";
#         exit 1;
#     }
# }
# '
