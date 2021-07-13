#!/bin/bash

# see https://stackoverflow.com/questions/8880603/loop-through-an-array-of-strings-in-bash
# see https://www.thegeekstuff.com/2010/06/bash-array-tutorial
# see https://unix.stackexchange.com/questions/4899/var-vs-var-and-to-quote-or-not-to-quote

declare -a FIG5=("OBST_SMALL" "OBST_LARGE" "OBST_LONG" "OBST_WIDE")
declare -a FIG6A=("HOTS_WEAK" "HOTS_MED" "HOTS_STRONG")
declare -a FIG6B=("HOTS_WEAK_SCALING")
declare -a FIG6C=("HOTS_CIRCLE_WEAK" "HOTS_CIRCLE_STRONG" "HOTS_LONG_WEAK" "HOTS_LONG_STRONG" "HOTS_WIDE_WEAK" "HOTS_WIDE_STRONG")
declare -a FIG7=("REFRACTION")
declare -a FIG7HELPER=("REFRACTION_DETSPEEDS")
declare -a NUMTEST=("OBST_WIDERCHANNEL" "HOTS_WIDERCHANNEL" "OBST_FINERDISC" "HOTS_FINERDISC")

declare -a allfigs=("${FIG5[@]}" "${FIG6A[@]}" "${FIG6B[@]}" "${FIG6C[@]}" "${FIG7[@]}" "${FIG7HELPER[@]}" "${NUMTEST[@]}")

# https://stackoverflow.com/questions/2427995/bash-no-arguments-warning-and-case-decisions
if [[ $# -eq 0 ]] ; then
    echo "Need argument, specifier for this run. Exiting."
    exit 0
fi

for scen in "${allfigs[@]}"
do
    # create new directory
    cp -rp ManyFeatures_template ManyFeatures_${scen}_$1
    cd ManyFeatures_${scen}_$1
    if [ ! -d LOGFILES ]; then
        mkdir LOGFILES
    fi
    # create appropriate parameter file
    cp param_${scen}.py param.py
    # set random number seeds
    for rhoid in $(seq 0 $(python -c "import param; print(param.rhos.size-1)"))
    do
        for gammaid in $(seq 0 $(python -c "import param; print(param.gammas.size-1)"))
        do
            ./setrandnumseed.sh 512 > fixedrandnumseeds_rhoid_${rhoid}_gammid_${gammaid}.txt
        done
    done
    # ./submit.sh
    cd ..
done
