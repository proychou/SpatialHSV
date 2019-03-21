#!/bin/bash
# script to run some number of parameter sets for a given model
# set numbering can start at a number other than one to add to existing
if [ "$1" != "" ]
then
model=$1
else
model=1
fi

if [ "$2" != "" ]
then
sets=$2
else
sets=10
fi

if [ "$3" != "" ]
then
starts=$3
else
starts=1
fi

# run all sims varying these parameters (via run_sets.pl)

count=$starts
while [ $sets -gt 0 ]
do
    thisdir="final_fit_"
    thisdir+=$count
    if [ ! -d $thisdir ]
    then
    mkdir $thisdir
    fi
    cd $thisdir
    cp ../fitting.in .
    sbatch -n 1 -t 1:00:00 ../../scripts/run_fits.pl fitting.in $model
    cd ..

    sets=`expr $sets - 1`
    count=`expr $count + 1`
done
