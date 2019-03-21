#!/bin/bash
# This script can be used to reanalyze runs (in case scoring changes)
# Counts on scores.csv files to NOT be in tar files.
if [ "$1" != "" ]
then
sets=$1
else
sets=1
fi

if [ "$2" != "" ]
then
model=$2
else
model=1
fi

#MODELS FOR FITTING / ADDITIONAL ASSUMPTIONS:
#Non-baseline (any model but model 1):
    # adds 2 variables: cytokine diffusion, cytokine uptake

#1) Cytokines decrease infected cell death rate
    #a. Death_cy = max death rate_cy * (1- (1 / (1 + (cy/IC50) ) ) ) + 1.25
    # adds 2 variables: inf cell max death rate, inf cell death ic50,

#2) Cytokines decrease infectivity
    #a. Beta_cy = beta / ( 1 + (cy/IC50) )
    # adds 1 variable: beta effect ic50

#3) Cytokines decrease viral replication rate
    #a. P_cy = p / ( 1 + (cy/IC50) )
    # adds 1 variable: prod effect ic50

#4) Cytokines activate all other Trm including bystanders & Trm specific T cells (cause them to proliferate & secrete cytokines)
    #a. Prolif_cy = theta * (1- (1 / (1 + (cy/IC50) ) ) )
    #b. Cyprod_cy = cyprod * (1- (1 / (1 + (cy/IC50) ) ) )
    #c. Theta & cyprod are the rates observed in Trm that recognize HSV infected cells
    # adds 1 variable: cyto/trm effect ic50

#16 COMPETING MODELS:
#1) BASELINE
#2) BASELINE + 1
#3) BASELINE + 2
#4) BASELINE + 3
#5) BASELINE + 4
#6) BASELINE + 1 + 2
#7) BASELINE + 1 + 3
#8) BASELINE + 1 + 4
#9) BASELINE + 2 + 3
#10) BASELINE + 2 + 4
#11) BASELINE + 3 + 4
#12) BASELINE + 1 + 2 + 3
#13) BASELINE + 1 + 2 + 4
#14) BASELINE + 1 + 3 + 4
#15) BASELINE + 2 + 3 + 4
#16) BASELINE + 1 + 2 + 3 + 4

rm -f model_fitting.log
while [ $sets -gt 0 ]
do
    echo -n "Extracting simulation data for model $model at " >> model_fitting.log
    tar xvzf ../fitting_dirs_model$model\.tar.gz
    date >> model_fitting.log
    #baseline
    if [ "$model" -eq "1" -o "$model" -eq "17"  ]
    then
        k=0
    fi
    if [ "$model" -eq "2" -o "$model" -eq "11" -o "$model" -eq "18" -o "$model" -eq "27" ]
    then
	k=3
    fi
    if [ "$model" -gt "2" -a "$model" -lt "6"  -o "$model" -gt "18" -a "$model" -lt "22"  ]
    then
	k=2
    fi
    if [ "$model" -ge "6" -a "$model" -le "11" -o "$model" -ge "22" -a "$model" -le "27" ]
    then
	k=4
    fi
    if [ "$model" -ge "12" -a "$model" -lt "16" -o "$model" -ge "28" -a "$model" -lt "32" ]
    then
	k=5
    fi
    if [ "$model" -eq "16" -o "$model" -eq "32" ]
    then
	k=6
    fi
    echo "Reanalyzing & scoring sim data for $model" >> model_fitting.log
    ../scripts/analyze_fits.sh 1 $model $k
    sleep 60
    missing=`../scripts/find_missing_scores.sh | wc -l`
    while [ $missing -gt 1 ]
    do
	echo -n "Waiting for $missing scores to complete at " >> model_fitting.log
	date >> model_fitting.log
	sleep 120
	missing=`../scripts/find_missing_scores.sh | wc -l`
    done
    echo -n "Waiting for final scoring to complete at " >> model_fitting.log
    date >> model_fitting.log
    sleep 60
    echo -n "Combining run scores at " >> model_fitting.log
    date >> model_fitting.log
    ../scripts/combo_all_cats.sh
    echo -n "Copying scoring files at " >> model_fitting.log
    date >> model_fitting.log
    top_run=`../scripts/get_top_run.sh`
    ../scripts/copy_overall_files.sh $model $top_run

    echo "Removing interim result directories for model $model" >> model_fitting.log
    rm -r final_fit_[1-9]*

    model=`expr $model + 1`
    sets=`expr $sets - 1`
done
