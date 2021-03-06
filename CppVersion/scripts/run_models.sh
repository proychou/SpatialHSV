#!/bin/bash
if [ "$1" != "" ]
then
runs=$1
else
runs=1000
fi

if [ "$2" != "" ]
then
sets=$2
else
sets=1
fi

if [ "$3" != "" ]
then
model=$3
else
model=1
fi

#MODELS FOR FITTING / ADDITIONAL ASSUMPTIONS:
#Baseline: varies viral diffusion

#Non-baseline (any model but model 1):
    # adds 1 variables: cytokine uptake (cytokine diffusion is fixed)

#1) Cytokines decrease infected cell death rate
    #a. Death_cy = max death rate_cy * (1- (1 / (1 + (cy/IC50) ) ) ) + 1.25
    # adds 2 variables: inf cell max death rate, inf cell death ic50,
    # CHANGE 10/11 - add just one variable (max death rate fixed)

#2) Cytokines decrease infectivity
    #a. Beta_cy = beta / ( 1 + (cy/IC50) )
    # adds 1 variable: beta effect ic50 (effect only since ic50 now fixed)

#3) Cytokines decrease viral replication rate
    #a. P_cy = p / ( 1 + (cy/IC50) )
    # adds 1 variable: prod effect ic50 (effect only since ic50 now fixed)

#4) Cytokines activate all other Trm including bystanders & Trm specific T cells (cause them to proliferate & secrete cytokines)
    #a. Prolif_cy = theta * (1- (1 / (1 + (cy/IC50) ) ) )
    #b. Cyprod_cy = cyprod * (1- (1 / (1 + (cy/IC50) ) ) )
    #c. Theta & cyprod are the rates observed in Trm that recognize HSV infected cells
    # adds 1 variable: cyto/trm effect ic50 (effect only since ic50 now fixed)
#5) Cytokine attraction gradient
    # does not add any variables; just changes Trm motion

#32 COMPETING MODELS:
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

rm -f model_fitting$model\.log
while [ $sets -gt 0 ]
do
    #cp base_fitting.in ic50_fitting.in
    echo -n "Launching simulations for model $model at " >> model_fitting$model\.log
    date >> model_fitting$model\.log
    #baseline
    if [ "$model" -eq "1" ]
    then
        k=1
    fi
    if [ "$model" -eq "2" -o "$model" -eq "11" ]
    then
	k=3
    fi
    if [ "$model" -gt "2" -a "$model" -lt "6"  ]
    then
	k=2
    fi
    if [ "$model" -ge "6" -a "$model" -le "11" ]
    then
	k=4
    fi
    if [ "$model" -ge "12" -a "$model" -lt "16"]
    then
	k=5
    fi
    if [ "$model" -ge "16" ]
    then
	k=6
    fi
    ../scripts/run_final.sh $model $runs
    sleep 60
    running=`squeue -u dswan -M beagle | wc -l`
    #running=`squeue -u dswan | wc -l`
    while [ $running -gt 2 ]
    do
	running=`expr $running - 1`
	echo -n "Waiting for $running runs to complete at " >> model_fitting$model\.log
	date >> model_fitting$model\.log
	sleep 300
	#running=`squeue -u dswan | wc -l`
	running=`squeue -u dswan -M beagle | wc -l`
    done
    echo -n "Analyzing result files at " >> model_fitting$model\.log
    date >> model_fitting$model\.log
    ./analyze_ic50_types.sh 1 $model $k
    sleep 60
    missing=`../scripts/find_missing_scores.sh | wc -l`
    while [ $missing -gt 1 ]
    do
	echo -n "Waiting for $missing scores to complete at " >> model_fitting$model\.log
	date >> model_fitting$model\.log
	sleep 30
	missing=`../scripts/find_missing_scores.sh | wc -l`
    done
    echo -n "Waiting for final scoring to complete at " >> model_fitting$model\.log
    date >> model_fitting$model\.log
    sleep 30
    echo -n "Combining run scores at " >> model_fitting$model\.log
    date >> model_fitting$model\.log
    ../scripts/combo_all_cats.sh $model
    echo -n "Copying scoring files at " >> model_fitting$model\.log
    date >> model_fitting$model\.log
    top_run=`../scripts/get_top_cat2.sh`
    ../scripts/copy_overall_files.sh $model $top_run

    echo "Saving interim result directories for model $model" >> model_fitting$model\.log
    tar -cvzf new_fitting_dirs_model$model\.tar.gz final_fit_[1-9]*/all_episodes.csv final_fit_[1-9]*/ic50_fitting.in final_fit_[1-9]*/results.csv

    echo "Removing interim result directories for model $model" >> model_fitting$model\.log
    rm -r final_fit_[1-9]*

    sets=`expr $sets - 1`
    model=`expr $model + 1`
done
