#!/bin/sh
# changed to drop duration in saving top scores (cat2 vs cat3)
if [ "$1" != "" ]
then
model=$1
else
model=1
fi
if [ "$2" != "" ]
then
the_run=$2
else
the_run=1
fi
cp final_fit_scores.csv final_scores_model$model\.csv
sort -n -k21 -t',' final_fit_scores.csv | grep -v model | head -n 50 > top_cat2_scores_model$model\.csv
cp final_fit_$the_run/all_episodes.csv cat2_epis_model$model\.csv 
cp final_fit_$the_run/results.csv cat2_results_model$model\.csv 
cp final_fit_$the_run/*.in cat2_model$model\.in 
