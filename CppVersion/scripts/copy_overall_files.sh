#!/bin/sh
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
sort -n -k23 -t',' final_fit_scores.csv | grep -v model | head -n 50 > top_final_scores_model$model\.csv
cp final_fit_$the_run/all_episodes.csv final_episodes_model$model\.csv 
cp final_fit_$the_run/results.csv final_results_model$model\.csv 
cp final_fit_$the_run/fitting.in final_inputs_model$model\.in 
