#!/bin/sh
top_score=`sort -n -k21 -t',' final_fit_scores.csv | grep -v model | head -n 1 | awk -F',' '{print $20;}'`
echo "Top cat2 score is $top_score (RSS)" >> model_fitting.log
top_file=`grep -l $top_score final_fit_*/scores.csv`
echo "Top file is $top_file" >> model_fitting.log
top_run=`echo "$top_file" | awk -F'/' '{split($1,a,"_");printf("%s\n",a[3]);}'`
echo "Top run is $top_run" >> model_fitting.log
echo "$top_run"
