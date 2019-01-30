#!/bin/tcsh
# This scripts makes the overall score file from those in all the subdirectories
#model,k,log beta,log viral prod,log viral diff,log cyto diff,log cyto uptake,log inf ic50,log beta ic50,log prod ic50,log cyt_trm_ic50,cyt_inf_death,hsv_fract,peak VL RSS,peak VL AIC,peak time RSS,peak time AIC,duration RSS,duration AIC,cat2 RSS,cat2 AIC,cat3 RSS,cat3 AIC
rm -f final_fit_scores.csv
if ( -f final_fit_1/scores.csv ) then
head -n 1 final_fit_1/scores.csv > final_fit_scores.csv
else if (-f final_fit_10/scores.csv ) then
head -n 1 final_fit_10/scores.csv > final_fit_scores.csv
endif
foreach file (final_fit_*/scores.csv)
tail -n 1 $file | grep -v NA >> final_fit_scores.csv
end
head -n 1 final_fit_1/scores.csv 
echo "Top score using peak RSS fit"
sort -n -k14 -t',' final_fit_scores.csv | grep -v model | head -n 1
echo
echo "Top score using peak AIC fit"
sort -n -k15 -t',' final_fit_scores.csv | grep -v model | head -n 1
echo
echo "Top 2 cat RSS score"
sort -n -k20 -t',' final_fit_scores.csv | grep -v model | head -n 1
echo
echo "Top 2 cat AIC score"
sort -n -k21 -t',' final_fit_scores.csv | grep -v model | head -n 1
echo
echo "Top 3 cat RSS score"
sort -n -k22 -t',' final_fit_scores.csv | grep -v model | head -n 1
echo
echo "Top 3 cat AIC score"
sort -n -k23 -t',' final_fit_scores.csv | grep -v model | head -n 1
