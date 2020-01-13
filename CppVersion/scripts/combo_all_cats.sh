#!/bin/tcsh
setenv model $1
if ( "$model" == "" ) then
    setenv logfile "model_fitting.log"
else
    setenv logfile model_fitting$model\.log
endif
#model,k,log beta,log viral prod,log viral diff,log cyto diff,log cyto uptake,log inf ic50,log beta ic50,log prod ic50,log cyt_trm_ic50,cyt_inf_death,hsv_fract,peak VL RSS,peak VL AIC,peak time RSS,peak time AIC,duration RSS,duration AIC,cat2 RSS,cat2 AIC,cat3 RSS,cat3 AIC
rm -f final_fit_scores.csv
setenv header 0
foreach dir (final_fit_[1-9]*)
    if ( -f $dir/scores.csv && -f $dir/all_episodes.csv ) then
	if ($header == 0) then
	    setenv header 1
	    head -n 1 $dir/scores.csv > final_fit_scores.csv
	    echo "Header from $dir." >> $logfile
	endif
	setenv epis `wc -l $dir/all_episodes.csv | awk '{print $1;}'`
	if ($epis == 84) then
	    tail -n 1 $dir/scores.csv >> final_fit_scores.csv
	else
	    echo "Skipped $dir.  Too few episodes!" >> $logfile
	endif
    endif
end
head -n 1 final_fit_scores.csv >> $logfile
echo "Top score using peak RSS fit" >> $logfile
sort -n -k14 -t',' final_fit_scores.csv | grep -v model | head -n 1 >> $logfile
echo >> $logfile
echo "Top score using peak AIC fit" >> $logfile
sort -n -k15 -t',' final_fit_scores.csv | grep -v model | head -n 1 >> $logfile
echo >> $logfile
echo "Top 2 cat RSS score" >> $logfile
sort -n -k20 -t',' final_fit_scores.csv | grep -v model | head -n 1 >> $logfile
echo >> $logfile
echo "Top 2 cat AIC score" >> $logfile
sort -n -k21 -t',' final_fit_scores.csv | grep -v model | head -n 1 >> $logfile
echo >> $logfile
echo "Top 3 cat RSS score" >> $logfile
sort -n -k22 -t',' final_fit_scores.csv | grep -v model | head -n 1 >> $logfile
echo >> $logfile
echo "Top 3 cat AIC score" >> $logfile
sort -n -k23 -t',' final_fit_scores.csv | grep -v model | head -n 1 >> $logfile
