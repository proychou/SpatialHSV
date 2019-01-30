#!/bin/tcsh
setenv scale $1
setenv model $2
setenv k $3
foreach thisdir ( final_fit_2 final_fit_[1-9]2 final_fit_[1-9]2[0-9] final_fit_[1-9][0-9]2[0-9])
    cd $thisdir
    echo $thisdir > /dev/stderr
    setenv ok_epis `cat results.csv | awk -F',' '{if (NR > 1 && $6 > 1400) {n++;}}END {printf("%d\n",n);}'`
    if ( $ok_epis > 0 && ! -f scores.csv ) then
	setenv beta `grep infectivity_free *.in|tail -n 1 | awk '{print $2;}'`
	setenv viral_prod `grep vir_prod_rate_mean *.in|tail -n 1 | awk '{print $2;}'`
	setenv viral_diff `grep diff_rate_virus *.in|tail -n 1 | awk '{print $2;}'`
	setenv trm_speed `grep trm_motility *.in|tail -n 1 | awk '{print $2;}'`
	setenv cyt_diff `grep cyt_diff_rate *.in|tail -n 1 | awk '{print $2;}'`
	setenv cyt_uptake `grep cyt_uptake *.in|tail -n 1 | awk '{print $2;}'`
	setenv beta_ic50 `grep beta_ic50 *.in|tail -n 1 | awk '{print $2;}'`
	setenv inf_ic50 `grep inf_ic50 *.in|tail -n 1 | awk '{print $2;}'`
	setenv prod_ic50 `grep prod_ic50 *.in|tail -n 1 | awk '{print $2;}'`
	setenv cyt_trm_ic50 `grep cyt_trm_ic50 *.in|tail -n 1 | awk '{print $2;}'`
	setenv cyt_inf_death `grep cyt_effect_infncelldeath *.in|tail -n 1 | awk '{print $2;}'`
	setenv hsv_fract `grep hsv_fract *.in|tail -n 1 | awk '{print $2;}'`
	Rscript ../char_final.R $beta $viral_prod $viral_diff $cyt_diff $cyt_uptake $beta_ic50 $inf_ic50 $prod_ic50 $scale $cyt_trm_ic50 $cyt_inf_death $hsv_fract
	head -n 84 all_episodes.csv > temp_episodes.csv
	mv temp_episodes.csv all_episodes.csv
	Rscript ../new_rescore_runs.R $beta $viral_prod $viral_diff $cyt_diff $cyt_uptake $beta_ic50 $inf_ic50 $prod_ic50 $scale $cyt_trm_ic50 $cyt_inf_death $model $k $hsv_fract
    else
	echo "Skipped" > /dev/stderr
	touch scores.csv
    endif
    cd ..
end
