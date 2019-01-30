The scripts directory contains scripts used when analyzing or plotting 
simulation data.  Each of the scripts are listed below with their purpose.

##################################################
1. Scripts for running various model (including parameter fitting)
##################################################

run_models.sh - This script is used to run and score one or more models. 
run_final.sh - script to run some number of parameter sets for a given model
run_fits.pl - perl script to generate a single parameter set from ranges

run_best_model.sh - used in fine tuning (make directories & call perl script)
run_best_sets.pl - used in fine tuning (parameter generation & sim launch)

##################################################
2. Scripts for analyzing & scoring run data (from a particular model)
##################################################

analyze_fits.sh - main analysis script (used on directories final_fit_<N>)
		  It launches analyze_fits_0-9.sh each of which does a portion
		  of the sub-directories (typically 100)

analyze_fits_0-4.sh - used to launch 1st 1/2 of full set of analyzation scripts

analyze_fits_0.sh - see notes above
analyze_fits_1.sh
analyze_fits_2.sh
analyze_fits_3.sh
analyze_fits_4.sh

analyze_fits_5-9.sh - used to launch 2nd 1/2 of full set of analyzation scripts
analyze_fits_5.sh - see notes above
analyze_fits_6.sh
analyze_fits_7.sh
analyze_fits_8.sh
analyze_fits_9.sh

char_final.R - used to extract episode data from result.csv files
new_rescore_runs.R - used to score parameter sets based on episode data

find_missing_scores.sh - used to check for absence of scores.csv files
combo_all_cats.sh - makes the overall score file from those in subdirectories
get_top_run.sh - used after analysis to find top AIC scoring parameter set
copy_overall_files.sh - used to copy files from best set (for plotting)

reanalyze.sh - This script can be used to reanalyze runs (for scoring changes)

##################################################
3. Scripts for making plots and/or movies...
##################################################

one_frame.R - converts state files (.csvs) to 2-4 panel plots (.png)
big_frame.R - same but for larger grids (425x425 vs 125x125)

et_ratio_plot.R - same as above, but for et_ratio data files

make_movie.sh - wrapper around call to menencoder to make avi from pngs
make_et_ratio_movie.sh - same as above, but for et_ratio pngs

