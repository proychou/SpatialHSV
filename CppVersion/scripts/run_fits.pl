#!/usr/bin/perl
# script to launch a single run with the fitting parameters drawn randomly from 
# given range of values (either log-based or as absolute values)
# comment out variables that are not currently being varied
# certain parameters are only relevant to certain models (see run_models.sh)
if ($#ARGV < 0) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 1 input arguments ($args given)\n";
}
if ($ARGV[1] =~ /^\d+/) {
    $model=$ARGV[1];
} else {
    die "2nd arg is the model (vs $ARGV[1])\n";
    exit (1);
}
my $infile;
my $outdir;

if ( -f $ARGV[0]) {
    $infile=$ARGV[0]
} else {
    die "1st arg is input file\n";
    exit (1);
}

# log scale for all variables

#my $low_viral_beta=-1;
#my $high_viral_beta=1;
#my $range_viral_beta= $high_viral_beta-$low_viral_beta;

#my $low_viral_diff=1;
#my $high_viral_diff=4;
#my $range_viral_diff= $high_viral_diff-$low_viral_diff;

#my $low_viral_prod=4;
#my $high_viral_prod=5;
#my $range_viral_prod= $high_viral_prod-$low_viral_prod;

#my $low_hsv_fract=0.05;
#my $high_hsv_fract=0.4;
#my $range_hsv_fract= $high_hsv_fract-$low_hsv_fract;

my $low_cyto_diff=1;
my $high_cyto_diff=4;
my $range_cyto_diff= $high_cyto_diff-$low_cyto_diff;

my $low_cyto_uptake=-2;
my $high_cyto_uptake=1;
my $range_cyto_uptake= $high_cyto_uptake-$low_cyto_uptake;

my $low_prod_ic50=-7;
my $high_prod_ic50=-4;
my $range_prod_ic50= $high_prod_ic50-$low_prod_ic50;

my $low_beta_ic50=-7;
my $high_beta_ic50=-4;
my $range_beta_ic50= $high_beta_ic50-$low_beta_ic50;

my $low_inf_ic50=-7;
my $high_inf_ic50=-4;
my $range_inf_ic50= $high_inf_ic50-$low_inf_ic50;

my $low_cyt_trm_ic50=-7;
my $high_cyt_trm_ic50=-4;
my $range_cyt_trm_ic50= $high_cyt_trm_ic50-$low_cyt_trm_ic50;

#my $low_cyt_inf_death=0;
#my $high_cyt_inf_death=2;
#my $range_cyt_inf_death= $high_cyt_inf_death-$low_cyt_inf_death;

#my $viral_beta= $range_viral_beta*rand() + $low_viral_beta;
#my $viral_diff= $range_viral_diff*rand() + $low_viral_diff;
#my $viral_prod= $range_viral_prod*rand() + $low_viral_prod;
#my $hsv_fract= $range_hsv_fract*rand() + $low_hsv_fract;
my $cyto_diff= $range_cyto_diff*rand() + $low_cyto_diff;
my $cyto_uptake= $range_cyto_uptake*rand() + $low_cyto_uptake;
my $beta_ic50= $range_beta_ic50*rand() + $low_beta_ic50;
my $prod_ic50= $range_prod_ic50*rand() + $low_prod_ic50;
my $inf_ic50= $range_inf_ic50*rand() + $low_inf_ic50;
#my $cyt_inf_death= $range_cyt_inf_death*rand() + $low_cyt_inf_death;
my $cyt_trm_ic50= $range_cyt_trm_ic50*rand() + $low_cyt_trm_ic50;

open (TSTFILE, ">>$infile");

if ( $model != 1 && $model != 17 ) {
    printf TSTFILE ("cyt_diff_rate %e\n",10**$cyto_diff);
    printf TSTFILE ("cyt_uptake %e\n",10**$cyto_uptake);
} else {
    printf TSTFILE ("cyt_diff_rate 0\n");
    printf TSTFILE ("cyt_uptake 0\n");
}
if ( $model == 2 || ($model >= 6 && $model <= 8) || 
    ($model >= 12 && $model <= 14) || $model == 16 ||
     $model == 18 || ($model >= 22 && $model < 25) || 
    ($model >= 28 && $model <= 30) || $model == 32 ) {
    printf TSTFILE ("cyt_effect_infpcelldeath 48\n");
    printf TSTFILE ("cyt_effect_infncelldeath 48\n");
    printf TSTFILE ("inf_ic50 %e\n",10**$inf_ic50);
} else {
    printf TSTFILE ("cyt_effect_infpcelldeath 0\n");
    printf TSTFILE ("cyt_effect_infncelldeath 0\n");
    printf TSTFILE ("inf_ic50 0\n");
}
if ( $model == 3 || $model == 6 || $model == 9 || $model == 10 ||
    $model == 12 || $model == 13 || $model == 15 || $model == 16||
     $model == 19 || $model == 22 || $model == 25 || $model == 26 ||
    $model == 28 || $model == 29 || $model == 31 || $model == 32) {
    printf TSTFILE ("cyt_effect_infectivity 1\n");
    printf TSTFILE ("beta_ic50 %e\n",10**$beta_ic50);
} else {
    printf TSTFILE ("cyt_effect_infectivity 0\n");
    printf TSTFILE ("beta_ic50 0\n");
}
if ( $model == 4 || $model == 7 || $model == 9 || $model == 11 ||
    $model == 12 || $model == 14 || $model == 15 || $model == 16||
     $model == 20 || $model == 23 || $model == 25 || $model == 27 ||
    $model == 28 || $model == 30 || $model == 31 || $model == 32) {
    printf TSTFILE ("cyt_effect_virprod 1\n");
    printf TSTFILE ("prod_ic50 %e\n",10**$prod_ic50);
} else {
    printf TSTFILE ("cyt_effect_virprod 0\n");
    printf TSTFILE ("prod_ic50 0\n");
}
if ( $model == 5 || $model == 8 || $model == 10 || $model == 11 ||
    ($model >= 13 &&  $model <= 16) ||
     $model == 21 || $model == 24 || $model == 26 || $model == 27 ||
    $model >= 29) {
    printf TSTFILE ("cyto_act 1\n");
    printf TSTFILE ("cyt_trm_ic50 %e\n",10**$cyt_trm_ic50);
} else {
    printf TSTFILE ("cyto_act 0\n");
    printf TSTFILE ("cyt_trm_ic50 0\n");
}
if ( $model >= 17 ) {
    printf TSTFILE ("trm_move_method 2\n");
}

# these parameters held the same for all sets
my $viral_prod=4.5;
printf TSTFILE ("infectivity_free 0.1\n");
printf TSTFILE ("diff_rate_virus 100\n");
printf TSTFILE ("hsv_fract 0.4\n");
printf TSTFILE ("vir_prod_rate_mean %e\n",10**$viral_prod);
close(TSTFILE);

system("../Spatial_HSV -r -f $infile -d");
