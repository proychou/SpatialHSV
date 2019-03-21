#!/usr/bin/perl
if ($#ARGV < 0) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 1 input arguments ($args given)\n";
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

my $low_viral_diff=1;
my $high_viral_diff=3;
my $range_viral_diff= $high_viral_diff-$low_viral_diff;

my $low_cyto_diff=1;
my $high_cyto_diff=3;
my $range_cyto_diff= $high_cyto_diff-$low_cyto_diff;

my $low_cyto_uptake=-2;
my $high_cyto_uptake=2;
my $range_cyto_uptake= $high_cyto_uptake-$low_cyto_uptake;

my $viral_diff= $range_viral_diff*rand() + $low_viral_diff;
my $cyto_diff= $range_cyto_diff*rand() + $low_cyto_diff;
my $cyto_uptake= $range_cyto_uptake*rand() + $low_cyto_uptake;

open (TSTFILE, ">>$infile");

printf TSTFILE ("cyt_diff_rate %e\n",10**$cyto_diff);
printf TSTFILE ("cyt_uptake %e\n",10**$cyto_uptake);
printf TSTFILE ("diff_rate_virus %e\n",10**$viral_diff);
close(TSTFILE);

system("../../Spatial_HSV -r -f $infile -d");
