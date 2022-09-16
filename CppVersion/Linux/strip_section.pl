#!/usr/bin/perl
if ($#ARGV < 1) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 2 input arguments ($args given)\n";
}
my $tag;
my $infile;

if ( "$ARGV[0]" ne "") {
    $tag=$ARGV[0]
} else {
    die "1st arg is tag\n";
    exit (1);
}
printf STDERR ("tag is %s\n",$tag);
if ( -f $ARGV[1]) {
    $infile=$ARGV[1]
} else {
    die "2nd arg is the input file\n";
    exit (1);
}

open(INF,"<$infile");

my $nesting=0;
my $in_ifdef=0;
my $skip_ifdef=0;
my $line_cnt=0;
my $out_cnt=0;
my $ifdef_cnt=0;
if ($line=<INF>) {
    chomp($line);              # remove the newline from $line.
    $line_cnt++;
    my @pieces = split(/ /,$line);

    if ($pieces[0] ne "//NOTE:") {
	printf ("//NOTE: This read-only file is generated from the file %s.\n//Any editing should be done in that file.\n\n",$infile);
    }
    printf("%s\n",$line);
}

while( $line=<INF>) {
    chomp($line);              # remove the newline from $line.

    $line_cnt++;
    my @pieces = split(/ /,$line);

    if ($pieces[0] eq "#ifdef" && $pieces[1] =~ /$tag/) {
	printf STDERR ("Encountered tag %s on line %d\n",$tag,$line_cnt);
	$nesting=1;
	$in_ifdef=1;
	$ifdef_cnt++;
    }
    elsif ($pieces[0] eq "#ifdef" && $in_ifdef==1) {
	printf STDERR ("Encountered nested ifdef (line %d,level=%d)\n",$line_cnt,$nesting);
	$nesting++;
    }
    elsif ($pieces[0] =~ /#else/ && $nesting==1) {
	printf STDERR ("Entered else for %s on line %d\n",$tag,$line_cnt);
	$in_ifdef=0;
	$skip_ifdef=1;
    }
    elsif ($pieces[0] =~ /#endif/ && $nesting > 0) {
	$nesting--;
	if ($nesting==0) {
	    printf STDERR ("Left tagged section %s on line %d\n",$tag,$line_cnt);
	    $in_ifdef=0;
	    $skip_ifdef = 0;
	} else {
	    printf STDERR ("Left nested ifdef (line %d,new level=%d)\n",$line_cnt,$nesting);
	}
    } elsif ($nesting<=0 || $skip_ifdef==1) {
	printf("%s\n",$line);
	$out_cnt++;
    }
}
close(INF);
printf STDERR ("Read %d lines and wrote %d lines (%d ifdefs)\n",$line_cnt,$out_cnt,$ifdef_cnt);
