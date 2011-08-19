#!/usr/bin/env perl
use strict;

my $ext="seq";
my $usage="
Write a separate file for each FASTA sequence found in file
Usage: splitfasta.pl infile [option]
Option:
-fam       : use family-based name (for SCOP/ASTRAL sequences
-name      : use sequence name as file name (default)
-ext <ext> : extension for sequence files (default=$ext)
\n";

if (@ARGV<1) {die $usage;;}

my $line;
my $infile=$ARGV[0];
my $outfile;
my $sequence="";
my $options="";
my $fam=0;             # option -fam?
my $famid="";
my %numfams=(); 
my $n=0;               # number of name lines read in so far

if (@ARGV>1) {
    $options.=join(" ",@ARGV[1..$#ARGV]);
}

# Set number of cpus to use
if ($options=~s/-fam//g)         {$fam=1;}
if ($options=~s/-name//g)        {$fam=0;}
if ($options=~s/-ext\s+(\S+)//g) {$ext=$1;}


open (INFILE,"<$infile") || die("ERROR: Can't open $infile: $!\n");

if ($fam) {

    while ($line=<INFILE>) {
	if ($line=~/^>(\S+)\s+(\S+)/) {
	    $famid=$2;
	    if ($n) {
		open (OUTFILE,">$outfile") || die("ERROR: Can't open $outfile: $!\n");
		print(OUTFILE $sequence);
		close(OUTFILE);
	    }
	    if (defined $numfams{$fam}) {$numfams{$fam}++;} else {$numfams{$fam}=1};
	    $outfile="$fam.".$numfams{$fam}.".seq";
	    $sequence=$line;
	    $n++;
	} else {
	    $sequence.=$line;
	}
    }
    if ($n) {
	open (OUTFILE,">$outfile") || die("ERROR: Can't open $outfile: $!\n");
	print(OUTFILE $sequence);
	close(OUTFILE);
    }

} else {

    my %exists=();
    while ($line=<INFILE>) {
	if ($line=~/^>(\S+)/) {
	    if ($n) {
		open (OUTFILE,">$outfile") || die("ERROR: Can't open $outfile: $!\n");
		print(OUTFILE $sequence);
		close(OUTFILE);
	    }
	    if ($exists{$1}) {print("\nWarning: id $1 appears more than once in $infile\n");}
	    $exists{$1}=1;
	    my $tmp = $1;
	    $tmp =~ s/\|/_/g;
	    $tmp =~ s/\./_/g;
	    $outfile="$tmp.seq";
	    $sequence=$line;
	    $n++;
	} else {
	    $sequence.=$line;
	}
    }
    if ($n) {
	open (OUTFILE,">$outfile") || die("ERROR: Can't open $outfile: $!\n");
	print(OUTFILE $sequence);
	close(OUTFILE);
    }
}


close(INFILE);
printf("Created %i sequence files\n",$n);



