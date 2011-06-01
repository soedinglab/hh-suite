#! /usr/bin/perl -w
#
# build_a3m.pl
# Build an A3M-file with HHblits
#

use strict;

# Default values:
my $v=2;            # verbose mode
my $cpu=1;
my $mact=0.5;
my $neff_thresh=5;
my $neff = 0;

my $addss = 1;      # default: add secondary structure

my $hh = "/cluster/bioprogs/hhblits";

my $infile = "";
my $outfile = "";

my $cmd = "";
my $line;

my $usage="
Build A3M-file with HHblits

Usage: build_a3m.pl -i <INFILE> [options] 

 -i    <FILE>   input file with query in FASTA-format

General options:

 -o    <FILE>   output A3M-file (default: infile.a3m)

 -mact <float>  mact-value for HHblits (default: $mact)
 -neff <float>  min. diversity of A3Ms (default: $neff_thresh)

 -noss          don't add secondary structure

 -cpu  <int>    number of cores (default: $cpu)
 -v    <int>    verbose mode (default: $v)
  
\n";

if (@ARGV<1) {die ($usage);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

# General options
if ($options=~s/ -v\s*(\d) / /g) {$v=$1;}
if ($options=~s/ -v / /g) {$v=2;}

if ($options=~s/ -o\s*(\S+)/ /g) {$outfile=$1;}
if ($options=~s/ -i\s*(\S+)/ /g) {$infile=$1;}

if ($options=~s/ -mact\s*(\S+)/ /g) {$mact=$1;}
if ($options=~s/ -neff\s*(\S+)/ /g) {$neff_thresh=$1;}

if ($options=~s/ -noss\s*/ /g) {$addss=0;}

if ($options=~s/ -cpu\s*(\S+)/ /g) {$cpu=$1;}

# Check options
if ($infile eq "") {
    die "ERROR! No infile given!\n";
}
if ($outfile eq "") {
    $infile =~ /^(\S+)\.\S+?$/;
    $outfile = "$1.a3m";
}

my $tmpdir="/dev/shm/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 


# MAIN part
############


# run 2 rounds HHblits
$cmd = "$hh/hhblits -i $infile -o /dev/null -oa3m $tmpdir/HHblits.a3m -cpu $cpu -mact $mact -n 2";
if (&System($cmd) != 0) {
    die ("ERROR with command $cmd!\n");
}

# check NEFF
$cmd = "$hh/hhmake -i $tmpdir/HHblits.a3m";
open (IN, "$cmd |");
while ($line = <IN>) {
    if ($line =~ /Effective number of sequences exp\(entropy\) =\s+(\S+)/) {
	$neff = $1;
	last;
    }
}
close IN;

# if diversity < neff_thresh, perform upto 4 additional rounds
if ($neff < $neff_thresh) {
    print "Diversity $neff is below threshold ($neff_thresh). Perform additional rounds...\n";
    $cmd = "$hh/hhblits -i $tmpdir/HHblits.a3m -o /dev/null -oa3m $outfile -cpu $cpu -mact $mact -n 4 -neffmax $neff_thresh";
    if (&System($cmd) != 0) {
	die ("ERROR with command $cmd!\n");
    }
} else {
    system("cp $tmpdir/HHblits.a3m $outfile");
}

if ($addss) {
    $cmd = "$hh/addss.pl $outfile");
    if (&System($cmd) != 0) {
	die ("ERROR with command $cmd!\n");
    }
}

if ($v < 4) {
    system("rm -rf $tmpdir");
}

exit;

sub System() {
    if ($v>=2) {print("Command: $_[0]\n");} 
    return system($_[0]);
}
