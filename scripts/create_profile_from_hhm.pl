#!/usr/bin/env perl
# 
# Create a profile (.prf) from a given HHM file

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode

my $help="
Create a profile (.prf) from a given HHM file

Usage: perl create_profile_from_hmmer.pl -i <infile> [-o <outfile>]

Options:
  -i <infile>   Input file in HHM format
  -o <outfile>  Output file in prf-format (default: infile.prf)

  -v [0-5]      verbose mode (default: $v)
\n";

# Variable declarations
my $line;
my $infile;
my $outfile;
my $i;
my $a;

my @counts;     # count profile
my @neffs;
my $name;       # name of HHM profile
my $len;        # length of HHM profile

                 #   A  C  D  E   F  G  H  I   K   L   M  N   P  Q  R   S   T   V   W   Y   
my @hhmaa2csaa = ( 0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18);
my @aminoacids   = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');
my %aa2i;
for ($a = 0; $a < 20; $a++) {
    $aa2i{$aminoacids[$a]} = $a;
}

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

if ($options=~s/ -i\s+(\S+) //) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}

if ($options=~s/ -v\s+(\S+) //) {$v=$1;}

if (!$infile) {print($help); print "ERROR! No input file!\n"; exit(1);}

if (!$outfile) {
    $infile =~ /^(\S+)\.\S+?$/;
    $outfile = "$1.prf";
}

##############################################################################################
# Main part
##############################################################################################

######################################
# Read HHM profile
######################################
open (IN, $infile);
while ($line = <IN>) {
    if ($line =~ /^NAME\s+(\S+)/) {
	$name = $1;
    } elsif ($line =~ /^LENG\s+(\d+)/) {
	$len = $1;
    } elsif ($line =~ /^HMM/) {
	last;
    }
}

$i = 0;
while ($line = <IN>) {
    if ($line =~ /^\/\//) { last; }

    if ($line =~ s/^\S \d+ //) {

	for ($a = 0; $a < 20; $a++) {
	    $line =~ s/^\s*(\S+)\s/ /;
	    $counts[$i][$hhmaa2csaa[$a]] = $1;
	    if ($counts[$i][$hhmaa2csaa[$a]] !~ /\*/ && $counts[$i][$hhmaa2csaa[$a]] == 0) { 
		$counts[$i][$hhmaa2csaa[$a]] = 1; 
	    }
	}

	$line = <IN>;
	$line =~ /^\s*\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+/;
	$neffs[$i] = $1;

	$i++;
    }
}

######################################
# write count_profile
######################################

open (OUT, ">$outfile");
# Write header
printf(OUT "CountProfile\n");
printf(OUT "NAME\t%s\n", $name);
printf(OUT "LENG\t%i\n", $len);
printf(OUT "ALPH\t20\n");              # 20 amino acid alphabet
printf(OUT "COUNTS");
for ($a = 0; $a < 20; $a++) {
    printf(OUT "\t%s", $aminoacids[$a]);
}
printf(OUT "\tNEFF\n");

# Write profile
for ($i = 0; $i < $len; $i++) {
    printf(OUT "%i", $i+1);
    for ($a = 0; $a < 20; $a++) {
	if ($counts[$i][$a] == '*') {
	    printf(OUT "\t*");
	} else {
	    printf(OUT "\t%i", $counts[$i][$a]);
	}
    }
    printf(OUT "\t%i\n", $neffs[$i]);
}

printf(OUT "//\n");
close OUT;

exit;


