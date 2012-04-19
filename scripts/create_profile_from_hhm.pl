#!/usr/bin/env perl
# 
# create_profile_from_hhm.pl 
# Create a profile (.prf) from a given HHM file

#     HHsuite version 2.0.15 (April 2012)
#
#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#     (C) Michael Remmert and Johannes Soeding, 2012

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use strict;

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode

my $help="
create_profile_from_hhm.pl from HHsuite $VERSION  
Create a profile (.prf) from a given HHM file

Usage: perl create_profile_from_hhm.pl -i <infile> [-o <outfile>]

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


