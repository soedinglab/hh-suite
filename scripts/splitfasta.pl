#! /usr/bin/env perl
# splitfasta.pl 
# Split a file with multiple, FASTA formatted sequences into many single-sequence FASTA files
#
# (C) Johannes Soeding, 2012
#
#     HHsuite version 2.0.14 (May 2012)
#
#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

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
use warnings;

my $ext="seq";
my $usage="
splitfasta.pl from HHsuite $VERSION  
Split a file with multiple, FASTA formatted sequences into multiple single-sequence FASTA files.
Write files into current directory and name each file by the first word after \">\" in the name line. 

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



