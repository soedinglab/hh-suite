#! /usr/bin/perl -w
# COG0488 and before -> delete

my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.

my $v=2;
my $usage="
Generates *.a3m and *.hhm files with 'NAME Name_of_domain description'
for each globbed file in cdd database (comprising COGs, KOGs, Pfam, SMART, and CD).
Calls buildali.pl for all files in fileglob.
Usage:   cdd.pl [options] 'fileglob'
Options (in addition to buildali.pl options): 
 -u       update: skip buildali.pl for infile if infile.a3m exists
 -f file  give location of cddid.tbl file (default: ./cddid.tbl)
Example: cdd.pl -n 0 -f cddid.tbl '*.fas'
\n";

if (@ARGV<1) {die("$usage");}

my @infiles=glob($ARGV[@ARGV-1]);
my $infile;
my $cddid="./cddid.tbl";
my $format;
my $update=0;
my $base;
my $root;
my $line;
my %consensusname;  # contains the name line of the consensus sequence of the domain

my $options="";
for (my $i=0; $i<@ARGV-1; $i++) {$options.=" $ARGV[$i] ";}
if ($options=~s/ -f\s+(\S+) //) {$cddid=$1;}
print("$options\n");

printf("Reading table of identifiers, names, and descriptions in $cddid\n");
my $COGs=0;
my $KOGs=0;
my $Pfams=0;
my $smarts=0;
my $cds=0;

# Read consensus name from cddid.tbl
open(CDDFILE,"<$cddid") || die("Error: cannot open $cddid for reading: $!\n");
while($line=<CDDFILE>) {
    if      ($line=~/^\S+\s+(COG\d+)\s+(\S+)\s+(.*?)\s*\d*$/) {
    } elsif ($line=~/^\S+\s+(KOG\d+)\s+(\S+)\s+\S+\s+(.*?)\s*\d*$/) {
    } elsif ($line=~/^\S+\s+(pfam\d+)\s+(\S+)\s+(.*?)\s*\d*$/) {
    } elsif ($line=~/^\S+\s+(smart\d+)\s+(\S+)\s+(.*?)\s*\d*$/) {
    } elsif ($line=~/^\S+\s+(cd\d+)\s+(\S+)\s+(.*?)\s*\d*$/) {
    } elsif ($line=~/^\S+\s+(\S+)\s+(\S+)\s+(.*?)\s*\d*$/) {
    }
    if ($1 eq $2) {
	$consensusname{$1}="$2 $3";
    } else {
	$consensusname{$1}="$1 $2 $3";
    }
#    print("id=$1 name=$consensusname{$1}\n");
}
close(CDDFILE);

# Replace name of first sequence with consensus name from cddid.tbl
printf("%i files to process\n",scalar(@infiles));
foreach $infile (@infiles) {
    if ($infile=~/^(.*)\..*$/) {$base=$1;} else {$base=$infile;}
    if ($base=~/^.*\/(.*)$/) {$root=$1;} else {$root=$base;}
    if ($update && -e $base.".a3m") {next;}

    # Replace name of consensus sequence with '>Name_of_domain description' line
    open(INFILE,"<$infile") || die("Error: cannot open $infile for reading: $!\n");
    my @lines=<INFILE>;
    close(INFILE);
    my $nameline=shift(@lines);
    if (defined $consensusname{$root}) {$nameline=sprintf(">%s\n",$consensusname{$root});}
    else {
	$nameline=~s/>lcl\|consensus\s*(KOG\d+),?/>$1/;
    }
    open(OUTFILE,">$base.tmp") || die("Error: cannot open $infile for writing: $!\n");
    printf(OUTFILE "%s",$nameline);
    foreach $line (@lines) {print(OUTFILE $line);}
    close(OUTFILE);

    # Call buildali.pl with ncbi-determined consensus sequence as query => has gaps for columns with >50% gaps!
    print("rootname=$root\n");
    if ($root=~/^COG/ || $root=~/^KOG/) {
	&System("$hh/buildali.pl -v 2 -n 0 -fas $options $base.tmp ");
    } else {
	&System("$hh/buildali.pl -v 2 -e 1e-3 -core -n 1 -fas $options $base.tmp ");
    }
    &System("$hh/hhmake -i $base.a3m -id 90");
    unlink("$base.tmp");
} # end foreach $infile
printf("Finished $0 %s\n",join(" ",@ARGV));
exit;

################################################################################################
### System command
################################################################################################
sub System()
{
    if ($v>=2) {printf("%s\n",$_[0]);} 
    return system($_[0])/256;
}

