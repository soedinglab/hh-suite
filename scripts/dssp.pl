#! /usr/bin/env perl
#
# For each sequence found in infile, read dssp file and write secondary structure
# assignments for each residue into a separate file in fasta format.
# Make sure the amino acids of the dssp file match those of the infile sequence and
# warn if an assignment with <=6 mismatches cannot be done.
# Usage: dssp.pl infile [dssp-dir]


use strict;
$|= 1; # Activate autoflushing on STDOUT

if (scalar(@ARGV)<1)
{
    print("\nFor each sequence found in infile, read dssp file and write secondary structure assignments \n");
    print("for each residue into a separate file in fasta format. Make sure the amino acids of the dssp file\n");
    print("match those of the in file and warn if an assignment with <=6 mismatches cannot be done.\n");
    print("Usage: dssp.pl infile [dssp-dir]\n\n");
    exit;
}

my $infile = $ARGV[0];
my $dsspdir = "/raid/db/dssp/data";
my $dssp="/raid/db/dssp/bin/dsspcmbi";
my $pdbdir="/raid/db/pdb";

if (scalar(@ARGV)==2) {$dsspdir = $ARGV[1];}

my $i;           #counts residues
my $line;        #input line
my $name;        #name of sequence in in file, e.g. d1g8ma1
my %longname;    #full name of sequence incl. description
my %range;       #chain and residue range of in sequence in pdb file
my %aas;         #amino acids from in file for each $name

my $dsspfile;
my $maxsproutfile;
my $nfile=0;     #number of files read in
my $pdbcode;     #pdb code for accessing dssp file; shortened from in code, e.g. 1g8m
my @ss_dssp;     #dssp states for residues (H,E,L)
my @aa_dssp;     #residues read from dssp file
my @aa_astr;     #residues in infile
my $aa_astr;     #residues from infile as string
my $ss_dssp;     #dssp states as string
my $aa_dssp;     #residues from dssp file as string
my $length;      #length of sequence

my $aa;          #current amino acid in dssp file
my $ss;          #current secondary structure state in dssp file
my $chain;       #two-letter code for chain in dssp file
my $range;       #range description for $name in dssp file
my $thisres;     #index of present residue in dssp file
my $lastres;     #index of last residue in dssp file
my $contained;   #is the current line in dssp file contained in the chain and residue ranges of $name?
my $gap;         #did last line in dssp file indicate missng residues? ->1 else 0
my $error;       #1: wrong assignment of residues dssp<->in
my $errinarow;   #how many errors have been made in a row?
my $assigned;    #number of correctly assigned residues
my $skipchain;   #found too many errors in this chain; try next chain with same 2nd letter (eg. PA SA A)
my $thischain;

# Read infile and prepare hashes %first, %last, %chain, $longname
open (INFILE, "<$infile") || die ("cannot open $infile: $!");
while ($line=<INFILE>)
{
    if ($line=~/>(\S+)\s+\S+\s+\((\S+)\)/)
    {
	$name=$1;
	chomp($line);
	$range{$name}=$2;
	$line=~/^>(.*)/;
	$longname{$name}=$1;
    }
    else
    {
	chomp($line);
	$line=~tr/a-zX \t/A-Z/d;
	$aas{$name}.=$line;
    }
}
close(INFILE);


# Read dssp files one after another
foreach $name (keys(%range))
{
    # prepare in sequence
    $aa_astr = $aas{$name};
    @aa_astr = split(//,$aa_astr);
    $length=length($aa_astr);

    # Read dssp file into @ss_dssp
    @ss_dssp={}; @aa_dssp={};
    for ($i=0; $i<$length; $i++) {$ss_dssp[$i]=" "; $aa_dssp[$i]=" ";}
    $name=~/[a-z](\S{4})/;
    $pdbcode=$1;
    $dsspfile="$dsspdir/$pdbcode.dssp";
    if (! open (DSSPFILE, "<$dsspfile"))
    {
	printf(STDERR "Cannot open $dsspfile: $!\n");
	$maxsproutfile="$pdbdir/$pdbcode.brk_ca_mod";
	if (! -e $maxsproutfile)
	{
	    printf(STDERR "ERROR: could open $maxsproutfile either: $!  Skipping $name\n");
	    next;
	}
	else
	{
	    &System("$dssp $maxsproutfile $dsspfile > dssp.out 2>&1");
	    if (! open (DSSPFILE, "<$dsspfile"))
	    {
		printf(STDERR "ERROR: dssp couldn't generate file from maxsprout-file. Skipping $name\n");
		next;
	    }
	}
    }
    while ($line=<DSSPFILE>) {if ($line=~/^\s*\#\s*RESIDUE\s+AA/) {last;}}

    $gap=0; $i=-1; $error=0; $assigned=0; $skipchain=""; $errinarow=0;
    while ($line=<DSSPFILE>)
    {
	if ($line=~/^.{5}(.{5})(.)(.)\s(.).\s(.)/)
	{
	    $thisres=$1;
	    $thisres=~tr/ //d;
	    if ($2.$3 eq $skipchain) {next;}
	    $skipchain="";
	    $thischain=$2.$3;
	    $chain=$3;
	    $chain=~tr/ //d;
	    $aa=$4;
	    $ss=$5;
	    if ($aa eq "!")  {$gap=1; next;}    # missing residues!
	    $range=$range{$name};  $contained=0;

	    # Check whether at least one range descriptor fits, e.g. (B:1-127,B:254-362)
	    do{
		if    ($range=~s/^(\S):(-?\d+)\w?-(\d+)(\w?)//      #syntax (A:56S-135S) or (R:56-135)
		    && $chain eq "$1" && $2<=$thisres && $thisres<=$3) {$contained=1;}
		elsif ($range=~s/^(-?\d+)\w?-(\d+)(\w?)//           #syntax (56-135)
		    && $chain eq "" && $1<=$thisres && $thisres<=$2) {$contained=1;}
		elsif ($range=~s/^(\S):// && $chain eq "$1") {$contained=1;} #syntax (A:) or (A:,2:)
		elsif ($range=~s/^-$// && $chain eq "") {$contained=1;}      #syntax (-)
		$range=~s/^,//;
	    } while($contained==0 && $range ne "");

	    # line not contained in specified range for $name?
	    if (!$contained) {$lastres=$thisres; $gap=0; next;}

	    if (!exists $aa_astr[$i+1]) {last;}
	    if($gap && $i!=-1)
	    {
		while ($aa_astr[$i+1] eq "X") {$i++;}
		if ((exists($aa_astr[$i+$thisres-$lastres]))
		    && $aa_astr[$i+$thisres-$lastres] eq $aa && $aa_astr[$i+1] ne $aa)
		{
		    $i+=$thisres-$lastres;
		    if (!exists $aa_astr[$i]) {$error=1000; $i=0;}
		}
		else
		{
		    $i++;
		    if (!exists $aa_astr[$i]) {$error=1000; $i=0;}
		}
	    }
	    else {$i++;}

	    $aa=~tr/a-z/CCCCCCCCCCCCCCCCCCCCCCCCCC/;

	    # Count wrong assignments
	    if ($error<1000 && $aa ne $aa_astr[$i] && $aa ne "X")
	    {
		if ($errinarow>=1)
		{
		    $error++; $errinarow++;
		    my $j=$i;
		    while (exists $aa_astr[$i] && $aa_astr[$i] ne $aa) {$i++;}
		    if (!exists $aa_astr[$i]) {$i=$j;}
		}
		else {$errinarow++; $error++;}

	    }
	    else {$errinarow=0;}


	    if ($error>6)
	    {
#		print("Errors:$error   i:$i   aa_astr:$aa_astr[$i]   aa:$aa\n");
#		print("$longname{$name}\n");
		$aa_dssp= join("",@aa_dssp);
		$ss_dssp= join("",@ss_dssp);
#		print("DSSP   '$aa_dssp'\n");
#		print("In '$aa_astr'\n");
		for ($i=0; $i<$length; $i++) {$ss_dssp[$i]=" "; $aa_dssp[$i]=" ";}
		$skipchain=$thischain;
		$error=0; $assigned=0; $i=-1; $gap=0; $errinarow=0;
		next;
	    }
	    $assigned++;
	    if    ($ss eq " ") {$ss="~";}
	    elsif ($ss eq "I") {$ss="~";}
#	    elsif ($ss eq "S") {$ss="~";}
#	    elsif ($ss eq "T") {$ss="~";}
#	    elsif ($ss eq "G") {$ss="~";}
#	    elsif ($ss eq "B") {$ss="~";}
	    $ss_dssp[$i]=$ss;
	    $aa_dssp[$i]=$aa;
	    $lastres=$thisres;
	    $gap=0;
	}
    }
    close(DSSPFILE);
    $nfile++;

    if ($skipchain || $error>6)
    {printf (STDERR "ERROR in $name ($nfile): too many wrong assingments of residues\n"); next;}
    elsif ($assigned<0.5*$length)
    {printf (STDERR "ERROR in $name ($nfile): too few correct assignments of residues\n"); next;}

#    printf(STDERR "%4i: %-100.100s  ",$nfile,$longname{$name});
#    printf(STDERR "%1i errors, %3i%% assigned\n",$error,int(100*$assigned/$length));

    $ss_dssp= join("",@ss_dssp);
    $aa_dssp= join("",@aa_dssp);
    $ss_dssp=~s/ \S /   /g;
    $ss_dssp=~s/ \S\S /    /g;
    $aa_dssp=~tr/ /-/;
    $ss_dssp=~tr/ /-/;

    if (!open (SSFILE,">$name.ss"))
    {printf(STDERR "ERROR: cannot open $name.ss: $!\n"); next;}
    printf(SSFILE ">ss_dssp $longname{$name}\n$ss_dssp\n");
    printf(SSFILE ">aa_dssp $longname{$name}\n$aa_dssp\n");
    printf(SSFILE ">aa_astral $longname{$name}\n$aa_astr\n\n");
    close (SSFILE);

}


exit;


sub System()
{
    my $command=$_[0];
    print("\nCalling '$command'\n");
    return system($command)/256; # && die("\nERROR: $!\n\n");
}
