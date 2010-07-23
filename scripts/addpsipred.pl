#!/usr/bin/perl
# addpsipred.pl version 1.1.4 (April 2005)
# Add PSIPRED secondary structure prediction to a FASTA of A3M alignment.
# Output format is A3M (see HHsearch README file).
#
# Please report bugs to johannes@soeding.com. Thank you.

# Delete the following 6 lines and set the variables in the next paragraph to your blast etc. directories
my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.

#my $bioprogs_dir="/home/soeding/programs";   # see next two lines
#my $ncbidir="$bioprogs_dir/blast";           # Put the directory path with the BLAST executables 
#my $perl="/home/soeding/perl";               # Put the directory path where reformat.pl is lying
#my $dummydb="/home/soeding/nr/do_no_delete"; # Put the name given to the dummy blast directory (or leave this name)

my $psipreddir="$bioprogs_dir/psipred/";     # Put the directory path with the PSIPRED executables 
my $execdir=$psipreddir."/bin";
my $datadir=$psipreddir."/data";
my $ss_cit="PSIPRED: Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. JMB 292:195-202.";

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=5;              # verbose mode

my $numres=100;        # number of residues per line for secondary  structure
my $informat="fas";    # input format

my $help="
Add PSIPRED secondary structure prediction to a multiple sequence alignment or HMMER file.
Input is a  multiple sequence alignment or a HMMER (multi-)model file. Allowed input formats are 
A2M/FASTA (default), A3M (-a3m), CLUSTAL (-clu), STOCKHOLM (-sto), HMMER (-hmm).
If the input file is an alignment, the output file is in A3M with default name <basename>.a3m.
If the input file is in HMMER format, the output is the same as the input, except that records SSPRD 
and SSCON are added to each model which contain predicted secondary structure and confidence values. 
In this case the output file name is obligatory and must be different from the input file name.
(( Remark: A3M looks misaligned but it is not. To reconvert to FASTA, type ))
((   'reformat.pl file.a3m file.fas'.                                      ))
(( For an explanation of the A3M format, see the HHsearch README file.     ))

Usage: perl addpsipred.pl <ali file> [<outfile>] [-fas|-a3m|-clu|sto]  
  or   perl addpsipred.pl <HMM file> <outfile> -hmm  
\n";

# Variable declarations
my $line;
my @seqs;              # sequences from infile (except >aa_ and >ss_pred sequences)
my $query_length;
my $qseq;              # residues of query sequence
my $name;              # query in fasta format: '>$name [^\n]*\n$qseq'
my $infile;
my $outfile;
my $ss_pred="";        # psipred ss states
my $ss_conf="";        # psipred confidence values

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

#Input format fasta?
if    ($options=~s/ -fas\s/ /g) {$informat="fas";}
elsif ($options=~s/ -a2m\s/ /g) {$informat="a2m";}
elsif ($options=~s/ -a3m\s/ /g) {$informat="a3m";}
elsif ($options=~s/ -clu\s/ /g) {$informat="clu";}
elsif ($options=~s/ -sto\s/ /g) {$informat="sto";}
elsif ($options=~s/ -hmm\s/ /g) {$informat="hmm";}

if ($options=~s/ -v\s+(\S+) //) {$v=$1;}

# Set input and output file
if ($options=~s/ -i\s+(\S+) //) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$infile=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$outfile=$1;}

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!$infile) {print($help); exit(1);}

my $v2 = $v-1;
if ($v2>2) {$v2-=2;}
if ($v2<0) {$v2=0;}


###############################################################################################
# Reformat input alignment to a3m and psiblast-readable format and generate file with query sequence
###############################################################################################

my $inbase; # $inbasename of infile: remove extension
my $inroot; # $inbasename of infile: remove path and extension
if ($infile=~/(.*)\..*/) {$inbase=$1;} else {$inbase=$infile;}  # remove extension
if ($inbase=~/.*\/(.*)/)  {$inroot=$1;} else {$inroot=$inbase;} # remove path 

############################################################################################
open (INFILE, "<$infile") || die("ERROR: cannot open '$infile': $!\n");
$line=<INFILE>;
if ($line=~/^HMMER/) {$informat="hmm";}
elsif ($line=~/^>/ || $line=~/^\#/) {}
else {die("\nError: wrong format! Only HMMER and FASTA/A2M/A3M format allowed\n");}
close(INFILE);
if ($informat eq "hmm" && !$outfile) {
    print("Error: no output file given. With the -hmm option an output file is obligatory\n"); exit(1);
}

if ($informat ne "hmm")
{
    if (!$outfile) {$outfile="$inbase.a3m";}

    # Use first sequence to define match states and reformat input file to a3m and psi
    &System("perl $perl/reformat.pl -v $v2 -M first $informat a3m $infile $inbase.in.a3m");
        
    # Read query sequence
    open (INFILE, "<$inbase.in.a3m") or die ("ERROR: cannot open $inbase.in.a3m: $!\n");
    $/=">"; # set input field separator
    my $i=0;
    $qseq="";
    while ($line=<INFILE>) {
	if ($line eq ">") {next;}
	$line=~s/>$//;
	if ($line!~/^ss_dssp/ && ($line=~/^ss_/ || $line=~/^aa_/)) {next;}
	$seqs[$i++]=">$line";
	if(!$qseq && $line!~/^ss_dssp/) {
	    $line=~s/^(\S*)[^\n]*//;
	    $name=$1;
	    $qseq=uc($line);
	    $qseq=~s/\n//g;
	}
    }
    
    close(INFILE);
    $query_length = ($qseq=~tr/A-Z/A-Z/);
    $qseq=~tr/a-zA-Z//cd;
    
    # If less than 26 match states => add sufficient number of Xs to the end of each sequence in $inbase.in.a3m
    my $q_match = ($qseq=~tr/A-Z/A-Z/); # count number of capital letters
    if ($q_match<=25) {                 # Psiblast needs at least 26 residues in query
	my $addedXs=('X' x (26-$q_match))."\n";
	$qseq.=$addedXs;     # add 'X' to query to make it at least 26 resiudes long
	for ($i=0; $i<@seqs; $i++) {	    
	    $seqs[$i]=~s/\n$//g;
	    $seqs[$i].=$addedXs;
	}
	open (INFILE,">$inbase.in.a3m");
	for ($i=0; $i<@seqs; $i++) {
	    printf(INFILE "%s",$seqs[$i]);
	}
	close INFILE;
    }
    $/="\n"; # set input field separator

    # Write query sequence file in FASTA format
    open (QFILE, ">$inbase.sq") or die("ERROR: can't open $inbase.sq: $!\n");
    printf(QFILE ">%s\n%s\n",$name,$qseq);
    close (QFILE);
    
    # Reformat into PSI-BLAST readable file for jumpstarting 
    &System("perl $perl/reformat.pl -v $v2 -r -noss a3m psi $inbase.in.a3m $inbase.psi");
    
    
    # Secondary structure prediction with psipred
    if ($v>=1) {printf ("\nPredicting secondary structure with PSIPRED ...\n");}
    &RunPsipred("$inbase.sq");
    
    
    open (ALIFILE, ">$outfile") || die("ERROR: cannot open $inbase.a3m: $!\n");
    if (open (PSIPREDFILE, "<$inbase.horiz")) {
	$ss_conf="";
	$ss_pred="";
	# Read Psipred file
	while ($line=<PSIPREDFILE>) {
	    if    ($line=~/^Conf:\s+(\S+)/) {$ss_conf.=$1;}
	    elsif ($line=~/^Pred:\s+(\S+)/) {$ss_pred.=$1;}
	}
	close(PSIPREDFILE);
	$ss_conf=~tr/0-9/0/c; # replace all non-numerical symbols with a 0
	$ss_pred=~s/(\S{$numres})/$1\n/g;
	$ss_conf=~s/(\S{$numres})/$1\n/g;
	print(ALIFILE ">ss_pred PSIPRED predicted secondary structure\n$ss_pred\n");
	print(ALIFILE ">ss_conf PSIPRED confidence values\n$ss_conf\n");
    }
    
    # Append alignment sequences to psipred sequences
    for ($i=0; $i<@seqs; $i++) {
	print(ALIFILE $seqs[$i]);
    }
    close(ALIFILE);
} 

############################################################################################
# $informat eq "hmm"
else 
{
    my @logoddsmat;
    my @lines;
    my $length;
    my $query;
    my $scale=0.13; # empirically determined scale factor between HMMER bit score and PSI-BLAST score, 0.3 for HMMER3
    my $acc;
    my $name;
    my $desc;
    my $nmodels=0;

    open (INFILE, "<$infile") || die("ERROR: cannot open $infile: $!\n");
    open (OUTFILE, ">$outfile") || die("ERROR: cannot open $outfile: $!\n");

    # Read HMMER file model by model
    while ($line=<INFILE>) {
	# Search for start of next model
	while ($line && $line!~/^HMMER/ && $line!~/^NAME /) {
	    $line=<INFILE>; 
	}
	if ($line=~/^HMMER3/) {

	    $scale = 0.3;
	    @logoddsmat=();
	    @lines=($line); 
	    
	    while ($line=<INFILE>) {push(@lines,$line); if ($line=~/^LENG/) {last;}}
	    $line=~/^LENG\s+(\d+)/;
	    $length=$1;  # number of match states in HMM
	    $query="";   # query residues from NULL emission lines
	    while ($line=<INFILE>) {push(@lines,$line); if ($line=~/^\s*m->m/) {last;}}
	    push(@lines,$line=<INFILE>);
	    if ($line !~ /^\s*COMPO/) {
		die("Error: need null-model probablities (Parameter COMPO)!\n");
	    }
	    $line=~s/^\s*COMPO\s+(\S.*\S)\s*$/$1/;
	    my @nullmodel = split(/\s+/,$line);
	    @nullmodel = map {$_ = exp(-1 * $_)} @nullmodel;  # Transform to probabilities
	    push(@lines,$line=<INFILE>); # state 0 insert emission
	    push(@lines,$line=<INFILE>); # transisitions from begin state
	    
	    while ($line=<INFILE>) {
		push(@lines,$line); 
		if ($line=~/^\/\//) {last;}
		$line=~s/^\s*\d+\s+(\S.*\S)\s+\d+\s+(\S)\s+\S\s*$/$1/;
		$query .= $2;
		my @probs = split(/\s+/,$line);
		@probs = map {$_ = exp(-1 * $_)} @probs;  # Transform to probabilities
		# calculate log-odds
		my @logodds = ();
		for (my $a = 0; $a < scalar(@probs); $a++) {
		    my $logodd = (log($probs[$a] / $nullmodel[$a]) / $log2) * 1000;
		    push(@logodds, $logodd);
		}
		
		push(@logoddsmat,\@logodds);
		push(@lines,$line=<INFILE>);
		push(@lines,$line=<INFILE>);
	    }

	} else {

	    $scale=0.13;
	    if ($line!~/^HMMER/ && $line!~/^NAME /) {last;}  # first line in each model must begin with 'HMMER...'
	    @logoddsmat=();
	    @lines=($line); 
	    
	    while ($line=<INFILE>) {push(@lines,$line); if ($line=~/^LENG/) {last;}}
	    $line=~/^LENG\s+(\d+)/;
	    $length=$1;  # number of match states in HMM
	    $query="";   # query residues from NULL emission lines
	    while ($line=<INFILE>) {push(@lines,$line); if ($line=~/^\s*m->m/) {last;}}
	    push(@lines,$line=<INFILE>);
	    while ($line=<INFILE>) {
		push(@lines,$line); 
		if ($line=~/^\/\//) {last;}
		$line=~s/^\s*\d+\s+(\S.*\S)\s*$/$1/;
		my @logodds = split(/\s+/,$line);
		push(@logoddsmat,\@logodds);
		push(@lines,$line=<INFILE>);
		$line=~/^\s*(\S)/;
		$query .= $1;
		push(@lines,$line=<INFILE>);
	    }
	}
	    
	# Write mtx matrix
	open (MTXFILE, ">$inbase.mtx") || die("ERROR: cannot open $inbase.mtx: $!\n");
	printf(MTXFILE "%i\n",$length);
	printf(MTXFILE "%s\n",$query);
	printf(MTXFILE "2.670000e-03\n4.100000e-02\n-3.194183e+00\n1.400000e-01\n2.670000e-03\n4.420198e-02\n-3.118986e+00\n1.400000e-01\n3.176060e-03\n1.339561e-01\n-2.010243e+00\n4.012145e-01\n");
	while (@logoddsmat) {
	    my @logodds = @{shift(@logoddsmat)};
	    print(MTXFILE "-32768 ");
	    splice(@logodds, 1,0,-32768/$scale);   # insert logodds value for B
	    splice(@logodds,20,0,  -100/$scale);   # insert logodds value for X
	    splice(@logodds,22,0,-32768/$scale);   # insert logodds value for Z
	    for (my $i=0; $i<23; $i++) {
		printf(MTXFILE "%4.0f ",$scale*$logodds[$i]);
	    }
	    print(MTXFILE "-32768 -400\n");
	}
	close(MTXFILE);
	
	# Call PSIPRED
	&System("$execdir/psipred $inbase.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $inbase.ss");
	
	# READ PSIPRED file
	if (open (PSIPRED, "$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $inbase.ss2 $inbase.ss |")) {
	    $ss_conf="";
	    $ss_pred="";
	    # Read Psipred file
	    while ($line=<PSIPRED>) {
		if    ($line=~/^Conf:\s+(\d+)/) {$ss_conf.=$1;}
		elsif ($line=~/^Pred:\s+(\S+)/) {$ss_pred.=$1;}
	    }
	    close(PSIPREDFILE);
	}
	
	# Add secondary structure to HMMER output file and print
	foreach $line (@lines) {
	    if ($line=~/^SSPRD/ || $line=~/^SSCON/|| $line=~/^SSCIT/) {next;}
	    if ($line=~/^HMM /) {
		$ss_pred=~s/(\S{$numres})/$1\nSSPRD /g;
		$ss_conf=~s/(\S{$numres})/$1\nSSCON /g;
		printf(OUTFILE "SSCIT HHsearch-readable PSIPRED secondary structure prediction (http://protevo.eb.tuebingen.mpg.de/hhpred/)\n");
		printf(OUTFILE "SSPRD %s\n",$ss_pred);
		printf(OUTFILE "SSCON %s\n",$ss_conf);
		printf(OUTFILE "SSCIT %s\n",$ss_cit);
	    }
	    printf(OUTFILE $line);
	}
	$nmodels++;
    }
	
    close(OUTFILE);
    close(INFILE);
    System("rm $inbase.pn $inbase.sq $inbase.sn $inbase.mn $inbase.chk $inbase.blalog $inbase.mtx $inbase.aux $inbase.ss $inbase.ss2");
    if ($v>=2) {printf("Added PSIPRED secondary structure to %i models\n",$nmodels);}
}

# Reformat output alignment to FASTA
#&System("perl $perl/reformat.pl -v 3 $inbase.a3m $inbase.fas");

if ($v<=2) {
    unlink("$inbase.in.a3m");
    unlink("$inbase.psi");
    unlink("$inbase.horiz");
} 
exit;

##############################################################################################
# Run SS prediction starting from alignment in $inbase.psi (called by BuildAlignment)
##############################################################################################
sub RunPsipred() {
    # This is a simple script which will carry out all of the basic steps
    # required to make a PSIPRED V2 prediction. Note that it assumes that the
    # following programs are in the appropriate directories:
    # blastpgp - PSIBLAST executable (from NCBI toolkit)
    # makemat - IMPALA utility (from NCBI toolkit)
    # psipred - PSIPRED V2 program
    # psipass2 - PSIPRED V2 program
    
    my $infile=$_[0];
    my $basename;  #file name without extension
    my $rootname;  #basename without directory path
    if ($infile =~/^(.*)\..*?$/)  {$basename=$1;} else {$basename=$infile;}
    if ($basename=~/^.*\/(.*?)$/) {$rootname=$1;} else {$rootname=$basename;}

    # Does dummy database exist?
    if (!-e "$dummydb.phr") {
	&System("cp $infile $dummydb");
	&System("$ncbidir/formatdb -i $dummydb");
    }

    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
    &System("$ncbidir/blastpgp -b 1 -j 1 -h 0.001 -d $dummydb -i $infile -B $inbase.psi -C $inbase.chk 1> $inbase.blalog 2> $inbase.blalog");
    
#    print("Predicting secondary structure...\n");
    
    system("echo $inroot.chk > $inbase.pn\n");
    system("echo $inroot.sq > $inbase.sn\n");
    system("$ncbidir/makemat -P $inbase");
    
    &System("$execdir/psipred $inbase.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $inbase.ss");

    &System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $inbase.ss2 $inbase.ss > $inbase.horiz");
    
    # Remove temporary files
    unlink(split ' ', "$inbase.pn $inbase.sn $inbase.mn $inbase.chk $inbase.blalog $inbase.mtx $inbase.aux $inbase.ss $inbase.ss2 $inbase.sq");
    return;
}


################################################################################################
### System command
################################################################################################
sub System()
{
    if ($v>=2) {printf("%s\n",$_[0]);} 
    return system($_[0])/256;
}


