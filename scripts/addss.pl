#!/usr/bin/env perl
# addss.pl version 1.0.0 (October 2009)
# Add DSSP states (if available) and PSIPRED secondary structure prediction to a FASTA or A3M alignment or HMMER file.
# Output format is A3M (for input alignments) or HMMER (see User Guide).

#########################################################################################################
# Delete the following 8 lines and set the variables in the next paragraph to your blast etc. directories
my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.
#########################################################################################################

#my $bioprogs_dir="/cluster/bioprogs";          # see next two lines
#my $ncbidir="$bioprogs_dir/blast";             # Put the directory path with the BLAST executables 
#my $hh="$bioprogs/hh";                         # Put the directory path with hhfilter and hhmake
#my $perl="$bioprogs/perl";                     # Put the directory path where reformat.pl is lying
#my $dummydb="/cluster/databases/do_no_delete"; # Put the name given to the dummy blast directory (or leave this name)

# The following variables have to set for using DSSP states (not neccessary, you can leave them empty!)
#my $dsspdir="";                                # Put the directory with dssp files 
#my $dssp="";                                   # Put the directory with dssp executable
#my $pdbdir="";                                 # Put the directory with PDB files

#########################################################################################################

my $psipreddir="$bioprogs_dir/psipred/";        # Put the directory path with the PSIPRED executables 
my $execdir=$psipreddir."/bin";
my $datadir=$psipreddir."/data";
my $ss_cit="PSIPRED: Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. JMB 292:195-202.";

use strict;

# Module needed for aligning DSSP-sequence
use Align;

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode

my $numres=100;        # number of residues per line for secondary structure
my $informat="a3m";    # input format
my $neff = 7;          # use alignment with this diversity for PSIPRED prediction

my $help="
Add DSSP states (if available) and PSIPRED secondary structure prediction to a multiple sequence alignment.
Input is a  multiple sequence alignment or a HMMER (multi-)model file. Allowed input formats are 
A2M/FASTA (default), A3M (-a3m), CLUSTAL (-clu), STOCKHOLM (-sto), HMMER (-hmm).
If the input file is an alignment, the output file is in A3M with default name <basename>.a3m.
If the input file is in HMMER format, the output is the same as the input, except that records SSPRD 
and SSCON are added to each model which contain predicted secondary structure and confidence values. 
In this case the output file name is obligatory and must be different from the input file name.
(( Remark: A3M looks misaligned but it is not. To reconvert to FASTA, type ))
((   'reformat.pl file.a3m file.fas'.                                      ))
(( For an explanation of the A3M format, see the User Guide.               ))

Usage: perl addss.pl <ali file> [<outfile>] [-fas|-a3m|-clu|-sto]  
  or   perl addss.pl <ali file> <outfile> -hmm  
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
my $ss_dssp;           # dssp states as string
my $sa_dssp;           # relative solvent accessibility from dssp as string {A,B,C,D,E} A:absolutely buried, B:buried, E:exposed
my $aa_dssp;           # residues from dssp file as string
my $aa_astr;           # residues from infile as string

my $xseq;              # sequence x returned from Align.pm
my $yseq;              # sequence y returned from Align.pm  
my $Sstr;              # match sequence returned from Align.pm

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

if ($informat eq "hmm" && !$outfile) {
    print("Error: no output file given. With the -hmm option an output file is obligatory\n"); exit(1);
}

###############################################################################################
# Reformat input alignment to a3m and psiblast-readable format and generate file with query sequence
###############################################################################################

my $inbase; # $inbasename of infile: remove extension
my $inroot; # $inbasename of infile: remove path and extension
if ($infile=~/(.*)\..*/) {$inbase=$1;} else {$inbase=$infile;}  # remove extension
if ($inbase=~/.*\/(.*)/)  {$inroot=$1;} else {$inroot=$inbase;} # remove path 

############################################################################################

if ($informat ne "hmm") {
    if (!$outfile) {$outfile="$inbase.a3m";}

    # Use first sequence to define match states and reformat input file to a3m and psi
    if ($informat ne "a3m") {
	&System("perl $perl/reformat.pl -v $v2 -M first $informat a3m $infile $inbase.in.a3m");
    } else {
	&System("cp $infile $inbase.in.a3m");
    }
    
    # Read query sequence
    open (INFILE, "<$inbase.in.a3m") or die ("ERROR: cannot open $inbase.in.a3m!\n");
    $/=">"; # set input field separator
    my $i=0;
    $qseq="";
    while ($line=<INFILE>) {
	if ($line eq ">") {next;}
	$line=~s/>$//;
	if ($line=~/^ss_/ || $line=~/^aa_/) {next;}
	$seqs[$i++]=">$line";
	if(!$qseq) {
	    $line=~s/^(.*)[^\n]*//;
	    $name=$1;
	    $qseq=uc($line);
	    $qseq=~s/\n//g;
	}
    }
    close(INFILE);
    
    if ($qseq =~ /\-/) {
	
	$/="\n"; # set input field separator
	
	# First sequence contains gaps => calculate consensus sequence
	&System("$hh/hhconsensus -i $inbase.in.a3m -s $inbase.sq -o $inbase.in.a3m > /dev/null");
	
    } else {
	
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
    }
    
    # Filter alignment to diversity $neff 
    if ($v>=1) {printf ("\nFiltering alignment to diversity $neff ...\n");}
    &System("$hh/hhfilter -v $v2 -neff $neff -i $inbase.in.a3m -o $inbase.in.a3m");
    
    # Reformat into PSI-BLAST readable file for jumpstarting 
    &System("perl $perl/reformat.pl -v $v2 -r -noss a3m psi $inbase.in.a3m $inbase.in.psi");
    
    open (ALIFILE, ">$outfile") || die("ERROR: cannot open $inbase.a3m: $!\n");
    
    # Add DSSP sequence (if available)
    if ($dssp ne "") {
        if ($v>=1) {printf ("\nRead DSSP state sequence (if available) ...\n");}
        if (!&AppendDsspSequences("$inbase.sq")) {
	    $ss_dssp=~s/(\S{$numres})/$1\n/g;
	    print(ALIFILE ">ss_dssp\n$ss_dssp\n");
        }
    }

    # Secondary structure prediction with psipred
    if ($v>=1) {printf ("\nPredicting secondary structure with PSIPRED ...\n");}
    &RunPsipred("$inbase.sq");
    
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
##############################################################
# HMMER format
else
{
    if (!$outfile) {$outfile="$inbase.hmm";}

    my $log2 = log(2);
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
	if (-e "$datadir/weights.dat4") { # Psipred version < 3.0
	    &System("$execdir/psipred $inbase.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $inbase.ss");
	} else {
	    &System("$execdir/psipred $inbase.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $inbase.ss");
	}

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
    System("rm $inbase.mtx $inbase.ss $inbase.ss2");
    if ($v>=2) {printf("Added PSIPRED secondary structure to %i models\n",$nmodels);}
}    

if ($v<=4) {
    unlink("$inbase.in.a3m");
    unlink("$inbase.in.psi");
    unlink("$inbase.horiz");
    unlink("$inbase.dssp");
} 

exit;
    
##############################################################################################
# Run SS prediction starting from alignment in $inbase.in.psi (called by BuildAlignment)
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
    &System("$ncbidir/blastpgp -b 1 -j 1 -h 0.001 -d $dummydb -i $infile -B $inbase.in.psi -C $inbase.chk 1> $inbase.blalog 2> $inbase.blalog");
    
    #print("Predicting secondary structure...\n");
    
    system("echo $inroot.chk > $inbase.pn\n");
    system("echo $inroot.sq > $inbase.sn\n");
    system("$ncbidir/makemat -P $inbase");

    if (-e "$datadir/weights.dat4") { # Psipred version < 3.0
	&System("$execdir/psipred $inbase.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $inbase.ss");
    } else {
	&System("$execdir/psipred $inbase.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $inbase.ss");
    }

    &System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $inbase.ss2 $inbase.ss > $inbase.horiz");
    
    # Remove temporary files
    unlink(split ' ', "$inbase.pn $inbase.sn $inbase.mn $inbase.chk $inbase.blalog $inbase.mtx $inbase.aux $inbase.ss $inbase.ss2 $inbase.sq");
    return;
}

##############################################################################################
# Read query sequence and extract dssp sequence
##############################################################################################
sub AppendDsspSequences() {
    my $qfile=$_[0];

    my $line;        #input line
    my $name;        #name of sequence in in file, e.g. d1g8ma1
    my $qrange;      #chain and residue range of query sequence
    my $aas="";      #amino acids from in file for each $name
    
    my $dsspfile;
    my $pdbfile;
    my $pdbcode;     #pdb code for accessing dssp file; shortened from in code, e.g. 1g8m 
    my @ss_dssp=();  #dssp states for residues (H,E,L)
    my @sa_dssp=();  #dssp states for residues (H,E,L)
    my @aa_dssp=();  #residues in dssp file
    my @aa_astr=();  #residues from infile
    my $length;      #length of sequence

    # Default parameters for Align.pm
    our $d=3;    # gap opening penatlty for Align.pm
    our $e=0.1;  # gap extension penatlty for Align.pm
    our $g=0.09; # endgap penatlty for Align.pm
    our $matrix="identity";

    # Read query sequence -> $name, $nameline, $range, $aas 
    open (QFILE, "<$qfile") || die ("cannot open $qfile: $!");
    while ($line=<QFILE>) {
	if ($line=~/>(\S+)/) {
	    $name=$1;

	    # SCOP ID? (d3lkfa_,d3grs_3,d3pmga1,g1m26.1)
	    if ($line=~/^>[defgh](\d[a-z0-9]{3})[a-z0-9_.][a-z0-9_]\s+[a-z]\.\d+\.\d+\.\d+\s+\((\S+)\)/) {
		$pdbcode=$1;
		$qrange=$2;
	    } 
	    
	    # PDB ID? (8fab_A, 1a0i)
	    elsif ($line=~/^>(\d[a-z0-9]{3})_?(\S?)\s/) {
		$pdbcode=$1;
		if ($2 ne "") {$qrange="$2:";} else {$qrange="-";}
	    }
	    
	    # DALI ID? (8fabA_0,1a0i_2)
	    elsif ($line=~/^>(\d[a-z0-9]{3})[A-Za-z0-9]?_\d+\s+\d+\.\d+.\d+.\d+.\d+.\d+\s+\((\S+)\)/) {
		$pdbcode=$1;
		$qrange=$2;
	    }
	    
	    else {
		if ($v>=3) {print("Warning: no pdb code found in sequence name '$name'\n");} 
		close(QFILE);
		return 1; # no astral/DALI/pdb sequence => no dssp states available
	    }
	    $aas="";

	}
	else
	{
	    chomp($line);
	    $line=~tr/a-z \t/A-Z/d;
	    $aas.=$line;
	}
    }
    close(QFILE);
    if ($v>=2) {printf("Searching DSSP state assignments...\nname=%s  range=%s\n",$name,$qrange);}

    # Try to open dssp file 
    $dsspfile="$dsspdir/$pdbcode.dssp";
    if (! open (DSSPFILE, "<$dsspfile")) {
	printf(STDOUT "WARNING: Cannot open $dsspfile!\n"); 
	$pdbfile="$pdbdir/pdb$pdbcode.ent";
	if (! -e $pdbfile) {
	    printf(STDOUT "WARNING Cannot open $pdbfile!\n"); 
	    return 1;
	} else  {
	    &System("$dssp $pdbfile $inbase.dssp > /dev/null");
	    &System("cp $inbase.dssp $dsspfile ");
	    $dsspfile="$inbase.dssp";
	    if (! open (DSSPFILE, "<$dsspfile")) {
		printf(STDERR "ERROR: dssp couldn't generate file from $pdbfile. Skipping $name\n");
		return 1;
	    } 
	}
    }

    #....+....1....+....2....+....3....+....4
    #  #  RESIDUE AA STRUCTURE BP1 BP2  ACC  etc.
    #  623  630 A R     <        0   0  280  etc. 
    #  624        !*             0   0    0  etc. 
    #  625    8 B A              0   0  105  etc. 
    #  626    9 B P    >>  -     0   0   71  etc. 
    #  292   28SA K  H  4 S+     0   0   71  etc.  (1qdm.dssp)
    #  293   29SA K  H  > S+     0   0   28  etc.    

    # Read in whole DSSP file
    for (my $try = 1; $try<=2; $try++) {
	$aa_dssp="";
	$sa_dssp="";
	$ss_dssp="";
	while ($line=<DSSPFILE>) {if ($line=~/^\s*\#\s*RESIDUE\s+AA/) {last;}}
	while ($line=<DSSPFILE>) 
	{
	    if ($line=~/^.{5}(.{5})(.)(.)\s(.).\s(.).{18}(...)/)
	    {
		my $thisres=$1;
		my $icode=$2;
		my $chain=$3;
		my $aa=$4;
		my $ss=$5;
		my $sa=$6;
		my $contained=0;
		my $range=$qrange;  
		if ($aa eq "!")  {next;}    # missing residues!
		$thisres=~tr/ //d;
		$chain=~tr/ //d;
		$icode=~tr/ //d;
		$sa=~tr/ //d;
		if ($try==1) {
		    do{
			if    ($range=~s/^(\S):(-?\d+)[A-Z]-(\d+)([A-Z])// && $chain eq $1 && $icode eq $4 && $2<=$thisres && $thisres<=$3) {
			    $contained=1; #syntax (A:56S-135S)
			}
			elsif ($range=~s/^(\S):(-?\d+)[A-Z]?-(\d+)[A-Z]?// && $chain eq $1 && $2<=$thisres && $thisres<=$3) {
			    $contained=1; #syntax (R:56-135)
			}
			elsif ($range=~s/^(-?\d+)[A-Z]-(\d+)([A-Z])// && $chain eq "" && $icode eq $3 && $1<=$thisres && $thisres<=$2) {
			    $contained=1; #syntax (56-135)
			}
			elsif ($range=~s/^(-?\d+)[A-Z]?-(\d+)[A-Z]?// && $chain eq "" && $1<=$thisres && $thisres<=$2) {
			    $contained=1; #syntax (56-135)
			}
			elsif ($range=~s/^(\S):// && $chain eq $1) {
			    $contained=1; #syntax (A:) or (A:,2:)
			} 
			elsif ($range=~s/^-$// && $chain eq "") {
			    $contained=1; #syntax (-) 
			}
			$range=~s/^,//;
#			print("qrange=$qrange  range='$range'  ires=$thisres  chain=$chain contained=$contained\n");
		    } while($contained==0 && $range ne "");
		    if ($contained==0) {next;}
		} # end if try==1
		$aa_dssp.=$aa;
		$ss_dssp.=$ss;
		$sa_dssp.=&sa2c($sa,$aa);
	    }
	}
	# if not enough residues were found: chain id is wrong => repeat extraction without checking chain id 
	if (length($aa_dssp)>=10) {last;} 
	close(DSSPFILE);
	open (DSSPFILE, "<$dsspfile");
   }
    close(DSSPFILE);

    if (length($aa_dssp)==0) {print("WARNING: no residues found in $dsspdir/$pdbcode.dssp\n"); return 1;} 
    if (length($aa_dssp)<=20) {printf("WARNING: only %i residues found in $dsspdir/$pdbcode.dssp\n",length($aa_dssp)); return 1;} 

    # Postprocess $aa_dssp etc
    $aa_dssp =~ tr/a-z/CCCCCCCCCCCCCCCCCCCCCCCCCC/;
    $ss_dssp =~ tr/ I/CC/;
    $ss_dssp =~ s/ \S /   /g;
    $ss_dssp =~ s/ \S\S /    /g;

    # Align query with dssp sequence
    $aa_astr = $aas;
    $xseq=$aas;
    $yseq=$aa_dssp;
    my ($imax,$imin,$jmax,$jmin);
    my (@i,@j);
    my $score=&AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr);  

    # Initialize strings (=arrays) for dssp states with "----...-"
    my @ss_dssp_ali=();   # $ss_dssp_ali[$i] is dssp state aligned to $aa_astr[$i] 
    my @sa_dssp_ali=();   # $sa_dssp_ali[$i] is solvent accessibility string
    my @aa_dssp_ali=();   # $aa_dssp_ali[$i] is dssp residue aligned to $aa_astr[$i] 
    for (my $i=0; $i<=length($aa_astr); $i++) { # sum up to len+1 
	                                        # because 0'th element in @ss_dssp and @aa_dssp is dummy "-" 
	$ss_dssp_ali[$i]="-";	
	$sa_dssp_ali[$i]="-";	
	$aa_dssp_ali[$i]="-";	
    }
    
    # To each residue (from i=0 to len-1) of input sequence $aa_astr assign aligned dssp state
    @ss_dssp = split(//,$ss_dssp);
    @sa_dssp = split(//,$sa_dssp);
    @aa_dssp = split(//,$aa_dssp);
    @aa_astr = split(//,$aa_astr);
    my $len = 0;
    unshift(@aa_dssp,"-"); #add a gap symbol at beginning -> first residue is at 1!
    unshift(@ss_dssp,"-"); #add a gap symbol at beginning -> first residue is at 1!
    unshift(@sa_dssp,"-"); #add a gap symbol at beginning -> first residue is at 1!
    unshift(@aa_astr,"-"); #add a gap symbol at beginning -> first residue is at 1!
    for (my $col=0; $col<@i; $col++) {
	if ($i[$col]>0) {
	    if ($j[$col]>0) {$len++;} # count match states (for score/len calculation)
	    $ss_dssp_ali[$i[$col]]=$ss_dssp[$j[$col]];
	    $sa_dssp_ali[$i[$col]]=$sa_dssp[$j[$col]];
	    $aa_dssp_ali[$i[$col]]=$aa_dssp[$j[$col]];
	}
	if ($v>=4) {
	    printf ("%s %3i   %s %3i\n",$aa_astr[$i[$col]],$i[$col],$aa_dssp[$j[$col]],$j[$col]);
	}
    }
    shift (@ss_dssp_ali);   # throw out first "-" 
    shift (@sa_dssp_ali);   # throw out first "-" 
    shift (@aa_dssp_ali);   # throw out first "-" 
    $aa_dssp=join("",@aa_dssp_ali);
    $ss_dssp=join("",@ss_dssp_ali);
    $sa_dssp=join("",@sa_dssp_ali);

    # Debugging output
    if ($v>=4) {printf(STDOUT "DSSP: %s: length=%-3i  score/len:%-5.3f\n",$name,$len,$score/$len);}
    if ($v>=4) {
	printf("IN:    %s\n",$xseq);
	printf("MATCH: %s\n",$Sstr);
	printf("DSSP:  %s\n",$yseq);
	printf("\n");
	printf(">ss_dssp $name\n$ss_dssp\n");
	printf(">sa_dssp $name\n$sa_dssp\n");
	printf(">aa_dssp $name\n$aa_dssp\n");
	printf(">aa_astra $name\n$aa_astr\n\n");
    }    
    if ($score/$len<0.5) {
	printf (STDOUT "\nWARNING: in $name: alignment score with dssp residues too low: Score/len=%f.\n\n",$score/$len);
	printf("IN:    %s\n",$xseq);
	printf("MATCH: %s\n",$Sstr);
	printf("DSSP:  %s\n",$yseq);
	return 1;
    }

    return 0;
}

################################################################################################
### Return solvent accessibility code
################################################################################################
sub sa2c ()
{
    my %maxsa = (A=>106, B=>160, C=>135, D=>163, E=>194, F=>197,  G=>84, H=>184, I=>169, K=>205, L=>164, M=>188, 
		 N=>157, P=>136, Q=>198, R=>248, S=>130, T=>142, V=>142, W=>227, X=>180, Y=>222, Z=>196); # maximum solvent accessiblity
    if ($_[1]=~/[a-z]/)  {return "F";}      # disulphide bridge
    if (!defined $maxsa{$_[1]}) {return "-";} # no amino acid
    my $rsa=$_[0]/$maxsa{$_[1]};
#    printf("aa=%s  sa=%5.1f  max_sa=%5.1f  rsa=%5.3f\n",$_[1],$_[0],$maxsa{$_[1]},$rsa);
    if    ($rsa<=0.02) {return "A";}
    elsif ($rsa<=0.14) {return "B";}
    elsif ($rsa<=0.33) {return "C";}
    elsif ($rsa<=0.55) {return "D";}
    else               {return "E";}
}

################################################################################################
### System command
################################################################################################
sub System()
{
    if ($v>2) {printf("%s\n",$_[0]);} 
    return system($_[0])/256;
}


