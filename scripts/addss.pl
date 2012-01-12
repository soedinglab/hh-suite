#!/usr/bin/env perl
# addss.pl version 1.0.0 (October 2009)
# Add PSIPRED secondary structure prediction (and DSSP annotation) to an MSA or HMMER file.
# Output format is A3M (for input alignments) or HMMER (see User Guide).

#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#     (C) Johannes Soeding and Michael Remmert, 2012

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http:#www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
my $ss_cit="PSIPRED: Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. JMB 292:195-202.";

use Align;     # Needleman-Wunsch and Smith-Waterman alignment functions
use File::Temp qw/ tempfile tempdir /;
use strict;

my $dummydb = $ENV{"HHLIB"}."/data/do_not_delete";
my $ss_cit="PSIPRED: Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. JMB 292:195-202.";


# Module needed for aligning DSSP-sequence

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode

my $numres=100;        # number of residues per line for secondary structure
my $informat="a3m";    # input format
my $neff = 7;          # use alignment with this diversity for PSIPRED prediction

my $help="
Add PSIPRED secondary structure prediction (and DSSP annotation) to a multiple sequence alignment (MSA) 
or HMMER (multi-)model file. 

If the input file is an MSA, the predicted secondary structure and confidence values are added as 
special annotation sequences with names >ss_pred, >ss_conf, and >ss_dssp to the top of the output 
A3M alignment. If no output file is given, the output file will have the same name as the input file, 
except for the extension being replaced by '.a3m'. Allowed input formats are A3m (default), 
A2M/FASTA (-fas, -a2m), CLUSTAL (-clu), STOCKHOLM (-sto), HMMER (-hmm).

If the input file contains HMMER models, records SSPRD and SSCON containing predicted secondary 
structure and confidence values are added to each model. In this case the output file name is 
obligatory and must be different from the input file name.

Usage: perl addss.pl <ali_file> [<outfile>] [-fas|-a3m|-clu|-sto]  
  or   perl addss.pl <hhm_file> <outfile> -hmm  
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

# Create tmpfile
my $tmpdir = tempdir( CLEANUP => 1 );
my ($tmpf, $tmpfile) = tempfile( DIR => $tmpdir );
my $tmpfile_no_dir;
if ($tmpfile=~/.*\/(.*)/)  {$tmpfile_no_dir=$1;} else {$tmpfile_no_dir=$tmpfile;} # remove path 



############################################################################################

if ($informat ne "hmm") {
    if (!$outfile) {$outfile="$inbase.a3m";}

    # Use first sequence to define match states and reformat input file to a3m and psi
    if ($informat ne "a3m") {
	&System("reformat.pl -v $v2 -M first $informat a3m $infile $tmpfile.in.a3m");
    } else {
	&System("cp $infile $tmpfile.in.a3m");
    }
    
    # Read query sequence
    open (INFILE, "<$tmpfile.in.a3m") or die ("ERROR: cannot open $tmpfile.in.a3m!\n");
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
	&System("hhconsensus -i $tmpfile.in.a3m -s $tmpfile.sq -o $tmpfile.in.a3m > /dev/null");
	
    } else {
	
	$query_length = ($qseq=~tr/A-Z/A-Z/);
	$qseq=~tr/a-zA-Z//cd;
	
	# If less than 26 match states => add sufficient number of Xs to the end of each sequence in $tmpfile.in.a3m
	my $q_match = ($qseq=~tr/A-Z/A-Z/); # count number of capital letters
	if ($q_match<=25) {                 # Psiblast needs at least 26 residues in query
	    my $addedXs=('X' x (26-$q_match))."\n";
	    $qseq.=$addedXs;     # add 'X' to query to make it at least 26 resiudes long
	    for ($i=0; $i<@seqs; $i++) {	    
		$seqs[$i]=~s/\n$//g;
		$seqs[$i].=$addedXs;
	    }
	    open (INFILE,">$tmpfile.in.a3m");
	    for ($i=0; $i<@seqs; $i++) {
		printf(INFILE "%s",$seqs[$i]);
	    }
	    close INFILE;
	}
	$/="\n"; # set input field separator
	
	# Write query sequence file in FASTA format
	open (QFILE, ">$tmpfile.sq") or die("ERROR: can't open $tmpfile.sq: $!\n");
	printf(QFILE ">%s\n%s\n",$name,$qseq);
	close (QFILE);
    }
    
    # Filter alignment to diversity $neff 
    if ($v>=1) {printf ("Filtering alignment to diversity $neff ...\n");}
    &System("hhfilter -v $v2 -neff $neff -i $tmpfile.in.a3m -o $tmpfile.in.a3m");
    
    # Reformat into PSI-BLAST readable file for jumpstarting 
    &System("reformat.pl -v $v2 -r -noss a3m psi $tmpfile.in.a3m $tmpfile.in.psi");
    
    open (ALIFILE, ">$outfile") || die("ERROR: cannot open $outfile: $!\n");
    
    # Add DSSP sequence (if available)
    if ($dssp ne "") {
        if (!&AppendDsspSequences("$tmpfile.sq")) {
	    $ss_dssp=~s/(\S{$numres})/$1\n/g;
	    print(ALIFILE ">ss_dssp\n$ss_dssp\n");
	    if ($v>=1) {printf ("\nAdding DSSP state sequence ...\n");}
        }
    }

    # Secondary structure prediction with psipred
    if ($v>=2) {printf ("Predicting secondary structure with PSIPRED ... ");}
    &RunPsipred("$tmpfile.sq");
    
    if (open (PSIPREDFILE, "<$tmpfile.horiz")) {
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
    if ($v>=2) {printf ("done \n");}
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
	open (MTXFILE, ">$tmpfile.mtx") || die("ERROR: cannot open $tmpfile.mtx: $!\n");
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
	&System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $tmpfile.ss");
	
	# READ PSIPRED file
	if (open (PSIPRED, "$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $tmpfile.ss2 $tmpfile.ss |")) {
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
		printf(OUTFILE "SSCIT HHsearch-readable PSIPRED secondary structure prediction:\n");
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
    System("rm $tmpfile.mtx $tmpfile.ss $tmpfile.ss2");
    if ($v>=2) {printf("Added PSIPRED secondary structure to %i models\n",$nmodels);}
}    

if ($v<=4) {
    unlink("$tmpfile.in.a3m");
    unlink("$tmpfile.in.psi");
    unlink("$tmpfile.horiz");
    unlink("$tmpfile.dssp");
} 

exit;
    
##############################################################################################
# Run SS prediction starting from alignment in $tmpfile.in.psi (called by BuildAlignment)
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
    &System("$ncbidir/blastpgp -b 1 -j 1 -h 0.001 -d $dummydb -i $infile -B $tmpfile.in.psi -C $tmpfile.chk 1> $tmpfile.blalog 2> $tmpfile.blalog");
    
    #print("Predicting secondary structure...\n");
    
    &System("echo "."$tmpfile_no_dir".".chk > $tmpfile.pn\n");
    &System("echo "."$tmpfile_no_dir".".sq  > $tmpfile.sn\n");
    &System("$ncbidir/makemat -P $tmpfile");
    
    &System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $tmpfile.ss");

    &System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $tmpfile.ss2 $tmpfile.ss > $tmpfile.horiz");
    
    # Remove temporary files
    unlink(split ' ', "$tmpfile.pn $tmpfile.sn $tmpfile.mn $tmpfile.chk $tmpfile.blalog $tmpfile.mtx $tmpfile.aux $tmpfile.ss $tmpfile.ss2 $tmpfile.sq");
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
	    &System("$dssp $pdbfile $tmpfile.dssp > /dev/null");
	    &System("cp $tmpfile.dssp $dsspfile ");
	    $dsspfile="$tmpfile.dssp";
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


#############################################################################
# Package Align
# Author: Johannes Soeding, 2006
# This is free software. You may use it under the terms of the GNU public license
# No warranty of any kind is given.
#############################################################################

#############################################################################
# Subroutine AlignSW
# Smith-Waterman local alignment
# usage: 
# 1. Use global variables of package Align.pm:
#    $score = &AlignSW();
#    printf("  XSEQ: $Align::xseq\n");
#    printf(" MATCH: $Align::Sstr\n");
#    printf("  YSEQ: $Align::yseq\n");
#    etc.
# 
# 2. Use references and/or global variables
#    $score = &AlignSW(\$xseq,\$yseq);
#    $score = &AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr,\@S);
#    printf("  XSEQ: $xseq\n");
#    printf(" MATCH: $Sstr\n");
#    printf("  YSEQ: $yseq\n");
#
# Input:  $xseq, $yseq   : sequences x and y as strings
# Param:  $main::d       : gap opening penalty
#         $main::e       : gap extension penalty 
# Output: return value   : bit score
#         $xseq, $yseq   : aligned residues of x and y (with - as gap)           
#         @i             : $i[$col],$j[$col] are aligned residues in column $col 
#         @j             :                   (first is 1 (NOT 0!), 0 means gap)
#         $imin          : first aligned residue of sequence x
#         $imax          : last  aligned residue of sequence x
#         $jmin          : first aligned residue of sequence y
#         $jmax          : last  aligned residue of sequence y
#         $Sstr          : string belonging to $xseq and $yseq showing quality of alignment
#         $S[$col]       : match score for aligning positions $i[$col] and $j[$col] 
#############################################################################

#############################################################################
# Subroutine AlignNW
# Needleman-Wunsch global alignment
# usage: $score = &AlignNW();
#        $score = &AlignNW(\$xseq,\$yseq);
#        $score = &AlignNW(\$xseq,\$yseq,\@i,\@j);
#        $score = &AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr,\@S);
#
# Input:  $xseq, $yseq   : sequences x and y as strings
# Param:  $main::d       : gap opening penalty
#         $main::e       : gap extension penalty 
#         $main::g       : end gap penalty 
# Output: return value   : bit score
#         $xseq, $yseq   : aligned residues of x and y (with - as gap)           
#         @i             : $i[$col],$j[$col] are aligned residues in column $col 
#         @j             :                   (first is 1 (NOT 0!), 0 means gap)
#         $imin          : first aligned residue of sequence x
#         $imax          : last  aligned residue of sequence x
#         $jmin          : first aligned residue of sequence y
#         $jmax          : last  aligned residue of sequence y
#         $Sstr          : string belonging to $xseq and $yseq showing quality of alingment
#         $S[$col]       : match score for aligning positions $i[$col] and $j[$col] 
#############################################################################

package Align;

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION=1.00;
our @ISA          = qw(Exporter);
our @EXPORT       = qw(&AlignSW &AlignNW $matrix);

our $xseq;      # first sequence
our $yseq;      # second sequence
our $ri;        # reference to input array: $i[$col] -> $ri->[$col] 
our $rj;        # reference to input array: $j[$col] -> $rj->[$col] 
our $imin;      # first aligned residue of sequence x
our $imax;      # last aligned residue of sequence x
our $jmax;      # first aligned residue of sequence y
our $jmin;      # last aligned residue of sequence y
our $Sstr;      # $Sstr annotates the match quality
our $rS;        # reference $rS->[$col] ->  $S[$col] = match score for aligning positions $i[$col] and $j[$col]  
our $matrix;

my $firstcall=1;
my @Sab;               # Substitution matrix in bit
#          A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
my @ch2i=( 0, 3, 4, 3, 6,13, 7, 8, 9,20,11,10,12, 2,20,14, 5, 1,15,16, 4,19,17,20,18, 6);
my @Gonnet = (
#  A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    X     
#   The Gonnet matrix is in units of 10*log10()
[ 2.4,-0.6,-0.3,-0.3, 0.5,-0.2, 0.0, 0.5,-0.8,-0.8,-1.2,-0.4,-0.7,-2.3, 0.3, 1.1, 0.6,-3.6,-2.2, 0.1,-1.0,-9.9], # A
[-0.6, 4.7, 0.3,-0.3,-2.2, 1.5, 0.4,-1.0, 0.6,-2.4,-2.2, 2.7,-1.7,-3.2,-0.9,-0.2,-0.2,-1.6,-1.8,-2.0,-1.0,-9.9], # R
[-0.3, 0.3, 3.8, 2.2,-1.8, 0.7, 0.9, 0.4, 1.2,-2.8,-3.0, 0.8,-2.2,-3.1,-0.9, 0.9, 0.5,-3.6,-1.4,-2.2,-1.0,-9.9], # N
[-0.3,-0.3, 2.2, 4.7,-3.2, 0.9, 2.7, 0.1, 0.4,-3.8,-4.0, 0.5,-3.0,-4.5,-0.7, 0.5, 0.0,-5.2,-2.8,-2.9,-1.0,-9.9], # D
[ 0.5,-2.2,-1.8,-3.2,11.5,-2.4,-3.0,-2.0,-1.3,-1.1,-1.5,-2.8,-0.9,-0.8,-3.1, 0.1,-0.5,-1.0,-0.5, 0.0,-1.0,-9.9], # C
[-0.2, 1.5, 0.7, 0.9,-2.4, 2.7, 1.7,-1.0, 1.2,-1.9,-1.6, 1.5,-1.0,-2.6,-0.2, 0.2, 0.0,-2.7,-1.7,-1.5,-1.0,-9.9], # Q
[ 0.0, 0.4, 0.9, 2.7,-3.0, 1.7, 3.6,-0.8, 0.4,-2.7,-2.8, 1.2,-2.0,-3.9,-0.5, 0.2,-0.1,-4.3,-2.7,-1.9,-1.0,-9.9], # E
[ 0.5,-1.0, 0.4, 0.1,-2.0,-1.0,-0.8, 6.6,-1.4,-4.5,-4.4,-1.1,-3.5,-5.2,-1.6, 0.4,-1.1,-4.0,-4.0,-3.3,-1.0,-9.9], # G
[-0.8, 0.6, 1.2, 0.4,-1.3, 1.2, 0.4,-1.4, 6.0,-2.2,-1.9, 0.6,-1.3,-0.1,-1.1,-0.2,-0.3,-0.8,-2.2,-2.0,-1.0,-9.9], # H
[-0.8,-2.4,-2.8,-3.8,-1.1,-1.9,-2.7,-4.5,-2.2, 4.0, 2.8,-2.1, 2.5, 1.0,-2.6,-1.8,-0.6,-1.8,-0.7, 3.1,-1.0,-9.9], # I
[-1.2,-2.2,-3.0,-4.0,-1.5,-1.6,-2.8,-4.4,-1.9, 2.8, 4.0,-2.1, 2.8, 2.0,-2.3,-2.1,-1.3,-0.7, 0.0, 1.8,-1.0,-9.9], # L
[-0.4, 2.7, 0.8, 0.5,-2.8, 1.5, 1.2,-1.1, 0.6,-2.1,-2.1, 3.2,-1.4,-3.3,-0.6, 0.1, 0.1,-3.5,-2.1,-1.7,-1.0,-9.9], # K
[-0.7,-1.7,-2.2,-3.0,-0.9,-1.0,-2.0,-3.5,-1.3, 2.5, 2.8,-1.4, 4.3, 1.6,-2.4,-1.4,-0.6,-1.0,-0.2, 1.6,-1.0,-9.9], # M
[-2.3,-3.2,-3.1,-4.5,-0.8,-2.6,-3.9,-5.2,-0.1, 1.0, 2.0,-3.3, 1.6, 7.0,-3.8,-2.8,-2.2, 3.6, 5.1, 0.1,-1.0,-9.9], # F
[ 0.3,-0.9,-0.9,-0.7,-3.1,-0.2,-0.5,-1.6,-1.1,-2.6,-2.3,-0.6,-2.4,-3.8, 7.6, 0.4, 0.1,-5.0,-3.1,-1.8,-1.0,-9.9], # P
[ 1.1,-0.2, 0.9, 0.5, 0.1, 0.2, 0.2, 0.4,-0.2,-1.8,-2.1, 0.1,-1.4,-2.8, 0.4, 2.2, 1.5,-3.3,-1.9,-1.0,-1.0,-9.9], # S
[ 0.6,-0.2, 0.5, 0.0,-0.5, 0.0,-0.1,-1.1,-0.3,-0.6,-1.3, 0.1,-0.6,-2.2, 0.1, 1.5, 2.5,-3.5,-1.9, 0.0,-1.0,-9.9], # T
[-3.6,-1.6,-3.6,-5.2,-1.0,-2.7,-4.3,-4.0,-0.8,-1.8,-0.7,-3.5,-1.0, 3.6,-5.0,-3.3,-3.5,14.2, 4.1,-2.6,-1.0,-9.9], # W
[-2.2,-1.8,-1.4,-2.8,-0.5,-1.7,-2.7,-4.0,-2.2,-0.7, 0.0,-2.1,-0.2, 5.1,-3.1,-1.9,-1.9, 4.1, 7.8,-1.1,-1.0,-9.9], # Y
[ 0.1,-2.0,-2.2,-2.9, 0.0,-1.5,-1.9,-3.3,-2.0, 3.1, 1.8,-1.7, 1.6, 0.1,-1.8,-1.0, 0.0,-2.6,-1.1, 3.4,-1.0,-9.9], # V
[-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-9.9], # X	  
[-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9]  # ~	  
      );

# A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
my @BLOSUM62 = (
[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0, 0,-9],
[-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1,-9],
[-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-1,-9],
[-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-1,-9],
[ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-2,-9],
[-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-1,-9],
[-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-1,-9],
[ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-9],
[-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-1,-9],
[-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-1,-9],
[-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-1,-9],
[-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-1,-9],
[-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-1,-9],
[-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-1,-9],
[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-9],
[ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0,-9],
[ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0,-9],
[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-2,-9],
[-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-1,-9],
[ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-1,-9],
[ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-9],
[-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9,-9]	  
    );

#    print("Substitution matrix:\n");
#    for ($a=0; $a<=20; $a++) {
#	for ($b=0; $b<=20; $b++) {
#	    printf("%6.1f ",$Sab[$a][$b]);
#	}
#	printf("\n");
#    }


# Set substitution matrix in bits (do only at first call of one of the alignment routines)
sub SetSubstitutionMatrix {
    if ($firstcall) {
	# Transform to bits;
	if (defined($matrix) && $matrix eq "Gonnet") {
	    for (my $a=0; $a<=20; ++$a) {
		for (my $b=0; $b<=20; ++$b) {
		    $Sab[$a][$b] = $Gonnet[$a][$b]*0.3322; # 1*log(10)/log(2);
		}
	    } 
	} elsif (defined($matrix) && $matrix eq "Blosum62") {
	    {printf("Using Blosum62 matrix...\n");}
	    for (my $a=0; $a<=20; $a++) {
		for (my $b=0; $b<=20; $b++) {
		    $Sab[$a][$b] = $BLOSUM62[$a][$b];
		}
	    }
	} else {
	    for (my $a=0; $a<20; ++$a) {
		for (my $b=0; $b<20; ++$b) {
		    $Sab[$a][$b] = -1;  
		}
		$Sab[$a][$a] = 2;
	    }
	    for (my $b=0; $b<=20; ++$b) {
		$Sab[20][$b] = $Sab[$b][20] = 0;  
		$Sab[21][$b] = $Sab[$b][21] = -10;  
	    }
	}

	$firstcall=0;
    }
}

# maxbt(val1,...,valx,\$bt) finds maximum of values and puts index of maximum into $bt
sub maxbt {
    my $rbt=pop @_; # last element of @_ is address of $bt
    my $max = shift;
    my $i=0;
    $$rbt = 0;
    foreach $_ (@_) {
	$i++;
	if ($_>$max) {$max=$_; $$rbt=$i;} 
    }
    return $max;
}

# max3bt(val1,val2,val3,\$bt) finds maximum of values and puts index of maximum into $bt
sub max3bt {
    if ($_[1] < $_[0]) {
	if ($_[2] < $_[0]) {
	    ${$_[3]}=0;
	    return $_[0];
	} else {
	    ${$_[3]}=2;
	    return $_[2];
	}
    } else {
	if ($_[2] < $_[1]) {
	    ${$_[3]}=1;
	    return $_[1];
	} else {
	    ${$_[3]}=2;
	    return $_[2];
	}
    }
}

# max2bt(val1,val2,\$bt) finds maximum of values and puts index of maximum into $bt
sub max2bt {
    if ($_[1] < $_[0]) {
	${$_[2]}=0;
	return $_[0];
    } else {
	${$_[2]}=1;
	return $_[1];
    }
}


#############################################################################
# Subroutien AlignSW
# Smith-Waterman local alignment
#############################################################################
sub AlignSW {
    if (@_>=1) {$xseq=$_[0];}
    if (@_>=2) {$yseq=$_[1];}
    if (@_>=3) {$ri=$_[2];}
    if (@_>=4) {$rj=$_[3];}
    if (@_>=5) {$imin=$_[4];}
    if (@_>=6) {$imax=$_[5];}
    if (@_>=7) {$jmin=$_[6];}
    if (@_>=8) {$jmax=$_[7];}
    if (@_>=9) {$Sstr=$_[8];}
    if (@_>=10) {$rS=$_[9];}

    if (length($$xseq)<1) {warn ("ERROR in Align.pm: sequence x is empty\n"); return 0;}
    if (length($$yseq)<1) {warn ("ERROR in Align.pm: sequence x is empty\n"); return 0;}

    my @xchr;            # ASCII characters of $xseq
    my @ychr;            # ASCII characters of $yseq
    my @xres;            # internal integer representation of residues of x
    my @yres;            # internal integer representation of residues of y

    $$xseq =~ s/\s//g;
    $$yseq =~ s/\s//g;
    @xchr = split(//,$$xseq);
    @ychr = split(//,$$yseq);

    my $Lx=@xchr;        # length of sequence x
    my $Ly=@ychr;        # length of sequence y
    my @M;               # $M[a][b] = score of best alignment of x[1..a] and y[1..b] ending in match state
    my @A;               # $A[a][b] = score of best alignment of x[1..a] and y[1..b] ending in gap in x
    my @B;               # $B[a][b] = score of best alignment of x[1..a] and y[1..b] ending in gap in y
    my @Mbt;             # $Mbt[a][b] = 0:STOP  1:M  2:A  3:B 
    my @Abt;             # $Abt[a][b] = 0:A     1:M
    my @Bbt;             # $Bbt[a][b] = 0:B     1:M
    my $score;           # bit score of alignment
    my $bt;              # backtracing variable set by &maxbt: which argument was largest? (first=0)
    my $state;           # STOP:0  M:1  A:2  B:3 
    my ($i, $j);    # indices for sequence x and y, respectively

    # Transform @xres and @yres to integer
    for ($i=0; $i<@xchr; $i++) {
	my $a=ord(uc($xchr[$i]));
	if ($a<65 || $a>90) {
	    if ($a!=ord(".") && $a!=ord("-") && $a!=ord("~")) {
		printf(STDERR "\nWARNING: invalid symbol '%s' in pos $i of first sequence to be aligned\n",$xchr[$i]);
	    }
	    $xres[$i]=21;
	} else {
	    $xres[$i]=$ch2i[$a-65];
	}
    }
    for ($j=0; $j<@ychr; $j++) {
	my $a=ord(uc($ychr[$j]));
	if ($a<65 || $a>90) {
	    if ($a!=ord(".") && $a!=ord("-") && $a!=ord("~")) {
		printf(STDERR "\nWARNING: invalid symbol '%s' in pos $j of second sequence to be aligned\n",$ychr[$j]);
	    }
	    $yres[$j]=21;
	} else {
	    $yres[$j]=$ch2i[$a-65];
	}
    }
    unshift (@xres,21); unshift (@xchr," "); # insert dummy 0'th element
    unshift (@yres,21); unshift (@ychr," "); # insert dummy 0'th element

    &SetSubstitutionMatrix;

    # Initialization
    for ($i=0; $i<=$Lx; $i++) {
	$M[$i][0]=-999;	$A[$i][0]=-999;	$B[$i][0]=-999;
    }
    for ($j=1; $j<=$Ly; $j++) {
	$M[0][$j]=-999;	$A[0][$j]=-999;	$B[0][$j]=-999;
    }
    
    # Iteration
    for ($i=1; $i<=$Lx; ++$i) {
	my $Mi =$M[$i];
	my $Mi1=$M[$i-1];
	my $Ai =$A[$i];
	my $Ai1=$A[$i-1];
	my $Bi =$B[$i];
	my $Bi1=$B[$i-1];
	my $Sabx=$Sab[$xres[$i]];
	my $j1=0;
	for ($j=1; $j<=$Ly; ++$j, ++$j1) {
	    ${$Mi}[$j] = max3bt(${$Mi1}[$j1],  ${$Ai1}[$j1],  ${$Bi1}[$j1], \$Mbt[$i][$j]) + ${$Sabx}[$yres[$j]];
	    ${$Ai}[$j] = max2bt(${$Ai}[$j1]-$main::e, ${$Mi}[$j1]-$main::d, \$Abt[$i][$j]);
	    ${$Bi}[$j] = max2bt(${$Bi1}[$j]-$main::e, ${$Mi1}[$j]-$main::d, \$Bbt[$i][$j]);
	}
    }

    # Finding maximum
    $score = -1000;
    for ($i=1; $i<=$Lx; $i++) {
	my $Mi =$M[$i];
	for ($j=1; $j<=$Ly; $j++) {
	    if (${$Mi}[$j]>$score) {$score=${$Mi}[$j]; $$imax=$i; $$jmax=$j;}
	}
    }

    # Backtracing
    @$ri=();
    @$rj=();
    @$rS=();
    $state=1; # last state is M
    $i=$$imax; $j=$$jmax;
    $$xseq=""; $$yseq="";
    while ($state) {
	if ($state==1) {        
	    # current state is M (match-match)
	    unshift(@$ri,$i);
	    unshift(@$rj,$j);
	    $state = $Mbt[$i][$j];
	    $$xseq=$xchr[$i].$$xseq;
	    $$yseq=$ychr[$j].$$yseq;
	    unshift(@$rS, $Sab[$xres[$i]][$yres[$j]]);
	    $$imin=$i; $$jmin=$j;
	    $i--; $j--;
	} elsif ($state==2) {
	    # current state is A (gap in x)
	    unshift(@$ri,0);
	    unshift(@$rj,$j);
	    $$xseq="-".$$xseq;
	    $$yseq=$ychr[$j].$$yseq;
	    $bt = $Abt[$i][$j--];
	    if ($bt) {
		# previous state was M
		unshift(@$rS,-$main::d);
		$state = 1;
	    } else {
		# previous state was A
		unshift(@$rS,-$main::e);
	    }
	} else {
	    # current state is B (gap in y)
	    unshift(@$ri,$i);
	    unshift(@$rj,0);
	    $$xseq=$xchr[$i].$$xseq;
	    $$yseq="-".$$yseq;
	    $bt = $Bbt[$i--][$j];
	    if ($bt) {
		# previous state was M
		unshift(@$rS,-$main::d);
		$state = 1;
	    } else {
		# previous state was B
		unshift(@$rS,-$main::e);
	    }
	}
    }

    # Set annotation string representing match quality
    $$Sstr="";
    for (my $col=0; $col<@$ri; $col++) {
	if ($xres[$ri->[$col]] eq $yres[$rj->[$col]]) {
	    $$Sstr.=uc($xchr[$ri->[$col]]);
	    } elsif ($rS->[$col] > 0 ) {
		$$Sstr.="+";
	    } else {
	       $$Sstr.=".";
	    }
    }
    return $score;
}


#############################################################################
# Subroutien AlignNW
# Needleman-Wunsch global alignment
#############################################################################
sub AlignNW {                             
    if (@_>=1) {$xseq=$_[0];}
    if (@_>=2) {$yseq=$_[1];}
    if (@_>=3) {$ri=$_[2];}
    if (@_>=4) {$rj=$_[3];}
    if (@_>=5) {$imin=$_[4];}
    if (@_>=6) {$imax=$_[5];}
    if (@_>=7) {$jmin=$_[6];}
    if (@_>=8) {$jmax=$_[7];}
    if (@_>=9) {$Sstr=$_[8];}
    if (@_>=10) {$rS=$_[9];}

    if (length($$xseq)<1) {warn ("ERROR in Align.pm: sequence x is empty\n"); return 0;}
    if (length($$yseq)<1) {warn ("ERROR in Align.pm: sequence x is empty\n"); return 0;}

    my @xchr;            # ASCII characters of $xseq
    my @ychr;            # ASCII characters of $yseq
    my @xres;            # internal integer representation of residues of x
    my @yres;            # internal integer representation of residues of y

    $$xseq =~ s/\s//g;
    $$yseq =~ s/\s//g;
    @xchr = split(//,$$xseq);
    @ychr = split(//,$$yseq);

    my $Lx=@xchr;        # length of sequence x
    my $Ly=@ychr;        # length of sequence y
    my @M;               # $M[a][b] = score of best alignment of x[1..a] and y[1..b] ending in match state
    my @A;               # $A[a][b] = score of best alignment of x[1..a] and y[1..b] ending in gap in x
    my @B;               # $B[a][b] = score of best alignment of x[1..a] and y[1..b] ending in gap in y
    my @Mbt;             # $Mbt[a][b] = 0:STOP  1:M  2:A  3:B 
    my @Abt;             # $Abt[a][b] = 0:A     1:M
    my @Bbt;             # $Bbt[a][b] = 0:B     1:M
    my $score;           # bit score of alignment
    my $bt;              # backtracing variable set by &maxbt: which argument was largest? (first=0)
    my $state;           # STOP:0  M:1  A:2  B:3 
    my ($i, $j);    # indices for sequence x and y, respectively

    # Transform @xres and @yres to integer
    for ($i=0; $i<@xchr; $i++) {
	my $a=ord(uc($xchr[$i]));
	if ($a<65 || $a>90) {
	    if ($a!=ord(".") && $a!=ord("-") && $a!=ord("~")) {
		printf(STDERR "\nWARNING: invalid symbol '%s' in pos $i of first sequence to be aligned\n",$xchr[$i]);
	    }
	    $xres[$i]=21;
	} else {
	    $xres[$i]=$ch2i[$a-65];
	}
    }
    for ($j=0; $j<@ychr; $j++) {
	my $a=ord(uc($ychr[$j]));
	if ($a<65 || $a>90) {
	    if ($a!=ord(".") && $a!=ord("-") && $a!=ord("~")) {
		printf(STDERR "\nWARNING: invalid symbol '%s' in pos $j of second sequence to be aligned\n",$ychr[$j]);
	    }
	    $yres[$j]=21;
	} else {
	    $yres[$j]=$ch2i[$a-65];
	}
    }
    unshift (@xres,21); unshift (@xchr," "); # insert dummy 0'th element
    unshift (@yres,21); unshift (@ychr," "); # insert dummy 0'th element
    
    &SetSubstitutionMatrix;

    # Initialization
    $M[0][0]=$A[0][0]=$B[0][0]=0;
    for ($i=1; $i<=$Lx; $i++) {
	$M[$i][0] = -999;	
	$A[$i][0] = -999;	
	$B[$i][0] = -$i*$main::g;
	$Bbt[$i][0] = 0; # previous state was B as well (gap in y)
    }
    for ($j=1; $j<=$Ly; $j++) {
	$M[0][$j] = -999;	
	$A[0][$j] = -$j*$main::g;	
	$B[0][$j] = -999;
	$Abt[0][$j] = 0; # previous state was A as well (gap in x)
    }
    
    # Iteration
    for ($i=1; $i<=$Lx; ++$i) {
	my $Mi =$M[$i];
	my $Mi1=$M[$i-1];
	my $Ai =$A[$i];
	my $Ai1=$A[$i-1];
	my $Bi =$B[$i];
	my $Bi1=$B[$i-1];
	my $Sabx=$Sab[$xres[$i]];
	my $j1=0;
	for ($j=1; $j<=$Ly; ++$j, ++$j1) {
	    ${$Mi}[$j] = max3bt(${$Mi1}[$j1],  ${$Ai1}[$j1],  ${$Bi1}[$j1], \$Mbt[$i][$j]) + ${$Sabx}[$yres[$j]];
	    ${$Ai}[$j] = max2bt(${$Ai}[$j1]-$main::e, ${$Mi}[$j1]-$main::d, \$Abt[$i][$j]);
	    ${$Bi}[$j] = max2bt(${$Bi1}[$j]-$main::e, ${$Mi1}[$j]-$main::d, \$Bbt[$i][$j]);
	}
    }

    # Finding maximum
    $score = -1000;
    for ($i=1; $i<=$Lx; $i++) {
	my $endgappenalty = ($Lx-$i)*$main::g;
	if ($M[$i][$Ly]-$endgappenalty > $score) {
	    $score=$M[$i][$Ly]-$endgappenalty; $$imax=$i; $$jmax=$Ly; $state = 1;
	}
	if ($A[$i][$Ly]-$endgappenalty > $score) {
	    $score=$A[$i][$Ly]-$endgappenalty; $$imax=$i; $$jmax=$Ly; $state = 2;
	}
	if ($B[$i][$Ly]-$endgappenalty > $score) {
	    $score=$B[$i][$Ly]-$endgappenalty; $$imax=$i; $$jmax=$Ly; $state = 3;
	}
    }
    for ($j=1; $j<$Ly; $j++) {
	my $endgappenalty = ($Ly-$j)*$main::g;
	if ($M[$Lx][$j]-$endgappenalty > $score) {
	    $score=$M[$Lx][$j]-$endgappenalty; $$imax=$Lx; $$jmax=$j; $state = 1;
	}
	if ($A[$Lx][$j]-$endgappenalty > $score) {
	    $score=$A[$Lx][$j]-$endgappenalty; $$imax=$Lx; $$jmax=$j; $state = 2;
	}
	if ($B[$Lx][$j]-$endgappenalty > $score) {
	    $score=$B[$Lx][$j]-$endgappenalty; $$imax=$Lx; $$jmax=$j; $state = 3;
	}
    }

    # Make sure the end gapped regions are also backtraced
    if ($$jmax<$Ly) {
	$Abt[$Lx][$$jmax+1] = $state;
	for ($j=$$jmax+2; $j<=$Ly; $j++) {$Abt[$Lx][$j] = 0;}
	$state = 2;
    } elsif ($$imax<$Lx) {
	$Bbt[$$imax+1][$Ly] = $state;
	for ($i=$$imax+2; $i<=$Lx; $i++) {$Bbt[$i][$Ly] = 0;}
	$state = 3;
    } else {
	$state = 1;
    }



    # Backtracing
    @$ri=();
    @$rj=();

    @$rS=();
    $i=$Lx; $j=$Ly;
    $$xseq=""; $$yseq="";
    while ($i || $j) {
	if ($state==1) {        
	    # current state is M (match-match)
	    unshift(@$ri,$i);
	    unshift(@$rj,$j);
	    $state = $Mbt[$i][$j]+1; # previous state 
	    $$xseq=$xchr[$i].$$xseq;
	    $$yseq=$ychr[$j].$$yseq;
	    unshift(@$rS, $Sab[$xres[$i]][$yres[$j]]);
	    $$imin=$i; $$jmin=$j;
	    $i--; $j--;
	} elsif ($state==2) {
	    # current state is A (gap in x)
	    unshift(@$ri,0);     # $ri->[$col]=0 for gap in $x
	    unshift(@$rj,$j);
	    $$xseq="-".$$xseq;
	    $$yseq=$ychr[$j].$$yseq;
	    $bt = $Abt[$i][$j--];
	    if ($bt) {
		# previous state was M
		if ($i==$Lx || $i==0) {
		    unshift(@$rS,-$main::g); # end gap
		} else {
		    unshift(@$rS,-$main::d); # gap opening
		}
		$state = 1;
	    } else {
		# previous state was A
		if ($i==$Lx || $i==0) {
		    unshift(@$rS,-$main::g); # end gap
		} else {
		    unshift(@$rS,-$main::e); # gap extension
		}
	    }
	} else {
	    # current state is B (gap in y)
	    unshift(@$ri,$i);
	    unshift(@$rj,0);     # $j[$col]=0 for gap in $y
	    $$xseq=$xchr[$i].$$xseq;
	    $$yseq="-".$$yseq;
	    $bt = $Bbt[$i--][$j];
	    if ($bt) {
		# previous state was M
		if ($j==$Ly || $j==0) {
		    unshift(@$rS,-$main::g); # end gap
		} else {
		    unshift(@$rS,-$main::d); # gap opening
		}
		$state = 1;
	    } else {
		# previous state was B
		if ($j==$Ly || $j==0) {
		    unshift(@$rS,-$main::g); # end gap
		} else {
		    unshift(@$rS,-$main::e); # gap extension
		}
	    }
	}
    }

    # Set annotation string representing match quality
    $$Sstr="";
    for (my $col=0; $col<@$ri; $col++) {
	if ($xres[$ri->[$col]] eq $yres[$rj->[$col]]) {
	    $$Sstr.=uc($xchr[$ri->[$col]]);
	    } elsif ($rS->[$col] > 0 ) {
		$$Sstr.="+";
	    } else {
	       $$Sstr.=".";
	    }
    }
    return $score;
}

1;
