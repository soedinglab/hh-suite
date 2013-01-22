#!/usr/bin/env perl
#
# addss.pl
# Add PSIPRED secondary structure prediction (and DSSP annotation) to an MSA or HMMER file.
# Output format is A3M (for input alignments) or HMMER (see User Guide).


#     HHsuite version 2.0.16 (January 2013)
#
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
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;     # Needleman-Wunsch and Smith-Waterman alignment functions
use File::Temp qw/ tempfile tempdir /;
use strict;

my $ss_cit="PSIPRED: Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. JMB 292:195-202.";


# Module needed for aligning DSSP-sequence

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode

my $numres=0;          # number of residues per line for secondary structure
my $informat="a3m";    # input format
my $neff = 7;          # use alignment with this diversity for PSIPRED prediction
my $program=$0;        # name of perl script
my $pdbfile;

my $help="
addss.pl from HHsuite $VERSION  
Add PSIPRED secondary structure prediction (and DSSP annotation) to a multiple sequence alignment (MSA) 
or HMMER (multi-)model file. 

If the input file is an MSA, the predicted secondary structure and confidence values are added as 
special annotation sequences with names >ss_pred, >ss_conf, and >ss_dssp to the top of the output 
A3M alignment. If no output file is given, the output file will have the same name as the input file, 
except for the extension being replaced by '.a3m'. Allowed input formats are A3M (default), 
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
my $header;            # header of MSA: everything before first '>'
my $name;              # query in fasta format: '>$name [^\n]*\n$qseq\n'
my $qseq;              # residues of query sequence
my $infile;
my $outfile;
my $ss_pred="";        # psipred ss states
my $ss_conf="";        # psipred confidence values
my $ss_dssp;           # dssp states as string
my $sa_dssp;           # relative solvent accessibility from dssp as string {A,B,C,D,E} A:absolutely buried, B:buried, E:exposed
my $aa_dssp;           # residues from dssp file as string
my $aa_astr;           # residues from infile as string
my $q_match;           # number of match states in query sequence
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

if ($options=~s/ -v\s+(\d+) / /g)  {$v=$1;}

# Set input and output file
if ($options=~s/ -i\s+(\S+) //) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$infile=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$outfile=$1;}

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!$infile) {print($help); exit(1);}

my $v2 = $v-1;
if ($v2>2) {$v2--;}
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
my $tmpdir;
if ($v<=3) {$tmpdir = tempdir( CLEANUP => 1);} else {$tmpdir = tempdir( CLEANUP => 0);}
my ($tmpf, $tmpfile) = tempfile( DIR => $tmpdir );
my $tmpfile_no_dir;
if ($tmpfile=~/.*\/(.*)/)  {$tmpfile_no_dir=$1;} else {$tmpfile_no_dir=$tmpfile;} # remove path 



############################################################################################

if ($informat ne "hmm") {
    if (!$outfile) {$outfile="$inbase.a3m";}

    # Use first sequence to define match states and reformat input file to a3m and psi
    if ($informat ne "a3m") {
	&HHPaths::System("$hhscripts/reformat.pl -v $v2 -M first $informat a3m $infile $tmpfile.in.a3m");
    } else {
	&HHPaths::System("cp $infile $tmpfile.in.a3m");
    }
    
    # Read query sequence
    open (INFILE, "<$tmpfile.in.a3m") or die ("ERROR: cannot open $tmpfile.in.a3m!\n");
    $/=">"; # set input field separator
    my $i=0;
    $qseq="";
    $header = <INFILE>;
    $header =~s />$//; 
    while ($line=<INFILE>) {
	$line=~s/>$//;
	if ($line=~/^ss_/ || $line=~/^aa_/) {next;}
	$seqs[$i++]=">$line";
	if(!$qseq) {
	    $line=~s/^(.*)[^\n]*//;
	    $name=$1;
	    $qseq=$line;
	    $qseq=~s/\n//g;
	}
    }
    close(INFILE);
    $/="\n"; # set input field separator

    if ($qseq =~ /\-/) {
	
	# First sequence contains gaps => calculate consensus sequence
	&HHPaths::System("hhconsensus -i $tmpfile.in.a3m -s $tmpfile.sq -o $tmpfile.in.a3m > /dev/null");
	
    } else {
	
	$query_length = ($qseq=~tr/A-Z/A-Z/);
	$qseq=~tr/A-Z//cd; # remove everything except capital letters
	
	# Write query sequence file in FASTA format
	open (QFILE, ">$tmpfile.sq") or die("ERROR: can't open $tmpfile.sq: $!\n");
	printf(QFILE ">%s\n%s\n",$name,$qseq);
	close (QFILE);
    }
    
    # Filter alignment to diversity $neff 
    if ($v>=1) {printf ("Filtering alignment to diversity $neff ...\n");}
    &HHPaths::System("hhfilter -v $v2 -neff $neff -i $tmpfile.in.a3m -o $tmpfile.in.a3m");
    
    # Reformat into PSI-BLAST readable file for jumpstarting 
    &HHPaths::System("$hhscripts/reformat.pl -v $v2 -r -noss a3m psi $tmpfile.in.a3m $tmpfile.in.psi");
    
    open (ALIFILE, ">$outfile") || die("ERROR: cannot open $outfile: $!\n");
    printf (ALIFILE "%s",$header);
    
    # Add DSSP sequence (if available)
    if ($dssp ne "") {
        if (!&AppendDsspSequences("$tmpfile.sq")) {
	    if ($numres) {
		$ss_dssp=~s/(\S{$numres})/$1\n/g;  # insert a line break every $numres residues
	    }
	    printf (ALIFILE ">ss_dssp\n%s\n",$ss_dssp);
	    if ($v>=1) {print("\nAdding DSSP state sequence ...\n");}
        }
    }

    # Secondary structure prediction with psipred
    if ($v>=2) {print("Predicting secondary structure with PSIPRED ... ");}
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
	    if ($numres) {
		$ss_pred=~s/(\S{$numres})/$1\n/g; # insert a line break every $numres residues
		$ss_conf=~s/(\S{$numres})/$1\n/g; # insert a line break every $numres residues
	    }
	printf(ALIFILE ">ss_pred PSIPRED predicted secondary structure\n%s\n",$ss_pred);
	printf(ALIFILE ">ss_conf PSIPRED confidence values\n%s\n",$ss_conf);
    }
    
    # Append alignment sequences to psipred sequences
    for ($i=0; $i<@seqs; $i++) {
	printf(ALIFILE "%s",$seqs[$i]);
    }
    close(ALIFILE);
    if ($v>=2) {print("done \n");}
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
	
	# Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
	if (-e "$datadir/weights.dat4") { # Psipred version < 3.0
	    &HHPaths::System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $tmpfile.ss");
	} else {
	    &HHPaths::System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $tmpfile.ss");
	}
	
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
		$ss_pred=~s/(\S{80})/$1\nSSPRD /g; # insert a line break every 80 residues
		$ss_conf=~s/(\S{80})/$1\nSSCON /g; # insert a line break every 80 residues
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
    &HHPaths::System("rm $tmpfile.mtx $tmpfile.ss $tmpfile.ss2");
    if ($v>=2) {printf("Added PSIPRED secondary structure to %i models\n",$nmodels);}
}    

if ($v<=3) {
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
	if (!-e "$dummydb") {die "Error in addss.pl: Could not find $dummydb\n";}

	&HHPaths::System("cp $infile $dummydb");
	&HHPaths::System("$ncbidir/formatdb -i $dummydb");
	if (!-e "$dummydb.phr") {die "Error in addss.pl: Could not find nor create index files for $dummydb\n";}
    }

    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
    &HHPaths::System("$ncbidir/blastpgp -b 1 -j 1 -h 0.001 -d $dummydb -i $infile -B $tmpfile.in.psi -C $tmpfile.chk 1> $tmpfile.blalog 2> $tmpfile.blalog");
    
    #print("Predicting secondary structure...\n");
    
    &HHPaths::System("echo "."$tmpfile_no_dir".".chk > $tmpfile.pn\n");
    &HHPaths::System("echo "."$tmpfile_no_dir".".sq  > $tmpfile.sn\n");
    &HHPaths::System("$ncbidir/makemat -P $tmpfile");
    
    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
    if (-e "$datadir/weights.dat4") { # Psipred version < 3.0
	&HHPaths::System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $tmpfile.ss");
    } else {
	&HHPaths::System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $tmpfile.ss");
    }

    &HHPaths::System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $tmpfile.ss2 $tmpfile.ss > $tmpfile.horiz");
    
    # Remove temporary files
    if ($v<=3) { unlink(split ' ', "$tmpfile.pn $tmpfile.sn $tmpfile.mn $tmpfile.chk $tmpfile.blalog $tmpfile.mtx $tmpfile.aux $tmpfile.ss $tmpfile.ss2 $tmpfile.sq");}
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
    if ($v>=3) {printf("Searching DSSP state assignments: name=%s  range=%s\n",$name,$qrange);}

    # Try to open dssp file 
    $dsspfile="$dsspdir/$pdbcode.dssp";
    if (! open (DSSPFILE, "<$dsspfile")) {
	if ($v>=3) {printf(STDERR "Warning in $program: Cannot open $dsspfile!\n");} 
	$pdbfile = &OpenPDBfile($pdbcode);
	if ($pdbfile eq "") {return 1;}

	system("$dssp $pdbfile $tmpfile.dssp 2> /dev/null");
	system("cp $tmpfile.dssp $dsspfile 2> /dev/null");
	$dsspfile="$tmpfile.dssp";
	if (! open (DSSPFILE, "<$dsspfile")) {
	    if ($v>=3) {printf(STDERR "Warning in $program: dssp couldn't generate file from $pdbfile. Skipping $name\n");}
	    return 1;
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

# Find the pdb file with $pdbcode in pdb directory 
sub OpenPDBfile() {
 
    my $pdbcode=lc($_[0]);
    if (! -e "$pdbdir") {
	if ($v>=3) {print(STDERR "Warning in $program: pdb directory '$pdbdir' does not exist!\n");} 
	return 1;
    }
    if (-e "$pdbdir/all") {$pdbfile="$pdbdir/all/";}
    elsif (-e "$pdbdir/divided") {$pdbfile="$pdbdir/divided/".substr($pdbcode,1,2)."/";}
    else {$pdbfile="$pdbdir/";}
    if ($pdbdir=~/divided.?$/) {$pdbfile.=substr($pdbcode,1,2)."/";}
    if    (-e $pdbfile."pdb$pdbcode.ent")   {$pdbfile.="pdb$pdbcode.ent";}
    elsif (-e $pdbfile."pdb$pdbcode.ent.gz") {$pdbfile="gunzip -c $pdbfile"."pdb$pdbcode.ent.gz |";}
    elsif (-e $pdbfile."pdb$pdbcode.ent.Z") {$pdbfile="gunzip -c $pdbfile"."pdb$pdbcode.ent.Z |";}
    elsif (-e $pdbfile."$pdbcode.pdb")      {$pdbfile."$pdbcode.pdb";}
    else {
	if ($v>=3) {printf(STDERR "Warning in $program: Cannot find pdb file $pdbfile"."pdb$pdbcode.ent!\n");}
	return "";
    }
    if (!open (PDBFILE, "$pdbfile")) {
	if ($v>=3) {printf(STDERR "Error in $program: Cannot open pdb file: $!\n");}
	return "";
    }
    return $pdbfile;
}
