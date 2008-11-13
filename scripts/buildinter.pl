#!/usr/bin/perl 
# Build child alignments around 'root' alignment[s]
# Usage: buildinter.pl [options] root-alignment[s]
#
# Flow diagram:
# For each input alignment file:
# * Use querydb.pl to jump-start PSI-BLAST with root alignment and pull out list of potential seed sequences
# * Push all seed sequences with seq-id to query (i.e. first sequence of root alignment) < qid_max onto @seeds
# * Store all seed sequences with seq-id to query >= qid_max in hash %neighborhood
# * While (@seeds not empty)  
#    - Remove first element of @seeds -> $seedname, $seedseq, $parent, $Evalue
#    - reject seed sequence if contained in %neighborhood or %falsepositives => next seed 
#    - build alignment $file.a3m and do SS prediction with PSIPRED
#    - locally align seed alignment ($tmp-n.a3m) to root superalignment ($tmp-X.a3m) using hhalign. 
#      Only columns matched to $tmp.a3m match columns will be upper case!
#    - determine E-value from Eseed and from P-value of HMM-HMM comparison (use hhcorr) 
#    - if E-value<=Einter
#       . accept seed alignment $file -> $tmp-n.a3m
#       . Use querydb.pl to jump-start PSI-BLAST with seed alignment $tmp-n.a3m and pull out list of potential seed sequences
#       . Push all seed sequences with seq-id to seed sequence < qid_max onto @seeds
#       . Store all seed sequences with seq-id to seed sequence >= qid_max in hash %neighborhood
#       . Add accepted alignment to superalignmetn $tmp-X.a3m
#    - if E-value>Ehh: reject seed alignment
#    - if P-value from HMM-HMM comparison > Pfalse: add all sequences in alignment to %falsepositives

my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;
use File::Basename;
use POSIX;
#$|= 1; # Activate autoflushing on STDOUT


# CUSTOMIZE
my $dbbase=$nre; 
#my $dbbase="/home/soeding/nr/nre"; # will use PSI-BLAST dbs $dbbase."90" and $dbbase."70" 
                                    # It is recommended to use nre = nr + env. 
                                    # 90: filtered with CD-HIT to 90% maximum pairwise sequence identity


# Sending signal 10 to the buildinter.pl process will terminate the search orderly, conserving sequences found so far.
setpgrp(0,0);
$SIG{'QUIT'}=\&Terminate;
$SIG{'KILL'}=\&Terminate;
$SIG{'INT'}=\&Terminate;
$SIG{'USR1'}=\&Terminate;
$SIG{10}=\&Terminate;

# Use nr+env_nr as standard database
my $tmpdir="/tmp";                    # directory where all temporary files are written: /tmp/UID/PID
if (-e "/mnt/spare0") {$tmpdir="/mnt/spare0";}
if (!-d "$tmpdir/$ENV{USER}") {mkdir("$tmpdir/$ENV{USER}",0777);}
$tmpdir.="/$ENV{USER}/$$";
if (!-d $tmpdir) {mkdir($tmpdir,0777);}
# my $tmpdir=".";                           

# Default values
our $v=2;                     # verbose mode
my $update=0;                 # 0:overwrite  1:do not overwrite
my $Einter=1E-3;              # maximum E-value in hhalign for HMM to be accepted as intermediary HMM
my $EinterY=1E-8;             # dito, but for the case that the PSI-BLAST E-value has been increased during PSI-BLAST search
my $Ey=1;                     # maximum E-value in hhalign for alignment to be added to root-Y.a3m
my $Emax=0.1;                 # maximum allowed psiblast E-value of seeds with parent alignment (can be dynamically adapted)
my $Emin=0.0;                 # minimum allowed psiblast E-value of seeds with parent alignment (can be dynamically adapted)
my $qid_max=0;                # maximum seq identity of seeds with first sequence of parent alignment (0: chose automatically)
my $bl=0.0;                   # lenient minimum per-residue bit score with query at ends of HSP
my $bs=0.167;                 # strict  minimum per-residue bit score with query at ends of HSP
my $bg=20;                    # maximum number of HSP residiues missing at ends of query to use lenient $bl
my $rmin=20;                  # minimum number of match residues in seed
my $extnd=20;                 # maximum number of residues by which seed HSPs will be extended on either side; 
my $nrejmax=20;               # maximum number of rejected seed sequences per parent (i.e. per seed file)
my $Nrejmax=500;              # maximum number of rejected seed sequences per root 
my $Naccmax=1000 ;            # maximum number of accepted seed sequences per root
my $Nalis1=10;                # up to this number of alignments (accepted or rejected), seeds are selected independent of their E-value; 
                              # idmax is adjusted to yield at least this number of seeds
my $Nalis2=20;                # up to this number of alignments, seeds are selected if E-value<1;
my $Nalis3=50;                # up to this number of alignments, seeds are selected if E-value<$Emax && E-value>1E-20;
my $Nalis4=100;               # up to this number of alignments, seeds are selected if E-value<$Emax && E-value>1E-10;
my $Nalis5=300;               # up to this number of alignments, seeds are selected if E-value<$Emax && E-value>1E-6;
my $Nmin_pfam=-1;             # if less than this number of sequences was found when @seeds=0, search PfamA/B
my $Xmax=100000;              # stop search for query when the -X alignment consitst of more than Xmax sequences
my $Ymax=100000;              # stop search for query when the -Y alignment consitst of more than Ymax sequences
my $Phh_max=1E-4;             # max allowable HMM-HMM P-value between PSI-BLAST iter's and for final seed alignment to be added to $tmp-Y.a3m
my $Phh_strict=1E-6;          # max allowable HMM-HMM P-value for alignment to be added to $tmp-X.a3m
my $Phh_false=0.1;            # add ALL sequences from alignment with Phh>Phh_false to %neighborhood

# Thresholds for inclusion into iterated PSI-BLAST alignments
my $maxiter=8;                # maximum number of psiblast iterations
my $E2=1E-3;                  # Psiblast Evalue threshold for building profile
my $Eult=1E-3;                # ultimate Evalue threshold for building output alignment
my $id=90;                    # maximum sequence identity in alignments
my $cov=10;                   # minimum coverage
my $min_hitlen=30;            # minimum number of residues in match states 
my $use_nr=0;                 # $use_nr=1: last search against nr (instead of nr90 or nr70)
my $bopt="";                  # this is NOT the default
my $sc="";                    # minimum score per column with query
my $pmax=1;                   # maximum p-value of HSP IN MATCH COLUMNS for inclusion into alignment
my $best="";                  # extract only the best HSP per sequence, except during last round
my $Nnr_eff=1.0E6;            # Effective number of sequences in $db90
my $qid="";                   # minimum sequence identity with query
my $mactX=0.5;                # column score offset in hhalign; controls greediness of HMM-HMM alignment for -X alignment
my $mact=0.2;                 # column score offset in hhalign; controls greediness of HMM-HMM alignment for -Y alignment
my $cpu=1;                    # number of CPUs to use
my $quick=0;                  # default: don't do quick search
my $skip_similar=1;           # Stop PSI-BLAST iteration if, after first search through nr70, no new sequences have been found


my $usage="
Build child alignments around root alignment[s] 
The program generates the following files for each input alignment root.a3m
 root-X.a3m   super-alignment of all PSI-BLAST alignments accepted as STRICTLY homologous to root.a3m with E-value < E 
 root-Y.a3m   sequences of all alignments accepted as homologous, with E-value < Ey or increased PSI-BLAST E-value threshold 
 root-Z.a3m   sequences of all alignments rejected with Ey < E-value but P-value<$Phh_false (not strongly rejected) 
 root-n.a3m   PSI-BLAST alignments found to be homologous to root.a3m for n=1,2,... (aligned to root.a3m)
The program may run for serveral hours or even days. To terminate it without loosing preliminary results, 
send the -10 signal to -PID (minus PID).

Usage:   buildinter.pl [options] root-alignment[s].a3m 

General options:
 -v  <int>      verbose mode (def=$v)
 -u             update: skip globbed alignment files for which a *.finished file exists already (def=off)
 -E  <float>    strict E-value for HMM to be accepted as strictly homologous (def=$Einter)
 -Ey <float>    loose E-value for HMM to be accepted as homologous (def=$Ey)
 -P  <float>    maximum P-value in HMM-HMM comparison for HMM to be accepted as intermediary HMM (def=$Phh_max)
 -d  <file>     HMM db file for calibration
 -mact <float>  maximum a-posteriory threshol in hhalign; controls greediness of HMM-HMM alignment for sequences in root-Y.a3m (default=$mact)
 -mactX <float> maximum a-posteriory threshol in hhalign; controls greediness of HMM-HMM alignment for sequences in root-X.a3m (default=$mactX)
 -cpu  <int>    number of CPUs to use when calling blastpgp
 -quick         do a quick half-hour hhsenser run in which seeds with *higher* E-values are processed first (Use -Emax 0.1 -tmax 0:30)

Options for seed sequences (i.e. sequences that \'seed\' new PSI-BLAST search)
 -Emax <float>  maximum allowed psiblast E-value of seeds with parent alignment (def=$Emax)
 -idmax <int>   MAX seq-id between seed sequences used to start new PSI-BLAST search (default=automatically adjust during search) 
 -rmin <int>    min coverage: min number of residues aligned with original query (default=$rmin)
 -extnd <int>   number of residues that are added to HSPs on either side to generate seed sequences (default=$extnd)
                This influences rmin: it will be increased to at least extnd+10 to suppress false positives
 -tmax <time>   maximum wall time in <hours>:<minutes> before premature termination (e.g. \"-tmax 0:45\")
 -Prob <int>    stop search when a match in the PDB database with at least this probability is found (default=OFF)

Options for PSI-BLAST search
 -n    <int>    maximum number of psiblast iterations  (def=$maxiter)
 -e    <float>  Evalue for last round (def=$Eult)
 -id   <int>    maximum sequence identity in PSI-BLAST alignments in % (def=$id)
 -cov  <int>    minimum coverage in % (def=$cov)
 -len  <int>    minimum number of template residues aligned with query match states (def=$min_hitlen)
 -bl  <float>   lenient minimum per-residue bit score for pruning ends of HSP (default=$bl)
 -bs  <float>   strict  minimum per-residue bit score for pruning ends of HSP (default=$bs)
 -bg  <float>   use the lenient bl below this number of end gaps and the strict bs otherwise (default=$bg)
 -nr            last iteration against nr (instead of nr70 or nr90).
 -db <basename>  basename of sequence database, e.g. /cluster/databases/nr_euk  ( =>  nr_euk90, nr_euk70)

Examples: 
 buildinter.pl AbrB.a3m 
 buildinter.pl -quick AbrB.a3m 
 buildinter.pl -accmax 200 -rejmax 100 -nr -id 99 -v 3 AbrB.a3m
 nohup nice -19 buildinter.pl *.a3m &> ./buildinter.log &
\n";

# Old options superseeded by -tmax:
# -rejmax <int>  stop search when rejmax seeds have been rejected (default=$Nrejmax)
# -accmax <int>  stop search when accmax seeds have been accepted (default=$Naccmax)
# -Xmax <int>    stop search when the -X alignment consists of more than Xmax sequences (default=$Xmax)
# -Ymax <int>    stop search when the -Y alignment consists of more than Ymax sequences (default=$Ymax)

# Variable declarations
my @files;             # input alignment files read in
my $curr;              # current seed alignment file (without extension)
our $file;             # input alignment file currently being processed (without extension)
our $tmp;              # input alignment file in $tempdir currently being processed (without extension)
our $dir;              # path of input alignment file currently being processed
my @seeds;             # array of pointers to seed sequences to process: @seed=($nameline,$seq,$parent,$Eseed)
my %inX;               # $inX{$nameid} contains index $i of record for sequence $nameid in @inX
my @inX;               # $inX[$i] contains all non-overlapping HSPs from sequence nameid found homologous to $tmp.a3m in FASTA format
my %inY;               # $inY{$nameid} contains index $i of record for sequence $nameid in @inY
my @inY;               # $inY[$i] contains all non-overlapping HSPs from sequence nameid found likely to be homologous to $tmp.a3m
my %inNeigh;           # $inNeigh{$nameid} contains index $i of record for sequence $nameid in @inNeigh
my @inNeigh;           # $inNeigh[$i] contains all non-overlapping HSPs from sequence nameid that have sequence identity > qid_max 
                       # with an already processed seed; sequences in %inNeigh are excluded as seeds
my %inSeeds;           # $inSeeds{$nameid} contains index of record for sequence $nameid in @inSeeds
my @inSeeds;           # $inSeeds[$i] contains all non-overlapping HSPs from sequence nameid that have been identified as seeds, 
                       # i.e. they have PSI-BLAST E-value < Emax with a previously accepted alignment

my $nameline;          # nameline of seed sequence to process
my $seq;               # residues of seed sequence to process
my $parent;            # parent of seed sequence, i.e. alignment by which seed was first detected
my $parent_index=0;    # indes $Nacc of parent of current seed sequence
my $name;              # $nameline="$name $description". 
my $nameid;            # First 14 letters of name (excluding a possible ':' for residue range) 
my $description;       # what comes after the name in the name line of the seed sequence
my $dssp="";           # dssp states of query sequence from root
my $tmax=time()+14*24*3600; # maximum wall time before termination of program
my $Nacc;              # number of accepted alignments for $tmp
my $Nrej;              # number of rejected alignments for $tmp
my %nrej=();           # $nrej{$parent} = number of rejected seed sequences for $parent 
my $Eseed;             # in the psiblast run started with $parent.a3m: E-value of seed sequence $tmp-n.seq 
my $Phh;               # P-value of putative intermediate HMM from HMM-HMM comparison with root alignment
my $Eval;              # total E-value of putative intermediate HMM, calculated from Phh, Eseed, and their correlation
my $Pmax;              # maximum of Phh and Pseed (=Eseed/Nnr_eff)
my $corr_seedHMM_seedSeq; # correlation coefficient between seed-HMM and seed-sequence
my $corr_X_parent;        # correlation coefficient between $tmp-X.hhm and parent HMM that found seed ($curr.hhm)
my $corr;              # total correlation between P(X,seedHMM) and P(parent,seedSeq) = $corr_X_parent x $corr_seedHMM_seedSeq
my $qid_auto;          # 0: take value from $qid_max  1: determine $qid_max automatically
my $line;              # line read in from file
my $nhits;             # number of hits found by alignhits
my $nhits_prev;        # number of nhits from previous psiblast round
my $cov0;              # effective minimum coverage (including $min_hitlen threshold)
my $pid=0;             # process id of currently running blastpgp child process
my $rank;              # rank in hit list of hhalign
my $cpus=0;            # no multiprocessing
my $pfamA=(glob "$newdbs/pfamA*")[0];
my $pfamB=(glob "$newdbs/pfamB*")[0];
my %hitname;           # hash contains names of all PfamA or B domains that have already been included
my $Nacc_last=1;       # used to ensure that querydb.pl with root-Y is only redone when new alignments were added to root-Y
my @path;              # search path for alignment indices, e.g. $path[$Nacc] = "#39 <= #7 <= #1 <= #0"
my $L;                 # length of query sequence
my $cov_rmin;          # = 100*$rmin/$L;
my $pdbdir;            # path to pdb.hhm
my $scalar_inX=0.5;    # value of scalar(@inX) at time of last search of PDB
my $Prob=0 ;           # stop search when a match in the PDB database with at least this probability is found (0=OFF)
my ($db90, $db70, $db);

# print("PfamA=$pfamA   PfamB=$pfamB\n");

$blastpgp.=" -I T";    # show gi's in defline

###############################################################################################
# Processing command line input
###############################################################################################

my $ARGC=scalar(@ARGV);
if ($ARGC<1) {die ($usage);}

my $options="";
for (my $i=0; $i<$ARGC; $i++) {$options.=" $ARGV[$i] ";}

# General options
if ($options=~s/ -v\s*(\d) / /) {$v=$1;}
if ($options=~s/ -v / /) {$v=2;}
if ($v>=2) {print("$0 $options\n");}
if ($options=~s/ -u / /) {$update=1;}
if ($options=~s/ -E\s+(\S+) / /)  {$Einter=$1;}
if ($options=~s/ -Ey\s+(\S+) / /) {$Ey=$1;}
if ($options=~s/ -P\s+(\S+) / /) {$Phh_max=$1;}
if ($options=~s/ -Px\s+(\S+) / /) {$Phh_strict=$1;}
if ($options=~s/ -d\s+(\S+) / /) {$calhhm=$1;}
if ($options=~s/ -mact \s+(\S+)/ /) {$mact=$1;}
if ($options=~s/ -mactX \s+(\S+)/ /) {$mactX=$1;}
if ($options=~s/ -lax / /) {
    $Emax=1; $Einter=$Eult=$E2=0.01; $bl=0.0; $bs=0.167; $mactX=0.5; $Phh_max=1E-4; $Phh_strict=1E-5; $Nrejmax=$Naccmax=1e5; $nrejmax=500;
}
if ($options=~s/ -quick\s+/ /g) {$quick=1; $Emax=0.1; $extnd=20; $tmax=time()+30*60;}
if ($options=~s/ -cpu \s+(\d+)/ /g) {$cpu=$1; $blastpgp.=" -a $1";}

# Options for seeds
if ($options=~s/ -Emax\s+(\S+) / /)  {$Emax=$1;}
if ($options=~s/ -rmin\s+(\d+) / /)  {$rmin=$1;}
if ($options=~s/ -idmax\s+(\S+) / /) {$qid_max=$1;}
if ($options=~s/ -extnd\s+(\S+) / /) {$extnd=$1;}
if ($options=~s/ -tmax\s+(\d+):(\d+) / /) {$tmax=time()+$1*3600+$2*60;}
if ($options=~s/ -rejmax\s+(\S+) / /) {$Nrejmax=$1;}
if ($options=~s/ -accmax\s+(\S+) / /) {$Naccmax=$1;}
if ($options=~s/ -Xmax\s+(\S+) / /)   {$Xmax=$1;}
if ($options=~s/ -Ymax\s+(\S+) / /)   {$Ymax=$1;}
if ($options=~s/ -Prob\s+(\S+) / /)   {$Prob=$1;}
if ($options=~s/ -pdb\s+(\S+) / /)    {$pdbdir=$1;}

# Options for building alignments
if ($options=~s/ -n\s+(\d+) / /) {$maxiter=$1;}
if ($options=~s/ -db\s+(\S+) / /) {$dbbase=$1;}
if ($options=~s/ -e2\s+(\S+) / /) {$E2=$1;}
if ($options=~s/ -e\s+(\S+) / /)  {$Eult=$1; $E2=$1;}
if ($options=~s/ -id\s+(\d+) / /) {$id=$1;}
if ($options=~s/ -cov\s+(\d+) / /) {$cov=$1;}
if ($options=~s/ -len\s+(\d+) / /) {$min_hitlen=$1;}
if ($options=~s/ -bl\s+(\S+) / /g)  {$bl=$1;}
if ($options=~s/ -bs\s+(\S+) / /g)  {$bs=$1;} 
if ($options=~s/ -bg\s+(\S+) / /g)  {$bg=$1;} 
if ($options=~s/ -best / /) {$best="-best";}
if ($options=~s/ -nr / /)   {$use_nr=1;}

# Warn if unknown options found
if ($options=~/-tmax/) {die("Error: wrong format for option -tmax <hours>:<minutes>\n");}
while ($options=~s/\s+(-\S+)//) {print("WARNING: unknown option $1\n");}

# Read remaining arguments into @files and check if they exist
$options=~s/^\s+//;
$options=~s/\s+$//;
@files=split(/\s+/,$options);
foreach $tmp (@files) {if (!-e $tmp) {die("\nError: input file '$tmp' does not exist\n");} }

my $v2 = $v-1;                      # verbose mode for sub-programs called
if ($v2>2) {$v2=2;}       
if ($v2<0) {$v2=0;}       
if ($qid_max==0) {$qid_auto=1;} else {$qid_auto=0;}
if ($bopt eq "") {$bopt="-bl 0 -bs 0.167";} # Default edge pruning paramenters
if ($rmin<$extnd+10) {$rmin=$extnd+10;}

# Set database versions filtered to 90% and 70% sequence identity. Use absolute paths to avoid errors during databas update
&UpdateDBLinks();

if ($v>=2) {
    if (@files==1) {
	print("Read in alignment $files[0]\n");
    } else {
	printf("Read in %i alignments\n",scalar(@files));
    }
}


###############################################################################################
####  Main loop: for each alignment $tmp.a3m that was read in ...
###############################################################################################
foreach $file (@files)
{
    if ($file=~/^(.*)\..*?$/)  {$file=$1;} # remove extension from input file 
    
    # Copy input file to $tmpdir
    if ($file=~/^(.*)\/(.*)?/)     {$dir=$1; $tmp="$tmpdir/$2";} else {$dir="."; $tmp="$tmpdir/$file";}  
    if ($tmp ne $file) {&System("cp $file.a3m $tmp.a3m");}
    
    # Update option: has $tmp already been done before?
    if ($update && -e "$file.finished") {
	print("Skipping $file.a3m\n");
	sleep(5);
	next;
    }

    if ($v>=2) {
	print("\n\n");
	print("*************************************************************************************\n");
	print("*************************************************************************************\n");
	print("*************************************************************************************\n");
	print("      Starting intermediate HMM search for $file.a3m \n");
	print("*************************************************************************************\n");
	print("*************************************************************************************\n");
	print("*************************************************************************************\n");
    } elsif ($v>=1) {
	$tmp=~/([^\/]*)$/;
	print("Starting intermediate HMM search for $1.a3m \n");
    }
    if ($tmax) {
	my $time=($tmax-time())/60;
	printf("Maximum time left before premature termination: %i days, %i hours, %i minutes\n",int($time/60/24),int($time/60)%24,$time%60);
    }    

    # Initialize 
    $parent=$tmp; # first parent of new seeds is file itself
    $path[0]="#0";
    @seeds=();
    %inX=(); @inX=();
    %inY=(); @inY=();
    %inNeigh=(); @inNeigh=();
    %inSeeds=(); @inSeeds=();
    &AddPath("$tmp.a3m",$path[0]);
    &AddSeqsToHash("$tmp.a3m",\%inX,\@inX);
    &WriteSeqHashToFile (\%inX,\@inX,"$tmp"."-X.a3m");
    &AddSeqsToHash("$tmp.a3m",\%inY,\@inY);
    &WriteSeqHashToFile (\%inY,\@inY,"$tmp"."-Y.a3m");

    unlink("$tmp"."-Z.a3m");
    &System("touch $tmp"."-Z.a3m");
    if ($v>=2) {print("\n");}

    # Determine length of query sequence (and make HMM)
    open(HHMAKE,"$hh/hhmake -v 2 -id $id -i $tmp"."-X.a3m -o $tmp"."-X.hhm |");
    <HHMAKE>; <HHMAKE>;$line=<HHMAKE>;
    $line=~/contains (\d+) match states/;
    $L=$1; # length of query sequence
    close(HHMAKE);
    if ($extnd>$L-40) {
	printf("WARNING: query length %i is too short for -extnd %i. ",$L,$extnd);
	printf("Reducing to -extnd %i\n",($extnd=($L-40>0? $L-40: 0)) );
	if ($rmin>$extnd+10) {$rmin=$extnd+10;}
    }
    if ($rmin>$L-30) {
	printf("WARNING: query length %i is too short for -rmin %i. ",$L,$rmin);
	printf("Reducing to -rmin %i\n", ($rmin=($L-30>0? $L-30: 0)) );
    }
    $cov_rmin = sprintf("%i",100*$rmin/$L);
    if ($min_hitlen<$rmin) {$min_hitlen=$rmin;}

    # Calibrate $tmp-X.hhm and $tmp-X-noSS.hhm
    &System("cp $tmp"."-X.hhm $tmp"."-X-noSS.hhm");
    &System("$hh/hhsearch -cpu $cpu -cal -norealign -i $tmp"."-X.hhm -d $calhhm",$v2);
    &System("$hh/hhsearch -cpu $cpu -cal -norealign -i $tmp"."-X-noSS.hhm -d $calhhm -ssm 0",$v2);

    &System("$perl/reformat.pl -noss a3m fas $tmp.a3m $tmp.fas",$v2);
    &System("$perl/reformat.pl -M first -r fas fas $tmp.fas $tmp.fas",$v2);

    # Find sequences around current seed alignment with Evalues <= Emax
    my $Emax0=$Emax;
    if ($quick==0 && $Emax<10) {$Emax0=10;}
    if ($quick==1 && $Emax<10)  {$Emax0=10;}
    $nhits= &System("$perl/querydb.pl -a $cpu -Emax $Emax0 -nmin $Nalis1 -rmin $rmin -extnd $extnd -sec -d $db90 -q $tmp.a3m $tmp.fas $tmp.seeds",$v2);

    if ($qid_auto) {
	# Determine optimal sequence identity threshold
	my @ids;
	$qid_max=90;
	open (SEEDS,"<$tmp.seeds") || die("Error: cannot open $tmp.seeds: $!\n"); 
	while ($line=<SEEDS>) {
	    if ($line=~/^>.* id=(\d+)/) {push(@ids,$1);} 
	}
	close(SEEDS);
	my @ids_sorted = sort {$a<=>$b} @ids; # make sure that at least $seeds_min1 seeds remain 
	if (@ids_sorted<$Nalis1)     {$qid_max=90;}
	elsif ($nhits<=10)   {$qid_max=&max(60,$ids_sorted[$Nalis1-1]); $Emin=0.0;}
	elsif ($nhits<=30)   {$qid_max=&max(50,$ids_sorted[$Nalis1-1]); $Emin=0.0;}
	elsif ($nhits<=100)  {$qid_max=&max(40,$ids_sorted[$Nalis1-1]); $Emin=1E-50;}
	elsif ($nhits<=300)  {$qid_max=&max(35,$ids_sorted[$Nalis1-1]); $Emin=1E-30;}
	elsif ($nhits<=1000) {$qid_max=&max(30,$ids_sorted[$Nalis1-1]); $Emin=1E-20;}
	elsif ($nhits<=3000) {$qid_max=&max(25,$ids_sorted[$Nalis1-1]); $Emin=1E-10;}
	elsif ($nhits<=10000) {$qid_max=&max(25,$ids_sorted[$Nalis1-1]); $Emin=1E-6;}
	else                 {$qid_max=&max(25,$ids_sorted[$Nalis1-1]); $Emin=1E-4;}
 	if ($v>=2) {printf("Choosing max sequence identity of parent seed with new seed idmax=%i%%\n",$qid_max);}
    }

    # Add all seqs with E<=Emax0 to @seeds and sequences with (qid>=qid_max OR E<=$Emin) to %neighborhood 
    &AddSeqsToHash("$tmp.seeds",\%inSeeds,\@inSeeds,100,$Emax0,"add seeds");
    &AddSeqsToHash("$tmp.seeds",\%inNeigh,\@inNeigh,$qid_max,$Emin,"add neigh");

    # Read dssp sequence in $tmp.a3m
    open (ROOT, "<$tmp.a3m") || die ("ERROR: cannot open $tmp.a3m: $!\n");
    $dssp="";
    while ($line=<ROOT>) { if($line=~/^>ss_dssp/) {last;} }
    while ($line=<ROOT>) {
	if ($line=~/^>/) {last;}
	chomp($line);
	$dssp.=$line;
    }
    close(ROOT);
	
    # Initialize for root
    $curr="$tmp.curr"; # $tmp.curr.a3m = name of current seed alignment
    $Nacc=0;
    $Nrej=0;

    # If no more seeds are left and fewer than $Nmin_pfam sequences found: search PfamA/B
#    if (@seeds==0 && scalar(@inY)<$Nmin_pfam) {&SearchPfam();}


    ###############################################################################################
    # While there are still seeds that need to be explored ...
#    printf("Seeds=%-3i Nrej=%-3i Nacc=%-3i\n",scalar(@seeds),$Nrej,$Nacc);
    while ($Nrej<$Nrejmax && $Nacc<$Naccmax && scalar(@inX)<=$Xmax && scalar(@inY)<=$Ymax && @seeds) 
    {
	
	# Take first element of seeds list 
	($nameline,$name,$description,$nameid,$seq,$parent,$Eseed)=@{shift(@seeds)};
#	printf("parent=$parent\nEseed=$Eseed\nname=$nameline"."seq=$seq\n");
	if ($parent=~/-(\d+)$/) {$parent_index=$1;} else {$parent_index=0;}

	# Check if E-value to high (seed likely not to be homolog)
	if ($Nacc+$Nrej>=$Nalis1 && ($Eseed>1 || $Nacc+$Nrej>=$Nalis2) && $Eseed>$Emax) {
	    if ($v>=4) {print("E-value $Eseed: , Alignments $Nacc+$Nrej => E-value of seed $name too high. Skip\n");} 
	    elsif ($v>=2) {print("H");}
	    next;
	}
	# Check if E-value too low (seed too similar to parent) 
	if ($Eseed<$Emin) {
	    if ($v>=4) {print("E-value $Eseed: , Alignments $Nacc+$Nrej => E-value of seed $name too low. Skip\n");} 
	    elsif ($v>=2) {print("E");}
	    next;
	}

	# Check if too many seeds from this parent have already been rejected
	if (!defined $nrej{$parent}) {
	    $nrej{$parent}=0;
	} elsif ($nrej{$parent}>=$nrejmax) {
	    if ($v>=3) {
		printf("Parent $parent has already $nrejmax rejected alignments. Skipping this one ...\n");
	    } elsif ($v>=2) {printf("R");}
	    next;
	} 

	# Check if this seed is contained in neigborhood of seed already used
	my $returned = &ContainedInHash($nameid,$seq,\%inNeigh,\@inNeigh);
	if ($returned>=2) {      # returns  2 if $nameid IS contained in %inNeigh AND $seq has overlap >=10 
	    if ($v>=2) {printf("O");}
	    next; # Skip this seed!
	} elsif ($returned==1) { # returns 1 if $nameid IS contained in %inNeigh BUT $seq has no overlap>=10 with any sequence in @inNeigh,
	    if ($v>=2) {printf("d");}
	} else {                 #  returned 0 if $nameid is NOT contained in %inNeigh
	    if ($v>=2) {printf("n");}
	}
	

	if ($tmax) {
	    my $time=int(($tmax+59-time())/60);
	    printf("\nMaximum time left until premature termination: %i days, %i hours, %i minutes\n",int($time/60/24),int($time/60)%24,$time%60);
	    if ($time<=0) {
		printf(" Maximum time until premature termination exceeded. Terminating search ...\n"); 
		last;
	    }
	}

	# Stop search when match in PDB found?
	if ($Prob>0 && scalar(@inX)>1.2*$scalar_inX) {
	    if (!-e $pdbdir) {
		$pdbdir = (glob $newdbs."/pdb*")[0]."/db/pdb.hhm";
	    }
	    my $prob=&System("$hh/hhsearch -cpu $cpu -i $tmp"."-X.hhm -d $pdbdir -o $tmp.hhr -B 20 -Z 20",$v2);
	    
#	    # Read hhr file
#	    open(RES,"<$tmp.hhr" ) || die("Error: Cannot open $tmp.hhr: $!\n");
#	    while ($line=<RES>) { if ($line=~/^\s*No Hit/) {last;} }
#	    $line =~ /^(.*) Prob/;
#	    my $cutres=length($1);
#	    $line=<RES>;
#	    close(RES);

#	    $line=~/.{$cutres}\s*(\S+)/;
	    if ($prob>=$Prob) {
		print("\nProbability of best hit = $prob > $Prob. Stopping HHsenser search\n");
		last;
	    } else {    
		print("\nProbability of best hit = $prob ...\n");
	    }
	    $scalar_inX = scalar(@inX);
	}

	# Change $seq to upper case AFTER comparing with %inNeigh
	# (=> loose information of seed-seq aligned to parent alignment!)
	$seq=uc($seq); 

	# Psiblast needs at least 26 residues in query. If query has less than that skip this sequence
	if (($seq=~tr/a-zA-Z/a-zA-Z/)<=25) {next}; # count number of residues
	
	# Build alignment around sequence in $curr
	my $ret=&BuildAlignmentWithSS($seq);
	if ($ret>=4) {next;} # if error occured skip this sequence
	if ($ret>=3) { # Alignment too similar (no new sequences after first search of nr70 => save time and skip
	    if ($v>=2) {print("No new sequences found after first search of $db70 => skip alignment to save time\n");}
	    elsif ($v>=1) {print("Too similar to previous alignment => skip\n");}
	    next;
	} 
	
#	if ($nhits<=1 && $Phh>0.1*$Phh_max) {
#	    print("Alignment in $curr consists of single sequence only. Skip it\n");
#	    unlink(glob("$curr.*"));
#	    next;
#	}
	
	# Test if PSIBLAST iterated to the end without Phh increasing above Phh_max (i.e. $ret<=1)
	if ($ret<=1) { 

	    my $added_seeds_neighbors=0;
	    $corr=-1;
	    for ($rank=1; $rank<=3; $rank++) {

		# Align new file with $tmp-X.hhm
		$Phh = &System("$hh/hhalign -vit -local -id $id -rank $rank -i $tmp"."-X.hhm -t $curr.a3m",$v2);
		if ($v>=2) {print ("returned $Phh\n");}
		if ($Phh>$Phh_max) {last;}

		# Omit sequences with fewer than $rmin residues aligned to master HMM
		unlink("$curr.tmp.a3m");
		if ($use_nr) {
		    &System("$hh/hhalign -local -id 100 -rank $rank -mact $mactX -i $tmp"."-X.hhm -t $curr.a3m -Aa3m $curr.tmp.a3m",$v2);
		    &System("$hh/hhfilter -diff 0 -id 100 -cov $cov_rmin -i $curr.tmp.a3m -o $curr.tmp.a3m",$v2);
		} else {
		    &System("$hh/hhalign -local -id $id -rank $rank -mact $mactX -i $tmp"."-X.hhm -t $curr.a3m -Aa3m $curr.tmp.a3m",$v2);
		    &System("$hh/hhfilter -diff 0 -id $id -cov $cov_rmin -i $curr.tmp.a3m -o $curr.tmp.a3m",$v2);
		}

		# Reject if E-value too high
		if ($Eseed>$Einter) {
		    if ($corr==-1) {
			# Make HMM with parameters as PSIBLAST profile, e.g. constant gap open/extend penalties of 0.5/5.5 bits
			&System("$hh/hhmake -id $id -i $curr.seq -o $curr.seq.hhm -pcm 2 -Blosum65 -gapb 10000 -gapd 2.6 -gape 1 -gapf 1 -gapg 1",$v2); 
			# Calculate correlation between single seed sequence and full alignment built from it 
			$corr_seedHMM_seedSeq = &System("$hh/hhcorr -id $id -i $curr.a3m $curr.seq.hhm -d $calhhm -psi2",$v2);
			$corr_X_parent        = &System("$hh/hhcorr -id $id -i $parent.a3m $tmp"."-X.hhm -d $calhhm",$v2);
			$corr=$corr_seedHMM_seedSeq*$corr_X_parent;
			if ($v>=3) {
			    printf("Correlation coeff. seed_HMM   <-> seed_seq:   %4.2f \n",$corr_seedHMM_seedSeq);
			    printf("Correlation coeff. root-X_HMM <-> parent_HMM: %4.2f \n",$corr_X_parent);
			}
			if ($corr==0) {die("Error: hhcorr returned $corr_seedHMM_seedSeq and $corr_X_parent. Something does not work. Stopping job\n");}
		    }

		    $Pmax=&max($Phh,$Eseed/$Nnr_eff);
		    $Eval = $Pmax**(-$corr) * $Phh*$Eseed;
		    if ($v>=2) {
			printf ("P(H2H)=%7.2G  E(seed)=%7.2G  corr=%4.2f x %4.2f => E-value=%7.2G  ",$Phh,$Eseed,$corr_seedHMM_seedSeq,$corr_X_parent,$Eval);
		    } elsif ($v>=1) {
			printf ("P(HMM-HMM)=%7.2G  E(seed)=%7.2G  corr=%4.2f x %4.2f => E-value=%7.2G  ",$Phh,$Eseed,$corr_seedHMM_seedSeq,$corr_X_parent,$Eval);
			
		    }
		} else {
		    $corr_seedHMM_seedSeq=$corr_X_parent=$corr=0;
		    $Eval=$Eseed;
		    if ($v>=2) {printf ("P(H2H)=%7.2G  E(seed)=%7.2G ",$Phh,$Eseed);}
		    elsif ($v>=1) {printf ("P(HMM-HMM)=%7.2G  E(seed)=%7.2G  ",$Phh,$Eseed);}
		}
	    
		# Accept new alignment?
		if ($Eval<$Ey) {
		    
		    $Nacc++;
		    my $new_parent="$tmp"."-$Nacc";
		    $path[$Nacc] = "#$Nacc <= ".$path[$parent_index];
		    &AddPath("$curr.tmp.a3m",$path[$Nacc]);
		    &AddPath("$curr.a3m",$path[$Nacc]);
		    
		    if ($v>=3) {print("Eval=$Eval  ret=$ret  Phh=$Phh\n");}
		    if (($ret==0 && $Eval<$Einter) && $Phh<=$Phh_strict) { # accept new alignment strictly?

			if ($v>=2) {
			    print("=>  HMM strictly accepted, alignment recorded in $new_parent, and sequences added to $tmp-X.a3m and $tmp-Y.a3m\n\n");
			} elsif ($v>=1) {
			    print("=>  HMM accepted under strict conditions\n");
			}

			
			# Create alignment for root-X.a3m file
			&AddSeqsToHash("$curr.tmp.a3m",\%inX,\@inX,0,1e8);
			&WriteSeqHashToFile (\%inX,\@inX,"$tmp"."-X.a3m");
			
			# Recalibrate $tmp-X.hhm and $tmp-X-noSS.hhm
			&System("$hh/hhmake -id $id -i $tmp"."-X.a3m -o $tmp"."-X.hhm",$v2);
			&System("cp $tmp"."-X.hhm $tmp"."-X-noSS.hhm");
			&System("$hh/hhsearch -cpu $cpu -cal -norealign -i $tmp"."-X.hhm -d $calhhm",$v2);
			&System("$hh/hhsearch -cpu $cpu -cal -norealign -i $tmp"."-X-noSS.hhm -d $calhhm -ssm 0",$v2);

		    } else {  # accept new alignment not so strictly?

			if ($v>=2) {print("=>  HMM accepted, alignment recorded in $new_parent, and sequences added to $tmp-Y.a3m\n\n");}
			elsif ($v>=1) {print("=>  HMM accepted under permissive conditions\n");}
#			# Add sequences with (E<=Emax OR E<=$Emin) to @seeds and add all seqs with qid>=qid_max to neighborhood 
# why this line??	&AddSeqsToHash("$curr.tmp.a3m",\%inNeigh,\@inNeigh,$qid_max,$Emin,"add neigh");
		    }

		    # Create alignment for root-Y.a3m file
		    unlink("$curr.tmp.a3m");
		    &System("$hh/hhalign -id $id -cov 0 -mact $mact -local -rank $rank -i $tmp"."-X.hhm -t $curr.a3m -Aa3m $curr.tmp.a3m",$v2);
		    &AddSeqsToHash("$curr.tmp.a3m",\%inY,\@inY,0,1e8);
		    &WriteSeqHashToFile (\%inY,\@inY,"$tmp"."-Y.a3m");

		    # Add remarks to alignment
		    open (FILE,"<$curr.a3m") or die ("Error: Couldn't open $curr.a3m: $!\n");
		    my @lines=<FILE>;
		    close (FILE);
		    open (FILE,">$new_parent.a3m") or die ("Error: Couldn't open $new_parent.a3m: $!\n");
		    printf(FILE "#ROOT   $file.a3m\n");
		    printf(FILE "#EVALUE %.2G\n",$Eval);
		    printf(FILE "#PATH   %s\n",$path[$Nacc]);
		    printf(FILE "#ESEED  %.2G\n",$Eseed);
		    printf(FILE "#PVALHH %.2G\n",$Phh);
		    printf(FILE "#CORR   %.2f\n",$corr);
		    print(FILE @lines);
		    close (FILE);
		    
		    if (!$added_seeds_neighbors) {
			# Find sequences around current seed alignment with Evalues <= Emax
			$nhits = &System("$perl/reformat.pl -noss a3m fas $new_parent.a3m $new_parent.fas",$v2);
			$nhits = &System("$perl/reformat.pl -M first -r fas fas $new_parent.fas $new_parent.fas",$v2);
			&System("$perl/querydb.pl -a $cpu -Emax $Emax -rmin $rmin -extnd $extnd -sec -d $db70 -q $curr.a3m $new_parent.fas $new_parent.seeds",$v2);
			
			# Add sequences with E<=Emax to @seeds and all seqs with (qid>=qid_max or E-value<=$Emin) to neighborhood
			&AddSeqsToHash("$new_parent.seeds",\%inSeeds,\@inSeeds, 100,$Emax,"add seeds");
			&AddSeqsToHash("$new_parent.seeds",\%inNeigh,\@inNeigh,$qid_max,$Emin,"add neigh");
			$added_seeds_neighbors=1;
			
			open(SEEDFILE,">$tmp.seeds") || die("Could not open $tmp.seeds: $!\n");
			foreach (my $i=0; $i<@seeds; $i++) {
			    printf(SEEDFILE ">%i %s parent=%s\n%s\n",$i,${$seeds[$i]}[3],${$seeds[$i]}[5],${$seeds[$i]}[4]);
			}
			close(SEEDFILE);

			# Adjust value of idmax to number of seeds to process?
			if ($qid_auto) {
			    $nhits=scalar(@seeds);
			    my $qid_tmp;
			    my $Emin_tmp;
			    if ($nhits<=5)       {$qid_tmp=80; $Emin_tmp=0.0;}
			    elsif ($nhits<=10)   {$qid_tmp=60; $Emin_tmp=0.0;}
			    elsif ($nhits<=30)   {$qid_tmp=50; $Emin_tmp=0.0;}
			    elsif ($nhits<=100)  {$qid_tmp=40; $Emin_tmp=1E-50;}
			    elsif ($nhits<=300)  {$qid_tmp=35; $Emin_tmp=1E-30;}
			    elsif ($nhits<=1000) {$qid_tmp=30; $Emin_tmp=1E-20;}
			    elsif ($nhits<=3000) {$qid_tmp=25; $Emin_tmp=1E-10;}
			    elsif ($nhits<=10000) {$qid_tmp=25; $Emin_tmp=1E-6;}
			    else                 {$qid_tmp=25; $Emin_tmp=1E-4;}
			    if ($v>=2) {printf("Choosing max sequence identity of parent seed with new seed idmax=%i%%\n",$qid_max);}
			    if ($qid_tmp<$qid_max) {$qid_max=$qid_tmp; $Emin=$Emin_tmp;}
			}
			
		    }

		} else {last;}
		if ($v>=2) {print("\n");}
	    } # end for($rank=1; $rank<=3; $rank++)
	    
	} 

	if ($ret>=2 || $rank==1) { # test rank, not Phh!

	    if ($v>=1) {printf ("P(H2H)=%7.2G  ",$Phh);}
	    # reject new alignment
	    $nrej{$parent}++;
	    $Nrej++;
	    if ($Phh>$Phh_false && scalar(@inX)>=100) {
		if ($v>=2)  {printf("=>  HMM strongly rejected! Phh=$Phh, |X|=%i\n",scalar(@inX));}
		elsif ($v>=1)  {print("=>  HMM strongly rejected\n");}
		# Add sequences with E<=Eult to neighborhood
		&AddSeqsToHash("$curr.a3m",\%inNeigh,\@inNeigh,100,$Eult,"add neigh");
	    } else {
		if ($v>=2) {print("=>  HMM rejected \n\n");}
		elsif ($v>=1)  {print("=>  HMM rejected\n");}
		# Add sequences with qid>=qid_max to neighborhood and
		&AddSeqsToHash("$curr.a3m",\%inNeigh,\@inNeigh,$qid_max,0,"add neigh");
		if ($use_nr) {
		    &System("$hh/hhfilter -diff 0 -id 100 -cov $cov_rmin -i $curr.a3m -o $curr.a3m",$v2);
		} else {
		    &System("$hh/hhfilter -diff 0 -id $id -cov $cov_rmin -i $curr.a3m -o $curr.a3m",$v2);
		}
		&System("cat $curr.a3m >> $tmp-Z.a3m");
	    }	    
	} 

	if ($v>=2) {print("\n");}
	unlink(glob("$curr.*"));

	# If no more seeds are left and fewer than $Nmin_pfam sequences found: search PfamA/B
	if (@seeds==0 && scalar(@inY)<$Nmin_pfam) {&SearchPfam();}

	# If no more seeds are left and fewer than $Nmin_pfam sequences found: search with whole root-Y for new seeds
	if (@seeds==0 && scalar(@inY)<$Nmin_pfam && $Nacc>=$Nacc_last) {
	    $Nacc_last=$Nacc;

	    # Find new seeds
	    &System("$perl/reformat.pl -noss a3m fas $tmp-Y.a3m $tmp-Y.fas",$v2);
	    &System("$perl/reformat.pl -M first -r fas fas $tmp-Y.fas $tmp-Y.fas",$v2);
	    &System("$perl/querydb.pl -a $cpu -Emax $Emax -rmin $rmin -extnd $extnd -sec -d $db70 -q $tmp-Y.a3m $tmp-Y.fas $tmp-Y.seeds",$v2);
	    
	    # Add sequences with (E<=Emax OR E<=$Emin) to @seeds and add all seqs with qid>=qid_max to neighborhood 
	    &AddSeqsToHash("$tmp-Y.seeds",\%inSeeds,\@inSeeds, 100,$Emax,"add seeds");
	    &AddSeqsToHash("$tmp-Y.seeds",\%inNeigh,\@inNeigh,$qid_max,$Emin,"add neigh");
	}

    } # end while(@seeds)
    
    if ($v>=2) {print("\n");}
    if ($Nacc>=1) {
	if ($use_nr) {
	    &System("$hh/hhfilter -M a3m -diff 0 -id 100 -i $tmp"."-X.a3m -o $tmp"."-X.a3m",$v2);
	    &System("$hh/hhfilter -M a3m -diff 0 -id 100 -i $tmp"."-Y.a3m -o $tmp"."-Y.a3m",$v2);
	} else {
	    &System("$hh/hhfilter -M a3m -diff 0 -id $id -i $tmp"."-X.a3m -o $tmp"."-X.a3m",$v2);
	    &System("$hh/hhfilter -M a3m -diff 0 -id $id -i $tmp"."-Y.a3m -o $tmp"."-Y.a3m",$v2);
	}
    }

    print("\n*************************************************************************************\n");
    print("*************************************************************************************\n");
    print("*************************************************************************************\n");
    printf(" Seeds left:              %i\n",scalar(@seeds)+1);
    printf(" Accepted alignments:     %i\n",$Nacc);
    printf(" Rejected alignments:     %i\n",$Nrej);
    printf(" Neighbor sequences:      %i\n",scalar(@inNeigh));
    printf(" Sequences in strict:     %i\n",scalar(@inX));
    printf(" Sequences in permissive: %i\n",scalar(@inY));
    print("*************************************************************************************\n");
    print("*************************************************************************************\n");
    print("*************************************************************************************\n");

    # Format $tmp-X 
#    &System("$perl/reformat.pl -num a3m fas $tmp"."-X.a3m $tmp"."-X.fas",$v2);
    
    &CleanUp();
    &System("touch $file.finished");
    
} # end for each input alignment

if ($v<3) {rmdir("$tmpdir");}
if ($v>=2) {print ("\nFinished  buildinter.pl @ARGV\n");}
elsif ($v>=1) {print ("\nFinished  buildinter.pl\n");}
exit 0;

# Removed temporary files and copy results files to $dir
sub CleanUp() 
{
    #unlink(glob("*.seeds"));
    &System("cp $tmp"."-X.a3m $dir/");
    &System("cp $tmp"."-Y.a3m $dir/");
    &System("cp $tmp"."-Z.a3m $dir/");
    &System("cp $tmp"."-X.hhm $dir/");
    &System("cp $tmp"."-*.fas $dir/");
    if ($v>=3) {
	print ("Temporary files are not removed\n");
    } else {
	unlink( glob "$tmp*");
    }
}

###############################################################################################
#### Build an alignment around sequence in $curr (called by main)	    
###############################################################################################
sub BuildAlignmentWithSS()
{
    my $qseq=$_[0];
    my $read=0;       # read current line with residues?
    my $iter;         # number of psiblast iterations done
    my @nam;          # names of sequences in alignment  
    my @seq;          # residues of sequences in alignment
    
    if ($v>=1) {
	print("\n\n");
	print("*************************************************************************************\n");
	printf(" Building alignment for %-80.80s\n",$nameline);
	printf(" E-seed:                  %.2G\n",$Eseed);
	printf(" Path:                    %s\n",$path[$parent_index]);
	printf(" Seeds to process:        %i\n",scalar(@seeds)+1);
	printf(" Accepted alignments:     %i\n",$Nacc);
	printf(" Rejected alignments:     %i\n",$Nrej);
	printf(" Rejected from parent:    %i\n",$nrej{$parent});
	printf(" Neighbor sequences:      %i\n",scalar(@inNeigh));
	printf(" Sequences in strict:     %i\n",scalar(@inX));
	printf(" Sequences in permissive: %i\n",scalar(@inY));
	print("*************************************************************************************\n");
    }
    
# The following lines are only needed if the dssp states of the pdb sequence should be included.
# This can only be done if the original capitalization of the seed sequence is conserved and  
# the seed sequence is NOT changed into upper case during reading.
#    # Remove all gaps "-" from $qseq AND corresponding states in $ss_dssp
#    my $qmatch=$qseq;
#    $qmatch=~tr/a-z.//d;
#    $dssp=~tr/a-z.//d;
#    my @qres=split(//,$qmatch);
#    my @qdssp=split(//,$dssp);
#    if ( @qres != @qdssp) {
#	die("Error in $curr.a3m: \n>$name\n$qmatch\n has different number of match states than \n>dssp\n$dssp\n");
#    }
#    $dssp="";
#    $query_length=0;
#    for (my $i=0; $i<@qres; $i++) {
#	if ($qres[$i] ne "-") {
#	    $dssp.=$dssp[$i];
#	}
#    }
#    $query_length = ($qseq=~tr/A-Za-z/A-Za-z/);
#    $qseq=~tr/-//d;; 
#    if ($v>=2) {
#	printf ("Residues not matched to pdb: %i   Residues matched to pdb: %i \n",($qseq=~tr/a-z/a-z/),($qseq=~tr/A-Z/A-Z/));
#    }

    # Write $curr.seq
    my $addwhat = "\#".($Nacc+1)." <= ".$path[$parent_index];
    if ($nameline!~s/(\s+E=\d+)/ ((seed: $addwhat))$1/) {$nameline.=" ((seed: $addwhat))";} 
    open (SEQFILE, ">$curr.seq")|| die("ERROR: can't open $curr.seq: $!\n");
    printf(SEQFILE "%s\n%s\n",$nameline,$qseq); 
    close (SEQFILE);

    # Minimum length of hits is $min_hitlen residues in match states
    my $q_match = ($qseq=~tr/A-Z/A-Z/); #count number of capital letters
    $cov0 = sprintf("%.0f",100*$min_hitlen/$q_match);
    if ($cov0>=80) {$cov0=80;}
    if ($cov>=$cov0) {$cov0=$cov;}
    

    # Set database versions filtered to 90% and 70% sequence identity. Use absolute paths to avoid errors during databas update
    &UpdateDBLinks();

    # Iterative psiblast search with query sequence
    if ($v>=2) {printf ("\nBuilding alignment for $curr with PSI-BLAST ...\n");}
    my $ret=&BuildAlignment("$curr");
    if ($use_nr) {
	&System("$hh/hhfilter -M a3m -diff 0 -id 100 -i $curr.tmp.a3m -o $curr.fil",$v2);
    } else {
	&System("$hh/hhfilter -M a3m -diff 0 -id $id -i $curr.tmp.a3m -o $curr.fil",$v2);
    }
    if ($ret>=2) {
	&System("cp $curr.fil $curr.a3m");
	return $ret;
    }
    &System("$perl/reformat.pl -r a3m psi $curr.fil $curr.psi",$v2);

    $nhits--; # subtract query sequence itself
    
    # Query needs no extension because the seeds around the parent sequence are already extended by querydb.pl

    # Secondary structure prediction with psipred
    if ($v>=2) {printf ("\nPredicting secondary structure with PSIPRED ...\n");}
    $v--;
    &RunPsipred("$curr.seq");
    $v++;
    @nam=();            #names of sequences in alignment  
    @seq=();            #residues of sequences in alignment
    if (open (PSIPREDFILE, "<$curr.horiz")) {
	# Read Psipred file
	my $in;           # input line
	my $aa_pred="";
	my $ss_pred="";
	my $ss_conf="";
	while ($in=<PSIPREDFILE>) {
	    if    ($in=~/^Conf:\s+(\d+)/) {$ss_conf.=$1;}
	    elsif ($in=~/^Pred:\s+(\S+)/) {$ss_pred.=$1;}
	    elsif ($in=~/^  AA:\s+(\S+)/) {$aa_pred.=$1;}
	}
	close(PSIPREDFILE);

	push(@nam,">aa_pred");
	push(@seq,"$aa_pred");
	push(@nam,">ss_pred");
	push(@seq,"$ss_pred");
	push(@nam,">ss_conf");
	push(@seq,"$ss_conf");
    }
    
    # Now write changed psipred states #and dssp states (unchanged)
    if (!open (ALIFILE, ">$curr.a3m")) {
	warn ("ERROR: cannot open $curr.a3m: $!\n");
	return 4;
    }
# The following lines are only needed if the dssp states of the pdb sequence should be included.
#    # Write dssp state sequence into $curr.tmp.a3m
#    print(ALIFILE ">ss_dssp $name $description\n$dssp\n"); # Write query name in first line of file

    # Then write to $curr.tmp.a3m again, setting to lower case those states that are lower case in query
    my @qres=unpack("C*",$qseq); # we need ALL query residues here, also the extended ones!
    for (my $m=0; $m<@seq; $m++) 	{
	my @res=split(//,$seq[$m]);
	for (my $i=0; $i<@res; $i++) {
	    if ($qres[$i]>=97) {	
		$res[$i]=~ tr/A-Z0123456789-/a-z.........../;
	    } else {
		$res[$i]=~ tr/a-z/A-Z/;
	    }
	}
	printf(ALIFILE "%s\n%s\n",$nam[$m],join("",@res));
    }
    close(ALIFILE);

    # Append alignment sequences to dssp- and psipred sequences
    &System("cat $curr.fil >> $curr.a3m");

    # The insert states at the beginning and the end are suppressed in order to conserve
    # the residue numbering when feeding the alignment to hhsearch. Otherwise, shifted 
    # 3D models would result
    open (ALIFILE, "<$curr.a3m");
    my @lines=<ALIFILE>;
    close (ALIFILE);
    open (ALIFILE, ">$curr.a3m"); 
    foreach $line (@lines) {
	chomp($line);
	$line=~s/^[a-z.]*([A-Z0-9-].*[A-Z0-9-])[a-z.]*$/$1/;
	print (ALIFILE "$line\n");
    }
    close (ALIFILE);

#    &System("$perl/reformat.pl -r -num a3m clu  $curr.a3m $curr.clu",$v2);
#    &System("$perl/reformat.pl -r -num a3m fas  $curr.a3m $curr.fas",$v2);
#    &System("$perl/reformat.pl -num a2m fas  $curr.a2m $curr.fas",$v2);

    if ($v<=2) {
	unlink(glob "*.ext.*");
	unlink(glob "*.core.*");
#       unlink("$curr.seq");   #no!
#       unlink("$curr.bla");
#       unlink("$curr.psi");  #no!
	unlink("$curr.tmp");
	unlink("$curr.fil");
	unlink("$curr.horiz");
	unlink("$curr.in.a3m");
	if ($v>=3) {print ("Deleting temporary files\n");}
    } 
    return $ret;
} 


###################################################################################################
# Build alignment for seqfile
###################################################################################################
sub BuildAlignment() {
    my $curr=$_[0];
    my $seqfile=$_[0].".seq";         # Query sequence (or extended query sequence) in fast format
    my $psifile=$_[0].".psi";        # current psiblast alignment in psiblast-readable format
    my $tmpa3mfile=$_[0].".tmp.a3m"; # final psiblast alignment in a3m format
    my $tmpfile=$_[0].".tmp";        # 
    my $blafile=$_[0].".bla";        # BLAST output file
    my $coreali=$_[0].".core.psi";   # current core alignment in psiblast-readable format 
    my $iter;                        # number of psiblast iterations done
    my $bopt="-bl $bl -bs $bs -bg $bg"; # score per col in bits for end pruning PSI-BLAST alignments
    my $bcore="-b 0.67";             # score per col in bits for end pruning of core alignment
    my $db=$db90;                    # do first iteration against nr90
    my $db_short=$db;                # do first iteration against nr90
    $db_short=~s/.*\///;
    my $pmaxopt;                     # in first iteration: "-pmax $pmax", then, if -P option given: "-pmax $pmax -P $psifile"
    my $ret=0;                       # return value of this method. 0: everythink ok. 1:used aggressive E-value threshold up to 1. 2:error
    my $Phh_prev;                    # P-value of previous hhalign command; recorded to see if much worse than $Phh
    my $checked_for_new_seqs=0;      # The first time a search against $db70 was done, check if any new seqs are found and if not reject
    $Phh=1000;

    if ($maxiter==0) {
	system("cp $seqfile $tmpa3mfile");
	return 0;
    } elsif ($maxiter==1) {
	
	&blastpgp("-e $Eult -d $db -i $seqfile","$blafile");
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt -a3m -q $seqfile $blafile $tmpa3mfile",$v2);
	$iter=1;
	if ($v>=1) {
	    printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$Eult,$cov0,$db_short);
	}
	return 0;
    }

    # Iterative psiblast search
    $bopt=~s/\s+-B .*//g;
    $bcore=~s/\s+-B .*//g;
    $iter=0; $nhits=-1;

    # First round with nr90
    $nhits_prev=$nhits; $iter++;
    &blastpgp("-e 1 -d $db -i $seqfile","$blafile");
    &System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bcore $pmaxopt -best -psi -q $seqfile $blafile $coreali",$v2);
    if ($bcore ne "") {$bcore.=" -B $coreali";} # From here on use $coreali for end pruning of $coreali
    if ($bopt ne "")  {$bopt.=" -B $coreali";}  # From here on use $coreali for end pruning of $psifile

    $nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bopt $pmaxopt $best -psi -q $seqfile $blafile $psifile",$v2);
    $nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpa3mfile",$v2);

    # Align new file with $tmp-X-noSS.hhm
#    $Phh_prev=$Phh;
    $Phh = &System("$hh/hhalign -vit -id $id -local -i $tmp"."-X-noSS.hhm -t $tmpa3mfile ",$v2);
    if ($v2>=2) {print ("returned $Phh\n");}
    if ($nhits>=10 && $Phh>$Phh_max) {return 2;} # return if Phh not sufficiently significant

    if ($pmaxopt ne "" && $pmaxopt!~/-P/) {$pmaxopt.=" -P $psifile";} # From here on use $psifile to score match states

    if ($v>=1) {
	printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$E2,$cov0,$db_short);}

   # Second to penultimate round
    while ($nhits<3000 && $nhits_prev<$nhits && $iter<$maxiter-1) {
	
	if ($db eq $db90 && $nhits>=50) { # switch from nr90 to nr70
	    $db=$db70; $db_short=$db; $db_short=~s/.*\///; $nhits=0;
	}


	$nhits_prev=$nhits; $iter++;
	&blastpgp("-e 1 -d $db -i $seqfile -B $psifile","$blafile");

	# Align new file with $tmp-X-noSS.hhm
#	$Phh_prev=$Phh;
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpfile",$v2);
	$Phh = &System("$hh/hhalign -vit -local -i $tmp"."-X-noSS.hhm -t $tmpfile ",$v2);
	if ($v2>=2) {print ("returned $Phh\n");}
	# Phh not sufficiently significant? 
	if ($nhits>=10 && $Phh>$Phh_max) {  # alignment not significant enough and no dearth of seqs?
	    if ($ret==1) {  # has E-value threshold ever been loosened during iterations?
		&System("mv $tmpa3mfile.bk $tmpa3mfile"); 
 		&System("mv $psifile.bk $psifile");       # restore alignment from time BEFORE loosening E-value 
		last;                                     
	    } else {
		return 2;   # => reject alignment
	    }
	}
	
	# If first search against nr70 but no new sequences found => reject for too high similarity
	if ($skip_similar) {
	    if (!$checked_for_new_seqs && $db eq $db70) {
		$/=">"; # set input field seperator
		open(TMP,"<$tmpfile") || die("Error: couldn't open $tmp: $!\n");
		my $returned;
		while (my $seq=<TMP>) {
		    if ($seq eq ">") {next;}
		    while (1) {if($seq=~s/(.)>/$1/) {$seq.=<TMP>;} else {last;}} # in the case that nameline contains a '>'
		    $seq=~s/^(.*)//;      # divide into nameline and residues; '.' matches anything except '\n'
		    my $nameline=">$1";      # don't move this line away from previous line $seq=~s/([^\n]*)//;
		    my $nameid="";
		    if ($nameline=~/^>(\S{1,14}):/ || $nameline=~/^>(\S{1,14})/) {$nameid=$1;} 
		    $seq=~tr/\n> .//d;  # remove all newlines, '.'
		    $returned = &ContainedInHash($nameid,$seq,\%inY,\@inY);
#		printf("%i ",$returned);
		    if ($returned<=1) {last;} # found sequence not yet contained in (%inY,@inY) 
		}
		close(TMP);
#	    print("\n");
		$/="\n"; # set input field seperator
		if ($returned>=2) {return 3;} # found no single sequence not yet contained in (%inY,@inY) => too similar, return
		$checked_for_new_seqs=1;
	    }
	}

	# If less then 10 sequences found and iteration would normally stop...
	my $E3=$E2;
	if ($nhits_prev>=$nhits && $nhits<=10) {
	    # ... aggressively include more sequences by increasing PSI-BLAST E-value threshold 
	    if ($ret==0) { # Is this first time in this BuildAlidngment-call that E-value threshold is increased? => make back-up of files
		&System("cp $tmpa3mfile $tmpa3mfile.bk");
		&System("cp $psifile $psifile.bk");
		$nhits_prev=$nhits=0;
		$ret=1;
	    }
	    my $E3_prev;
	    while ($nhits_prev>=$nhits) { 
		$E3_prev=$E3;
		if (0.01>$E3) {$E3*=10;}
		elsif (0.02>$E3) {$E3=0.02;}
		elsif (0.05>$E3) {$E3=0.05;}
		elsif (0.1>$E3) {$E3=0.1;}
		elsif (0.2>$E3) {$E3=0.2;}
		elsif (0.5>$E3) {$E3=0.5;}
		else {last;}
		# Align new file with $tmp-X-noSS.hhm
		$Phh_prev=$Phh; 
		$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E3 $qid $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpfile",$v2);
		$Phh = &System("$hh/hhalign -vit -local -i $tmp"."-X-noSS.hhm -t $tmpfile",$v2);
		if ($v2>=2) {print ("returned $Phh\n");}
		if ($Phh>$Phh_prev) {$E3=$E3_prev; last;}
	    }	    
	} 

	if ($nhits_prev>$nhits) {last;}
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E3 $qid $sc $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpfile",$v2);
	&System("cat $tmpfile >> $tmpa3mfile");
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E3 $qid $sc $bopt $pmaxopt $best -psi -q $seqfile $blafile $tmpfile",$v2);
	&System("cat $tmpfile >> $psifile");
	
	if ($v>=1) {
	    printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$E3,$cov0,$db_short);
	}
	
	
    } 
    
    # Last round 
    # If last search against nr OR if new sequences were found
    if ($use_nr || $nhits_prev<$nhits || $nhits==0) { 
	if ($use_nr) {$db=$nr; $db_short=$db; $db_short=~s/.*\///;}
	elsif ($db eq $db90 && $nhits>=50) { # switch from nr90 to nr70
	    $db=$db70; $db_short=$db; $db_short=~s/.*\///; $nhits=0;
	}
	
	$nhits_prev=$nhits; $iter++;
	&blastpgp("-e 1 -d $db -i $seqfile -B $psifile","$blafile");
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpfile",$v2);
	if ($nhits_prev>$nhits) {
	    $iter--;
	} else {
	    &System("cat $tmpfile >> $tmpa3mfile");
	    $nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt $best -psi -q $seqfile $blafile $tmpfile",$v2);
	    &System("cat $tmpfile >> $psifile");
	}
    } 

    if ($v>=1) {
	    printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$Eult,$cov0,$db_short);
    }
	
    return $ret;
}


##############################################################################################
# Run SS prediction starting from alignment in $infile.psi (called by BuildAlignment)
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
    my $basename;  # file name without extension
    my $tmpname;   # basename without directory path
    
    
    if ($infile =~/^(.*)\..*?$/)  {$basename=$1;} else {$basename=$infile;}
    if ($basename=~/^.*\/(.*?)$/) {$tmpname=$1;} else {$tmpname=$basename;}
#    print("infile=$infile\nbasename=$basename\nrootname=$tmpname\n");

    # Does dummy database exist?
    if (!-e "$dummydb.phr") {print("WARNING: did not find $dummydb.phr... using $db70\n"); $dummydb=$db70; } 

    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
    if (&SupSystem("$blastpgp -b 0 -j 1 -h 0.001 -d $dummydb -i $infile -B $basename.psi -C $basename.chk &> $basename.bla")) 
	{return;};
    

#    print("Predicting secondary structure...\n");
    
    system("echo $tmpname.chk > $basename.pn\n");
    system("echo $tmpname.seq > $basename.sn\n");
    system("$ncbidir/makemat -P $basename");
    
    &System("$execdir/psipred $basename.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $basename.ss");

    &System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $basename.ss2 $basename.ss > $basename.horiz");
    
    # Remove temporary files
    system ("rm -f $basename.pn $basename.sn $basename.mn $basename.chk $basename.bla $basename.mtx $basename.aux $basename.ss $basename.ss2");
    return;
}

##############################################################################################
# Search PfamA+B for intermediate alignments
##############################################################################################
sub SearchPfam()
{
    if ($v>=2) {
	print("\n***************************************************************************\n");
	print("*******           Searching PfamA+B for intermediate HMMs           *******\n");
	print("***************************************************************************\n\n");
    }
    
    my @hitname=(); # names of all hits in PfamA+B with E-value < Einter
    
    # Search PfamA+B with -X file
    &System("$hh/hhsearch -i $tmp"."-X.hhm -d '$pfamA/db/pfam.hhm' -o $tmp"."-PfamAB.hhr");
    
    # Read all database hits with E-value<$Einter
    open(HHR,"<$tmp"."-PfamAB.hhr") || die("Error: could not open $tmp"."-PfamAB.hhr: $!\n");
    while ($line=<HHR>) {if ($line=~/^\s*No Hit/) {last;}}  # advance until title line
    $line=~/^( No Hit.*Prob)/;
    my $cutres=length($1);
    my $n=0;
    while ($line=<HHR>) {
	if ($line=~/^\s*$/) {last;}
	$line=~/^.{$cutres}\s*(\S+)/;
	if ($1>=$Einter) {next;}
	$line=~/\s*\d+\s+(\S+)/;
	$hitname[$n++]=$1;
    }
    close(HHR);
    if ($v>=2) {print("Found $n hits below E-value threshold\n");}
    
    # Add alignment of each significant hit to $tmp-X.a3m
    for ($n=0; $n<@hitname; $n++) {
	if (defined $hitname{$hitname[$n]}) {next;} # has $hitname already been included in alignment?
	$hitname{$hitname[$n]}=1;
	
	$Nacc++;
	my $new_parent="$tmp"."-$Nacc";
	my $pfamfile;
	if ($hitname[$n]=~/^PF/) {$pfamfile="$pfamA/$hitname[$n]";} else {$pfamfile="$pfamB/$hitname[$n]";} 
	
	# Align pfam file to root-X.a3m file and add sequences to root-X
	unlink("$curr.tmp.a3m");
	&System("$hh/hhalign -id $id -cov 0 -mact $mactX -local -i $tmp"."-X.hhm -t $pfamfile.a3m -Aa3m $curr.tmp.a3m",$v2);
	&AddSeqsToHash("$curr.tmp.a3m",\%inX,\@inX,0,1E8,"found in Pfam");
	&WriteSeqHashToFile (\%inX,\@inX,"$tmp"."-X.a3m");
	
	# Add sequences to root-Y
	&AddSeqsToHash("$curr.tmp.a3m",\%inY,\@inY,0,1E8,"found in Pfam");
	&WriteSeqHashToFile (\%inY,\@inY,"$tmp"."-Y.a3m",);
	
	# Find new seeds
	&System("$perl/reformat.pl -noss a3m fas $pfamfile.a3m $new_parent.fas",$v2);
	&System("$perl/reformat.pl -M first -r fas fas $new_parent.fas $new_parent.fas",$v2);
	&System("$perl/querydb.pl -a $cpu -Emax $Emax -rmin $rmin -extnd $extnd -sec -d $db70 -q $pfamfile.a3m $new_parent.fas $new_parent.seeds",$v2);
	
	# Add sequences with (E<=Emax OR E<=$Emin) to @seeds and add all seqs with qid>=qid_max to neighborhood 
	&AddSeqsToHash("$new_parent.seeds",\%inSeeds,\@inSeeds, 100,$Emax,"add seeds");
	&AddSeqsToHash("$new_parent.seeds",\%inNeigh,\@inNeigh,$qid_max,$Emin,"add neigh");
    }
}
    








###############################################################################################
# Read one sequence after the other from file and store names in hash if
#   qid>qid_min && E<Emax. 
# If sequence name already exists in hash, check overlap of 5-mers between new sequence and 
# sequences stored under this name. If overlap is <10, add new sequence.
# If overlap is >=10 store only longer of the two.
###############################################################################################
sub AddSeqsToHash()
{
    my $seqfile = $_[0];      # add sequences contained in this file to hash
    my $hashref = $_[1];      # ${$hashref}{$nameid} = index of the record for sequences with $nameid in @{\arrayref}
    my $arrayref= $_[2];      # ${$arrayref}{$index} = ">$nameid range_1\nSeq1\n>$nameid range_2\nSeq2\n..."
    my $qid_min = $_[3];      # seq-id with query of sequence to be added must be >=qid_min
    my $Emax    = $_[4];      # E-value of sequence to be added must be <=Emax
    my $addwhat = $_[5];      # "add seeds" if seed seqs are added, "add neigh" if neighbourhod seqs are added, undefined otherwise
    my $nameline="";          # name of sequence
    my $nameid;               # first 14 letters of name
    my $seq;                  # sequence record (name and residues), then only residues
    my $sequence_added;       # Flag that says if the current sequence has been added to the hash 
    my $Evalue;               
    my $parent=$seqfile;      # parent name = basename of seeds file $seqfile
    if ($seqfile =~/(.*)\..*?/) {$parent=$1;}
    $parent=~/\d+$/;
    $parent_index=$1;

    my $i;
    $/=">"; # set input field seperator
    open (FILE, "<$seqfile") || die "ERROR: Couldn't open $seqfile: $!\n";

    while ($seq=<FILE>) {
	if ($seq eq ">") {next;}
	while (1) {if($seq=~s/(.)>/$1/) {$seq.=<FILE>;} else {last;}} # in the case that nameline contains a '>': '.' matches anything except '\n'
	$seq=~s/^(.*)//;      # divide into nameline and residues; '.' matches anything except '\n'
	$nameline=">$1";      # don't move this line away from previous line $seq=~s/([^\n]*)//;
	if ($nameline=~/^>aa_/) {next;}
	if ($nameline=~/^>ss_/ && defined $addwhat) {next;}
	if ($nameline=~/^>(\S{1,14}):/ || $nameline=~/^>(\S{1,14})/) {$nameid=$1;} 
	else {
	    print("Error in AddSeqsToHash():\nfile=$seqfile\nnameline=$nameline\nseq=$seq\ninput line=$.  Hit Return to continue\n"); 
	    <STDIN>;
	    next;
	}
	$seq=~tr/\n> .//d;  # remove all newlines, '.'
	
	# Does the sequence fulfill the E-value and seq.id criteria?
	if ($nameline=~/^>.* E=(\S+)\s+.* id=(\d+)/) {
	    $Evalue=$1;
	    if ($2<=$qid_min && $Evalue>=$Emax) {next;}
	} else {
	    $Evalue=0; # no E-value is given when the original a3m file is read into %inX, @inX and %inY, @inY
	}
	
	
	if ($seq!~/^[a-z]*([A-Z0-9-]*\S*[A-Z0-9-])[a-z]*$/) {
	    print("Error in AddSeqsToHash():\nfile=$seqfile\nnameline=$nameline\nseq=$seq\ninput line=$.  Hit Return to continue\n"); 
	    <STDIN>;
	    next;
	}
	my $center=$1;

	# Skip sequence in new alignment if too low coverage (too few match chars)
	my $matches=($center=~tr/A-Z0-9/A-Z0-9/);
	my $tmplen=($center=~tr/A-Z0-9-/A-Z0-9-/);
	if ($matches/$tmplen<0.01*$cov) {next;}
	
	if (defined $addwhat) {
	    # If neighborhood or seed seqs are stored, shorten nameline and set name and description
	    if ($addwhat eq "add seeds" || $addwhat eq "add neigh") { 
		$nameline=~/^(>\S+)\s*(.*)/;
		$name=$1;	
		$description=$2;
		if ($addwhat eq "add seeds") {
		    $nameline=sprintf("%-80.80s",$nameline);
		} else {
		    $nameline=$1;
		}
	    } 
	} 

	$sequence_added=0;
	my $index=${$hashref}{$nameid};
	if (defined $index) {
	    if ($nameline=~/^>ss_/) {next;}
#	    print("domains:\n${$arrayref}[$index]\n");
	    my @domains=split(/>/,${$arrayref}[$index]);
#	    for($i=1; $i<scalar(@domains); $i++) {printf("domain[%i]=%s\n",$i,$domains[$i]);} # DEBUG
	    for($i=1; $i<scalar(@domains); $i++) {$domains[$i]=">".$domains[$i];}
	    for($i=1; $i<scalar(@domains); $i++) {
		$domains[$i]=~/^.*\n[a-z]*([A-Z0-9-]*\S*[A-Z0-9-])[a-z]*$/ || die("Error: domain=$domains[$i]\n");
		my $centerdom=$1;
		if (&Overlap($center,$centerdom)>=10) {
		    if ($matches>($centerdom=~tr/A-Z0-9/A-Z0-9/)) {
			if ($v>=4) {printf("Replacing \n%swith\n%s\n%s\n",$domains[$i],$nameline,$seq);}
			$domains[$i]=sprintf("%s\n%s\n",$nameline,$seq);
		    }
		    last;
		}
	    }
	    if ($i>=scalar(@domains)) {
		push(@domains,sprintf("%s\n%s\n",$nameline,$seq)); 
		$sequence_added=1;
	    }
	    ${$arrayref}[$index]=join("",@domains);
#	    for($i=1; $i<scalar(@domains); $i++) {printf("domain_after[%i]=%s\n",$i,$domains[$i]);} # DEBUG
	} else {
	    ${$hashref}{$nameid}=scalar(@{$arrayref});
	    push(@{$arrayref},sprintf("%s\n%s\n",$nameline,$seq));
	    $sequence_added=1;
	}

	# Sequences in $seqfile are seeds? Then store: nameline, name, desc, seq, parent id and Evalue of seed
	if (defined $addwhat && $addwhat eq "add seeds" && $sequence_added==1) {
	    # Use $seq as seed sequence
	    $seq=~tr/Uu.-/Cc/d;   # replace selenocysteine by cysteine (to prevent blastpgp -B from crashing)
	    my @seed=($nameline,$name,$description,$nameid,$seq,$parent,$Evalue);
	    if ($quick && $Evalue<0.1 && $Nacc+$Nrej==0) {
		unshift(@seeds,\@seed); # seeds with E-values up to 0.1 will be processed first, in order of DEcreasing E-values
	    } else {
		push(@seeds,\@seed);    # seeds will be processed in order of INcreasing E-values (decreasing significance)
	    }
	}


    }
    close FILE;
    $/="\n"; # reset input field seperator
    return;
}

###############################################################################################
##  Returns  0 if $nameid is NOT contained in %{$hashred}
##  Returns  1 if $nameid IS contained in %{$hashred} BUT $seq has no overlap>=10 with any sequence in @{$arrayref}, 
##  returns  2 if $nameid IS contained in %{$hashred} AND $seq has overlap >=10 
###############################################################################################
sub ContainedInHash()
{
    my $nameid  = $_[0];      # nameid to be tested
    my $seq     = $_[1];      # sequences to be tested
    my $hashref = $_[2];      # ${$hashref}{$nameid} = index of the record for sequences with $nameid in @{\arrayref}
    my $arrayref= $_[3];      # ${$arrayref}{$index} = ">$nameid range_1\nSeq1\n>$nameid range_2\nSeq2\n..."
    my $i;

    $seq=~/^[a-z]*([A-Z0-9-]*\S*[A-Z0-9-])[a-z]*$/ 
	|| die("Error in ContainedInHash():\nnameid=$nameid\nseq=$seq\n");
    my $center=$1;
    my $index=${$hashref}{$nameid};
    if (defined $index) {
	my @domains=split(/>/,${$arrayref}[$index]);
	for($i=1; $i<scalar(@domains); $i++) {$domains[$i]=">".$domains[$i];}
	for($i=1; $i<scalar(@domains); $i++) {
	    $domains[$i]=~/^.*\n[a-z]*([A-Z0-9-]*\S*[A-Z0-9-])[a-z]*$/ || die("Error: domain=$domains[$i]\n");
	    my $centerdom=$1;
	    if (&Overlap($center,$centerdom)>=10) {
		return 2;
	    }
	}
	if ($i<scalar(@domains)) {return 2;}
	return 1;
    } 
    return 0;
}

###############################################################################################
# Write all sequences in hash to $tmpX
###############################################################################################
sub WriteSeqHashToFile() {
    my $hashref = $_[0];
    my $arrayref = $_[1];
    my $file = $_[2];
    my $nameid;               # first 14 letters of name
    open (FILE, ">$file") || die "ERROR: Couldn't open $file: $!\n";
#    foreach $nameid (keys(%{$hashref})) {print("$nameid  ");} print("\n");
    for (my $index=0; $index<@{$arrayref}; $index++) {
	printf(FILE "%s",${$arrayref}[$index]);
    }
    close FILE;
}

###############################################################################################
# Write all sequences in hash to $tmpX
###############################################################################################
sub AddPath() {
    my $file = $_[0];      # add path to sequences contained in this file 
    my $path = $_[1];
    open (FILE, "<$file") || die "ERROR: Couldn't open $file: $!\n";
    my @lines=<FILE>;
    close FILE;
    open (FILE, ">$file") || die "ERROR: Couldn't open $file: $!\n";
    foreach my $line (@lines) {
	if ($line=~/^>/) {
	    # Remove previous ((path)), but not ((seed: ...))
	    if ($line!~/\(\(.*\)\)/) {
		# Add new ((path))
		if ($line!~s/(\s+E=\d+)/ (($path))$1/) {
		    chomp($line); $line.=" (($path))\n"; 
		}
	    }		
	}
	print(FILE $line);
    }
    close FILE;
}


# Minimum
sub max {
    my $max = shift @_;
    foreach (@_) {
	if ($_>$max) {$max=$_} 
    }
    return $max;
}

################################################################################################
### Execute blastpgp
###########################################################################################
#####
sub blastpgp() {
    if ($v>=2) {printf("\$ $blastpgp -a $cpu -b 20000 -v 1 %s &> $_[1]\n",$_[0]);} 
    if (system("$blastpgp -a $cpu -b 20000 -v 1 $_[0] &> $_[1]")) {
	die("Error: PSI-BLAST returned with error. Consults blast output in $_[1]\n\n");
    }
    return;
}

################################################################################################
### System command with return value parsed from output
################################################################################################
sub System()
{
    my %regexs=(
	"alignhits.pl"=>"^\(\\d+\) sequences extracted",
	"reformat.pl" =>"^Reformatted \\S+ with \(\\d+\)",
	"querydb.pl"  =>"^Database file \\S+ with \(\\d+\)",
	"hhfilter"    =>"^\(\\d+\) out of \\d+ sequences passed filter",
	"hhmake"      =>"^\(\\d+\) out of \\d+ sequences passed filter",
	"hhsearch"    =>"^  1 \\S.\{29\}\\s+\(\\d+\.\\d+\)",
	"hhalign"     =>"^Aligned .* P-value = \(\\S+\)",
	"hhcorr"      =>"^Correlation coefficient <~ \(\\S+\)",
	);
    my $v2;
    if ($_[1]) {$v2=$_[1];} else {$v2=0;}
    $_[0]=~/^(\S+)/;
    my $program=$1;
    $program=~s/^\S+\///; # remove path
#    printf("program=$program  regex=%s\n",$regexs{$program});
    if(defined $regexs{$program}) {
	my $line;
	if ($v>=2) {printf("\$ %s",$_[0]); if ($v2>=2) {print("\n");}} 
	my $regex=$regexs{$program};
	open(SYS,"$_[0] |") || die("\nError in system(\"$_[0] |\"): $!\n");
	while ($line=<SYS>) {
	    if ($v2>=2) {print($line);}
	    if ($line=~/$regex/) {last;}
	}
	close(SYS);
	if (!$line) {die("\nError in $0: could not parse the result of $_[0]: regex=$regex\n");}
	$line=~/$regex/;
	if ($v>=2 && $v2<2) {printf(" -> $1\n");} 
 	return $1;
    } else {
	if ($v>=2) {printf("\$ %s\n",$_[0]);} 
	return system($_[0])/256;
    }
}


################################################################################################
### Supervised system command where process is killed after maximum time threshold 
################################################################################################
sub SupSystem()
{
    my $command=$_[0];
    my $FIRSTTRY=900;  
    my $sys=0;
    $SIG{'CHLD'}='IGNORE'; # remove child from process table after termination

    if ($v>=2) {print("$command\n");} 
    if ($pid=fork()) {
	# Parent process for supervision of child process
	my $time=0;
	while (kill(0,$pid) && $time<$FIRSTTRY) { # kill($pid,0) returns true while process is running 
	    sleep(1);
	    if ($v>=2) {print(".");}
	    $time++;
	}
	if ($time>=$FIRSTTRY) { 
	    kill(-9,$pid); # send child process the KILL signal
	    if ($v>=1) {print(" killed after $time seconds !!!!!!!!!!!!!!!!!!!!!!\nStarting second try ... \n");}
	    $sys=1;
	} else {
	    kill(-9,$pid); # send child process the KILL signal
	    if ($v>=2) {print("\n");}
	}
    } else {
	# Child process for system command
	system($command);
	exit 0;
    }
    
    $SIG{'CHLD'}='DEFAULT'; # reset signal handler; otherwise sytem() will not return Perl script's exit value!!)
    return $sys;
}


################################################################################################
### Set database versions filtered to 90% and 70% sequence identity. Use absolute paths to avoid errors during databas update
################################################################################################
sub UpdateDBLinks() 
{
    $db90 = $dbbase."90";  
    $db70 = $dbbase."70";  
    $db   = $dbbase;  
    if    (-l $db90  && defined readlink($db90)) {$db90=readlink($db90);}
    elsif (-l "$db90.pal"  && defined readlink("$db90.pal")) {$db90=readlink("$db90.pal"); $db90=~s/\.pal$//;}
    if    (-l $db70 && defined readlink($db70)) {$db70=readlink($db70);}
    elsif (-l "$db70.pal"  && defined readlink("$db70.pal")) {$db70=readlink("$db70.pal"); $db70=~s/\.pal$//;}
    if    (-l $db  && defined readlink($db))  {$db =readlink($db);}
    elsif (-l "$db.pal"  && defined readlink("$db.pal")) {$db=readlink("$db.pal"); $db=~s/\.pal$//;}
}


sub ChildTerminated()
{
    if ($pid) {kill(-9,$pid);}
    return;
}


# Calculate the number of identical 5-tupels in $x and $y; remove extended ends from $x
sub Overlap() 
{
    my $x=$_[0];      # string 1
    my $y=$_[1];      # string 2
    my $lentup=5;     # length of tupels
    my $ovlap=0;      # number of 5-tupels occuring in both $x and $y

    $x=~s/^[a-z]*//;     # remove all lower case residues at the beginning of $x
    $x=~s/[a-z]*$//;     # remove all lower case residues at the end of $x
    $x=~tr/a-z .-/A-Z/d; # remove all gaps and turn into upper case
    $y=~s/^[a-z]*//;     # remove all lower case residues at the beginning of $y
    $y=~s/[a-z]*$//;     # remove all lower case residues at the end of $y
    $y=~tr/a-z .-/A-Z/d; # remove all gaps and turn into upper case
    for (my $pos=0; $pos<length($x)-$lentup+1; $pos++) {
	my $substr=substr($x,$pos,$lentup);
	if ($y=~/$substr/) {$ovlap++};
    }
    return $ovlap;
}

sub Terminate()
{
    &CleanUp();
    die("\nHHsearch was forced to terminate. Preliminary results files copied to $dir/\n");
}
