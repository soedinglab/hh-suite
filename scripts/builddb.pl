#!/usr/bin/perl -w
# Build a database from a fasta sequence file
# Usage: builddb.pl infile outdir
#
# $db90 must have been formatted with -o option for fastacmd to work: 'formatdb -i $db90 -o T' 
my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;


# CUSTOMIZE
my $dbbase=$nre; 
#my $dbbase="/home/soeding/nr/nre"; # will use PSI-BLAST dbs $dbbase."90" and $dbbase."70" 
                                    # It is recommended to use nre = nr + env. 
                                    # 90: filtered with CD-HIT to 90% maximum pairwise sequence identity

# Default values:
our $v=2;            # verbose mode

# Default values for building alignments
my $maxiter=8;       # maximum number of psiblast iterations
my $forceext=0;      # search with extended sequence and with original
my $E2=1E-3;         # Psiblast Evalue threshold for building profile
my $Eult=1E-3;       # ultimate Evalue threshold for building output alignment
my $id=90;           # maximum sequence identity
my $cov=20;          # minimum coverage
my $min_hitlen=0;    # minimum number of residues in match states 
my $Ndiff=0;         # number of maximally different sequences
my $bl=0.0;          # lenient minimum per-residue bit score with query at ends of HSP
my $bs=0.167;        # strict  minimum per-residue bit score with query at ends of HSP
my $bg=20;           # maximum number of HSP residiues missing at ends of query to use lenient $bl
my $qid="";          # minimum sequence identity with query
my $sc="";           # minimum score per column with query
my $pmax=1e-7;       # maximum p-value of HSP IN MATCH COLUMNS for inclusion into alignment
my $best="";         # extract only the best HSP per sequence, except during last round
my $names=0;         # Derive file names from 0:sequence names, 1: family ids, 2:input file name 3:input file name + seqnames
my $rmin=26;         # below this number of match state residues, sequence will be removed from PSI-BLAST jumpstart file   
my $ihs_min=-1;      # if less sequences found, hhsenser (buildinter.pl) is called
my $do_dssp=1;       # DO determine DSSP states by default
my $noss=0;          # default= include predicted and DSSP secondary structure
my $cpu=1;           # no multithreading
my $outdir=".";      # output directory where a3m alignment is written
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID

# Create directory path
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 


my $usage="
Build a database of alignments in a3m format from the sequences in the fasta input file.
Read sequences from file one by one and 
* build alignment from original query sequence
* if query alignment contains <20 sequqnces and query has <50 residues and query can be extended 
  by more than 10\% of residues then do search with extended query and merge alignments
* include dssp states
* include psipred secondary structure prediction  

Usage:   builddb.pl infile outdir [options] 

General options:
 -v   <int>   verbose mode (def=$v)
 -u           update: do not overwrite *.a3m files already existing (def=off)
 -cpu <int>   number of CPUs to use when calling blastpgp
 -fam         derive file names from family ids (for SCOP sequences) (default: derive from sequence names)
 -index       derive file names by appending a file index to the name of the input file

Options for building alignments:
 -n   <int>   maximum number of psiblast iterations  (def=$maxiter)
 -e   <float> E-value for inclusion in PSI-BLAST profile (def=$Eult)
 -id  <int>   maximum sequence identity in % (def=$id)
 -cov <int>   minimum coverage in % (Coverage = length of HSP / length of query) (def=$cov)
 -len <int>   minimum number of residues in HSP (def=$min_hitlen)
 -bl  <float> lenient minimum per-residue bit score for pruning ends of HSP (default=$bl)
 -bs  <float> strict  minimum per-residue bit score for pruning ends of HSP (default=$bs)
 -bg  <float> use the lenient bl below this number of end gaps and the strict bs otherwise (default=$bg)
 -p   <float> only for extended search: maximum p-value of HSP IN MATCH COLUMNS for inclusion into alignment (def=$pmax)
 -extend      force query extension: psiblast with original AND extended sequence and merge alignments (def=off)
 -lc           filter out low complexity regions in query sequence (only for PSI-BLAST search)
 -ihs <int>   run quick intermediate HMM search (HHsenser) if less than <int> sequences found (def=off)

Example: nohup nice -19   builddb.pl scop70.1.63 . > ./builddb.log &
\n";

# Variable declarations
my $tmp="delme.tmp";
my $outbase;           # basename (no extension) for output file name
my $line;              # input line
my %nfiles=();         # how many files already exist for a certain scop classification?
my $q_match;           # number of match states in query
my $qseq;              # residues of query sequence
my $nameline;          # query in fasta format: '>$nameline\n$qseq'
my $name;              # $nameline="$name $description". 
my $description;       # For astral sequences: $name="$scopid $famid $range"
my $root;              # root name of current file
my $nhits_last=0;
my $famid="";          # for SCOP sequences: the family id of current sequence
my $update=0;          # 0:overwrite  1:do not overwrite
my $nfile=0;           # number of alignent file currently being generated
my $xseq;              # sequence x returned from Align.pm
my $yseq;              # sequence y returned from Align.pm  
my $Sstr;              # match sequence returned from Align.pm
my $tseq="";           # residues of template sequence (for ExtendSequence() etc.)
my $nhits;             # number of hits found by alignhits
my $nhits_prev;        # number of nhits from previous psiblast round
my $cov0;              # effective minimum coverage (including $mi_-hitlen threshold)
my $infile="";
my $nseqin=1;          # Needed to use same subroutine as for buildali.pl
my $ss_dssp;           # dssp states as string
my $aa_dssp;           # residues from dssp file as string
my $aa_astr;           # residues from infile as string
my $ss_pred="";        # psipred ss states
my $ss_conf="";        # psipred confidence values
my $aa_pred="";        # residues from psipred file
my $foundX=0;          # found an X that stands for 'inserted domain' in SCOP/ASTRAL sequences?
my $blastpgp=$ncbidir."/blastpgp -I T -s T"; # show gi's in defline; use Smith-Waterman
my $rfilt=$bioprogs_dir."/rfilt"; # rfilt is a inhouse extended version of David Jones' pfilt
my $lc_filter=0;       # default: low-complexity filter on
my ($db90, $db70, $db);

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($usage);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

# General options
if ($options=~s/ -v\s*(\d) / /g) {$v=$1;}
if ($options=~s/ -v / /g) {$v=2;}
if ($options=~s/ -u / /g) {$update=1;}
if ($options=~s/ -cpu\s+(\d+) / /g) {$cpu=$1; $blastpgp.=" -a $1";}
if ($options=~s/ -names\s*(\S+) / /g)   {$names=$1;}
if ($options=~s/ -fam / /g)   {$names=1;}
if ($options=~s/ -index / /g) {$names=2;}

# Options for building alignments:
if ($options=~s/ -n\s+(\d+) / /g) {$maxiter=$1;}
if ($options=~s/ -extend / /g) {$forceext=1;}
if ($options=~s/ -slow / /g) {$forceext=1;}
if ($options=~s/ -fast / /g) {;}
if ($options=~s/ -e\s+(\S+) / /g)   {$Eult=$1; $E2=$1;}
if ($options=~s/ -e2\s+(\S+) / /g)  {$E2=$1;}
if ($options=~s/ -id\s+(\d+) / /g) {$id=$1;}
if ($options=~s/( -qid\s+\S+) / /g) {$qid=$1;}
if ($options=~s/ -cov\s+(\d+) / /g) {$cov=$1;}
if ($options=~s/ -len\s+(\d+) / /g) {$min_hitlen=$1;}
if ($options=~s/ -best / /g) {$best="-best";}
if ($options=~s/( -s\/c\s+\S+) / /g) {$sc=$1;}
if ($options=~s/ -b\s+(\S+) / /g)   {$bl=$bs=$1;}
if ($options=~s/ -bl\s+(\S+) / /g)  {$bl=$1;}
if ($options=~s/ -bs\s+(\S+) / /g)  {$bs=$1;} 
if ($options=~s/ -bg\s+(\S+) / /g)  {$bg=$1;} 
if ($options=~s/ -pmax\s+(\S+) / /g){$pmax=$1;}
if ($options=~s/ -lc / /g){$lc_filter=1;}
if ($options=~s/ -ihs\s+(\S+) / /g){$ihs_min=$1;}
if ($options=~s/ -p\s+(\S+) / /g){$pmax=$1;}
if ($options=~s/ -dssp / /g){$do_dssp=1;}
if ($options=~s/ -noss / /g) {$noss=1;}

# Set input and output file and tmp dir
if ($options=~s/ -i\s+(\S+) / /g) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) / /g) {$outdir=$1;}
if ($options=~s/ -tmp\s+(\S+) / /g) {$tmpdir=$1;}
if (!$infile  && $options=~s/^\s*([^- ]\S*)\s+/ /) {$infile=$1;} 
if ($options=~s/^\s*([^-]\S*)\s* / /) {$outdir=$1;} 

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!$infile) {print($usage); exit(1);}

# Warn if unknown options found
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; print("WARNING: unknown options '$options'\n");}

my $v2 = $v-1;
if ($v2>2) {$v2=2;}
if ($v2<0) {$v2=0;}

if ($v>=3) {$|= 1;} # Activate autoflushing on STDOUT

# Set database versions filtered to 90% and 70% sequence identity. Use absolute paths to avoid errors during database update
&UpdateDBLinks();

if ($v>=2) {printf("Hostname=%s",`hostname`);}

my $nfiless=0;
###############################################################################################
# Main loop: write each sequence in infile into its file $tmp.seq and BuildAlignment for it
###############################################################################################
my $inroot; # root name of infile: remove path and extension
if ($infile=~/(.*)\..*?/) {$inroot=$1;} else {$inroot=$infile;} # remove extension
if ($inroot=~/.*\/(.*)?/)  {$inroot=$1;}                        # remove path 
my $nseq=0;


open (INFILE, "<$infile") || die ("ERROR: cannot open $infile: $!\n");
while (1) #assign next line of file to $inline while there is a next line
{
    $line=<INFILE>;
#    print("Reading $line");
    if( (!$line) || ($line=~/^>/) )
    {
	if ($tmp ne "delme.tmp") {

	    $q_match = ($qseq=~tr/A-Z/A-Z/);    # count number of capital letters
	    if ($q_match<=25) {                 # Psiblast needs at least 26 residues in query
		$qseq.='X' x (26-$q_match);     # add 'X' to query to make it at least 26 resiudes long
		print("Resiudes='$qseq'\n"); 
		printf(ALIFILE "%s\n",'X' x (26-$q_match));
	    }
	    close ALIFILE;
	    
	    # Filter query sequence file (if not SCOP sequence) (only own short-repeat filtering -r)
	    if (!$famid) {
		&System("cp $tmp.seq $tmp.org.seq");
		&System("$rfilt -t -c -x -r $tmp.org.seq > $tmp.seq");
		open(FILTFILE,"<$tmp.seq") || die("Error: could not execute $rfilt -t -c -x -r $tmp.seq: $!\n");
		$nameline=<FILTFILE>; # read first line
		undef $/;
		$qseq = <FILTFILE>; # slurp in all rest of sequence
		$/="\n";
		close(FILTFILE);
		$qseq=~tr/\n//d;
	    }

	    
	    &BuildAlignmentWithSS();
	    
	    # If not enough sequences in alignment, do intermediate HMM search
	    if ($nhits<$ihs_min) {
		print("\n**************************************************************************************\n");
		print("********  Using intermeditate HHM search (HHsenser) to broaden alignment  ************\n");
		print("**************************************************************************************\n");
		&System("$hh/buildinter.pl -cpu $cpu -quick $outbase.a3m",$v);
		&System("mv $outbase.a3m $outbase.a3m.orig");
		&System("mv $outbase"."-X.a3m $outbase.a3m");
		print("\n********               Finished intermeditate HHM search                  ************\n");
		print("**************************************************************************************\n");
	    }
	} 
	if(!$line) {last;} # if no more new line => stop reading in

	# Read nameline and extract SCOP family code
	if ($line=~/^>(\S+)\s+(.*)/) {
	    $name=$1;
	    $description=$2;
	    $nameline=$line;
	    if ($line=~/^>[a-z]\S+\s+([a-k]\.\d+\.\d+\.\d+)\s+\S+\s+/){$famid=$1;} else {$famid="";}

	    if ($names==1 && $famid ne "") {
		# name file after family id
		# determine number $nfiles{$famid} of sequences already found with same scop-id	    
		if (exists $nfiles{$famid}) {$nfiles{$famid}++;} else {$nfiles{$famid}=1;}
		$root="$famid.$nfiles{$famid}";
	    } elsif ($names==2) {
		# name file after input file plus index $nseq
		if ($nseq) {$root="$inroot"."_$nseq";} else { $root="$inroot";}
		$nseq++;
	    } elsif ($names==3) {
		# name file after input file plus index sequence name
	        $root=$inroot."-".$name;
		$name=$inroot.":".$name;
		$nseq++;
	    } else {
		# name file after sequence name
		$root=$name;
		if ($root!~s/^gi\|(\d+)\|.*/gi_$1/) {$root=~tr/a-zA-Z0-9+_.-//cd;}
	    }
	    $tmp="$tmpdir/$root";
	    $outbase="$outdir/$root";
	}

	if (-s "$outbase.a3m" && $update==1) { # does the a3m file already exist?
	    if ($v>=3) {print("$outbase.a3m already exists... skipping...\n");}
	    $root="delme.tmp"; # do not overwrite existing file
	    $tmp="$tmpdir/$root";
	    open(ALIFILE,"$tmp");
	} else {	
	    # Open new file to hold query sequence. Later used to append aligned hits 
	    open (ALIFILE, ">$tmp.seq")|| die("ERROR: can't open $tmp.seq: $!\n");
	    print(ALIFILE "$nameline"); # Write query name in first line of file
	    $qseq="";
	}
    }
    elsif ($root ne "delme.tmp") {
#	print("Residues: $line");
#	$line=~tr/a-z/A-Z/; # transform into upper case
	chomp($line);
	if ($line=~tr/a-zA-Z//cd) {print(STDOUT "WARNING: removed invalid symbols from sequence $name: $line\n");}
#	$line=~tr/Uu/Cc/;        # BLAST bug removed; not necessary anymore
	if ($famid) {            # If query is scop sequence replace ahdiemXdgavfsd by AHDIEMxDGAVFSD: 
	    $foundX+=($line=~tr/X/*/);     # Replace all domain insertion codes X by *,
	    $line=uc($line);     # ... turn sequence into upper case letters,
	    $line=~tr/*/x/;      # ... and transform * to x. In ExtendSequence() split on /x/ pattern
	    $qseq.=$line; 
	    $line=~tr/x//d;      # remove small x from $qseq (otherwise PSIPRED will contain one letter too much)
	    printf(ALIFILE "%s\n",$line);
	} else {
	    $qseq.=uc($line);
	    printf(ALIFILE "%s\n",uc($line));
	}
    }
}
close ALIFILE;
close INFILE;

if ($v<3 && $tmpdir ne $outdir) {rmdir("$tmpdir");}
if ($v>=2) {print ("\nFinished  builddb.pl @ARGV\n");}
elsif ($v>=1) {print ("\nFinished  builddb.pl\n");}
exit 0;


###############################################################################################
#### Build an alignment around sequence in $tmp (called by main)	    
###############################################################################################
sub BuildAlignmentWithSS()
{
    my $read=0;         # read current line with residues?
    my $length=0;       # length of last sequence printed into file
    my $iter;           # number of psiblast iterations done
    my @nam;            # names of sequences in alignment  
    my @seq;            # residues of sequences in alignment
    my $nseq=0;         # number of sequences read in from alignment
    
    $nfile++;
    if ($v>=2) {
	print("\n");
	print("************************************************\n");
	print(" Building alignment for $tmp ($name)\n");
	print("************************************************\n");
    }

    # Minimum length of hits is $min_hitlen residues in match states
    $cov0 = sprintf("%.0f",100*$min_hitlen/$q_match);
    if ($cov0>=80) {$cov0=80;}
    if ($cov>=$cov0) {$cov0=$cov;}
    
    if (($qseq=~tr/A-WYZ/A-WYZ/)>=20) {
	# Start from single sequence
	
	# Iterative psiblast search with ORIGINAL sequence
	if ($v>=1) {printf ("\nBuilding alignment for query with PSI-BLAST ...\n");}
	&BuildAlignment("$tmp",0);
	
	# Extend query? (extended residues are lower case)
	if ($v>=3) {print("forceext=$forceext  foundX=$foundX  nhits=<$nhits  q_match=$q_match\n");}
	if ($forceext || ($foundX && $nhits<150) || ($nhits<20 && $q_match<=150)) { # if fewer than ??? hits found     
	    if ($v>=1) {printf ("\nExtending query sequence ...\n");}
	    my $eseq=&ExtendSequence($db90);
	    if ($v>=2) {printf ("Extended residues=%i  original residues=%i\n",($eseq=~tr/a-z/a-z/),($eseq=~tr/A-Z/A-Z/));}
	    
	    # Iterative psiblast search with EXTENDED query sequence? (extensions in lower-case)
	    if (($eseq=~tr/a-z/a-z/) > 0.1*($eseq=~tr/A-Z/A-Z/)) {
		if ($v>=1) {printf ("\nBuilding alignment for extended query with PSI-BLAST ...\n");}
		&BuildAlignment("$tmp.ext",1);
		&System("$hh/hhfilter -M a3m -id $id $qid -cov $cov0 -diff $Ndiff -i $tmp.ext.a3m -o $tmp.fil",$v2);
		
		# Remove all insert states from $tmp.seq to be able to call psipred
		&System("$perl/reformat.pl -r a3m a3m $tmp.seq $tmp.seq",$v2);
		
		# Append results for extended query to results for original query
		&System("cat $tmp.fil  >> $tmp.a3m");
		if ($v>=2) {print(STDOUT "\n");}
	    }
	}

    } else {
	# Too few non-X residues left after filtering -> skip PSI-BLAST search
	if (!$famid) {
	    &System("cp $tmp.org.seq $tmp.seq");
	}
	&System("cp $tmp.seq $tmp.a3m");
    }

#    $qseq=~tr/a-z//d; # remove 'x' symbols that signify domain insertions in scop/astral sequences

    # Filter for maximum sequence identity
    $nhits=&System("$hh/hhfilter -M a3m -id $id $qid -cov $cov0 -diff $Ndiff -i $tmp.a3m -o $tmp.fil",$v);
    &System("$perl/reformat.pl -r a3m psi $tmp.fil $tmp.psi",$v2);

    # Count sequences 
    open (FILE, "<$tmp.fil") or die ("ERROR: cannot open $tmp.fil: $!\n");
    my @lines=<FILE>;
    $nhits=scalar(@lines)/2;
    close(FILE);
    @lines=();

    @nam=();            #names of sequences in alignment  
    @seq=();            #residues of sequences in alignment

    unlink ("$tmp.a3m");
    if (!$noss) { 
	
	# Read query sequence and write dssp state sequence (plus dssp residues for verif) into $tmp.a3m
	if ($do_dssp && $dsspdir ne "") {
	    if (!&AppendDsspSequences("$tmp.seq")) {
		push(@nam,">ss_dssp"); 
		push(@seq,"$ss_dssp");
#		push(@nam,">aa_dssp");
#		push(@seq,"$aa_dssp");
	    }
	}
	
	# Secondary structure prediction with psipred
	if ($v>=1) {printf ("\nPredicting secondary structure with PSIPRED ...\n");}
	&RunPsipred("$tmp.seq");
	if (open (PSIPREDFILE, "<$tmp.horiz")) {
	    my $in;             # input line
	    $ss_conf="";
	    $ss_pred="";
	    $aa_pred="";
	    # Read Psipred file
	    while ($in=<PSIPREDFILE>) {
		if    ($in=~/^Conf:\s+(\S+)/) {$ss_conf.=$1;}
		elsif ($in=~/^Pred:\s+(\S+)/) {$ss_pred.=$1;}
		elsif ($in=~/^  AA:\s+(\S+)/) {$aa_pred.=$1;}
	    }
	    close(PSIPREDFILE);
	    $ss_conf=~tr/0-9/0/c; # replace all non-numerical symbols with a 0
	    
	    push(@nam,">ss_pred");
	    push(@seq,"$ss_pred");
	    push(@nam,">ss_conf");
	    push(@seq,"$ss_conf");
#	    push(@nam,">aa_pred");
#	    push(@seq,"$aa_pred");
	}
	
	# Write psipred seqs to $tmp.a3m, setting to lower case those states that are lower case in query
	if (!open (ALIFILE, ">$tmp.a3m")) {
	    warn ("ERROR: cannot open $tmp.a3m: $!\n");
	    return 1;
	}
	for (my $i=0; $i<@nam; $i++) {print(ALIFILE "$nam[$i]\n$seq[$i]\n");}
	
	close(ALIFILE);

    } # end if (!$noss)
    
    # Append alignment sequences to dssp- and psipred sequences
    &System("cat $tmp.a3m $tmp.org.seq $tmp.fil  > $outbase.a3m");


#    &System("$perl/reformat.pl -v $v2 -r -num a3m clu  $outbase.a3m $tmp.clu",$v);
#    &System("$perl/reformat.pl -v $v2 -r -num a3m fas  $outbase.a3m $tmp.fas",$v);
#    &System("$perl/reformat.pl -v $v2 -num a2m fas  $tmp.a2m $tmp.fas",$v);

    if ($v>=3) {
	print ("Temporary files are not removed\n");
    } else {
	unlink( glob "$tmp*");
    }
    return;
} 

##############################################################################################
### Extend query sequence with longer matching sequences found in nr (called by BuildAlignment)
##############################################################################################
sub ExtendSequence {
    my $dbfile=$_[0];

    my $line;
    my $nseqs=4;  # number of full-length db sequences to be stored per query
    my $n=0;      # counts number of full-length db sequences read in
    my $id;       # key (i.e. gi number) of template (full length db sequence)
    my @ids=();   # keys (i.e.gi numbers) for best templates found in nr for query
    my $tname=""; # template name
    my $maxlen=0; # length of longest of 5 best hits
    my $next;     # returned by AlignWithTemplate(): 1 means try next alignment (current alignment unsuccessful)
    my $ndbfile;

    # Fast psiblast search for very similar sequences
    &System("$blastpgp -I T -b 10 -v 10 -e 1e-6 -A 10 -f 15 -d $dbfile -i $tmp.seq &> $tmp.bla");

    # Search names of longest, highest-scoring template sequences in blast results 
    open (BLASTFILE, "<$tmp.bla") or die "ERROR: Couldn't open $tmp.bla: $!\n";
    while ($line=<BLASTFILE>) {if ($line=~/Sequences producing significant alignments:/) {last;}}  
    while ($line=<BLASTFILE>) {
	if ($line=~/^>\s*gi\|(\d+)/ || $line=~/^>\s*(\S+)/) {
	    $id=$1;

	    do {$line=<BLASTFILE>;} while ($line && $line!~/Length =\s+(\d+)/);
	    my $len=$1;
	    do {$line=<BLASTFILE>;} while ($line && $line!~/Identities =\s+\S+\s+\((\d+)/);
	    # make sure there are more than 50% positives
	    if ($1>50) {
		# make sure the LONGEST sequence is first
		if ($len>$maxlen && $maxlen<$q_match+400) {
		    $maxlen=$len;
		    unshift(@ids,$id);
		} else {
		    push(@ids,$id);
		}
		if (++$n>=$nseqs) {last;}
	    }
	}
    }
    close(BLASTFILE);


    # Extract full-length sequences from database and align to query until good alignment found
    if ($v>=3) {printf(STDOUT "Extracting full-length sequences from database ... (%i)\n",scalar(@ids));}
    for ($n=0; $n<@ids; $n++) {
	$id = $ids[$n];  # id of n'th template 
	
	if ($v>=2) {print("$ncbidir/fastacmd -d $dbfile -s '$id'\n");} 
	if (open(FASTACMD,"$ncbidir/fastacmd -d $dbfile -s '$id' |")) { # db must have been formatted with -o T !
	    $tname=<FASTACMD>;
	    if ($tname!~/^>(\S+)/) {
		print("\nError in fastacmd! Skipping ExtendSequence\n\n");
		return $qseq;
	    }
	    $tname=$1;
	    $tseq="";
	    while ($line=<FASTACMD>) {chomp($line); $tseq.=$line;};
	    close(FASTACMD);
	} else {
	    print("WARNING: Could not find sequence gi|$id\n");
	    next;
	}

	# Try to align n'th template sequence $tseq retrieved in the db with the query
	if ($v>=3) {printf(STDOUT "\nAligning %s with sequence %-20.20s ...  ",$tmp,$tname); }
	$next = &AlignWithTemplate($qseq,$tseq);
#	print ("id=$id  next=$next\n");
	if (defined $next && $next==0) {last;}
    } #end for $n

    # Write extended query sequence into $outfile
    open (QFILE, ">$tmp.ext.seq") or die "ERROR: Couldn't open $tmp.ext: $!\n";
    if (! defined $next) {
	if ($n>0) {warn ("WARNING: could not find extension sequences!!\n\n");}
#	$qseq=~tr/*/X/d;
#	$qseq=uc($qseq);
	printf(QFILE ">%s %s\n",$name,$description);
	printf(QFILE "%s\n",uc($qseq));
    } elsif ($next==0) {
	# Leave a maximum of 100 lower-case residues at either end of sequence
	$tseq=~/([a-z]{0,100}[A-Z].*[A-Z][a-z]{0,100})/;
	$qseq=$1;
#	$qseq=~tr/uU/cC/;        # BLAST bug removed; not necessary anymore
	if ($v>=3) {
	    # Good alignment found
	    printf(STDOUT "\n");
	    printf(STDOUT " COMPL: $qseq\n");
	}
	printf(QFILE ">%s (%s) %s\n",$name,$tname,$description);
	printf(QFILE "%s\n",$qseq);
    } else {
	if ($v>=3) {printf(STDOUT "Using original sequence from $infile\n");}
#	$qseq=~tr/*/X/d;
#	$qseq=uc($qseq);
	printf(QFILE ">%s %s\n",$name,$description);
	printf(QFILE "%s\n",uc($qseq));
    }
    close(QFILE);
    if ($v>=3) {print("\n");}

    return $qseq;
}


###################################################################################################
# Align sequence $tname with residues $tseq with query with residues $qseq (called by ExtendSequence)
###################################################################################################
sub AlignWithTemplate() {
    # Default parameters
    our $d=5;    # gap opening penalty for Align.pm
    our $e=0.5;  # gap extension penatlty for Align.pm
    our $g=0;    # endgap penatlty for Align.pm
    our $matrix="Gonnet";

    # Align full-length template sequence to query
    my @subseqs;
    my $subseq;
    my $next;
    my $score;  
    my $len;
    my (@i,@j,@S);
    my ($imin,$imax,$jmin,$jmax);

    # Split query by x's (signifying inserted domains in astral sequences)
    if ($v>=3) {printf("qsseq = $qseq subseqs=%i\n",scalar(@subseqs));}
   @subseqs = split(/x/,$qseq); # 'x' marks domain or residue insertions (replaced the X symbol in scop/astral)
#   printf(STDOUT ">tseq_before\n%s\n",$tseq);
    
    # Find $subseq in $tseq by NW-alignment and transform aligned residues of $tseq into upper case 
    $tseq = lc($tseq);
    $next = 0; 
    foreach $subseq (@subseqs) {
	# Align $subseq with $tseq and return alignment in @i, @j
	$xseq=$subseq;
	$yseq=$tseq;
	$len = length($subseq);
	$score=&AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr);  
	if ($v>=3) {printf(STDOUT " score/length = %6.2f  ",$score/$len);}
	if ($score/$len<1.1) {
	    $next=1;
	    last;
	} elsif ($imin>1 && $imax<$len) {
	    $next=2;
	    last;
	} elsif (substr($tseq,$jmin-1,$jmax-$jmin+1)=~tr/A-Z/A-Z/) {
	    $next=3;
	    last;
	} else {
	    # replace the substring of the template sequence that is aligned to the query with the query subsequence
	    substr($tseq,$jmin-1,$jmax-$jmin+1) = uc($subseq);	# first=0 (instead of first=1)
	}
    } # end foreach $subseq (@subseqs)
    
    if ($v>=4 && $next==1) {
	printf(STDOUT "\nbad match for $tmp: score/length=%6.2f  length=%i\n",$score/$len,$len);
	printf(STDOUT "  XSEQ: $xseq\n");
	printf(STDOUT " MATCH: $Sstr\n");
	printf(STDOUT "  YSEQ: $yseq\n");
    } elsif ($v>=4 && $next==2) {
	printf(STDOUT "\ncould not match whole length of sequence from $tmp: imin=%-3i imax=%-3i len=%-3i\n",$imin,$imax,$len);
	printf(STDOUT "  XSEQ: $xseq\n");
	printf(STDOUT " MATCH: $Sstr\n");
	printf(STDOUT "  YSEQ: $yseq\n");
    } elsif ($v>=4 && $next==3) {
	printf(STDOUT "\nQuery sequence seqments from $tmp overlap on template sequence: imin=%-3i imax=%-3i len=%-3i\n",$imin,$imax,$len);
	printf(STDOUT "  XSEQ: $xseq\n");
	printf(STDOUT " MATCH: $Sstr\n");
	printf(STDOUT "  YSEQ: $yseq\n");
    }

    return $next;
}




###################################################################################################
# Build alignment for seqfile
###################################################################################################
sub BuildAlignment() {
    my $tmp=$_[0];
    my $seqfile=$_[0].".seq";        # query sequence (or extended query sequence) in fasta format
    my $psifile=$_[0].".psi";        # current psiblast alignment in psiblast-readable format
    my $a3mfile=$_[0].".a3m";        # final psiblast alignment in a3m format
    my $tmpfile=$_[0].".tmp";        # 
    my $blafile=$_[0].".bla";        # BLAST output file
    my $coreali=$_[0].".core.psi";   # current core alignment in psiblast-readable format 
    my $extend=$_[1];                # 1: use extended query sequence  0: use original query
    my $iter;                        # number of psiblast iterations done
    my $bopt="-bl $bl -bs $bs -bg $bg"; # score per col in bits for end pruning PSI-BLAST alignments
    my $bcore="-b 0.67";                # score per col in bits for end pruning of core alignment
#    my $bcore="-bs 0.67 -bl 0.33";  # score per col in bits for end pruning of core alignment
    my $db=$db90;                    # do first iteration against nr90
    my $db_short=$db;                # do first iteration against nr90
    $db_short=~s/.*\///;
    my $iter_remove;                 # Remove leading and trailing ends after $iter_remove iterations
    my $pmaxopt;                     # in first iteration: "-pmax $pmax", then, if -P option given: "-pmax $pmax -P $psifile"

    # Do search with extended sequence?
    if ($extend) {
	$iter_remove=8;  # trim off extensions after 3rd round
	$pmaxopt="-pmax $pmax";
     } else {
	$iter_remove=0;
	$pmaxopt="";
    }

    if ($maxiter==0) {
	system("cp $seqfile $a3mfile");
	&System("$perl/reformat.pl a3m psi $a3mfile $psifile");
	$nhits=$nseqin;
	return;
    } elsif ($maxiter==1) {

	if ($nseqin<=1) {
	    &blastpgp("-e $Eult -d $db -i $seqfile","$blafile");
	} else {
	    &blastpgp("-e $Eult -d $db -i $seqfile -B $psifile","$blafile");
	}
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt -a3m -q $seqfile $blafile $a3mfile",$v2);
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt -psi -q $seqfile $blafile $psifile",$v2);
	$iter=1;
	if ($v>=1) {
	    printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$Eult,$cov0,$db_short);
	}
	if ($v>=2) {printf ("\n");}
	return;
    }

    # Iterative psiblast search
    $iter=0; $nhits=-1;

    # First round with nr90
    $nhits_prev=$nhits; $iter++;
    if ($nseqin<=1) {
	&blastpgp("-e $Eult -d $db -i $seqfile","$blafile");
	&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bcore $pmaxopt -best -psi -q $seqfile $blafile $coreali",$v2)
    } else {
	&blastpgp("-e $Eult -d $db -i $seqfile -B $psifile","$blafile");
	&System("cp $psifile $coreali");
    }

    if ($bcore ne "") {$bcore.=" -B $coreali";} # From here on use $coreali for end pruning of $coreali
    if ($bopt ne "")  {$bopt.=" -B $coreali";}  # From here on use $coreali for end pruning of $psifile
#    &System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bopt $pmaxopt -a3m -q $seqfile $blafile test1.a3m",$v2);
#    &System("$perl/reformat.pl -num -r a3m fas test1.a3m test1.fas",$v2);
    $nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bopt $pmaxopt -best -psi -q $seqfile $blafile $psifile",$v2);
    $nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bopt $pmaxopt -best -a3m -q $seqfile $blafile $a3mfile",$v2);
    if ($pmaxopt ne "" && $pmaxopt!~/-P/) {$pmaxopt.=" -P $psifile";} # From here on use $psifile to score match states

    if ($v>=1) {
	printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$E2,$cov0,$db_short);}
    if ($v>=2) {printf ("\n");}

   # Second to penultimate round
    while ($nhits<3000 && $nhits_prev<$nhits && $iter<$maxiter-1) {
	
	if($iter==$iter_remove) { # Use trimmed sequence from now on? (Only applies if searching with extended query)
	    if ($v>=3) {print("Trimming off query extensions ...\n");}
	    &RemoveLeadingAndTrailingInsertsPsi($seqfile,$psifile);
	    &RemoveLeadingAndTrailingInsertsPsi($seqfile,$coreali);
	    &RemoveLeadingAndTrailingInsertsA3M($a3mfile);
	    &RemoveLeadingAndTrailingInsertsQuery($seqfile);
	    $nhits=0;    # Make sure that prev blast result is not used even if present result has less hits
	    $iter_remove=0;
	}
	if ($db eq $db90 && $nhits>=50) { # switch from nr90 to nr70
	    $db=$db70; $db_short=$db; $db_short=~s/.*\///; $nhits=0;
	}


	$nhits_prev=$nhits; $iter++;
	if ($extend) {&System("$perl/reformat.pl -M first -uc $psifile $psifile");} # make all residues in psifile uppercase!
	&blastpgp("-e $Eult -d $db -i $seqfile -B $psifile","$blafile");

#       The following line iterates the core alignment with the strict end pruning threshold 
#       (instead of using the results form the first round). This has been found to lead to an impairment	
#	if($iter<=2) {&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $bcore $pmaxopt -best -psi -q $seqfile $blafile $coreali",$v2);}
#	&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $sc $bopt $pmaxopt $best -a3m -q $seqfile $blafile test".$iter.".a3m",$v2);
#	&System("$perl/reformat.pl -num -r a3m fas test".$iter.".a3m test".$iter.".fas",$v2);

	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $sc $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpfile",$v2);
	if ($nhits_prev>$nhits) {$iter--; last;}
	&System("cat $tmpfile >> $a3mfile");
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $E2 $qid $sc $bopt $pmaxopt $best -psi -q $seqfile $blafile $tmpfile",$v2);
	&System("cat $tmpfile >> $psifile");

	if ($v>=1) {
	    printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$E2,$cov0,$db_short);}
	if ($v>=2) {printf ("\n");}
    } 
    
    if($iter<=$iter_remove) { # Trim sequences?
	if ($v>=3) {print("Trimming off extensions ...\n");}
	&RemoveLeadingAndTrailingInsertsPsi($seqfile,$psifile);
	&RemoveLeadingAndTrailingInsertsPsi($seqfile,$coreali);
	&RemoveLeadingAndTrailingInsertsQuery($seqfile);
	&RemoveLeadingAndTrailingInsertsA3M($a3mfile);
	$nhits=0;    # Make sure that prev blast result is not used even if present result has less hits
	$iter_remove=0;
    }
    if ($db eq $db90 && $nhits>=50) { # switch from nr90 to nr70
	$db=$db70; $db_short=$db; $db_short=~s/.*\///; $nhits=0;
    }
     
    # Last round 
    # If new sequences were found
    if ($nhits_prev<$nhits || $nhits==0) { 

	$nhits_prev=$nhits; $iter++;
	if ($extend) {&System("$perl/reformat.pl -M first -uc $psifile $psifile");} # make all residues in psifile uppercase!
	&blastpgp("-e $Eult -d $db -i $seqfile -B $psifile","$blafile");
	$nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt $best -a3m -q $seqfile $blafile $tmpfile",$v2);
	if ($v>=1) {
	    printf ("Found %s sequences in PSI-BLAST round $iter (E-value<%5.0E, coverage>%2i%%, db=%s)\n",$nhits,$Eult,$cov0,$db_short);}
	if ($v>=2) {printf ("\n");}
	if ($nhits_prev>$nhits) {
	    $iter--;
	} else {
	    &System("cat $tmpfile >> $a3mfile");
	    $nhits=&System("$perl/alignhits.pl -cov $cov0 -e $Eult $qid $sc $bopt $pmaxopt $best -psi -q $seqfile $blafile $tmpfile",$v2);
	    &System("cat $tmpfile >> $psifile");
	}
    } 

	
    return;
}



##############################################################################################
# Remove Leading and trailing lower case residues from infile
##############################################################################################
sub RemoveLeadingAndTrailingInsertsPsi() {
    my $qfile=$_[0];
    my $infile=$_[1];
    my $line;
    my $seq="";
    open (QFILE, "<$qfile") || die ("cannot open $qfile: $!");
    while ($line=<QFILE>){if($line=~/^>/) {last;} }
    while ($line=<QFILE>){
	if( $line=~/^>/) {last;}
	chomp($line);
	$line =~ tr/ //d; # delete blanks
	$seq.=$line;
    } 
    close(QFILE);
    $seq=~/^\s*([a-z.]*)([A-Z0-9-]\S*[A-Z0-9-])[a-z.]*\s*$/ || die("Error: wrong format in $infile at seq '$seq'\n");
    my $offset=length($1);
    my $length=length($2);

    open (ALIFILE, "<$infile") or die ("Error: Couldn't open $infile: $!\n");
    my @lines=<ALIFILE>;
    close (ALIFILE);
    open (ALIFILE, ">$infile") or die ("Error: Couldn't open $infile: $!\n");; 
    my $n=0;
    foreach $line (@lines) {
	if ($line=~/^\s*$/) {next;}
	chomp($line);
	$line=~/^(.*\s+)(\S+)\s*$/ || die("Error: wrong format in $infile at line '$line'\n");
	$name=$1;
	$line=substr($2,$offset,$length);
	if (($line=~tr/A-Z/A-Z/)<$rmin && $n>0) {next;}
	print (ALIFILE $name.$line."\n");
	$n++;
    }
    close (ALIFILE);
}

sub RemoveLeadingAndTrailingInsertsA3M() {
    my $infile=$_[0];
    my $line;
    my $seq="";
    my $name;
    open (ALIFILE, "<$infile") or die ("Error: Couldn't open $infile: $!\n");
    my @lines=<ALIFILE>;
    close (ALIFILE);
    open (ALIFILE, ">$infile") or die ("Error: Couldn't open $infile: $!\n");; 
    foreach $line (@lines) {
	if ($line=~/^\s*$/) {next;}
	if ($line=~/^>(.*)/) {
	    if ($seq) {
		$seq=~/^\s*[a-z.]*([A-Z0-9-]\S*[A-Z0-9-])[a-z.]*\s*$/ || die("Error: wrong format in $infile at line '$line'\n");
		print (ALIFILE ">$name\n$1\n");
	    }
	    $name=$1;
	    $seq="";
	} else {
	    chomp($line);
	    $seq.=$line;
	}
    }

    if ($seq) {
	$seq=~/[a-z.]*([A-Z0-9-]\S*[A-Z0-9-])[a-z.]*\s*$/ || die("Error: wrong format in $infile at line '$line'\n");
	print (ALIFILE ">$name\n$1\n");
    }
    close (ALIFILE);
}


sub RemoveLeadingAndTrailingInsertsQuery() {
    my $qfile=$_[0];
    my $line;
    my $name;
    my $seq="";
    open (QFILE, "<$qfile") || die ("cannot open $qfile: $!");
    while ($line=<QFILE>){if($line=~/^>/) {last;} }
    $line=~/(>.*)/;
    $name=$1;
    while ($line=<QFILE>){
	if( $line=~/^>/) {last;}
	chomp($line);
	$line =~ tr/ //d; # delete blanks
	$seq.=$line;
    } 
    close(QFILE);
    $seq=~/[a-z.]*([A-Z0-9-]\S*[A-Z0-9-])[a-z.]*\s*$/ || die("Error: wrong format in $qfile at line '$line'\n");
    open (QFILE, ">$qfile") || die ("cannot open $qfile: $!");
    print(QFILE "$name\n$1\n");
    close(QFILE);
}



##############################################################################################
# Read query sequence and write dssp state sequence into $outbase.a3m (called by BuildAlignment)
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
	    $line.=" ";
		
	    # SCOP ID? (d3lkfa_,d3grs_3,d3pmga1,g1m26.1)
	    if ($line=~/^>[defgh](\d[a-z0-9]{3})[a-z0-9_.][a-z0-9_]\s+[a-z]\.\d+\.\d+\.\d+\s+\((\S+)\)/) {
		$pdbcode=$1;
		$qrange=$2;
	    } 

	    # DALI ID? (8fabA_0,1a0i_2)
	    elsif ($line=~/^>(\d[A-Za-z0-9]{3})[A-Za-z0-9]?_\d+\s+\d+\.\d+.\d+.\d+.\d+.\d+\s+\((\S+)\)/) {
		$pdbcode=lc($1);
		$qrange=$2;
	    }

	    # PDB ID? (8fab_A, 1a0i)
	    elsif ($line=~/^>(\d[A-Za-z0-9]{3})_?(\S?)\s/) {
		$pdbcode=lc($1);
		if ($2 ne "") {$qrange="$2:";} else {$qrange="-";}
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
    if ($v>=2) {printf("\nSearching DSSP state assignments...\nname=%s  range=%s\n",$name,$qrange);}

    # Try to open dssp file 
    $dsspfile="$dsspdir/$pdbcode.dssp";
    if (! open (DSSPFILE, "<$dsspfile")) {
	printf(STDOUT "WARNING: Cannot open $dsspfile: $!\n"); 
	$pdbfile="$pdbdir/pdb$pdbcode.ent";
	if (! -e $pdbfile) {
	    printf(STDOUT "WARNING Cannot open $pdbfile: $!\n"); 
	    return 1;
	} else  {
	    if (! open (DSSPFILE, "<$dsspfile")) {
		&System("$dssp $pdbfile $tmp.dssp &> $tmp.log");
#	        &System("cp $tmp.dssp $dsspfile "); # dssp directory is on a read-only filesystem
	    }
	    $dsspfile="$tmp.dssp";
	    if (! open (DSSPFILE, "<$dsspfile")) {
		printf(STDERR "ERROR: dssp couldn't generate file from $pdbfile. Skipping $name\n");
		return 1;
	    } 
	}
    }

    #....+....1....+....2....+....3....+....4
    #  623  630 A R     <        0   0  280  etc. 
    #  624        !*             0   0    0  etc. 
    #  625    8 B A              0   0  105  etc. 
    #  626    9 B P    >>  -     0   0   71  etc. 
    #  292   28SA K  H  4 S+     0   0   71  etc.  (1qdm.dssp)
    #  293   29SA K  H  > S+     0   0   28  etc.    

    # Read in whole DSSP file
    for (my $try = 1; $try<=2; $try++) {
	$aa_dssp="";
	$ss_dssp="";
	while ($line=<DSSPFILE>) {if ($line=~/^\s*\#\s*RESIDUE\s+AA/) {last;}}
	while ($line=<DSSPFILE>) 
	{
	    if ($line=~/^.{5}(.{5})(.)(.)\s(.).\s(.)/)
	    {
		my $thisres=$1;
		my $icode=$2;
		my $chain=$3;
		my $aa=$4;
		my $ss=$5;
		my $contained=0;
		my $range=$qrange;  
		if ($aa eq "!")  {next;}    # missing residues!
		$thisres=~tr/ //d;
		$chain=~tr/ //d;
		$icode=~tr/ //d;
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
    my @aa_dssp_ali=();   # $aa_dssp_ali[$i] is dssp residue aligned to $aa_astr[$i] 
    for (my $i=0; $i<=length($aa_astr); $i++) { # sum up to len+1 
	                                        # because 0'th element in @ss_dssp and @aa_dssp is dummy "-" 
	$ss_dssp_ali[$i]="-";	
	$aa_dssp_ali[$i]="-";	
    }
    
    # To each residue (from i=0 to len-1) of input sequence $aa_astr assign aligned dssp state
    @ss_dssp = split(//,$ss_dssp);
    @aa_dssp = split(//,$aa_dssp);
    @aa_astr = split(//,$aa_astr);
    my $len = 0;
    unshift(@aa_dssp,"-"); #add a gap symbol at beginning -> first residue is at 1!
    unshift(@ss_dssp,"-"); #add a gap symbol at beginning -> first residue is at 1!
    unshift(@aa_astr,"-"); #add a gap symbol at beginning -> first residue is at 1!
    for (my $col=0; $col<@i; $col++) {
	if ($i[$col]>0) {
	    if ($j[$col]>0) {$len++;} # count match states (for score/len calculation)
	    $ss_dssp_ali[$i[$col]]=$ss_dssp[$j[$col]];
	    $aa_dssp_ali[$i[$col]]=$aa_dssp[$j[$col]];
	}
	if ($v>=4) {
	    printf ("%s %3i   %s %3i\n",$aa_astr[$i[$col]],$i[$col],$aa_dssp[$j[$col]],$j[$col]);
	}
    }
    shift (@ss_dssp_ali);   # throw out first "-" 
    shift (@aa_dssp_ali);   # throw out first "-" 
    $aa_dssp=join("",@aa_dssp_ali);
    $ss_dssp=join("",@ss_dssp_ali);

    # Debugging output
    if ($v>=3) {printf(STDOUT "DSSP: %4i  %s: length=%-3i  score/len:%-5.3f\n",$nfile,$name,$len,$score/$len);}
    if ($v>=4) {
	printf("IN:    %s\n",$xseq);
	printf("MATCH: %s\n",$Sstr);
	printf("DSSP:  %s\n",$yseq);
	printf("\n");
	printf(">ss_dssp $name\n$ss_dssp\n");
	printf(">aa_dssp $name\n$aa_dssp\n");
	printf(">aa_astra $name\n$aa_astr\n\n");
    }    
    if ($score/$len<0.5) {
	printf (STDOUT "\nWARNING: in $name ($nfile): alignment score with dssp residues too low: Score/len=%f.\n\n",$score/$len);
	printf("IN:    %s\n",$xseq);
	printf("MATCH: %s\n",$Sstr);
	printf("DSSP:  %s\n",$yseq);
	return 1;
    }

#    printf(">ss_dssp\n$ss_dssp\n");
#    printf(">aa_dssp\n$aa_dssp\n");
    return 0;
}


##############################################################################################
# Run SS prediction starting from alignment in $tmp.psi (called by BuildAlignment)
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
    if (!-e "$dummydb.phr") {print("WARNING: did not find $dummydb.phr... using $db70\n"); $dummydb=$db70; } 

    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
    &System("$blastpgp -b 0 -j 1 -h 0.001 -d $dummydb -i $infile -B $basename.psi -C $basename.chk &> $basename.bla");
    
#    print("Predicting secondary structure...\n");
    
    system("echo $rootname.chk > $basename.pn\n");
    system("echo $rootname.seq > $basename.sn\n");
    system("$ncbidir/makemat -P $basename");
    
    &System("$execdir/psipred $basename.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $basename.ss");

    &System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $basename.ss2 $basename.ss > $basename.horiz");
    
    # Remove temporary files
    system ("rm -f $basename.pn $basename.sn $basename.mn $basename.bla $basename.mtx $basename.aux $basename.ss $basename.ss2");
     system ("rm -f $basename.chk");
    return;
}

################################################################################################
### Run blastpgp
################################################################################################
sub blastpgp() {
    if (&System("$blastpgp -b 20000 -v 1 $_[0] &> $_[1]")) {
	print("\nError: PSI-BLAST returned with error:\n");
	&System("cat $_[1]"); die();
    }
    return;
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
    $_[0]=~s/^perl\s+//; # remove "perl "
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



