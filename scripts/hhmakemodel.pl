#! /usr/bin/env perl
#
# hhmakemodel.pl 
# Generate a model from an output alignment of hhsearch. 
# Usage: hhmakemodel.pl -i file.out (-ts file.pdb|-al file.al) [-m int|-m name|-m auto] [-pdb pdbdir] 

#     HHsuite version 2.0.3 (January 2012)
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
use strict;
use Align;


$|=1;  # force flush after each print

# Default parameters
our $d=7;     # gap opening penalty for Align.pm; more than 2 mismatches - 2 matches    ## previously: 1
our $e=0.01;  # gap extension penatlty for Align.pm; allow to leave large gaps bridging uncrystallized regions  ## previously: 0.1
our $g=0.1;   # endgap penalty for Align.pm; allow to shift SEQRES residues for uncrystallized aas to ends of alignment  ## previously: 0.9
my $v=2;      # 3: DEBUG

my $formatting="CASP";     # CASP or LIVEBENCH
my $servername="HHpred2";  # HHpred2 or HHpred3
my $MINRES=30;     # minumum number of new query residues required for a hit to be used as additional parent
my $infile="";
my $outfile="";
my $outformat="fas";
my $pickhits="1 "; # default: build one model from best hit 
my $Pthr=0;
my $Ethr=0;
my $Prob=0;
my $shift=0;       # ATTENTION: set to 0 as default!  
my $NLEN=14;       # length of the name field in alignments of hhsearch-output 
my $NUMRES=100;    # number of residues per line in FASTA, A2M, PIR format
my $program="hhmakemodel.pl";
my $usage="
hhmakemodel.pl from HHsuite $VERSION  
From the top hits in an hhsearch output file (hhr), you can  
  * generate a MSA (multiple sequence alignment) containing all representative 
    template sequences from all selected alignments (options -fas, -a2m, -a3m, -pir) 
  * generate several concatenated pairwise alignments in AL format (option -al)
  * generate several concatenated coarse 3D models in PDB format (option -ts) 
In PIR, PDB and AL format, the pdb files are required in order to read the pdb residue numbers 
and ATOM records.
The PIR formatted file can be used directly as input to the MODELLER homology modelling package.
Usage: $program [-i] file.hhr [options]

 Options:
 -i   <file.hhr>        results file from hhsearch with hit list and alignments
 -fas <file.fas>        write a FASTA-formatted multiple alignment to file.fas
 -a2m <file.a2m>        write an A2M-formatted multiple alignment to file.a2m
 -a3m <file.a3m>        write an A3M-formatted multiple alignment to file.a3m
 -m   <int> [<int> ...] pick hits with specified indices  (default='-m 1')
 -p   <probability>     minimum probability threshold (default=$Pthr)        
 -e   <E-value>         maximum E-value threshold (default=$Ethr)        
 -q   <query_ali>       use the full-length query sequence in the alignment 
                        (not only the aligned part);
                        the query alignment file must be in HHM, FASTA, A2M,  
                        or A3M format.
 -N                     use query name from hhr filename (default: use same  
                        name as in hhr file)
 -first                 include only first Q or T sequence of each hit in MSA
 -v                     verbose mode (default=$v)

Options when database matches in hhr file are PDB or SCOP sequences
 -pir <file.pir>        write a PIR-formatted multiple alignment to file.pir 
 -ts  <file.pdb>        write the PDB-formatted models based on *pairwise*  
                        alignments into file.pdb
 -al  <file.al>         write the AL-formatted *pairwise* alignments into file.al
 -d   <pdbdirs>         directories containing the pdb files (for PDB, SCOP, or DALI  
                        sequences) (default=$pdbdir)
 -s   <int>             shift the residue indices up/down by an integer (default=$shift);           
 -CASP                  formatting for CASP (for -ts, -al options) (default: LIVEBENCH  
                        formatting)
 Options when query is compared to itself (for repeat detection) 
 -conj                  include also conjugate alignments in MSA (with query and  
                        template exchanged) 
 -conjs                 include conjugate alignments and sort by ascending diagonal  
                        value (i.e. i0-j0)
\n"; 

# Options to help extract repeats from self-alignments:
# 1. default   2. -conj    3. -conj_diag   4. -conj_compact
#     ABCD         ABCD         ---A            ABCD      
#     BCD-         BCD-         --AB            BCDA
#     D---         CD--         -ABC            CDAB
#     CD--         D---         ABCD            DABC
#                  ---A         BCD-
#                  --AB         CD--
#                  -ABC         D---


# Variable declarations
my $line;        # input line
my $score=-1;    # score of the best model; at the moment: Probability  
my $qname="";    # name of query from hhsearch output file (infile)
my $tname="";    # name of template (hit) from hhsearch output file (infile)
my $qnameline=""; # nameline of query
my $tnameline;   # nameline of template
my $pdbfile;     # name of pdbfile to read
my $pdbcode;     # four-letter pdb code in lower case and _A if chain A (e.g. 1h7w_A)
my $aaq;         # query amino acids from hhsearch output
my @aaq;         # query amino acids from hhsearch output
my @qname;       # query names in present alignment as returned from ReadAlignment()
my @qfirst;      # indices of first residues in present alignmet as returned from ReadAlignment()
my @qlast;       # indices of last residues in present alignmet as returned from ReadAlignment()
my @qseq;        # sequences of querys in present alignment as returned from ReadAlignment()
my @tname;       # template names in present alignment as returned from ReadAlignment()
my @tfirst;      # indices of first residues in present alignmet as returned from ReadAlignment()
my @tlast;       # indices of last residues in present alignmet as returned from ReadAlignment()
my @tseq;        # sequences of templates in present alignment as returned from ReadAlignment()
my $aat;         # template amino acids from hhsearch output
my @aat;         # template amino acids from hhsearch output
my $aapdb;       # template amino acids from pdb file
my @aapdb;       # template amino acids from pdb file
my $qfirst=0;    # first residue of query
my $qlast=0;     # last residue of query
my $qlength;     # length of query sequence
my $tfirst=0;    # first residue of template
my $tlast=0;     # first residue of template
my $tlength;     # length of template sequence
my $l=1;         # counts template residues from pdb file (first=1, like for i[col2] and j[col2]
my $col1=0;      # counts columns from hhsearch alignment
my $col2=0;      # counts columns from alignment (by function &AlignNW) of $aat versus $aapdb 
my @i1;          # $i1[$col1] = index of query  residue in column $col1 of hhsearch-alignment
my @j1;          # $j1[$col1] = index of template residue in column $col1 of hhsearch-alignment
my @j2;          # $j2[$col2] = index of hhsearch template seq in $col2 of alignment against pdb template sequence
my @l2;          # $l2[$col2] = index of pdb template seq in $col2 of alignment against hhsearch template sequence
my @l1;          # $l1[$col1] = $l2[$col2]
my $res;         # residue name
my $chain;       # pdb chain from template name
my $qfile;       # name of query sequence file (for -q option)
my $qmatch;      # number of match states in alignment
my $hit;         # index of hit in hit list
my $k;           # index of hit sorted by position in alignment with query (k=1,...,k=@first-2)
my %picked=();   # $picked{$hit} is defined and =$k for hits that will be transformed into model 
my @remarks;
my @printblock;  # block 0: header  block k: k'th hit
my $keyword="";  # either METHOD for CASP format or REMARK for LIVEBENCH format
my $conj=0;      # include conjugate sequences? Sort in which way?
my $conjugate=0; # when query is compared to itself: do not include conjugate alignments
my $onlyfirst=0; # include only first representative sequence of each Q/T alignment
my $dummy;       # dummy
my $addchain=1;  # 1: PDB files contain chain-id as in 1abc_A.pdb (not 1abc.pdb or pdb1abc.pdb etc.) 
my $pdbdirs=$pdbdir; # default pdb directory with *.pdb files
my $options="";

# Processing command line options
if (@ARGV<1) {die $usage;}
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}


# Set options
if ($options=~s/ -i\s+(\S+) / /g)    {$infile=$1;}
if ($options=~s/ -q\s+(\S+) / /g)    {$qfile=$1;}
if ($options=~s/ -ts\s+(\S+) / /ig)  {$outfile=$1; $outformat="TS";}
if ($options=~s/ -pdb\s+(\S+) / /ig) {$outfile=$1; $outformat="TS";}
if ($options=~s/ -al\s+(\S+) / /ig)  {$outfile=$1; $outformat="AL";}
if ($options=~s/ -pir\s+(\S+) / /ig) {$outfile=$1; $outformat="PIR";}
if ($options=~s/ -fas\s+(\S+) / /ig) {$outfile=$1; $outformat="FASTA";}
if ($options=~s/ -a2m\s+(\S+) / /ig) {$outfile=$1; $outformat="A2M";}
if ($options=~s/ -a3m\s+(\S+) / /ig) {$outfile=$1; $outformat="A3M";}
if ($options=~s/ -p\s+(\S+) / /g)    {$Pthr=$1;}
if ($options=~s/ -e\s+(\S+) / /g)    {$Ethr=$1;}
if ($options=~s/ -s\s+(\S+) / /g)    {$shift=$1;}
if ($options=~s/ -d\s+(([^-\s]\S*\s+)*)/ /g) {$pdbdirs=$1;}
if ($options=~s/ -m\s+((\d+\s+)+)/ /g) {$pickhits=$1; }
if ($options=~s/ -first\s+/ /ig)     {$onlyfirst=1;}

# Self-alignment options
if ($options=~s/ -conj\s+/ /ig) {$conj=1;}
if ($options=~s/ -conjs\s+/ /ig) {$conj=2;}

# Switch formatting and method description
if ($options=~s/ -CASP\s+/ /ig) {$formatting="CASP";}
if ($options=~s/ -LIVEBENCH\s+/ /ig) {$formatting="LIVEBENCH";}
if ($options=~s/ -server\s+(\S+)/ /g) {$servername=$1;}

# Set verbose mode?
if ($options=~s/ -v\s+(\d+) / /g)  {$v=$1;}
elsif ($options=~s/ -v\s+/ /g)    {$v=1;}

# Read infile and outfile 
if (!$infile  && $options=~s/^\s*([^-]\S+)\s*/ /) {$infile=$1;} 
if (!$outfile && $options=~s/^\s*([^-]\S+)\s*/ /) {$outfile=$1;} 
if ($options=~s/ -N / /ig)  {
    $qname=$infile; 
    $qname=~s/^.*?([^\/]+)$/$1/;    # remove path
    $qname=~s/^(.*)\.[^\.]*$/$1/;   # remove extension
    $qnameline=$qname;
}

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if ($infile eq "")  {die("$usage\nError in $program: input file missing: $!\n");}
if ($outfile eq "") {die("$usage\nError in $program: output file missing: $!\n");}

my @pdbdirs = split(/\s+/,$pdbdirs);

# Find query name in input file
open (INFILE, "<$infile") || die "Error in $program: Couldn't open $infile: $!\n";
while ($line=<INFILE>) {
    if ($v>=3) {print("$line");}
    if ($qname eq "" && $line=~/^Query:?\s*(\S+)(.*)/)     {$qname=$1; $qnameline=$1.$2;}
    if ($line=~/^Match_columns:?\s*(\S+)/) {$qmatch=$1; last;}
}
if (!($line=<INFILE>)) {die ("Error in $program: wrong format in $infile: $!\n");}


# Prepare hash %pick with indices of hits that will be transformed into model
# No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
#  1 153l   Lysozyme (E.C.3.2.1.17) 100.0       0       0  381.0  19.4  185    1-185     1-185 (185)
#  2 1qsa_A Soluble lytic transglyc 100.0 2.1E-39 2.5E-43  225.8   8.3  149   21-182   423-600 (618)
#  3 1ltm   36 kDa soluble lytic tr  95.9 3.5E-06 4.1E-10   50.3  11.0   95   28-122    76-221 (320)
  
# option '-m m1 m2 m3': pick models manually
my @pickhits = split(/\s+/,$pickhits);
$k=1;
foreach $hit (@pickhits) {
    if (!defined $picked{$hit}) {$picked{$hit}=$k;}
    $k++;
}

if ($outformat eq "AL" || $outformat eq "TS") {
    &MakePairwiseAlignments();
} else {
    &MakeMultipleAlignment();
}
exit;
    

##################################################################################
# Construct AL or TS formatted alignment as a list of pairwise alignments
##################################################################################
sub MakePairwiseAlignments()
{   
    # Scan through query-vs-template-alignments from infile and create first (combination) model
    $hit=0;   # counts hits in hit list
    my $models=0;
    while ($line=<INFILE>) {
	if ($line=~/^>(\S+)/) {
	    $hit++; 
	    if ($Pthr || $Ethr || defined $picked{$hit}) {
		# Found right alignment (hit)
		if (defined $picked{$hit}) {$k=$picked{$hit};} else {$k=$hit;}

		if ($line=~/^>(.*?)\s+E=.*$/) {
		    $line=$1;         # remove E=1.3E-30 etc. at the end
		} else {
		    $line=~/^>(.*)/;
		    $line=$1;
		}
		my $nameline=$line;
		my $evalue;
		$line=<INFILE>;
		if ($line=~/Probab\s*=\s*(\S+).*E-value\s*=\s*(\S+)/) {$score=$1; $evalue=$2} 
		else {$score=0; warn("WARNING: could not print score $line");}
		if ($line=~/Aligned_columns=\s*(\S+)/) {;} else {warn("WARNING: could not find COLS\n");}
		if ($Pthr && $score<$Pthr) {last;}  # Probability too low -> finished
		if ($Ethr && $evalue>$Ethr) {last;} # Evalue too high > finished
		
		
		# Commented out in CASP format
		if ($formatting eq "LIVEBENCH") {
		    $printblock[$k] ="PFRMAT $outformat\n";
		    $printblock[$k].="TARGET $qname\n";
		}
		
		$remarks[$k]="REMARK $k: $nameline\n";
		$remarks[$k].="REMARK    $line";

		&ReadAlignment();
		$qfirst = $qfirst[0];
		$qlast  = $qlast[0];
		$aaq    = $qseq[0];
		$tfirst = $tfirst[0];
		$aat    = $tseq[0];
		$tname  = $tname[0];

		if ($v>=3) {
		    for (my $i=0; $i<@qfirst; $i++) {
			printf("Q %-14.14s %s\n",$qname[$i],$qseq[$i]);
		    }
		    printf("\n");
		    for (my $i=0; $i<@tfirst; $i++) {
			printf("T %-14.14s %s\n",$tname[$i],$tseq[$i]);
		    }
		    printf("\n");
		}	    
  
		# Extract pdbcode and construct name of pdbfile and return in global variables $pdbid and $chain
		if (&ExtractPdbcodeAndChain($tname[0])) {next;}
		if ($chain eq "[A ]") {$pdbcode.="_A";} elsif ($chain eq ".") {;} else {$pdbcode.="_$chain";}
	
		# Read score (=probability)
		$printblock[$k].="REMARK $nameline\n";
		$printblock[$k].="REMARK $line";
		$printblock[$k].="SCORE  $score\n";
		$printblock[$k].="PARENT $pdbcode\n";
		$printblock[$k].="MODEL  $k\n";

		&WritePairwiseAlignments();
		$printblock[$k].="END\n";
		$models++;
 
	    }
	}
    }
    $k=$#printblock;  # set $k to last index in @printblock
    if ($k<0) {
	$printblock[1]="PARENT NONE\nTER\n";
	$printblock[1].="END\n";
	if ($v>=1) {print("WARNING: no hits found for model!\n");}
    }
    close (INFILE);
    
    if ($v>=2) {
	printf("$models models built\n");
    }
    
    
    
    # Write model file header
    
    #---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
    
    # Print header
    my $date = scalar(localtime);
    if ($formatting eq "CASP") {    
	$printblock[0]="PFRMAT $outformat\n";
	$printblock[0].="TARGET $qname\n";
    }
    $printblock[0].="REMARK AUTHOR $servername\n";
    $printblock[0].="REMARK $date\n";
    $printblock[0].="REMARK J. Soeding (johannes.soeding\@tuebingen.mpg.de)\n";
#    $printblock[0].="REMARK J. Soeding \n";
    
# Add method description
#----+----1----+----2----+----3----+----4----+----5----|----6----+----7----+----8
    if ($formatting eq "CASP") {$keyword="METHOD";} else {$keyword="REMARK";}

    if ($servername eq "HHpred2") {
	$printblock[0].=
"$keyword Fold recognition by HMM-HMM comparison and secondary structure scoring
$keyword A multiple alignment is built from the template sequence with PSI-BLAST
$keyword (up to 8 rounds with E-value threshold 1E-3). PSIPRED is used for 
$keyword secondary structure prediction. The alignment is converted to an HMM and 
$keyword compared with a database of HMMs derived from representative sequences
$keyword in the PDB and SCOP (70% maximum sequence identity, updated weekly).   
$keyword Their secondary structure states are determined by DSSP. A fast and 
$keyword sensitive column score is used, as well as a secondary structure scoring 
$keyword scheme which is analogous to amino acid scoring with substitution
$keyword matrices.
$keyword Soding, J. (2005) Protein homology detection by HMM-HMM comparison.
$keyword Bioinformatics 21, 951-960.
REMARK Templates used:
";

} elsif ($servername eq "HHpred3") {
    $printblock[0].=
"$keyword Fold recognition by HMM-HMM comparison using intermediate HMM search
$keyword A multiple alignment is built from the template sequence with PSI-BLAST
$keyword (up to 8 rounds with E-value threshold 1E-3). PSIPRED is used for 
$keyword secondary structure prediction. The alignment is converted to an HMM and
$keyword compared with a database of HMMs derived from representative sequences
$keyword in the PDB and SCOP (70% maximum sequence identity, updated weekly).  
$keyword Their secondary structure states are determined by DSSP. A fast and 
$keyword sensitive column score is used, as well as a secondary structure scoring
$keyword scheme which is analogous to amino acid scoring with substitution 
$keyword matrices.
$keyword When no significant match can be found, our intermediate HMM search
$keyword method HHsenser is used to add neighboring, homologous alignments to the 
$keyword original query alignment by HMM-HMM comparison. The resulting super-
$keyword alignment is converted to an HMM and PDB and SCOP are searched again.
$keyword Soding, J. (2005) Protein homology detection by HMM-HMM comparison.
$keyword Bioinformatics 21, 951-960.
REMARK Templates used:
";
}  
 
    # Add remarks
    for ($k=0; $k<@remarks; $k++) {
	if (defined $remarks[$k]) {
	    $printblock[0].=$remarks[$k];
	}
    }
    $printblock[0].="REMARK \n";

   
    # Print @printblock into outfile
    
    open (OUTFILE, ">$outfile") || die "Error in $program: Couldn't open $outfile: $!\n";
    foreach my $printstr (@printblock) {
	my @printarr=split(/\n/,$printstr);
	if ($outformat eq "TS") {
	    foreach $printstr (@printarr) {
		printf(OUTFILE "%-80.80s\n",$printstr);
	    }
	} else {
	    foreach $printstr (@printarr) {
		printf(OUTFILE "%s\n",$printstr);
	    }
	}
    }
    close (OUTFILE);
        
    if ($outformat eq "TS") {
	# Call MaxSprout to generate sidechains
    }
    return;    
}


##################################################################################
# Construct multiple alignment in FASTA, A2M, or PIR format
##################################################################################
sub MakeMultipleAlignment()
{
    my @hitnames=();  # $hitnames[$k] is the nameline of the ihit'th hit
    my @hitseqs=();   # $hitseqs[$k] contains the residues of the ihit'th hit
    my @hitdiag=();   # $hitdiag[$k] = $qfirst[0]-$tfirst[0] 
    my @conjnames=(); # $hitnames[$k] is the nameline of the ihit'th conjugate hit
    my @conjseqs=();  # $hitseqs[$k] contains the residues of the ihit'th conjugate hit
    my @conjdiag=();  # $hitdiag[$k] = $qfirst[0]-$tfirst[0] for conjugate alignments
    
    my $new_hit;           # residues of new hit
    my $i;                 # residue index
    my $j;                 # residue index
    my $k;                 # sequence index

    $hitnames[0]="";
    $hitseqs[0]="";
    $hitdiag[0]=0;
    $conjnames[0]="";
    $conjseqs[0]="";
    $conjdiag[0]=0;

    open (INFILE, "<$infile") || die "Error in $program: Couldn't open $infile: $!\n";
    $hit=0;  # counts hits in hit list
    
    # Read one alignment after the other
    while ($line=<INFILE>) {
	
	# Found new aligment
	if ($line=~/^>(\S+)/) {
	    $hit++;
	    
	    # Is alignment selected by user?
	    if ($Pthr || $Ethr || defined $picked{$hit}) {

		if ($line=~/^>(\S+)(.*)/) {$tname=$1; $tnameline=$1.$2;}
		else {die("\nError: bad format in $infile, line $.: code 1\n");}
		
		$line = <INFILE>;
		if ($line=~/Probab\s*=\s*(\S+).*E-value\s*=\s*(\S+)/) {
		    if ($Pthr && $1<$Pthr) {last;}  # Probability too low -> finished
		    if ($Ethr && $2>$Ethr) {last;} # Evalue too high > finished
		} else { die("\nError: bad format in $infile, line $.: code 2\n"); }

		# Read next alignment with $aaq, $qfirst, @tseq, @first, and @tname
		&ReadAlignment();
		chomp($tnameline);
		if ($tnameline=~/\S+\s+(.*)/) {$tname[0].=" $1";}  # template seed sequence gets its description

		# Format sequences into @hitseqs and @hitnames
		&FormatSequences(\@hitnames,\@hitseqs,\@hitdiag,\@qname,\@qseq,\@qfirst,\@qlast,\$qlength,\@tname,\@tseq,\@tfirst,\@tlast,\$tlength);

		# Use conjugate alignments?
		if ($conj>0) {
		    &FormatSequences(\@conjnames,\@conjseqs,\@conjdiag,\@tname,\@tseq,\@tfirst,\@tlast,\$tlength,\@qname,\@qseq,\@qfirst,\@qlast,\$qlength);
		}

	    } # end: if ($Pthr>0 || defined $picked{$hit})

	} # end: if ($line=~/^>(\S+)/)  # found new alignment

    } # end while
    close (INFILE);

    # Insert full-length query sequence?
    if ($qfile) {
	$hitseqs[0]="";	
	open (QFILE, "<$qfile") || die "Error in $program: Couldn't open $qfile: $!\n";
	while ($line=<QFILE>) {
	    if ($line=~/^>/ && $line!~/^>ss_/ && $line!~/^>sa_/ && $line!~/^>aa_/ && $line!~/^>Consensus/) {last;}
	}
	while ($line=<QFILE>) {
	    if ($line=~/^>/ || $line=~/^\#/) {last;}
	    $line=~tr/\n\.-//d; 
	    $line=~tr/a-z/A-Z/;
	    $hitseqs[0].=$line;
	}
	close(QFILE);
	if ($v>=2) {printf("\nQ(full) %-14.14s %s\n",$qname,$hitseqs[0]);}
    }
    
    
    # DEBUG
    if ($v>=3) {
	printf("\nQuery    %-14.14s %s\n",$qname,$hitseqs[0]);
	for ($k=1; $k<@hitnames; $k++) {
	    printf("T hit %3i  %-14.14s %s\n",$k,$hitnames[$k],$hitseqs[$k]);
	}
	printf("\n");
	printf("\nQuery    %-14.14s %s\n",$qname,$conjseqs[0]);
	for ($k=1; $k<@conjnames; $k++) {
	    printf("T conj %3i %-14.14s %s\n",$k,$conjnames[$k],$conjseqs[$k]);
	}
	printf("\n");
    }


    # Include conjugate sequences?
    if ($conj>0) {
	shift(@conjseqs);   # delete zeroth ("query") sequence of @conjseqs
	shift(@conjnames);  # 
	shift(@conjdiag);   #
	
	# Sort by diagonals $hitdiag[], $conjdiag[]
	&Sort(\@hitdiag,\@hitseqs,\@hitnames);
	&Sort(\@conjdiag,\@conjseqs,\@conjnames);

	# Append conjugate sequences to hitseqs
	splice(@hitseqs,scalar(@hitseqs),0,@conjseqs);    
	splice(@hitnames,scalar(@hitnames),0,@conjnames); 
	
	if ($v>=3) {
	    printf("\nQuery    %-14.14s %s\n",$qname,$hitseqs[0]);
	    for ($k=1; $k<@hitnames; $k++) {
		chomp($hitnames[$k]);
		printf("T tot %3i  %-14.14s %s\n",$k,$hitnames[$k],$hitseqs[$k]);
		$hitnames[$k].="\n";
	    }
	}
  }

    # Insert gaps:
    
    my @len_ins; # $len_ins[$j] will count the maximum number of inserted residues after match state $j.
    my @inserts; # $inserts[$j] contains the insert (in small case) of sequence $k after the $j'th match state
    my $insert;
    my $ngap;
    
    # For each match state determine length of LONGEST insert after this match state and store in @len_ins
    for ($k=0; $k<@hitnames; $k++) {
	# split into list of single match states and variable-length inserts
	# ([A-Z]|-) is the split pattern. The parenthesis indicate that split patterns are to be included as list elements
	# The '#' symbol is prepended to get rid of a perl bug in split
	$j=0;
	@inserts = split(/([A-Z]|-)/,"#".$hitseqs[$k]."#");
#	    printf("Sequence $k: $hitseqs[$k]\n");
#	    printf("Sequence $k: @inserts\n");
	foreach $insert (@inserts) {
	    if( !defined $len_ins[$j] || length($insert)>$len_ins[$j]) {
		$len_ins[$j]=length($insert);
	    }
	    $j++;
#	    printf("$insert|");
	}
#	printf("\n");
    }
    
    
    
    # After each match state insert residues and fill up with gaps to $len_ins[$i] characters
    for ($k=0; $k<@hitnames; $k++) {
	# split into list of single match states and variable-length inserts
	@inserts = split(/([A-Z]|-)/,"#".$hitseqs[$k]."#");
	$j=0;
	
	# append the missing number of gaps after each match state
	foreach $insert (@inserts) {
	    if($outformat eq "FASTA") {
		for ($i=length($insert); $i<$len_ins[$j]; $i++) {$insert.="-";}
	    }
	    else {
		for ($i=length($insert); $i<$len_ins[$j]; $i++) {$insert.=".";}
	    }
	    $j++;
	}
	$hitseqs[$k] = join("",@inserts);
	$hitseqs[$k] =~ tr/\#//d; # remove the '#' symbols inserted at the beginning and end
    }
    
    # Remove columns at beginning and end with gaps in all sequences
    my $remove_start;
    my $remove_end;
    my $len;
    $hitseqs[0]=~/^(-*)/;
    $remove_start=length($1);
    $hitseqs[0]=~/(-*)$/;
    $remove_end=length($1);
    for ($k=0; $k<@hitnames; $k++) {
	$hitseqs[$k]=~s/^.{$remove_start}(.*).{$remove_end}/$1/;
    }
    $len=($hitseqs[0]=~tr/a-zA-Z/a-zA-Z/);
    
    # Prepare name line of query
    if ($outformat eq "PIR") {
	my $qnametmp=$qname;
	$qnametmp=~tr/:/;/;
	$qnameline=~/^\S+\s*(.*)/;
	my $qnamelinetmp=$1;
	$qnamelinetmp=~tr/:/;/;
	$hitnames[0] = sprintf(">P1;%s\nsequence:%s:%4i: :%4i: :%s: : 0.00: 0.00\n",$qnametmp,$qnametmp,$remove_start+1,$len+$remove_start,$qnamelinetmp);
    } else {
	# outformat is "FASTA" or "A2M" or "A3M" or ...
	$hitnames[0] = ">$qnameline\n";
    } 
    
    # If pretty diagonally sorted order is wanted...
    if ($conj>0) {
	if ($conj==2) {
	    my $center = 0.5*(scalar(@hitseqs)-1);
	    @conjseqs = splice(@hitseqs,$center+1,$center);
	    splice(@hitseqs,0,0,@conjseqs);
	    @hitseqs = reverse(@hitseqs);
	    
	    @conjnames = splice(@hitnames,$center+1,$center);
	    splice(@hitnames,0,0,@conjnames);
	    @hitnames = reverse(@hitnames);
	    # Shorten namelines of all but first sequence
	    my %count;
	    for ($k=0; $k<@hitnames; $k++) {
		if ($k==$center) {$k++;}
		$hitnames[$k]=~/(\S{1,14})/;
		if (!defined $count{$1}) {$count{$1}=0;}
		my $count = ++$count{$1};
#		printf("vorher: %s   ",$hitnames[$k]);
		$hitnames[$k]=~s/^(\S{1,14}).*/$1:$count/;
#		printf("nachher: %s\n",$hitnames[$k]);
	    }
	} else {
	    for ($k=0; $k<@hitnames; $k++) {$hitnames[$k]=">$qname\n";}
	}
    }
    
    # Remove gaps? Captialize?
    if ($outformat eq "PIR") {
	for ($k=0; $k<@hitnames; $k++) {
	    $hitseqs[$k].="*";;                    # Transform to upper case
	    $hitseqs[$k]=~tr/a-z./A-Z-/;           # Transform to upper case
	    $hitseqs[$k]=~s/(.{1,$NUMRES})/$1\n/g; # insert newlines every NUMRES positions
	}
    } elsif ($outformat eq "FASTA") {
	for ($k=0; $k<@hitnames; $k++) {
	    $hitseqs[$k]=~tr/a-z./A-Z-/;           # Transform to upper case
	    $hitseqs[$k]=~s/(.{1,$NUMRES})/$1\n/g; # insert newlines every NUMRES positions
	}
    } elsif ($outformat eq "A2M") {
	for ($k=0; $k<@hitnames; $k++) {$hitseqs[$k]=~s/(.{1,$NUMRES})/$1\n/g;}    # insert newlines every NUMRES positions
    } elsif ($outformat eq "A3M") {
	for ($k=0; $k<@hitnames; $k++) {$hitseqs[$k]=~tr/.//d;$hitseqs[$k].="\n";} # Remove gaps aligned to inserts
    }    
    
    # Write sequences into output file
    open (OUTFILE, ">$outfile") || die ("cannot open $outfile:$!");
    for ($k=0; $k<@hitnames; $k++) {
	print(OUTFILE "$hitnames[$k]$hitseqs[$k]");
    }
    close OUTFILE;
    
    
    if ($v>=2) {
	printf("%i sequences written to $outfile\n",scalar(@hitnames));
    }
}


# Format sequences into @hitseqs and @hitnames
# & Call with FormatSequences(\@hitnames,\@hitseqs,\@qname,\@qseq,\@qfirst,\@qlast,\$qlength,\@tname,\@tseq,\@tfirst,\@tlast,\$tlength);
sub FormatSequences()
{
    my $p_hitnames = $_[0];     # typeglob to $hitname
    my $p_hitseqs  = $_[1];     # ...
    my $p_hitdiag  = $_[2];     # ...

    my $p_qname   = $_[3];      # 
    my $p_qseq    = $_[4];      # 
    my $p_qfirst  = $_[5];      # 
    my $p_qlast   = $_[6];      # 
    my $p_qlength = $_[7];      # 

    my $p_tname   = $_[8];      # 
    my $p_tseq    = $_[9];      # 
    my $p_tfirst  = $_[10];     # 
    my $p_tlast   = $_[11];     # 
    my $p_tlength = $_[12];     # 
    my $i;
    
    if ($v>=2) {
	if (defined $picked{$hit}) {
	    print("hit=$hit  picked=$picked{$hit} tname=$tname[0]");
	} else {
	    print("hit=$hit  picked=evalue<$Ethr tname=$tname[0]");
	}
	for (my $i=1; $i<@{$p_tname}; $i++) {
	    print(", $tname[$i]");
	}
	print("\n");
    }
    
    
    my $qfirst = ${$p_qfirst}[0];
    my $qlast  = ${$p_qlast}[0];
    my $qlength = ${$p_qlength};
    my $aaq = ${$p_qseq}[0];

    @aaq = unpack("C*",$aaq); # needed for transforming template sequences into a3m based on query residues (NOT HMM match states!)
    $aaq=~tr/.-//d;           # remove all gaps from query sequence
    
    # For all template sequences in the present alignment
    for (my $k=0; $k<@{$p_tname}; $k++) {
	
	$tname =${$p_tname}[$k];
	$tfirst=${$p_tfirst}[$k];
	$aat   =${$p_tseq}[$k];
	
	# Transform template residues into a3m format: 
	# match states are those where query has residue (NOT where HMM has match state!)
	# This makes sense since we want to build a model for the query sequence.
	@aat = unpack("C*",$aat);
	$aat="";
	
	# Transform all columns with residue in query into match/delete states, all others to inserts
	for ($i=0; $i<scalar(@aaq); $i++) {
	    if ($aaq[$i]!=45 && $aaq[$i]!=46) { # no gap in query
		if($aat[$i]==46) {
		    $aat.="-";                  # transform '.' to '-' if aligned with a query residue
		} else {
		    $aat .= uc(chr($aat[$i]));  # UPPER case if aligned with a query residue (match state)
		}
	    } else {
		if($aat[$i]!=45 && $aat[$i]!=46) { # no gap in template?
		    $aat.=lc(chr($aat[$i]));       # lower case if aligned with a gap in the query (insert state) 
		}
	    }
	}
	
	if ($v>=2) {	
	    printf("\nQ %-14.14s %s\n",$qname,$aaq);
	    printf("T %-14.14s %s\n",$tname,$aat);
	}
	
	# Outformat is PIR? => read residues and indices from PDB ATOM records
	if ($outformat eq "PIR") {
	    
	    # Extract pdbcode and construct name of pdbfile and return in global variables $pdbid and $chain
	    if (&ExtractPdbcodeAndChain($tname)) {next;}
	    
	    # Read sequence from pdb file
	    if (!open (PDBFILE, "<$pdbfile")) {
		die ("Error in $program: Couldn't open $pdbfile: $!\n");
	    }
	    $aapdb="";
	    $l=0;
	    my @nres;        # $nres[$l] = pdb residue index for residue $aapdb[$l]
	    my $nres=-1e6;
	    my $resolution=-1.00;
	    my $rvalue=-1.00;
	    while ($line=<PDBFILE>) {
		if ($line=~/^REMARK.*RESOLUTION\.\s+(\d+\.?\d*)/) {$resolution=$1;}
		if ($line=~/^REMARK.*R VALUE\s+\(WORKING SET\)\s+:\s+(\d+\.?\d*)/) {$rvalue=$1;}
		if ($line=~/^ENDMDL/) {last;} # if file contains NMR models read only first one
		if (($line=~/^ATOM\s+\d+  .. [ A](\w{3}) $chain\s*(-?\d+.)/ ||
		    ($line=~/^HETATM\s+\d+  .. [ A](\w{3}) $chain\s*(-?\d+.)/ && &Three2OneLetter($1) ne "X") ) && 
		     $2 ne $nres ) { 
		    $res=$1;
		    $nres=$2;
		    $nres[$l]=$2;
		    $res=&Three2OneLetter($res);
		    $aapdb[$l++]=$res;
		    $aapdb.=$res;
		}
	    }
	    close (PDBFILE);
	    if (length($aapdb)<=0) {die("Error: chain $chain not found in pdb file $pdbfile\n");}

	    # Align template in hh-alignment ($aat) with template sequence in pdb ($aapdb)
	    my $xseq=$aat;
	    $xseq=~tr/-/~/; # transform Deletes to '~' to distinguish them from gaps '-' inserted by Align.pm   
	    my $yseq=$aapdb;
	    my ($jmin,$jmax,$lmin,$lmax);
	    my $Sstr;
	    my $score;  
	    # The aligned characters are returend in $j2[$col2] and $l2[$col2]
	    $score=&AlignNW(\$xseq,\$yseq,\@j2,\@l2,\$jmin,\$jmax,\$lmin,\$lmax,\$Sstr);  

	    # DEBUG
	    if ($v>=3) {
		printf("Template (hh)  $xseq\n");
		printf("Identities     $Sstr\n");
		printf("Template (pdb) $yseq\n");
		printf("\n");
		if ($v>=4) {
		    for ($col2=0; $col2<@l2 && $col2<1000; $col2++) {
			printf("%3i  %3i:%s  %3i:%s -> %i\n",$col2,$j2[$col2],substr($aat,$j2[$col2]-1,1),$l2[$col2],substr($aapdb,$l2[$col2]-1,1),$nres[$l2[$col2]-1]);
		    }
		}
	    }	
	    
            # check for reasonable alignment
	    my $num_match = 0;
	    for ($i=0; $i<@l2; $i++) {
		if ($j2[$i] > 0 && $l2[$i] > 0) {
		    $num_match++;
		}
	    }
	    if (($score/$num_match) < 1) {
		print "WARNING! Match score with PDBfile (score: $score   num: $num_match   score/num:".($score/$num_match).") to low => $pdbfile not included!\n";
		next;
	    }
   
	    # Assign a3m-formatted amino acid sequence from pdb file to $aapdb
	    $aapdb="";
	    my @xseq=unpack("C*",$xseq);
	    my @yseq=unpack("C*",$yseq);
	    for ($i=0; $i<@yseq; $i++) {
		if(($xseq[$i]>=65 && $xseq[$i]<=90) || $xseq[$i]==ord('~')) { # if $aat has upper case residue or Delete state
		    # Match state
		    $aapdb.=uc(chr($yseq[$i]));
		} else {
		    # Insert state
		    if ($yseq[$i]!=45) {$aapdb.=lc(chr($yseq[$i]));} # add only if not a gap '-'
		}
	    }
	    
	    # Remove overlapping ends of $aapdb
	    $aapdb=~s/^[a-z]*(.*?)[a-z]*$/$1/;

	    # Glue gaps at beginning and end of aligned pdb sequence and add sequence to alignment
	    push (@{$p_hitseqs}, ("-" x ($qfirst-1)).$aapdb.("-" x ($qlength-$qlast)) ); # use ATOM record residues $aapdb!
	    
	    # Write hitname in PIR format into @hitnames
	    my $descr;
	    my $organism;
	    my $struc=$pdbcode;
	    if ($tnameline=~/^(\S+)\s+(.*)/) {$descr=$2; $descr=~tr/://d;} else {$descr=" ";}
	    if ($tnameline=~/^(\S+)\s+.*\s+\{(.*)\}/) {$organism=$2;} else {$organism=" ";}
	    if (length($chain)>1 || $chain eq ".") { # MODELLER's special symbol for 'chain unspecified'
		$chain="."; 
	    } elsif ($addchain && $chain ne " ") {
		$struc.="_$chain";
	    } 
#	    push (@{$p_hitnames}, sprintf(">P1;%s\nstructureX:%4s:%4i:%1s:%4i:%1s:%s:%s:%-.2f:%-.2f\n",$struc,$struc,$nres[$lmin-1],$chain,$nres[$lmax-1],$chain,$descr,$organism,$resolution,$rvalue) );
	    
	    my $descrtmp=$descr;
	    $descrtmp=~tr/:/;/;
	    $organism=~tr/://d;
	    push (@{$p_hitnames}, sprintf(">P1;%s\nstructureX:%4s: :%1s: :%1s:%s:%s:%-.2f:%-.2f\n",$struc,$struc,$chain,$chain,$descrtmp,$organism,$resolution,$rvalue) );
	    push (@{$p_hitdiag}, $tfirst-$qfirst);
	} else {
	    # outformat is "FASTA" or "A2M" or "A3M" or ...
	    # Write hitname in FASTA format into @hitnames
	    push (@{$p_hitseqs}, ("-" x ($qfirst-1)).$aat.("-" x ($qlength-$qlast)) ); 
	    push (@{$p_hitnames}, ">$tname\n" );
	    push (@{$p_hitdiag}, $tfirst-$qfirst);
	} 
		
	if ($onlyfirst>0) {last;}  # extract only first (seed?) sequence in each alignment
	
    } # end:  for (my $k=0; $k<@{$tname}; $k++)
    
    # Paste aligned subsequence of query over $hitseqs[0] 
    if (${$p_hitseqs}[0] eq "") {${$p_hitseqs}[0] = "-" x $qlength;}
    if (!$qfile) {substr(${$p_hitseqs}[0],$qfirst-1,length($aaq),$aaq);}
    
    return;
}


##################################################################################
# Read Alignment from infile (*.hhr file)
# Results: 
# $aaq:    query residues in present alignment 
# $qfirst: index of first query residue in present alignment  
# @tname:  template names in present alignmen
# @tfirst: indices of first residues in present alignmet
# @tseq:   sequences of templates in present alignment
##################################################################################
sub ReadAlignment() {

    @qname=();     # name of $it'th query in this alignment
    @qfirst=();    # index of first residue in $it'th query in this alignment
    @qlast=();     # index of last residue in $it'th query in this alignment
    @qseq=();      # residues of $it'th query in this alignment

    @tname=();     # name of $it'th template in this alignment
    @tfirst=();    # index of first residue in $it'th template in this alignment
    @tlast=();     # index of last residue in $it'th template in this alignment
    @tseq=();      # residues of $it'th template in this alignment

    if ($v>=3) {printf("Searching for Q $qname vs T $tname\n");}
    $line=<INFILE>;    

    # Search for first line beginning with Q ot T and not followed by aa_, ss_pred, ss_conf, or Consensus
    while (1) {
	my $i; # index for query sequence in this alignment
	# Scan up to first line starting with Q; stop when line 'No\s+\d+' or 'Done' is found
	while (defined $line && $line!~/^Q\s(\S+)/) {
	    if ($line=~/^No\s+\d/ || $line=~/^Done/) {last;}
	    $line=<INFILE>; next;
	} 
	if (!defined $line || $line=~/^No\s+\d/ || $line=~/^Done/) {last;}

	# Scan up to first line that is not secondary structure line or consensus line	
	while (defined $line && $line=~/^Q\s+(ss_|sa_|aa_|Consens|Cons-)/)  {$line=<INFILE>;} 

	# Read next block of query sequences
	$i=0;
        while ($line=~/^Q\s+/) {
	    if ($line!~/^Q\s+(ss_|sa_|aa_|Consens|Cons-)/ && $line=~/^Q\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\((\d+)/) {
		$qname[$i]=$1;
		if (!$qfirst[$i]) {$qfirst[$i]=$2;} # if $qfirst is undefined then this is the first alignment block -> set $qfirst to $1
		if (!$qseq[$i]) {$qseq[$i]=$3;} else {$qseq[$i].=$3;}
		$qlast[$i]=$4; 
		if ($i==0) {$qlength=$5}
		$i++;
	    }
	    $line=<INFILE>;
	} 
	if ($i==0) {
	    die("\nError in $program: bad format in $infile, line $.: query block\n");
	}

	# Scan up to first line starting with T	
	while (defined $line && $line!~/^T\s+(\S+)/) {$line=<INFILE>;} 

	# Scan up to first line that is not secondary structure line or consensus line	
	while (defined $line && $line=~/^T\s+(ss_|sa_|aa_|Consens|Cons-)/)  {$line=<INFILE>;} 

	# Read next block of template sequences
	$i=0;
        while ($line=~/^T\s+/) {
	    if ($line!~/^T\s+(ss_|sa_|aa_|Consens|Cons-)/ && $line=~/^T\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\((\d+)/){
		$tname[$i]=$1;
		if (!$tfirst[$i]) {$tfirst[$i]=$2;} # if $tfirst is undefined then this is the first alignment block -> set $tfirst to $1
		if (!$tseq[$i]) {$tseq[$i]=$3;} else {$tseq[$i].=$3;}
		$tlast[$i]=$4; 
		if ($i==0) {$tlength=$5}
		$i++;
	    }
	    $line=<INFILE>;
	} 
	if ($i==0) {
	    die("\nError in $program: bad format in $infile, line $.: template block\n");
	}
	
    } # end while ($line=<INFILE>)  
    
#    if (!$qfirst)    {$qfirst=1;}  # if still $qfirst==0 then set $qfirst to 1
#    for (my $i=0; $i<@tfirst; $i++) {
#	if (!$tfirst[$i]) {$tfirst[$i]=1;}  # if still $tfirst[$i]==0 then set $tfirst to 1
#    }    
    
    # Check lengths
    if (length($qseq[0])!=length($tseq[0])) {
	print("\nError: query and template lines do not have the same length in $infile, line $.\n");
	for (my $i=0; $i<@qfirst; $i++) {
	    printf("Q %-14.14s %s\n",$qname[$i],$qseq[$i]);
	}
	printf("\n");
	for (my $i=0; $i<@tfirst; $i++) {
	    printf("T %-14.14s %s\n",$tname[$i],$tseq[$i]);
	}
	printf("\n");
	exit 1;
    }
    
    if ($v>=3) {
	for (my $i=0; $i<@qfirst; $i++) {
	    printf("Q %-14.14s %s\n",$qname[$i],$qseq[$i]);
	}
	printf("\n");
	for (my $i=0; $i<@tfirst; $i++) {
	    printf("T %-14.14s %s\n",$tname[$i],$tseq[$i]);
	}
	printf("\n");
    }	    
    return;
}


##################################################################################
# Write Alignment to $printblock[$k]
##################################################################################
sub WritePairwiseAlignments() {

    #Delete columns with gaps in both sequences
    $aaq=uc($aaq);
    $aat=uc($aat);
    @aaq=split(//,$aaq);
    @aat=split(//,$aat);
    my $col=0;
    for ($col1=0; $col1<@aaq; $col1++) {
	if ($aaq[$col1]=~tr/a-zA-Z/a-zA-Z/ || $aat[$col1]=~tr/a-zA-Z/a-zA-Z/) {
	    $aaq[$col]=$aaq[$col1];
	    $aat[$col]=$aat[$col1];
	    $col++;
	}
    }
    splice(@aaq,$col); # delete end of @aaq;
    splice(@aat,$col);
    $aaq=join("",@aaq);
    $aat=join("",@aat);
    
    
    # Count query and template residues into @i1 and @j1 
    for ($col1=0; $col1<@aaq; $col1++) {
	if ($aaq[$col1]=~tr/a-zA-Z/a-zA-Z/) {
	    $i1[$col1]=$qfirst++;  #found query residue in $col1
	} else {
	    $i1[$col1]=0;     #found gap in $col1
	}
	if ($aat[$col1]=~tr/a-zA-Z/a-zA-Z/) {
	    $j1[$col1]=$tfirst++;  #found template residue in $col1
	} else {
	    $j1[$col1]=0;     #found gap in $col1
	}
    }
    

    # DEBUG
    if ($v>=3) {
	printf ("col    Q  i1     T  j1\n");
	for ($col1=0; $col1<@aaq; $col1++) {
	    printf ("%3i    %s %3i     %s %3i\n",$col1,$aaq[$col1],$i1[$col1],$aat[$col1],$j1[$col1]);
	}
	printf ("\n");
    }


    # Read protein chain from pdb file
# ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
# ATOM      1  N   SER A  27      38.637  79.034  59.693  1.00 79.70          # ATOM   2083  CD1 LEU A  22S     15.343 -12.020  43.761  1.00  5.00           C
    
    # Extract pdbcode and construct name of pdbfile and return in global variables $pdbid and $chain
    if (&ExtractPdbcodeAndChain($tname)) {next;}

    # Read sequence from pdb file
    if (! defined $pdbfile) {die ("Error in $program: Couldn't find pdb code in $tname\n");}
    open (PDBFILE, "$pdbfile") || die ("Error in $program: Couldn't open $pdbfile: $!\n");
    if ($chain eq "[A ]") {$pdbcode.="_A";} elsif ($chain eq ".") {;} else {$pdbcode.="_$chain";}
    $aapdb=""; $l=1;
    $line=<PDBFILE>; 
    while ($line) {if ($line=~/^ATOM/) {last;} $line=<PDBFILE>;} # advance to ATOM records
    my @nres;        # $nres[$l] = pdb residue index for residue $aapdb[$l]
    my @coord;       # $coord[$l] = coordinates of CA atom of residue $aapdb[$l]
    while ($line) {
	if ($line=~/^ATOM\s+\d+  CA [ A](\w{3}) $chain\s*(-?\d+.)   (\s*\S+\s+\S+\s+\S+)/ ||
	    ($line=~/^HETATM\s+\d+  CA [ A](\w{3}) $chain\s*(-?\d+.)   (\s*\S+\s+\S+\s+\S+)/ && &Three2OneLetter($1) ne "X") ) {
	    $res=$1;
	    $nres[$l]=$2;
	    $coord[$l]=$3."  1.00";
	    $res=&Three2OneLetter($res);
	    $aapdb[$l]=$res;
	    $aapdb.=$res;
	    $l++;
	}
	elsif ($l>10 && $line=~/^ATOM\s+\d+  CA/) {last;}
	elsif ($line=~/^ENDMDL/) {last;} # if file contains NMR models read only first one
	$line=<PDBFILE>;
    }
    close (PDBFILE);
    
    # Align template in hh-alignment ($aat) with template sequence in pdb ($aapdb)
    
    my $xseq=$aat;
    my $yseq=$aapdb;
    my ($jmin,$jmax,$lmin,$lmax);
    my $Sstr;
    my $score;  
    $xseq=~tr/-/~/d; # transform Deletes to '~' to distinguish them from gaps inserted by Align.pm   
    #the aligned characters are returend in $j2[$col2] and $l2[$col2]
    $score=&AlignNW(\$xseq,\$yseq,\@j2,\@l2,\$jmin,\$jmax,\$lmin,\$lmax,\$Sstr);  
    
    # DEBUG
    if ($v>=3) {
	printf("Template (hh)  $xseq\n");
	printf("Identities     $Sstr\n");
	printf("Template (pdb) $yseq\n");
	printf("\n");
	if ($v>=4) {
	    for ($col2=0; $col2<@l2 && $col2<200; $col2++) {
		printf("%3i  %3i  %3i\n",$col2,$j2[$col2],$l2[$col2]);
	    }
	}
    }	
    
    # DEBUG
    
    # Construct alignment of $aaq <-> $aapdb via alignments $aaq <-> $aat and $aat <-> $aapdb:  
    # Find $l1[$col1] = line of pdb file corresponding to residue $aat[$col1] and $aaq[$col1]
    $col2=0;
    for ($col1=0; $col1<@aaq; $col1++) {
	if ($j1[$col1]==0 || $i1[$col1]==0) {$l1[$col1]=0; next;} # skip gaps in query and gaps in template
	while ($j2[$col2]<$col1+1) {$col2++;} # in $j2[col2] first index is 1, in $col1 first column is 0
	$l1[$col1] = $l2[$col2];
	if ($v>=4) {printf("l1[%i]=%i  l2[%i]=%i\n",$col1,$l1[$col1],$col2,$l2[$col2]);}
    }
    
    
    if ($pdbcode ne "NONE") {
	if ($outformat eq "TS") {
	    for ($col1=0; $col1<@aat; $col1++) {
		if ($i1[$col1]==0) {next;} # skip gaps in query
		if ($j1[$col1]==0) {next;} # skip gaps in template sequence
		if ($l1[$col1]==0) {next;} # skip if corresponding residue was skipped in pdb file
		
		$printblock[$k].=sprintf("ATOM  %5i  CA  %3s  %4i    %-50.50s\n",$i1[$col1],&One2ThreeLetter($aaq[$col1]),$i1[$col1]+$shift,$coord[$l1[$col1]]);
		if ($v>=4) {
		    printf("ATOM  %5i  CA  %3s  %4i    %-50.50s\n",$i1[$col1],&One2ThreeLetter($aaq[$col1]),$i1[$col1]+$shift,$coord[$l1[$col1]]);
		}
	    }
	} else {
	    for ($col1=0; $col1<@aat; $col1++) {
		if ($i1[$col1]==0) {next;} # skip gaps in query
		if ($j1[$col1]==0) {next;} # skip gaps in template sequence
		if ($l1[$col1]==0) {next;} # skip if corresponding residue was skipped in pdb file
		$printblock[$k].=sprintf("%1s %3i    %1s %s\n",$aaq[$col1],$i1[$col1],$aat[$col1],$nres[$l1[$col1]]);
		if ($v>=4) {printf("%1s %3i    %1s %s\n",$aaq[$col1],$i1[$col1],$aat[$col1],$nres[$l1[$col1]]);}
	    }
	}
    }
    $printblock[$k].=sprintf("TER\n");
    return;
}



# Extract pdbcode and construct name of pdbfile and return in global variables $pdbid and $chain
sub ExtractPdbcodeAndChain() 
{
    my $name=$_[0];
    $name=~/^(\S+)/;
    $name=$1;

    # SCOP ID? (d3lkfa_,d3grs_3,d3pmga1,g1m26.1)
    if ($name=~/^[defgh](\d[a-z0-9]{3})([a-z0-9_.])[a-z0-9_]$/) {
	$pdbcode=$1;
	if ($2 eq "_") {$chain="[A ]";} else {$chain=uc($2);}
    } 
    
    # PDB ID? (8fab, 1a0i)
    elsif ($name=~/^(\d[a-z0-9]{3})$/) {
	$pdbcode=$1;
	$chain="[A ]";
    }

    # PDB ID? (8fab_A)
    elsif ($name=~/^(\d[a-z0-9]{3})_(\S)$/) {
	$pdbcode=$1;
	$chain=$2;
    }
    
    # PDB ID? (1u1z_ABC)
    elsif ($name=~/^(\d[a-z0-9]{3})_(\S\S+)$/) {
	$pdbcode=$1;
	$chain="[$2]";
    }
    
    # DALI ID? (8fabA_0,1a0i_2)
    elsif ($name=~/^(\d[a-z0-9]{3})([A-Za-z0-9]?)_\d+$/) {
	$pdbcode=$1;
	$chain=$2;
    }
    
    else {
	$pdbcode=$name;
	$chain="A";
#	return 1; # no SCOP/DALI/pdb sequence 
    }

    my $div=substr($pdbcode,1,2);
    $pdbfile = "";

    foreach $pdbdir (@pdbdirs) {
	#print "PDB-dir: $pdbdir\n";
#	if (-e "$pdbdivdir/$div/pdb$pdbcode.ent")   {$pdbfile="$pdbdivdir/$div/pdb$pdbcode.ent"; last;}
#	if (-e "$pdbdivdir/$div/pdb$pdbcode.ent.Z") {$pdbfile="gunzip -c $pdbdivdir/$div/pdb$pdbcode.ent.Z |"; last;}
	if (-e "$pdbdir/pdb$pdbcode.ent") {$pdbfile="$pdbdir/pdb$pdbcode.ent"; last;}
	if (-e "$pdbdir/$pdbcode.pdb")    {$pdbfile="$pdbdir/$pdbcode.pdb"; last;}
	if (-e "$pdbdir/$name.pdb")       {$pdbfile="$pdbdir/$name.pdb"; last;}
	if (-e "$pdbdir/$pdbcode"."_$chain.pdb")    {$pdbfile="$pdbdir/$pdbcode"."_$chain.pdb"; last;}
    }

    if ($pdbfile eq "") {
	if ($v>=2) {print("Warning: no pdb file found for sequence name '$name'\n");} 
	return 1;
    }

    return 0;
}


# Resort arrays according to sorting array0:
# Resort(\@array0,\@array1,...,\@arrayN)
sub Sort() 
{
    my $p_array0 = $_[0];
    my @index=();
    for (my $i=0; $i<@{$p_array0}; $i++) {$index[$i]=$i;}
    @index = sort { ${$p_array0}[$a] <=> ${$p_array0}[$b] } @index;
    foreach my $p_array (@_) {
	my @dummy = @{$p_array};
	@{$p_array}=();
	foreach my $i (@index) {
	    push(@{$p_array}, $dummy[$i]);
	}
    }
}

##################################################################################
# Convert three-letter amino acid code into one-letter code
##################################################################################
sub Three2OneLetter {
    my $res=uc($_[0]);
    if    ($res eq "GLY") {return "G";}
    elsif ($res eq "ALA") {return "A";}
    elsif ($res eq "VAL") {return "V";}
    elsif ($res eq "LEU") {return "L";}
    elsif ($res eq "ILE") {return "I";}
    elsif ($res eq "MET") {return "M";}
    elsif ($res eq "PHE") {return "F";}
    elsif ($res eq "TYR") {return "Y";}
    elsif ($res eq "TRP") {return "W";}
    elsif ($res eq "ASN") {return "N";}
    elsif ($res eq "ASP") {return "D";}
    elsif ($res eq "GLN") {return "Q";}
    elsif ($res eq "GLU") {return "E";}
    elsif ($res eq "CYS") {return "C";}
    elsif ($res eq "PRO") {return "P";}
    elsif ($res eq "SER") {return "S";}
    elsif ($res eq "THR") {return "T";}
    elsif ($res eq "LYS") {return "K";}
    elsif ($res eq "HIS") {return "H";}
    elsif ($res eq "ARG") {return "R";}

    # The HETATM selenomethionine is read by MODELLER like a normal MET in both its HETATM_IO=off and on mode!
    elsif ($res eq "MSE") {return "M";} # SELENOMETHIONINE 
    elsif ($res eq "ASX") {return "B";}
    elsif ($res eq "GLX") {return "Z";}
    else                  {return "X";}

    # The following post-translationally modified residues are ignored by MODELLER in its default SET HETATM_IO=off mode
#    elsif ($res eq "SEC") {return "C";} # SELENOCYSTEINE
#    elsif ($res eq "SEP") {return "S";} # PHOSPHOSERINE 
#    elsif ($res eq "TPO") {return "T";} # PHOSPHOTHREONINE 
#    elsif ($res eq "TYS") {return "Y";} # SULFONATED TYROSINE 
#    elsif ($res eq "KCX") {return "K";} # LYSINE NZ-CARBOXYLIC ACID
}

##################################################################################
# Convert one-letter amino acid code into three-letter code
##################################################################################
sub One2ThreeLetter {
    my $res=uc($_[0]);
    if    ($res eq "G") {return "GLY";}
    elsif ($res eq "A") {return "ALA";}
    elsif ($res eq "V") {return "VAL";}
    elsif ($res eq "L") {return "LEU";}
    elsif ($res eq "I") {return "ILE";}
    elsif ($res eq "M") {return "MET";}
    elsif ($res eq "F") {return "PHE";}
    elsif ($res eq "Y") {return "TYR";}
    elsif ($res eq "W") {return "TRP";}
    elsif ($res eq "N") {return "ASN";}
    elsif ($res eq "D") {return "ASP";}
    elsif ($res eq "Q") {return "GLN";}
    elsif ($res eq "E") {return "GLU";}
    elsif ($res eq "C") {return "CYS";}
    elsif ($res eq "P") {return "PRO";}
    elsif ($res eq "S") {return "SER";}
    elsif ($res eq "T") {return "THR";}
    elsif ($res eq "K") {return "LYS";}
    elsif ($res eq "H") {return "HIS";}
    elsif ($res eq "R") {return "ARG";}
    elsif ($res eq "U") {return "SEC";}
    elsif ($res eq "B") {return "ASX";}
    elsif ($res eq "Z") {return "GLX";}
    else                {return "UNK";}
}
