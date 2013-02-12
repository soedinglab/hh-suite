#! /usr/bin/env perl
# pdbfilter.pl - Read pdb or SCOP sequences from infile and write representative set of sequences to outfile

#
#     HHsuite version 2.0.16 (April 2013)
#
#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#     (C) Johannes Soeding, 2012

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

$|= 1; # Activate autoflushing on STDOUT

# Default options
my $idmax=90;      # maximum sequence identity threshold 
my $Evalmin=0.01;  # minimum BLAST E-value (should be <0.01 to ensure that sequences being filtered out are homologous to some representative sequence
my $covmin=90;     # minimum coverage threshold (should be large enough to ensure that no structural domain is lost from the filtered set just because it occurs in a sequence with, say, 2 other domains similar to some representative sequence) 
my $v=2;
my $blastpgp="$ncbidir/blastpgp";

my $help="
 pdbfilter.pl - Read pdb or SCOP sequences from infile and write representative set of 
 sequences to outfile
 Compares each sequence with all others using BLAST (blastpgp). If two sequences A and B are
 sufficiently similar, the sequence with lower resolution will be removed. 
 The exact criterion for removal of sequence B is: 
 IF more than \$covmin\% residues of sequence B are aligned to sequence A
 AND the sequence identity of the alignment is larger than \$idmin 
 AND the E-value is better than \$Evalmin
 AND sequence B has lower resolution than A. 
 The input file must have been prepared with pdb2fasta.pl or scop2fasta.pl. 
 Sequences with fewer than 15 residues are suppressed.

 Usage:   pdbfilter.pl infile filtered_file [-u old_filtered_file] [-id int] [-cov int]
 Options
  -id <int>  maximum sequence identity between representative sequences (default=$idmax)
  -e <float> minimum E-value between representative sequences (default=$Evalmin)
  -cov <int> minimum coverage of a sequence with a similar one to get thrown out (default=$covmin)
  -u <file>  update the old filtered file; this saves a lot of execution time 
  -v <int>   verbose mode
 Example: pdbfilter.pl pdb.fas pdb70.fas -u pdb70.fas -id 70 -cov 90\n
";
if (@ARGV<2) {print($help); exit;}

my $infile;
my $outfile;
my $oldfile;
my $root="";       # $outfile without extension
my $line;
my $pdbid="";      # e.g. 1ahs_C
my $qpdbid;        # pdbid of query sequence
my $seq="";        # sequence record (name+residues) in FASTA
my @seqs;          # sequence records as they were read
my $len=0;         # length of sequence to be read in
my %len;           # $len{$pdbid} is length of sequence 
my $resol;         # experimental resolution in Angstroem 
my %resol;         # experimental resolution in Angstroem 
my %excluded;      # representative sequences are all those not in this hash
my %accepted;      # representative sequences accpeted so far
my %similar;       # $similar{$pdbid} contains list of all pdbids (SCOPIDs) that are represented by (i.e. similar to) $qpdbid
my %het;           # $het{$pdbid} is "*" if sequence $pdbid has a HET: record (i.e. contains HET group with >=10 atoms), "" otherwise
my $id;            # sequence identity to query
my $cov;           # coverage 
my $Evalue;        # BLAST E-value
my $nold=0;        # number of sequences in oldfile
my $ntot=0;        # total number of sequences in oldfile and infile
my $k=0;           # counts sequences read in so far
my $sort_by_family=1; # set to 0 if at least one non-SCOP sequence found 


# Read command line options
my $options="";
for (my $i=0; $i<=$#ARGV; $i++) {$options.=" $ARGV[$i]";}
if ($options=~s/ -id\s+(\S+)//) {$idmax=$1;}
if ($options=~s/ -cov\s+(\S+)//) {$covmin=$1;}
if ($options=~s/ -e\s+(\S+)//)  {$Evalue=$1;}
if ($options=~s/ -u\s+(\S+)//) {$oldfile=$1;} 
if ($options=~s/ -v\s*(\d+)//) {$v=$1;} 
if ($options=~s/ -v//) {$v=2;} 
if (!$infile  && $options=~s/^\s*([^- ]\S+)\s*//) {$infile=$1;} 
if (!$outfile && $options=~s/^\s*([^- ]\S+)\s*//) {$outfile=$1;} 

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!$infile)  {print($help); print("Error: no input file given\n"); exit;}
if (!$outfile) {print($help); print("Error: no output file given\n"); exit;}

$covmin*=0.01;
if ($outfile=~/^(\S*?)\.(.*)$/) {$root=$1;} else {$root=$outfile;}

# Read sequences from oldfile
if ($oldfile) {$nold=&ReadSequences($oldfile,0);}

# Read NEW sequences from infile (the ones that are not yet in oldfile)
$ntot=&ReadSequences($infile,$nold);

# Sort by resolution
@seqs=sort SortByResolution @seqs;
#foreach $seq (@seqs) {print("$seq");}

# Record resolution and length of all sequences in hashes
for ($k=0; $k<@seqs; $k++) {
    $seqs[$k]=~s/^(\S+?)!(\S+?)!>(\S+)(.*)/>$3$4/o;
    $len{$3}=$1;
    $resol{$3}=$2;
    my $pdbid=$3;
    if ($seqs[$k]=~s/ PDB:\s*(( \d\S{3,5}\*?)+)//) {
	$similar{$pdbid}=$1;
    } elsif ($seqs[$k]=~s/ SCOPID:\s*(( [a-z]\d\S{5}\*?)+)//) {
	$similar{$pdbid}=$1;
    } else {
	$similar{$pdbid}="";
    }
    if ($seqs[$k]=~/ HET:/) {
	$het{$pdbid}="*";
    } else {
	$het{$pdbid}="";
    }
}


# Format BLAST database, initialize
if ($v>=2) {printf("Doing PSI-BLAST runs\n");}
&PrintSequences("$root.db",$nold);
&System("$ncbidir/formatdb -i $root.db");
my $nacc=0;                # number of accepted sequences = number of BLAST searches already performed
my $done=0;                # number of sequences already processed
my $nexcl=0;               # number of already excluded chains
my $step=($ntot>1000)*0.05+($ntot<=1000)*0.1;
my $reformat=$step;        # when X% of the total number of sequences have been excluded, reformat database

####################################################################################################
# For each sequence in infile, do BLAST search and exclude hits that are too similar to query
foreach $seq (@seqs) {
    $seq=~/^>(\S+)/o;
    $done++;
    if($excluded{$1}) {next;}
    $qpdbid=$1;
    $resol=$resol{$1};
    $len=$len{$1};

    if ($v>=2 && !($nacc%100)) {
	printf("\nSearches:%-4i Done:%3i%%  DB-size:%3i%% ",$nacc,100*($nexcl+$nacc)/$ntot+0.5,100*($ntot-$nexcl)/$ntot+0.5);
    }

    # When a substantial number of sequences have been excluded, reformat database WITHOUT excluded sequences
    if (!$oldfile && ($nexcl+$nacc)/$ntot>$reformat) {
	if ($v>=2) {printf("\b\b\b%2i%%",100*$reformat+0.5);}
	elsif ($v>=3) {
	    printf("\nReformatting search database containing %i out of %i sequences\n",$ntot-$nexcl,$ntot);}
	# Write updated database with all seqs with index >=$done that have not yet been excluded
	# (Don't have to compare A->B and B->A, hence exclude seqs with index <$done)
	&PrintSequences("$root.db",$done);  
	&System("$ncbidir/formatdb -i $root.db");
	$reformat+=$step;
    }

    open (TMPFILE,">$root.tmp") || die ("ERROR: cannot open $root.tmp for writing: $!\n");
    print(TMPFILE $seq);
    close(TMPFILE);

    # Find hits that are too similar
    if ($v>=3) {print("\n$blastpgp -i $root.tmp -d $root.db -v 1 -b 1000 -s T -z $ntot");}
    open(BLAST,"$blastpgp -i $root.tmp -d $root.db -v 1 -b 1000 -s T -z $ntot|");
    while ($line=<BLAST>){
	
	if ($line=~/^>(\S+)/o && $1 ne $qpdbid && !$excluded{$1} && !$accepted{$1}) {
	    $pdbid=$1;
	    while ($line=<BLAST>){if ($line=~/^ Score = /o) {last;}}
	    $line=~/ Expect =\s*(\S+)/ or die("Error: format error in '$blastpgp -i $root.tmp -d $root.db -v 1 -b 1000 -s T -z $ntot|', line $.\n");
	    $Evalue=$1;
	    $Evalue=~s/^e/1e/;
	    $Evalue=~s/,$//;
	    $line=<BLAST>;
	    $line=~/^ Identities =\s*\d+\/(\d+)\s+\((\d+)\%\)/o or die("Error: format error in '$blastpgp -i $root.tmp -d $root.db -v 1 -b 1000 -s T -z $ntot|', line $.\n");
	    $len=$1;
	    $id=$2;
	    # Coverage = (length of whole alignment (including gaps) - gaps in query or HSP) / total length of matched sequence
	    if ($line=~/Gaps =\s*(\d+)\/\d+/) {$cov=($len-$1)/$len{$pdbid};} else {$cov=$len/$len{$pdbid};} 

	    ## Main filtering criterion: remove sequence from representative set if... 
	    if ($id>=$idmax && $Evalue<=$Evalmin && $cov>=$covmin && $resol<=$resol{$pdbid}) {
	
		$excluded{$pdbid}=1; 
#		if ($similar{$qpdbid} ne "") {$similar{$qpdbid}.=",";} # separate pdbids with identical sequences with ','
		if (!$similar{$qpdbid}) {$similar{$qpdbid}="";}
		$similar{$qpdbid}.=" ".$pdbid.$het{$pdbid}.$similar{$pdbid};
		$nexcl++;
		if ($v>=4) {
		    print($line);
		    printf("Rep: $qpdbid  excl: $pdbid  (len=%3i cov=%3.0f  id=$1)\n",$len{$pdbid},100*$cov);
		}
	    }
	}
    }
    close(BLAST);
    $nacc++;
    $accepted{$qpdbid}=1;
    if ($v>=2) {print(".");}
}
if ($v>=2) {print("\n");}
####################################################################################################

if ($sort_by_family) {@seqs=sort SortByFamily @seqs;}

# Print out all representative sequences to outfile
&PrintSequences($outfile,0);
if ($v>=2) {
    printf("Written %i out of %i sequences to $outfile\n",$ntot-$nexcl,$ntot);
}
unlink("$root.tmp");
unlink(glob("$root.db*"));
exit;




# Read sequences in infile beginning at index $k
sub ReadSequences() 
{
    my $infile=$_[0];
    my $k=$_[1];
    my $k0=$k;
    my $nres=0; # skip sequences with fewer than 15 real residues
    if ($v>=3) {printf("Reading $infile ... \n");}
    open (INFILE,"<$infile") || die ("ERROR: cannot open $infile for reading: $!\n");
    while ($line=<INFILE>) {
	if ($line=~/^>/) {
	    if ($pdbid && !$len{$pdbid} && $nres>=15) {
		if ($len<26) {chomp $seq; $seq.=('x'x(26-$len))."\n";} # add 'x' to make $seq at least 26 resiudes long
		$seqs[$k++]="$len!$resol!$seq";	    
	    }
	    # Read pdbid (or SCOP id) and resolution from name line
	    if ($line=~/^>(\S{4,6}) .*; (\d+\.?\d*|NMR)A? \{/o) {
		# PDB sequence
		$pdbid=$1;
		if ($2 eq "NMR") {$resol=10;} else {$resol=$2;}
		$seq=$line;
		$len=$nres=0;
		$sort_by_family=0;
	    } elsif ($line=~/^>([a-z]\S{6}) .* (\d+\.?\d*)\s*$/) {
		# SCOP sequence
		$pdbid=$1;
		if ($2 eq "" or $2 eq "NMR") {$resol=10;} else {$resol=$2;}
		$resol=$2;
		$line=~s/^(>.*) (\d+\.?\d*)\s*$/$1\n/; # remove resolution at the end of the name line
		$seq=$line;
		$len=0;
	    } else {print("WARNING: nameline format not recognized: $line");}
	    $nres=0;
	} else {
	    $seq.=$line;
	    $len+=length($line)-1;
	    $nres+=($line=~tr/a-wyA-WY/a-wyA-WY/);
	}
    }
    if ($pdbid && !$len{$pdbid}) {
	if ($len<26) {$seq.=('X'x(26-$len))."\n";} # add 'X' to make $seq at least 26 resiudes long
	$seqs[$k++]="$len!$resol!$seq";	    
    }
    close(INFILE);
    if ($v>=2) {printf("Read %i sequences from $infile ... \n",$k-$k0);}
    return $k;
}


# Print all sequences that have not been excluded so far to file, beginning at index $k
sub PrintSequences() 
{
    open (OUTFILE,">$_[0]") || die("Error: could not write to $_[0]: $!\n");
    for (my $k=$_[1]; $k<@seqs; $k++) {
	$seqs[$k]=~/>(\S+)/o;
	if (!$excluded{$1}) {
	    my $seq=$seqs[$k];
	    if ($similar{$1} ne "") {
		my $pdbs=$similar{$1};
		# Limit number of represented pdbs to 30
		if ($pdbs =~/^(\s+\S+){31,}/) {$pdbs =~/^((?:\s+\S+){30})/; $pdbs=$1." ...";}
		if ($pdbs=~/^\s*\d\S{3,5}/) {
		    $seq=~s/^(.*)/$1 PDB:$pdbs/;
		} elsif ($pdbs=~/^\s*[a-z]\d\S{5}/) {
		    $seq=~s/^(.*)/$1 SCOPIDS:$pdbs/;
		}
	    }
	    print(OUTFILE $seq);
	}
    }
    close(OUTFILE);
}

sub SortByResolution() {
    my $aa;
    my $bb;
    if ($a!~/^\d+!(\S+)!>/o) {
	$a=~/^(\S+)/o;
	printf("Error: sequence %s does not contain resolution in right format\n",$1);
	return 0;
    }
    $aa=$1;
    if ($b!~/^\d+!(\S+)!>/o) {
	$b=~/^(\S+)/o;
	printf("Error: sequence %s does not contain resolution in right format\n",$1);
	return 0;
    }
    $bb=$1;
    return ($aa<=>$bb);
}

sub SortByFamily() {
    my $aa;
    my $bb;
    if ($a!~/^>\S+ (\S+)/) {
	$a=~/^(\S+)/o;
	printf("Error: sequence %s does not contain resolution in right format\n",$1);
	return 0;
    }
    $aa=$1;
    if ($b!~/^>\S+ (\S+)/) {
	$b=~/^(\S+)/o;
	printf("Error: sequence %s does not contain resolution in right format\n",$1);
	return 0;
    }
    $bb=$1;
    return ($aa cmp $bb);
}

sub System() {
    $v--;
    &HHPaths::System($_[0]);
    $v++;
}
