#! /usr/bin/perl -w
# Read a SCOP or PDB sequence file (with SCOP or PDB sequence identifier in the 
# query sequence) and generate a pdb file for the query sequence. The pdb file 
# will contain only the residues present in the query sequence. The residue numbers
# in columns 23-26 refer to the input sequence.
# Usage: makepdbfile.pl infile [outfile]

my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;
$|=1; # autoflush on

# Default parameters for Align.pm
our $d=3;    # gap opening penalty for Align.pm
our $e=0.1;  # gap extension penalty for Align.pm
our $g=0.09; # endgap penalty for Align.pm
our $v=2;    # verbose mode
our $matrix="identity";

# Global variables
my $infile;            # input file
my $outfile;           # output file
my $line;              # input line
my $nameline="";       # characters following '>' in input file
my $aaq="";            # amino acid sequence read from input file
my $outdir;   # directory path of input file
my $program=$0;        # name of perl script
my $replaceMSE=1;      # replace HETATM MSE records by ATOM MET
$program=~s/.*\///;    # remove path
my $dopt=0;
my $altloc="A";

if (scalar(@ARGV)<1) {
    print("
 $program - generate a PDB file for a SCOP or PDB sequence with residue numbers 
 corresponding to this sequence
 The program reads a sequence which must have a SCOP domain, PDB chain, or DALI domain identifier 
 (e.g. d1hz4a_, 1hz4_A, or 1hz4A_1) from a FASTA file. It generates a PDB formatted file by 
 aligning the sequence to the sequence extracted from the ATOM records of the corresponding pdb 
 chain. HETATM records for MSE (selenomethionine) will be interpreted as ATOM records for MET in 
 the alignment. MSE will be changed to MET in the output file.
 
 Usage:   $program [options] infile [outfile] 
 Example: $program d1hz4a_.a3m /cluster/tmp/d1hz4a_.pdb 

 Options: 
   -d <pdbdir>  give directory of pdb files (default=$pdbdir)
   -a [A|B]     filter alternative location identifier (e.g. A or B)
   -o <file>    output file
\n"); 
    exit;
}


my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

# General options
if ($options=~s/ -d\s*(\S+) / /) {$pdbdir=$1; $dopt=1;}
if ($options=~s/ -a\s*(\S+) / /) {$altloc=$1;}
if ($options=~s/ -i\s*(\S*) //) {$infile=$1;}
if ($options=~s/ -o\s*(\S*) //) {$outfile=$1;}
if (!$infile  && $options=~s/^\s*([^-]\S*)\s* / /) {$infile=$1;} 
if (!$outfile && $options=~s/^\s*([^-]\S*)\s* / /) {$outdir=$1;} 
if (!$outfile) {
    $outfile=$infile;
    $outfile=~s/^(.*)\..*?$/$1.pdb/;
}

# Read query sequence in infile
open(INFILE,"<$infile") or die("Error: can't open $infile: $!\n");
while ($line=<INFILE>) {  # find query sequence
    if ( $line=~/^>/ && $line!~/^>(aa|ss)_/){last;}
}
$line=~/^>(.*)/;
$nameline=$1;
# Read query sequence
while ($line=<INFILE>) {
    if($line=~/^>/) {last;}
    chomp($line);
    $aaq.=$line;
}

if (&MakePdbFile($nameline,$aaq,$outfile) !=0) {exit(1);}
close(INFILE);
if ($v>=2) {print("Done\n");}
exit(0);



##################################################################################
# Make a pdb file $base.pdb with correct resdiue numbers from query sequence in $base.hhm 
##################################################################################
sub MakePdbFile() 
{


# >g1avo.1 a.24.8.1 (A:,B:) Proteasome activator reg(alpha) {Human (Homo sapiens)}    
    my $nameline=$_[0]; # everything in line following '>' in infile 
    my $aaq=$_[1];      # query amino acids
    my $outfile=$_[2];  # query amino acids
    my $pdbid;          # e.g. 1hz4
    my $chain;          # chain id of query sequence for hmmfile 
    my $line;           # currently read input line
    my $header;         # PDB file header
    my @remarks;        # PDB file remarks
    my $date=`date`;
    chomp($date);
    $nameline=~/^(\S+)/;
    my $name=$1;

    # SCOP ID? (d3lkfa_,d3grs_3,d3pmga1,g1m26.1)
    if ($nameline=~/^[d-u](\d[a-z0-9]{3})([a-z0-9_.])[a-z0-9_]\b/) {
	$pdbid=lc($1);
	$chain=uc($2);
	if ($chain eq "_") {$chain="[A ]";} else {$chain=uc($2);} 
	$header="HEADER    SCOP domain $name                    $date";
    }
    
    # DALI ID? (8fabA_0,1a0i_2)
    elsif ($nameline=~/^(\d[a-z0-9]{3})([A-Za-z0-9])?_\d+\s+\d+\.\d+.\d+.\d+.\d+.\d+/) {
	$pdbid=lc($1);
	$chain=$2;
	if ($chain eq "") {$chain="A";}
	$header="HEADER    DALI domain $name                    $date";
    }
    
    # PDB ID? (8fab_A, 1a0i, 1a0i_2)
    elsif ($nameline=~/^(\d[a-z0-9]{3})_?(\S?)\b/) {
	$pdbid=lc($1);
	$chain=$2;
	if ($chain eq "") {$chain="[A ]";}
	$header="HEADER    PDB chain $name                      $date";
    }

    # T0289_B or similar
    elsif ($dopt && $nameline=~/^([^_\s]+)/) {
	$pdbid=$1;
	if ($nameline=~/^[^_\s]+_(\S)/) {$chain=$1;} else {$chain=".";}
	$header="HEADER    $pdbid                             $date";
    }    
    else {
	die("Error: no pdb code found in sequence '$nameline'\n");
    }
    my $pdbfile="$pdbdir/pdb$pdbid.ent";
    if ($v>=2) {print("pdb file: $pdbfile  chain: $chain\n")}

    my $i;          # index for residues in scop sequence
    my $l;          # index for residue record in pdb file: l-1'st character in $aapdb <=> l'th record of $pdbline
    my @pdbline;    # $pdbline[$l][$natom] = line in pdb file for l'th residue in $aapdb and atom N, CB, or O 
    my $natom;      # runs from 0 up to 2 (N, CB, O)
    my $res;        # residue in pdb line
    my $atom;       # atom code in pdb file (N, CA, CB, O, ...)
    my $aapdb="";   # template amino acids from ATOM record of pdb file
    my $col;        # column of alignment query (from hhm file) versus pdb-residues
    my $nres;       # residue number in pdb file

    # Read the pdb file
    if (!open (PDBFILE,"<$pdbfile")) {
	$pdbfile="$pdbdir/$pdbid.pdb";
	if (!open (PDBFILE,"<$pdbfile")) {
	    if (!open (PDBFILE,"</raid/db/pdb_old/pdb$pdbid.ent")) {
		warn("Error: can't open $pdbfile: $!\n"); return 2;
	    }
	}
    }
    $l=0;
    $nres=-1e6;
    while ($line=<PDBFILE>) {
# ATOM      1  N   GLY A   1     -19.559   8.872   4.925  1.00 16.44           N
# ATOM      2  CA  GLY A   1     -19.004   8.179   6.112  1.00 14.30           C
	if ($line=~/^ENDMDL/) {last;} # if file contains NMR models read only first

	if ($line=~/^ATOM  \s*\d+ (....)[ $altloc](\w{3}) $chain\s*(-?\d+\w?).*/ ||
	    $line=~/^HETATM\s*\d+ (....)[ $altloc](MSE) $chain\s*(-?\d+).*/) {
	    $atom=$1;
	    $res=$2;
	    # New residue?
	    if ($3 ne $nres) {  # $3<$nres if new chain (A:,B:)
 		$nres=$3;
		$l++;
		$res=&Three2OneLetter($res);
		$aapdb.=$res;
		$natom=0;
	    }
	    if ($replaceMSE) {
		$line=~s/^(HETATM\s*\d+ )SE(...MSE)/$1 S$2/ && $line=~s/SE(\s*)$/ S$1/;
		$line=~s/^HETATM(\s*\d+ .....)MSE/ATOM  $1MET/;
	    }
	    $pdbline[$l][$natom++]=$line;
	}
    }
    close (PDBFILE);

    # Align scop query sequence ($aaq) with query sequence in pdb ($aapdb)	
    my $xseq=$aaq;
    my $yseq=$aapdb;
    my ($imin,$imax,$lmin,$lmax);
    my $Sstr;
    my $score;  
    my (@i,@l);    # The aligned characters are returend in $j[$col] and $l[$col]
    $score=&AlignNW(\$xseq,\$yseq,\@i,\@l,\$imin,\$imax,\$lmin,\$lmax,\$Sstr);  
    
    # DEBUG
    if ($v>=3) {
    }	
    
    # Set remarks etc.
    my $author ="AUTHOR    J. Soeding    johannes\@soeding.com";
    push(@remarks,"REMARK    This file was generated by $program. Its residue numbers refer");
    push(@remarks,"REMARK    to the following SEQRES sequence, for which $program was called:");
    push(@remarks,"REMARK    >$nameline");
    $line=uc($aaq);
    while ($line) {
	$line=~s/^(\S{1,60})//;
	push(@remarks,"REMARK    $1");
    }

    # Print scop-pdb file
    if (!open (OUTFILE,">$outfile")) {warn("Error: can't open $outfile: $!\n"); return 3;}
    printf(OUTFILE "%-80.80s\n",$header);
    printf(OUTFILE "%-80.80s\n",$author);
    foreach my $remark (@remarks) {printf(OUTFILE "%-80.80s\n",$remark);}
    my $match=0; 
    my $len=length($aaq);
    for ($col=0; $col<@i; $col++) {
	if ($i[$col] && $l[$col]) {
	    $match++;
	    for ($natom=0; $natom<scalar(@{$pdbline[$l[$col]]}); $natom++) {
		$pdbline[$l[$col]][$natom]=~/(^......\s*\d+ .....\w{3} $chain).{5}(.*)/;
		if ($v>=4) {printf("%s%4i %s\n",$1,$i[$col],$2);}
		$line=sprintf("%s%4i %s\n",$1,$i[$col],$2);
#		$pdbline[$l[$col]][$natom]=~/(^......\s*\d+ .....\w{3} $chain).{4}(.*)/;
#		if ($v>=4) {printf("%s%4i%s\n",$1,$i[$col],$2);}
#		$line=sprintf("%s%4i%s\n",$1,$i[$col],$2);
		print(OUTFILE $line);
	    }
	}
    }
    $line=~/^ATOM\s+(\d+)......(.........)/;
    printf(OUTFILE "TER   %5i      %s                                                    \n",1+$1,$2);
    print(OUTFILE "END                                                                             \n");
    close (OUTFILE);

    # Warn?
    if ($v>=2 || ($v>=1 && $len-$match>5)) {
	if ($v>=1 && $len-$match>1) {
	    printf("\nWARNING: could not find coordinates for %i query residues:\n",$len-$match);
	} else { printf("\n"); }
	$nameline=~/^(\S+)/;
	printf("%-14.14s $xseq\n","Q $1:");
	printf("%-14.14s $Sstr\n","Identities:");
	printf("%-14.14s $yseq\n","T $pdbid"."_$chain:");
	printf("\n");
	if ($v>=4) {
	    for ($col=0; $col<@l && $col<200; $col++) {
		printf("%3i  %3i  %3i\n",$col,$i[$col],$l[$col]);
	    }
	}
    }

    return 0;
}
# End MakePdbFile()


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
    elsif ($res eq "ASX") {return "D";}
    elsif ($res eq "GLX") {return "E";}
    elsif ($res eq "MSE") {return "M";} # SELENOMETHIONINE 
    elsif ($res eq "SEP") {return "S";} # PHOSPHOSERINE 
    elsif ($res eq "SEC") {return "C";} # SELENOCYSTEINE
    elsif ($res eq "TPO") {return "T";} # PHOSPHOTHREONINE 
    elsif ($res eq "TYS") {return "Y";} # SULFONATED TYROSINE 
    elsif ($res eq "KCX") {return "K";} # LYSINE NZ-CARBOXYLIC ACID
    else                  {return "X";}
}

sub System {
    if ($v>=2) {print($_[0]."\n");} 
    return system($_[0]); 
}
