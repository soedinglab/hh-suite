#! /usr/bin/env perl
# Build astral-style fasta file from scop parseable files and pdb files

my $rootdir;
BEGIN {
   if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} 
  else {$rootdir="/cluster";}
};

use strict;
use lib "/home/soeding/perl";        # for soeding's computer
use lib "/cluster/user/soeding/perl"; # for soeding's computer
use lib $rootdir."/bioprogs/hhpred"; # for toolkit
use lib $rootdir."/perl";            # for soeding
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;
$|=1; # autoflush on


my $help="
 Build ASTRAL-style FASTA file from scop parseable files and pdb SEQRES records
 (You may filter the output file by sequence identity with pdbfilter.pl)

 Usage:   scop2fasta.pl dir.des.scop.txt_<version> <options>
 Example: scop2fasta.pl dir.des.scop.txt_masterfile-2006-12-01-3.out -d /cluster/databases/pdb/all -o scop100_2006-04

 Options: 
   -d <pdbdir>  give directory of pdb files (default=$pdbdir)
   -o <file>    output file (default: scop100_<version>))
\n";

if (@ARGV<2) {print($help); exit;}

# Default parameters for Align.pm
our $d=3;    # gap opening penalty for Align.pm
our $e=0.1;  # gap extension penalty for Align.pm
our $g=0.09; # endgap penalty for Align.pm
our $v=2;    # verbose mode
our $matrix="identity";

my $outfile;
my $scopfile=$ARGV[0];
my %seq_SEQRES;      # $seq_SEQRES{1hz4A} contains the SEQRES sequence of that pdb id and chain
my %seq_ATOM;        # $seq_ATOM{1hz4A} contains the ATOM sequence of that pdb id and chain
my %idx_ATOM;        # $idx_ATOM{1hz4A} contains the array with ATOM index numbers for that pdb id and chain
my %ins_ATOM;        # $ins_ATOM{1hz4A} contains the array with insertion codes for that pdb id and chain
my %read;            # $read{1hz4} = defined if pdb record 1hz4 has already been read
my %pdbc2i;          # $pdbc2i{$pdbc} = aligned characters returned by AlignNW()
my %pdbc2j;          # $pdbc2j{$pdbc} = aligned characters returned by AlignNW()
my %resolution;      # $resolution{$pdb} = resolution in Angstrom
my %seq2resolution;  # $seq2resolution{scop_sequence} =  resolution of the corresponding pdb structure in A
my %seq2pdbc;        # $seq2pdbc{scop_sequence} =  pdbc code for that domain 
my %seq2i;           # $seq2i{scop_sequence} = index of @seqs where scop_sequence record is stored
my @seqs;            # $seqs[index] = ">d1gcja_ a.118.1.1 (A:) Importin beta {Mouse (Mus musculus)}\npdyas...."
my @chains=();       # chain A in residue range A:372-450
my @first=();        # first residue 372 in residue range A:372-450
my @last=();         # last residue 450 in residue range A:372-450
my $line;            # dummy for reading in
my %three2one=(
	       "ALA"=>"A","VAL"=>"V","PHE"=>"F","PRO"=>"P","MET"=>"M","ILE"=>"I","LEU"=>"L","ASP"=>"D","GLU"=>"E","LYS"=>"K",
	       "ARG"=>"R","SER"=>"S","THR"=>"T","TYR"=>"Y","HIS"=>"H","CYS"=>"C","ASN"=>"N","GLN"=>"Q","TRP"=>"W","GLY"=>"G",
	       "2AS"=>"D","3AH"=>"H","5HP"=>"E","ACL"=>"R","AIB"=>"A","ALM"=>"A","ALO"=>"T","ALY"=>"K","ARM"=>"R","ASA"=>"D",
	       "ASB"=>"D","ASK"=>"D","ASL"=>"D","ASQ"=>"D","AYA"=>"A","BCS"=>"C","BHD"=>"D","BMT"=>"T","BNN"=>"A","BUC"=>"C",
	       "BUG"=>"L","C5C"=>"C","C6C"=>"C","CCS"=>"C","CEA"=>"C","CHG"=>"A","CLE"=>"L","CME"=>"C","CSD"=>"A","CSO"=>"C",
	       "CSP"=>"C","CSS"=>"C","CSW"=>"C","CXM"=>"M","CY1"=>"C","CY3"=>"C","CYG"=>"C","CYM"=>"C","CYQ"=>"C","DAH"=>"F",
	       "DAL"=>"A","DAR"=>"R","DAS"=>"D","DCY"=>"C","DGL"=>"E","DGN"=>"Q","DHA"=>"A","DHI"=>"H","DIL"=>"I","DIV"=>"V",
	       "DLE"=>"L","DLY"=>"K","DNP"=>"A","DPN"=>"F","DPR"=>"P","DSN"=>"S","DSP"=>"D","DTH"=>"T","DTR"=>"W","DTY"=>"Y",
	       "DVA"=>"V","EFC"=>"C","FLA"=>"A","FME"=>"M","GGL"=>"E","GLZ"=>"G","GMA"=>"E","GSC"=>"G","HAC"=>"A","HAR"=>"R",
	       "HIC"=>"H","HIP"=>"H","HMR"=>"R","HPQ"=>"F","HTR"=>"W","HYP"=>"P","IIL"=>"I","IYR"=>"Y","KCX"=>"K","LLP"=>"K",
	       "LLY"=>"K","LTR"=>"W","LYM"=>"K","LYZ"=>"K","MAA"=>"A","MEN"=>"N","MHS"=>"H","MIS"=>"S","MLE"=>"L","MPQ"=>"G",
	       "MSA"=>"G","MSE"=>"M","MVA"=>"V","NEM"=>"H","NEP"=>"H","NLE"=>"L","NLN"=>"L","NLP"=>"L","NMC"=>"G","OAS"=>"S",
	       "OCS"=>"C","OMT"=>"M","PAQ"=>"Y","PCA"=>"E","PEC"=>"C","PHI"=>"F","PHL"=>"F","PR3"=>"C","PRR"=>"A","PTR"=>"Y",
	       "SAC"=>"S","SAR"=>"G","SCH"=>"C","SCS"=>"C","SCY"=>"C","SEL"=>"S","SEP"=>"S","SET"=>"S","SHC"=>"C","SHR"=>"K",
	       "SOC"=>"C","STY"=>"Y","SVA"=>"S","TIH"=>"A","TPL"=>"W","TPO"=>"T","TPQ"=>"A","TRG"=>"K","TRO"=>"W","TYB"=>"Y",
	       "TYQ"=>"Y","TYS"=>"Y","TYY"=>"Y","AGM"=>"R","GL3"=>"G","SMC"=>"C","ASX"=>"B","CGU"=>"E","CSX"=>"C","GLX"=>"Z",
	       "LED"=>"L"	       );


my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

# General options
if ($options=~s/ -d\s*(\S+) / /) {$pdbdir=$1;}
if ($options=~s/ -o\s*(\S+) //)  {$outfile=$1;}
if ($options=~s/ -v\s*(\S+) //)  {$v=$1;}

if (! $outfile) {
    $scopfile=~/dir\....\.scop\.txt_(.*)/;
    $outfile="scop100_$1";
    $outfile=~s/masterfile-//;
    $outfile=~s/\.out//;
    printf("Output will be written to %s in the end\n",$outfile);
}

# Read descriptions from dir.des.scop.txt_...
my %desc;   # hash containing all the descriptions, accessed by node-ID
my %cfsf;   # $cfsf{46458}="a.1.1" Class Fold Superfamily Family-code of node-ID 
my %count;  # $count{46458} gives the index of the last NEW fold, superfamily, or family assigned here (e.g. 99) 
open(DES,"<$scopfile") || die("Error: can not open $scopfile for reading: $!\n");
while ($line=<DES>) {
    if ($line=~/^\#/) {next;}
    $line=~/(\d+)\t(\S+)\t(\S+)\t\S+\t(.*)/;
    $desc{$1}=$4;
    if ($3 ne "unassigned-sccs") {$cfsf{$1}=$3;}
}
close(DES);


# Read protein records with ids for protein and species descriptions from dir.cla.scop.txt_...
$scopfile=~s/dir\.des\.scop/dir\.cla\.scop/;
my %pdbc2dom; # $pdb2dom{1hz4a}=1 for first domain of 1hz4A
# d6prch1	6prc	H:37-258	b.41.1.1	25457	cl=48724,cf=50345,sf=50346,fa=50347,dm=50348,sp=50349,px=25457
open(CLA,"<$scopfile") || die("Error: can not open $scopfile for reading: $!\n");
while ($line=<CLA>) {
    if ($line=~/^\#/) {next;}
    $line=~/(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\tcl=(\d+),cf=(\d+),sf=(\d+),fa=(\d+),dm=(\d+),sp=(\d+)/;
    my $sid=$1;
    my $pdb=$2;
    my $region=$3;
    my $fam=$4;
    my $cl=$5;
    my $fd=$6;
    my $sf=$7;
    my $fa=$8;
    my $dm=$9;
    my $sp=$10;
    my $chain; # chain symbol in scop-id (sid)
    my $dom="_";
    @chains=();  # chain A in residue range A:372-450
    @first=();   # first residue 372 in residue range A:372-450
    @last=();    # last residue 450 in residue range A:372-450

    # Read chains / residue ranges
    # -  223-314  A:118-211,A:372-450   A:,B:16-130
    my $region1=$region;
    while ($region1=~/-|:|\d+/) {
#	printf("-> '%s' ",$region1);
	if ($region1=~s/^-$//) {                          # -
	    push(@chains," "); push(@first,-1e5); push(@last,1e5);
	    $chain="A"; last;
	} elsif ($region1=~s/^(-?\d+)-(\d+)//) {          # 223-314
	    push(@chains," "); push(@first,$1);   push(@last,$2);
	    $chain="A";
	} elsif ($region1=~s/^(-?\d+\w)-(\d+\w)//) {      # 223S-314S
	    push(@chains," "); push(@first,$1);   push(@last,$2);
	    $chain="A";
	} elsif ($region1=~s/^(\w):(-?\d+)-(\d+)//) {     # A:223-314
	    push(@chains,$1);  push(@first,$2);   push(@last,$3); 
	    if (! defined $chain) {$chain=$1;} elsif ($chain ne $1) {$chain=".";} 
	} elsif ($region1=~s/^(\w):(-?\d+\w)-(\d+\w)//) { # A:1S-104S
	    push(@chains,$1);  push(@first,$2);   push(@last,$3);
	    if (! defined $chain) {$chain=$1;} elsif ($chain ne $1) {$chain=".";} 
	} elsif ($region1=~s/^(\S)://) {                  # A:
	    push(@chains,$1);  push(@first,-1e5); push(@last,1e5);
	    if (! defined $chain) {$chain=$1;} elsif ($chain ne $1) {$chain=".";} 
	}
	$region1=~s/^,//; # remove leading colon
	$region1=~s/^\s//; # remove leading colon
    }     my $pdbc=$pdb.lc($chain);
    
    # SCOP-ID (d1hz4a_) unassigned?
    if ($sid=~/^unassigned-sid/) {
	# Assign domain number (for unassigned sid's)
	if ($first[0]>-1000) {
	    if (defined $pdbc2dom{$pdbc}) {
		$dom = ++$pdbc2dom{$pdbc};
		if ($dom>9) { $dom=chr($dom+87); } # 10'th domain will have a, 11'th b, and so on
	    } else {
		$dom = $pdbc2dom{$pdbc} = 1;
	    }
	}
	$sid="e".$pdbc.$dom;
    } else {
	if ($sid=~/......(\d)/) {$pdbc2dom{$pdbc}=$1;}
    }

    # Family code unassigned? => assign preliminary code (e.g. a.999.98.99, or c.12.1.99)
    if ($fam=~/unassigned/) {

	if (defined $cfsf{$fa}) {
	    $fam=$cfsf{$fa};

	} elsif (defined $cfsf{$sf}) {
	    $fam=$cfsf{$sf};
	    if (defined $count{$sf}) {$count{$sf}--;} else {$count{$sf}=99;} 
	    $fam.=".".$count{$sf};
	    $cfsf{$fa}=$fam;

	} elsif (defined $cfsf{$fd}) {
	    $fam=$cfsf{$fd};
	    if (defined $count{$fd}) {$count{$fd}--;} else {$count{$fd}=99;} 
	    $fam.=".".$count{$fd};
	    $cfsf{$sf}=$fam;
	    $count{$sf}=99;
	    $fam.=".99";
	    $cfsf{$fa}=$fam;

	} elsif (defined $cfsf{$cl}) {
	    $fam=$cfsf{$cl};
	    if (defined $count{cl}) {$count{$cl}--;} else {$count{$cl}=999;} 
	    $fam.=".".$count{$cl};
	    $cfsf{$fd}=$fam;
	    $count{$fd}=99;
	    $fam.=".99";
	    $cfsf{$sf}=$fam;
	    $count{$sf}=99;
	    $fam.=".99";
	    $cfsf{$fa}=$fam;

	} else {
	    $fam="unassigned";
	}
    }

    # Print SCOP sequence
    my $nameline = sprintf(">%s %s (%s) %s {%s}",$sid,$fam,$region,$desc{$dm},$desc{$sp});
    &ExtractSCOPSeqFromPDBFile($pdb,$nameline);
}
close(CLA);
 
# Print SCOP sequences
open(OUT,">$outfile") || die("Error: can not open $outfile for writing: $!\n");
foreach my $seq (@seqs) {
    print(OUT $seq);
}
close(OUT);
exit;


######################################################################################################
# Extract ASTRAL-type SEQRES sequence for pdb code $pdb and list of chains and residues given in
# @chains, @first, @last. 
# Sequences are in lower case and segments are linked with a "X character.
# The residue numbering in @first and @last refers to the numbering in the pdb ATOM records. 
# Hence an alignment of the ATOM to the SEQRES sequence is necessary if not entire chain extracted. 
######################################################################################################
sub ExtractSCOPSeqFromPDBFile()
{
    my $pdb = $_[0]; # pdb code 
    my $nameline = $_[1];
    my $pdbfile = "$pdbdir/pdb$pdb.ent";
    my $seq;        # SCOP sequence to be returned
    my $chain;      # e.g. 1hz4
    my $pdbc;       # e.g. 1hz4A
    my $nres=-100;  # residue number
    my $j;          # 

    # Has PDB file not yet been read? => read it 
    if (! defined $read{$pdb}) {
	# Open the pdb file
	if (!open (PDBFILE,"<$pdbfile")) {printf("Error: can't open $pdbfile: $!\n\n"); return ;}
	$read{$pdb}=1;
	
	# Read resolution
	$resolution{$pdb}=99; # if not found, set to 100 for NMR
	# REMARK   2 RESOLUTION. 2.00 ANGSTROMS.                                          
	$line=<PDBFILE>;
	while ($line && $line!~/^REMARK   2 /o && $line!~/^SEQRES /o) {$line=<PDBFILE>;}			     
	if (!$line) {if ($v>=2) {print("Error: wrong format in $pdbfile. Skipping file ...\n");} next;} 
	if ($v>=2 && $line!~/^REMARK   2 /) {print("\n\nWarning: no REMARK   2 line found in $pdbfile\n");}
	while ($line=~/^REMARK   2 /o && $line!~/^SEQRES /o) {
	    if ($line=~/^REMARK   2\s+RESOLUTION\.\s+(\d+\.?\d*)/) {$resolution{$pdb}=$1; last;}
	    $line=<PDBFILE>;
	}

	# Read SEQRES records
	my $residues;   # residues in SEQRES record
	while ($line && $line!~/^SEQRES/) {$line=<PDBFILE>; }
	while ($line && $line=~/^SEQRES/) {
	    # SEQRES   1 A  478  VAL ARG ILE ALA LEU LYS LYS ARG PRO ILE ASP ARG ASN  
	    $line=~/^SEQRES\s+\d+\s(.)\s+\d+\s+[+-]?(\w.*)/;
	    $chain=$1;
	    if ($chain eq " ") {$chain="A";}
	    $residues=$2;
	    $pdbc=$pdb.$chain;
	    if (! defined $seq_SEQRES{$pdbc}) {$seq_SEQRES{$pdbc}="";}
	    while ($residues=~s/^(\w\w\w)\s+//) {$seq_SEQRES{$pdbc}.=&Three2OneLetter($1);}
	    $line=<PDBFILE>;
	};
	
	
	# Read ATOM records
	while ($line=<PDBFILE>) {
	    # ATOM      1  N   GLY A   1     -19.559   8.872   4.925  1.00 16.44           N
	    # ATOM      2  CA  GLY A   1     -19.004   8.179   6.112  1.00 14.30           C
	    if ($line=~/^ENDMDL/) {last;} # if file contains NMR models read only first
	    
	    if ( ($line=~/^ATOM  \s*\d+ ....[ A](\w{3}) (.)\s*(-?\d+)(.)/ ) || 
		 ($line=~/^HETATM\s*\d+ ....[ A](\w{3}) (.)\s*(-?\d+)(.)/ && defined $three2one{$1} ) ) {
		
		if ($nres==$3) {next;} 

		# Append residue to chain and record ATOM index of residue
		$nres=$3;
		$chain=$2;
		if ($chain eq " ") {$chain="A";}
		$pdbc=$pdb.$chain;
		if (! defined $seq_ATOM{$pdbc}) {
		    $seq_ATOM{$pdbc}="";
		    my @dummy1=(); 
		    my @dummy2=(); 
		    $idx_ATOM{$pdbc}=\@dummy1;
		    $ins_ATOM{$pdbc}=\@dummy2;
		}
		$seq_ATOM{$pdbc} .= &Three2OneLetter($1);
		push(@{$idx_ATOM{$pdbc}},$nres); 
		push(@{$ins_ATOM{$pdbc}},$4); 
	    }
	}
	close (PDBFILE);
    }    

    # Extract sequence segments and append to $seq
    while (scalar(@chains)>=1) {

	my $chain = shift(@chains);
	my $first = shift(@first);
	my $last  = shift(@last);
	my $inscode;
	my ($jfirst,$jlast);
	my ($ifirst,$ilast);
	my ($imin,$imax,$jmin,$jmax,$Sstr);
	
	# Extract only a region of chain
	my $col;     # column of alignment query (from hhm file) versus pdb-residues
	my ($iref,$jref);
	$pdbc = $pdb.$chain;
	if (! defined $seq_SEQRES{$pdbc}) {$pdbc = $pdb.lc($chain);}
	if (! defined $seq_SEQRES{$pdbc}) {printf("Error: could not find chain $chain in $pdbfile in \n\n"); next;}
	
	# Extract region, not entire chain?
	if ($first=~/^\d+[A-Z]/ || $first>-1000) {
   
	    # Find residue number $first in $idx_ATOM{$pdbc}
	    my $idx=\@{$idx_ATOM{$pdbc}};
	    if ($first=~s/([A-Z])$//) {$inscode=$1;} else {$inscode=" ";}
	    for ($j=0; $j<scalar(@{$idx}); $j++) {
		if (${$idx}[$j]==$first) {
		    if ($inscode eq ${$ins_ATOM{$pdbc}}[$j]) {
			$jfirst=$j+1; 
			last;
		    }
		}
	    }
	    if (! defined $jfirst) {printf("Error: can not find residue number $first in $pdbfile\n\n"); next;}

	    # Find residue number $last in $idx_ATOM{$pdbc}
	    if ($last=~s/([A-Z])$//) {$inscode=$1;} else {$inscode=" ";}
	    for ($j=$jfirst+1; $j<scalar(@{$idx}); $j++) {
		if (${$idx}[$j]==$last) {
		    if ($inscode eq ${$ins_ATOM{$pdbc}}[$j]) {
			$jlast=$j+1; 
			last;
		    }
		}
	    }
	    if (! defined $jlast) {printf("Error: can not find residue number $last in $pdbfile\n\n"); next;}
	    if ($jfirst>$jlast) {printf("Error: resdidue $first appears before $last ($jfirst, $jlast) in $pdbfile\n\n"); next;}
	    if ($v>=3) {printf("first PDB:%-3i -> %-3i  last PDB:%-3i -> %-3i\n",$first,$jfirst,$last,$jlast);}

	}	    
	
	# Alignment already done?
	if (! defined $pdbc2i{$pdbc}) {
	    
	    # Must align $seq_SEQRES{$pdbc} to $seq_ATOM{$pdbc}	    
	    my $xseq=$seq_SEQRES{$pdbc};
	    my $yseq=$seq_ATOM{$pdbc};
	    if (!defined $xseq) {print("\nWARNING: $pdbc does not contain any SEQRES residues\n"); $xseq=$yseq;}
	    if (!defined $yseq) {print("\nWARNING: $pdbc does not contain any ATOM residues. Skipping chain\n"); next;}
	    my $len=length($yseq);
	    my (@i,@j);  # The aligned characters are returend in $i[$col] and $j[$col]
	    my $score=&AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr);  
	    
	    if ($v>=3 || ($v>=1 && ($jmin>1 || $jmax-$jmin+1<$len))) {
		if ($v>=1 && ($jmin>1 || $jmax-$jmin+1<$len)) {
		    printf("\nWARNING: $pdbc ATOM sequence not fully covered by SEQRES sequence: imin=%-i   imax=%-i   jmin=%-i:  jmax=%-i  len=%i\n",$imin,$imax,$jmin,$jmax,length($yseq));
		} else { printf("\n"); }
		printf("%-14.14s $xseq\n","Q $pdbc"."_SEQRES");
		printf("%-14.14s $Sstr\n","Identities:");
		printf("%-14.14s $yseq\n","T $pdbc"."_ATOM");
		printf("\n");
		if ($v>=4) {
		    for ($col=0; $col<@j && $col<200; $col++) {
			printf("%3i  %3i  %3i\n",$col,$i[$col],$j[$col]);
		    }
		}
	    }
	    $iref = \@i;
	    $jref = \@j;
	    if ($first=~/^\d+[A-Z]/ || $first>-1000) { # do we need to store this alignment?
		$pdbc2i{$pdbc} = \@i;
		$pdbc2j{$pdbc} = \@j;
	    }
	} else {
	    $iref = $pdbc2i{$pdbc};
	    $jref = $pdbc2j{$pdbc};
	}
	
	# Extract region, not entire chain?
	if ($first=~/^\d+[A-Z]/ || $first>-1000) {

	    # Find position in SEQRES sequence corresponding to $jfirst 
	    for ($col=0; $col<@{$jref}; $col++) {
		if (${$jref}[$col]==$jfirst) {$ifirst=${$iref}[$col]; last}
	    }
	    if (! defined $ifirst) {printf("Error: Could not find residue jfirst=$jfirst in \@j=@{$jref}\n\n"); next;}
	    
	    # Find position in SEQRES sequence corresponding to $jlast 
	    for (; $col<@{$jref}; $col++) {
		if (${$jref}[$col]==$jlast) {$ilast=${$iref}[$col]; last}
	    }
	    if (! defined $ilast) {printf("Error: Could not find residue jlast=$jlast in \@j=@{$jref}\n\n"); next;}
	    if ($v>=3) {printf("jfirst:%-3i -> ifirst:%-3i  jlast:%-3i -> ilast:%-3i\n",$jfirst,$ifirst,$jlast,$ilast);}

	} else {

	    $ifirst = $imin;
	    $ilast = $imax;	    
	}	

	$seq .= "X".lc( substr($seq_SEQRES{$pdbc},$ifirst-1,$ilast-$ifirst+1) );
	
	if ($first!~/^\d+[A-Z]/ && $first<=-1000) {@{$iref}=(); @{$jref}=();}
	
    }
    
    if (defined $seq && $seq ne "") {
	
	# Remove leading X and format to 60 residues per line
	$seq=~s/^X//;
	$seq=~s/(.{60})/$1\n/g;
	$seq=~s/\n$//;
	
	# SCOP sequence already defined?
	if (defined $seq2resolution{$seq}) {
	    my $res_prev = $seq2resolution{$seq};
	    if ($res_prev<=$resolution{$pdb}) {return;} # prev domain has same or better resolution? => skip
	    $seqs[$seq2i{$seq}] = sprintf("%s %-.2f\n%s\n",$nameline,$resolution{$pdb},$seq);
	    printf("***** Replacing %s (%.2fA) with  %s  (%.2fA) *****\n",$seq2pdbc{$seq},$res_prev,$pdbc,$resolution{$pdb});
	    
	} else {
	    # SCOP sequence not yet defined
	    $seq2i{$seq} = scalar(@seqs);
	    push (@seqs, sprintf("%s %-.2f\n%s\n",$nameline,$resolution{$pdb},$seq) );
#	    push (@seqs, sprintf("%s\n%s\n",$nameline,$seq) );
	}
	$seq2resolution{$seq} = $resolution{$pdb};
	$seq2pdbc{$seq} = $pdbc;	
	
	if ($v==2) { printf("%s %-.2f\n",$nameline,$resolution{$pdb}); }
	if ($v>=3) { printf("%s %-.2f\n%s\n",$nameline,$resolution{$pdb},$seq); }
	
    }

    return;
}


##################################################################################
# Convert three-letter amino acid code into one-letter code
##################################################################################
sub Three2OneLetter {
    my $three=uc($_[0]);
    if (defined $three2one{$three}) {
	return $three2one{$three};
    } else {
	return "X";
    }
}

