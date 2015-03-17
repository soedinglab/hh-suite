#! /usr/bin/env perl
# mergeali.pl- Merge multiple alignments in A3M format via a multiple alignment of their seed sequences.
# Usage: mergeali.pl [-i] infile.fas [-o] outfile.a3m [options]

# (C) Johannes Soeding, 2012

#     HHsuite version 3.0.0 (15-03-2015)
#
#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

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

#     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;
use File::Temp "tempfile";
use strict;
$|= 1; # Activate autoflushing on STDOUT

# Default parameters
our $d=1;    # gap opening penatlty for Align.pm
our $e=0.1;  # gap extension penalty for Align.pm
our $g=0.0;  # endgap penalty for Align.pm
my $program="mergeali.pl";
my $v=2;     # 3: DEBUG
my $help="
 mergeali.pl - Merge multiple sequence alignments (MSAs) in a3m format via a FASTA-formatted 
 master MSA of (parts of) their seed sequences. 

 The file names of the 'slave' MSAs must be of the form <seqname>.a3m with <seqname> the name of 
 the corresponding sequence in the master MSA. For example, if the name line in the master MSA 
 is '>F6Y4V2_MONDO Uncharacterized protein...', the file containing the MSA for this sequence 
 must be called F6Y4V2_MONDO.a3m and be found in the working directory (or some other 
 directories specified using '-d <dirs>'. 
 The sequences in the master MSA must be subsequences of the seed sequences in the A3M-formatted 
 slave MSAs. They also must have a match state assignment by their seed (=first) sequence (as 
 produced by hhblits or by reformat.pl with  option -M first).
 DSSP state sequences will be included by default, and a consensus DSSP state sequence is 
 appended at the top. The output file is generated in A3M format. 

 Usage: $program <infile_master.fas> <outfile.a3m> [options]

 Options:
 -all         match states are all columns with a residue in either of the seed sequences (default)
 -first       match states are all columns with a residue in the first of the seed sequences
 -d <dirs>    directories (separated by spaces) where to find the *.a3m alignments (default='.')
 -mark        mark the seed sequences of slave MSAs in the master alignment by inserting a \@ 
              before their name: '>\@name' 
 -diff <int>  use only the <int> most different sequences from each alignment for merging 
              (calls hhfilter) 
 -name <name> write a name line '#name' as first line into the output file
 -full        for the first seed sequence in the master alignment, use the full-length sequence from the alignment 
              (including leading and trailing residues)
 -v <int>     verbose mode (0: no messages; 1: only warnings; 2: verbose)                                                                           
 Example: $program 1mkaA_templ.fas 1mkaA_full.a3m -d '. /home/soeding/pdb_8Jun05 /home/soeding/scop70_1.67'
\n";

# Processing command line options
if (@ARGV<1) {die $help;}

# Variable declarations
my $infile="";
my $outfile="";
my @indirs=(".",""); # array containing all input directories where alignment files can be found 
                 # (default: current directory or absolute path given)
my $indir;       # variable looping through @indirs
my $tmpfile;     # temporary file used for filtering
my $baseout;     # outfile without extension
my $alifile;     # one of the alignment files to be merged
my $nseed=0;     # counts seed sequences in infile
my $seed_in;     # sequence read from input file
my @seed_in;     # residues of seed sequenences in $infile
my @names;       # seed names
my $seed_ali;    # sequence read from alignment files and to be reformatted by inserting the gaps from $seed_in
my $i;           # counts match states of $seed_in from $infile
my $j;           # counts match states of $seed from $alifile  
my $col;         # columns in alignment of $seed_in (from $infile) against $seed (from $alifile)
my ($imin,$imax,$jmin,$jmax);
my @jj;	         # $jj[$i] is residue from $seed_ali aligned with residue $i of $seed_in
my $mark=0;      # do not mark seed sequences by introducing a @ before their name
my $diff="";     # use only the X most different sequences per alignment?
my $nseq=0;      # total number of sequences in mega-alignment
my $all=1;       # 0: match states according to first sequence
my @match=();    # index array for match states of first seed sequence
my $options="";
my (@H, @E, @C, @G, @B, @S, @T); # counts with DSSP states for consensus calculation
my @seqs=();     # sequence records to print into outfile
my $aliname;     # name of alignment specified by -name option
my $full=0;      # don't add leading and trailing residues from seed_ali

# Set options
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}
if ($options=~s/ -i\s+(\S+) //) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}
if ($options=~/ -d\s+(([^-]\S*\s+)+)/) {@indirs=split(/\s+/,$1); $options=~s/ -d\s+(([^-]\S*\s+)+)//g;}
if ($options=~s/ -mark //) {$mark=1;}
if ($options=~s/ -first //) {$all=0;}
if ($options=~s/ -all //) {$all=1;}
if ($options=~s/ -diff\s+(\S+) //) {$diff=$1;}
if ($options=~s/ -name\s+(.*?) -/-/) {$aliname=$1;}
elsif ($options=~s/ -name\s+(.*?)\s*$//) {$aliname=$1;}
if ($options=~s/ -full //) {$full=1;}
if ($options=~s/ -v\s+(\d+) //) {$v=$1;}
elsif ($options=~s/ -v //)      {$v=1;}

# Read infile and outfile 
if (!$infile  && $options=~s/^\s*([^-]\S+)\s*//) {$infile=$1;} 
if (!$outfile && $options=~s/^\s*([^-]\S+)\s*//) {$outfile=$1;} 

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if ($infile eq "")  {die("$help\nError in $program: input file missing: $!\n");}
if ($outfile eq "") {die("$help\nError in $program: output file missing: $!\n");}

# Determine output format (A3M or FASTA) and basename of outfile
if ($outfile=~/(\S*)\.(\S*?)$/) {
    $baseout=$1;
} else {
    $baseout=$outfile;
}

if ($diff ne "") {
    if (!-d "/tmp/$ENV{USER}") {mkdir("/tmp/$ENV{USER}",0777);}
    my $handle;
    ($handle,$tmpfile)=tempfile("XXXXXXXXXX",SUFFIX => ".a3m",DIR=>"/tmp/$ENV{USER}");
    close $handle or die "Error: unable to close $tmpfile\n";
}

$aliname=~s/\s+/ /g;
$aliname=~s/\s+/ /g;
$aliname=~s/\s+$//g;

# Read input file into @names, @seqs
$/=">";  # set input field seperator
open (INFILE,"<$infile") || die ("Error: could not open $infile for reading: $!\n");
$seed_in=<INFILE>;                                  # skip first line
while ($seed_in=~s/(.)>$/$1/) {$seed_in.=<INFILE>;} # skip first line
while ($seed_in=<INFILE>) {
    if ($seed_in eq ">") {next;}
    while ($seed_in=~s/(.)>$/$1/) {$seed_in.=<INFILE>;} # if nameline contains a '>'
    $seed_in=~s/(\S+)[^\n]*//;
    $names[$nseed]=$1;
    $seed_in=~tr/a-zA-Z.-//cd;     # remove all invalid symbols from sequence, including ">" at the end
    $seed_in=~tr/a-z./A-Z-/;       # transform to upper case 
#    $seed_in=~tr/a-z.-/A-Z~~/;     # transform to upper case 
    $seed_in[$nseed]=$seed_in;
    $nseed++;
} 
close INFILE;

# Add leading and trailing residues from first seed sequence?
if ($full) {

    $nseed=0;
    $seed_in=$seed_in[$nseed];
    $seed_in=~tr/.-//d;      # remove all gaps

    # Open alignment file
    foreach $indir (@indirs) {
	$alifile=$indir."/".$names[$nseed].".a3m";
	my $alifilefas=$indir."/".$names[$nseed].".fas";
	if (-e $alifile) {last;}
	if (-e $alifilefas) {
	    &HHPaths::System("perl $hhscripts/reformat.pl -M first $alifilefas $alifile");
	    last;
	}   
	$alifilefas="";
    }
    if ($alifile eq "") {die("Error: could not find $alifile in input directory/-ies.\n");}
    if ($diff ne "") {
	&HHPaths::System("hhfilter -v $v -i $alifile -o $tmpfile -diff $diff -cov 5");
	$alifile=$tmpfile;
    } 
    # Read seed sequence in alifile
    open (ALIFILE,"<$alifile") || die ("Error: could not open $alifile for reading: $!\n");
    $seed_ali=<ALIFILE>;                                   # skip first line
    while ($seed_ali=~s/(.)>$/$1/) {$seed_ali.=<ALIFILE>;} # skip first line
    while ($seed_ali=<ALIFILE>) {
	if ($seed_ali eq ">") {next;}
	while ($seed_ali=~s/(.)>$/$1/) {$seed_ali.=<ALIFILE>;} # if nameline contains a '>'
	$seed_ali=~s/^([^\n]+)//;
	my $nameline=$1;
	if ($nameline=~/^(ss_|aa_|sa_)/) {next;} # skip secondary structure sequences 
	$seed_ali=~s/[^a-zA-Z.-]//g;  # remove all invalid symbols from sequence, including ">" at the end
	if(($seed_ali=~tr/.-/.-/)>0) {die("Error: seed sequence in $alifile must not contain gaps!\n");}
	last;
    }
    close (ALIFILE);

    # Align $seed_in to $seed_ali 
    my $xseq=$seed_in;
    my $yseq=uc($seed_ali);
    my @i;
    my @j;
    my $Sstr;
    my $score;  
    # The aligned characters are returned in $i[$col] and $j[$col]
    $score = &AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr);  
    
   # DEBUG
    if ($v>=3) {
	printf("SEED_IN:   %s\n",$xseq);
	printf("MATCH:     %s\n",$Sstr);
	printf("SEED_ALI:  %s\n",$yseq);
	printf("score=%6.2f   imax=%i imin=%i   jmax=%i jmin=%i\n",$score,$imax,$imin,$jmax,$jmin);
	printf("\n");
	if ($v>=4) {
	    for ($col=0; $col<@j && $col<200; $col++) {
		printf("%3i  %3i %s %3i %s %3i  %3i\n",$col,$i[$col],substr($seed_in,$i[$col],1),$j[$col],substr($seed_ali,$j[$col],1),$j[$col],$jj[$i[$col]]);
	    }
	}
    }

    if ($jmin>1) {
	$seed_in[0] = substr($seed_ali,0,$jmin-1).$seed_in[0];
	my $leftgaps  = ('-' x ($jmin-1));
	for ($nseed=1; $nseed<scalar(@seed_in); $nseed++) {
	    $seed_in[$nseed] = $leftgaps.$seed_in[$nseed];
	}
    }
    
    my $L=length($seed_ali);
    if ($L-$jmax>0) {
	$seed_in[0] .= substr($seed_ali,$jmax,$L-$jmax);
	my $rightgaps = ('-') x ($L-$jmax);
	for ($nseed=1; $nseed<scalar(@seed_in); $nseed++) {
	    $seed_in[$nseed] .= $rightgaps;
	}
    } 
}

#########################################################################################################################
# For each seed sequence from infile: open alignment file and reformat its sequences by inserting gaps read from infile

for ($nseed=0; $nseed<scalar(@seed_in); $nseed++) {
    $seed_in=$seed_in[$nseed];
    my @seed_gaps=split(/[A-Z]/,$seed_in,-1); # $seed_gaps[$i] contains the gaps after the $i'th residue in $seed_in
    $seed_in=~tr/.-//d;      # remove all gaps

    # Open alignment file
    foreach $indir (@indirs) {
	$alifile=$indir."/".$names[$nseed].".a3m";
	my $alifilefas=$indir."/".$names[$nseed].".fas";
	if (-e $alifile) {last;}
	if (-e $alifilefas) {
	    &HHPaths::System("perl $hhscripts/reformat.pl -M first $alifilefas $alifile");
	    last;
	}   
	$alifilefas="";
    }
    if ($alifile eq "") {die("Error: could not find $alifile in input directory/-ies.\n");}
    if ($diff ne "") {
	&HHPaths::System("hhfilter -v $v -i $alifile -o $tmpfile -diff $diff -cov 5");
	$alifile=$tmpfile;
    } 


    open (ALIFILE,"<$alifile") || die ("Error: could not open $alifile for reading: $!\n");

    if ($v>=3) { #DEBUG
	printf("seed_in=$seed_in\n");
	printf("seed_gaps="); 
	for ($i=0; $i<@seed_gaps; $i++) {printf("%i%s ",$i,$seed_gaps[$i]);} 
	print("\n");
	print("Name=".$names[$nseed]."\n");
    }

    # Read seed sequence in alifile
    $seed_ali=<ALIFILE>;                                   # skip first line
    while ($seed_ali=~s/(.)>$/$1/) {$seed_ali.=<ALIFILE>;} # skip first line
    while ($seed_ali=<ALIFILE>) {
	if ($seed_ali eq ">") {next;}
	while ($seed_ali=~s/(.)>$/$1/) {$seed_ali.=<ALIFILE>;} # if nameline contains a '>'
	$seed_ali=~s/^([^\n]+)//;
	my $nameline=$1;
	if ($nameline=~/^(ss_|aa_|sa_)/) {next;} # skip secondary structure sequences 
	$seed_ali=~s/[^a-zA-Z.-]//g;  # remove all invalid symbols from sequence, including ">" at the end
	if(($seed_ali=~tr/.-/.-/)>0) {die("Error: seed sequence in $alifile must not contain gaps!\n");}
	last;
    }
    close (ALIFILE);

    # Align $seed_in to $seed_ali 
    my $xseq=$seed_in;
    my $yseq=uc($seed_ali);
    my @i;
    my @j;
    my $Sstr;
    my $score;  
    # The aligned characters are returned in $i[$col] and $j[$col]
    $score = &AlignNW(\$xseq,\$yseq,\@i,\@j,\$imin,\$imax,\$jmin,\$jmax,\$Sstr);  
    
    
   # DEBUG
    if ($v>=3) {
	printf("SEED_IN:   %s\n",$xseq);
	printf("MATCH:     %s\n",$Sstr);
	printf("SEED_ALI:  %s\n",$yseq);
	printf("score=%6.2f   imax=%i imin=%i   jmax=%i jmin=%i\n",$score,$imax,$imin,$jmax,$jmin);
	printf("\n");
	if ($v>=4) {
	    for ($col=0; $col<@j && $col<200; $col++) {
		printf("%3i  %3i %s %3i %s %3i  %3i\n",$col,$i[$col],substr($seed_in,$i[$col],1),$j[$col],substr($seed_ali,$j[$col],1),$j[$col],$jj[$i[$col]]);
	    }
	}
    }

    # Does alignment look reasonable?
    my $sc=$score/&max(1,&max($imax-$imin,$jmax-$jmin));
    if ($sc<0.5) {
	printf("\nWARNING: alignment for %s has score per column of only %6.2f   Skipping alignment\n",$names[$nseed],$sc,);
	if ($v>=2) {
	    printf("SEED_IN:   %s\n",$xseq);
	    printf("MATCH:     %s\n",$Sstr);
	    printf("SEED_ALI:  %s\n",$yseq);
	    printf("score=%6.2f   imax=%i imin=%i   jmax=%i jmin=%i\n",$score,$imax,$imin,$jmax,$jmin);
	}
	next;
    }
     # Calculate $jj[$i]: residue number from $seed_ali aligned with residue $i of $seed_in (first=1)
    @jj=(); $jj[0]=0;
    for ($col=0; $col<@i; $col++) {
	$i=$i[$col]; 
	$j=$j[$col]; 
	if ($i>0) {
	    if ($j>0) {
		$jj[$i]=$j;
	    } else {
		$jj[$i]=0;
	    }
	}
    }


    # For each sequence in alifile: reformat sequence by inserting gaps read from infile after each match state
    my $first=1; # seed sequences not yet marked?
    my $res;     # sequence from alignment file
    open (ALIFILE,"<$alifile") || die ("Error: could not open $alifile for reading: $!\n");
    $res=<ALIFILE>;                                                  # skip first line 
    while (1) {if($res=~s/(.)>/$1/) {$res.=<ALIFILE>;} else {last;}} # skip first line 
    while ($res=<ALIFILE>) {
	if ($res eq ">") {next;}
	while (1) {if($res=~s/(.)>/$1/) {$res.=<ALIFILE>;} else {last;}} # if nameline contains a '>'
	$res=~s/^([^\n]+)//;
	my $nameline=$1;
	# Skip secondary structure sequences (ss_dssp is needed for consensus determination)
	if ($nameline=~/^(ss_|sa_|aa_)/ && $nameline!~/^ss_dssp/) {next;} 
	$res=~s/[^a-zA-Z.-]//g;  # remove all invalid symbols from sequence, including ">" at the end
	
	# Split next $res sequence from alignment file into 'Match plus inserts' strings
	my @res=split(/(?=[A-Z-][a-z.]*)/,"X".$res); # $res[$j] contains $j'th match state plus inserts
#	printf("alifile: %-30.30s   seq = \n%s\n",$alifile,$res);
	if (@res<$jmax) {
	    printf("Error: sequence %-.20s in %s has fewer match states (%i) than there are residues \n",$nameline,scalar(@res),$alifile);
	    printf("in the seed sequence of this alignment in $infile (%i).\n",$jmax);
	    print("Remember alignment files must have match states exactly in those columns where the seed has a residue.\n");
	    die();
	}

	# THIS IS WHERE THE MERGING REALLY HAPPENS:
	# For each residue in $seed_in: append corresponding Match+inserts string(s) from $res,  
	# plus the gaps between residues in $seed_in
	$res=$seed_gaps[0];
	$j=$jmin;
	for ($i=1; $i<=$#jj; $i++) {
	    if ($jj[$i]>0) { 
		# Add Match+insert strings of $res that are NOT aligned to residues in $seed_in -> inserts
		for (; $j<$jj[$i]; $j++) {
		    $res[$j]=~tr/A-Z-/a-z/d; # transform to inserts
		    $res.=$res[$j];
		}
		$res.=$res[$j++].$seed_gaps[$i]; # add gaps from input alignment found after $i'th residue in $seed_in
	    } else {
		# Residue in $seed_in aligned to gap in $seed_ali => add '-' to represent missing match state
		$res.="-".$seed_gaps[$i];
	    }
	}
	
	
	if ($nameline=~/^ss_dssp/) {
	    $res=~tr/a-z//d;  # remove inserts from dssp sequence
	    @res=unpack("C*",$res);
	    for ($i=0; $i<@res; $i++) {
		if    (! defined $H[$i]) {$H[$i]=$E[$i]=$C[$i]=$G[$i]=$B[$i]=$S[$i]=$T[$i]=0;}
		if    ($res[$i]==72) {$H[$i]++;} # H
		elsif ($res[$i]==69) {$E[$i]++;} # E
		elsif ($res[$i]==67) {$C[$i]++;} # C
		elsif ($res[$i]==71) {$G[$i]++; $H[$i]+=0.35; $C[$i]+=0.4;} # G
		elsif ($res[$i]==66) {$B[$i]++; $E[$i]+=0.35; $C[$i]+=0.4;} # B
		elsif ($res[$i]==83) {$S[$i]++; $C[$i]+=0.4;} # S
		elsif ($res[$i]==84) {$T[$i]++; $C[$i]+=0.4;} # T 
		elsif ($res[$i]==73) {$C[$i]++;} # I
	    }
	    next; # Skip dssp sequence
	    
	} elsif ($first) {
	    # If seed sequences are to be marked, place a @ before the name of the seed sequence
	    if ($mark) {$nameline="\@"."$nameline";} 
	    $first=0;

	    if ($all==0 && $nseed==0) {
		# Determine match state index array
		@res=unpack("C*",$res);
		for ($i=$j=0; $i<@res; $i++) {
		    my $c=$res[$i];
		    if ($c>=65 && $c<=90) {
			$match[$j++]=1;
		    } elsif ($c==45) {
			$match[$j++]=0;
		    }
		}
	    }
	}
	

	# Turn all match states to inserts that don't have a residue in the first seed sequence
	&AssignMatchStates(\$res);

	push(@seqs,">$nameline\n$res\n");
	$nseq++;

	# DEBUG
	if ($v>=3) {  
	    printf(">%-80.80s\n%s\n",$nameline,$res);
	    if ($v>=4) { 
		printf("res=   "); 
		for ($i=$imin; $i<=$imax; $i++) {
		    printf("i=%-2i jj=%-2i  res[jj]=%-4.4s seed_gaps[i]=%s\n",$i,$jj[$i],$res[$jj[$i]],$seed_gaps[$i]);
		}
	    }
	    print("\n");
	} 

    } 
    close ALIFILE;
}
# For each seed sequence from infile: open alignment file and reformat its sequences by inserting gaps read from infile
#########################################################################################################################

# Add DSSP consensus sequence
if (@H>10) {
    my $res="";
    my $c;
    my $max;
    for ($i=0; $i<@H; $i++) {
	$max=0; $c="-";
	if ($C[$i]>$max) {$max=$C[$i]; $c="C";}
	if ($H[$i]>$max) {$max=$H[$i]; $c="H";}
	if ($E[$i]>$max) {$max=$E[$i]; $c="E";}
	if ($G[$i]>$max) {$max=$G[$i]; $c="G";}
	if ($B[$i]>$max) {$max=$B[$i]; $c="B";}
	if ($S[$i]>$max) {$max=$S[$i]; $c="S";}
	if ($T[$i]>$max) {$max=$B[$i]; $c="T";}
	$res.=$c;
    }
    &AssignMatchStates(\$res);
    $res=~tr/a-z.//d;
    unshift(@seqs,">ss_dssp consensus\n$res\n");
}

# Print outfile
open (OUTFILE,">$outfile") || die ("Error: could not open $outfile for writing: $!\n");
if ($aliname) {print(OUTFILE "#$aliname\n");}
foreach my $seq (@seqs) {
    print(OUTFILE "$seq");
}
close(OUTFILE);

if ($diff ne "") {unlink($tmpfile);}

if ($v>=2) {
    print("Created alignment $outfile with $nseq sequences\n");
}
#&HHPaths::System("perl $hhscripts/reformat.pl $baseout.a3m $baseout.fas -v 3");
# Reformat $baseout.a3m such that match columns are defined by first seed sequence?
#&HHPaths::System("perl $hhscripts/reformat.pl -M first $baseout.fas $baseout"."_first.a3m");

exit;

# Turn all match states to inserts that don't have a residue in the first seed sequence
sub AssignMatchStates() 
{
    my $pres=$_[0];
    if ($all==0) {
	my @res=unpack("C*",${$pres});
	for ($i=$j=0; $i<@res; $i++) {
	    my $c=$res[$i];
	    if ($c>=65 && $c<=90) {
		if (!$match[$j++]) {$res[$i]=$c+32;} # make character lower-case
	    } elsif ($c==45) {
		if (!$match[$j++]) {$res[$i]=46;}    # transform '-' to '.'
	    }
	}
	${$pres} = pack("C*",@res);
	${$pres} =~tr/.//d;
    }
}


# Maximum
sub max {
    my $max = shift @_;
    foreach (@_) {
	if ($_>$max) {$max=$_} 
    }
    return $max;
}
# Minimum
sub min {
    my $min = shift @_;
    foreach (@_) {
	if ($_<$min) {$min=$_} 
    }
    return $min;
}

