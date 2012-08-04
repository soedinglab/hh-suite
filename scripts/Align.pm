# Package Align.pl
# (c) Johannes Soeding, 2006
# Perl functions for Smith-Waterman and Needleman-Wunsch sequence alignment

#     HHsuite version 2.0
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
[-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,+1.0,-9.9], # X	  
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
[ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,+1,-9],
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
	    $Sab[20][20] = $Sab[20][20] = +1;# if in doubt, match X with X
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

    my $dx = $main::dx;
    my $dy = $main::dy;
    if (! defined $dx) {$dx = $main::d;}
    if (! defined $dy) {$dy = $main::d;}
    
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
	    ${$Ai}[$j] = max2bt(${$Ai}[$j1]-$main::e, ${$Mi}[$j1]-$dx, \$Abt[$i][$j]);
	    ${$Bi}[$j] = max2bt(${$Bi1}[$j]-$main::e, ${$Mi1}[$j]-$dy, \$Bbt[$i][$j]);
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
		unshift(@$rS,-$dx);
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
		unshift(@$rS,-$dy);
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

    my $dx = $main::dx;
    my $dy = $main::dy;
    if (! defined $dx) {$dx = $main::d;}
    if (! defined $dy) {$dy = $main::d;}
    printf("dx=%f  dy=%f\n",$dx,$dy); ##############DEBUG#############
    
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
	    ${$Ai}[$j] = max2bt(${$Ai}[$j1]-$main::e, ${$Mi}[$j1]-$dx, \$Abt[$i][$j]);
	    ${$Bi}[$j] = max2bt(${$Bi1}[$j]-$main::e, ${$Mi1}[$j]-$dy, \$Bbt[$i][$j]);
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
		    unshift(@$rS,-$dx); # gap opening
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
		    unshift(@$rS,-$dy); # gap opening
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
