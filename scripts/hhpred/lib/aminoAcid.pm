package aminoAcid;

require Exporter;

use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT = qw(@SIM);

my @GONNET = (
##  A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V 
 10227, 3430, 2875, 3869, 1625, 2393, 4590, 6500, 2352, 3225, 5819, 4172, 1435, 1579, 3728, 4610, 6264,  418, 1824, 5709, # A
  3430, 7780, 2209, 2589,  584, 2369, 3368, 3080, 2173, 1493, 3093, 5701,  763,  859, 1893, 2287, 3487,  444, 1338, 2356, # R 
  2875, 2209, 3868, 3601,  501, 1541, 2956, 3325, 1951, 1065, 2012, 2879,  532,  688, 1480, 2304, 3204,  219, 1148, 1759, # N
  3869, 2589, 3601, 8618,  488, 2172, 6021, 4176, 2184, 1139, 2151, 3616,  595,  670, 2086, 2828, 3843,  204, 1119, 2015, # D
  1625,  584,  501,  488, 5034,  355,  566,  900,  516,  741, 1336,  591,  337,  549,  419,  901, 1197,  187,  664, 1373, # C
  2393, 2369, 1541, 2172,  355, 1987, 2891, 1959, 1587, 1066, 2260, 2751,  570,  628, 1415, 1595, 2323,  219,  871, 1682, # Q
  4590, 3368, 2956, 6021,  566, 2891, 8201, 3758, 2418, 1624, 3140, 4704,  830,  852, 2418, 2923, 4159,  278, 1268, 2809, # E
  6500, 3080, 3325, 4176,  900, 1959, 3758,26066, 2016, 1354, 2741, 3496,  741,  797, 2369, 3863, 4169,  375, 1186, 2569, # G
  2352, 2173, 1951, 2184,  516, 1587, 2418, 2016, 5409, 1123, 2380, 2524,  600, 1259, 1298, 1642, 2446,  383,  876, 1691, # H
  3225, 1493, 1065, 1139,  741, 1066, 1624, 1354, 1123, 6417, 9630, 1858, 1975, 2225, 1260, 1558, 3131,  417, 1697, 7504, # I
  5819, 3093, 2012, 2151, 1336, 2260, 3140, 2741, 2380, 9630,25113, 3677, 4187, 5540, 2670, 2876, 5272, 1063, 3945,11005, # L
  4172, 5701, 2879, 3616,  591, 2751, 4704, 3496, 2524, 1858, 3677, 7430,  949,  975, 2355, 2847, 4340,  333, 1451, 2932, # K
  1435,  763,  532,  595,  337,  570,  830,  741,  600, 1975, 4187,  949, 1300, 1111,  573,  743, 1361,  218,  828, 2310, # M
  1579,  859,  688,  670,  549,  628,  852,  797, 1259, 2225, 5540,  975, 1111, 6126,  661,  856, 1498, 1000, 4464, 2602, # F
  3728, 1893, 1480, 2086,  419, 1415, 2418, 2369, 1298, 1260, 2670, 2355,  573,  661,11834, 2320, 3300,  179,  876, 2179, # P
  4610, 2287, 2304, 2828,  901, 1595, 2923, 3863, 1642, 1558, 2876, 2847,  743,  856, 2320, 3611, 4686,  272, 1188, 2695, # S
  6264, 3487, 3204, 3843, 1197, 2323, 4159, 4169, 2446, 3131, 5272, 4340, 1361, 1498, 3300, 4686, 8995,  397, 1812, 5172, # T
   418,  444,  219,  204,  187,  219,  278,  375,  383,  417, 1063,  333,  218, 1000,  179,  272,  397, 4101, 1266,  499, # W
  1824, 1338, 1148, 1119,  664,  871, 1268, 1186,  876, 1697, 3945, 1451,  828, 4464,  876, 1188, 1812, 1266, 9380, 2227, # Y
  5709, 2356, 1759, 2015, 1373, 1682, 2809, 2569, 1691, 7504,11005, 2932, 2310, 2602, 2179, 2695, 5172,  499, 2227,11569);# V


my @BLOSUM62 = (
  0.0222,
  0.0022,0.0181,
  0.0019,0.0019,0.0148,
  0.0021,0.0015,0.0037,0.0225,
  0.0016,0.0004,0.0004,0.0004,0.0127,
  0.0018,0.0024,0.0015,0.0016,0.0003,0.0076,
  0.0029,0.0025,0.0021,0.0049,0.0003,0.0034,0.0168,
  0.0057,0.0016,0.0027,0.0024,0.0007,0.0013,0.0018,0.0396,
  0.0010,0.0013,0.0014,0.0009,0.0002,0.0010,0.0013,0.0009,0.0096,
  0.0031,0.0012,0.0009,0.0011,0.0012,0.0008,0.0012,0.0013,0.0006,0.0191,
  0.0043,0.0023,0.0013,0.0014,0.0016,0.0016,0.0019,0.0020,0.0009,0.0115,0.0388,
  0.0032,0.0062,0.0024,0.0024,0.0005,0.0030,0.0040,0.0024,0.0012,0.0015,0.0024,0.0166,
  0.0013,0.0007,0.0005,0.0004,0.0004,0.0007,0.0006,0.0007,0.0003,0.0025,0.0052,0.0008,0.0045,
  0.0016,0.0009,0.0007,0.0007,0.0005,0.0005,0.0008,0.0011,0.0008,0.0030,0.0056,0.0009,0.0012,0.0186,
  0.0021,0.0009,0.0008,0.0011,0.0003,0.0008,0.0014,0.0013,0.0005,0.0009,0.0013,0.0015,0.0004,0.0005,0.0195,
  0.0065,0.0022,0.0030,0.0027,0.0011,0.0018,0.0029,0.0037,0.0011,0.0017,0.0024,0.0030,0.0008,0.0012,0.0016,0.0137,
  0.0037,0.0017,0.0022,0.0019,0.0009,0.0013,0.0021,0.0022,0.0007,0.0027,0.0033,0.0023,0.0010,0.0012,0.0013,0.0048,0.0133,
  0.0004,0.0003,0.0002,0.0001,0.0002,0.0002,0.0002,0.0004,0.0002,0.0004,0.0007,0.0003,0.0002,0.0009,0.0001,0.0003,0.0003,0.0074,
  0.0013,0.0009,0.0007,0.0006,0.0004,0.0006,0.0008,0.0008,0.0016,0.0015,0.0023,0.0010,0.0006,0.0043,0.0004,0.0010,0.0009,0.0010,0.0113,
    0.0049,0.0015,0.0011,0.0012,0.0014,0.0011,0.0016,0.0017,0.0006,0.0120,0.0094,0.0019,0.0023,0.0025,0.0012,0.0023,0.0035,0.0004,0.0015,0.0206);

use vars qw(@SIM);
@SIM = ();

sub new {
    my $class = shift;

    my $self = {};

    bless($self, $class);

    $self->setupMatrices();

    return $self;
}


sub setupMatrices {
    my $self = shift;

    $self->setupSubstitutionMatrix();
}


sub aa2i {
    my $self = shift;
    my $aa = shift;

    $aa = uc($aa);

    if ($aa eq 'A') {    return 0; } 
    elsif ($aa eq 'R') { return 1; } 
    elsif ($aa eq 'N') { return 2; } 
    elsif ($aa eq 'D') { return 3; } 
    elsif ($aa eq 'C') { return 4; } 
    elsif ($aa eq 'Q') { return 5; } 
    elsif ($aa eq 'E') { return 6; } 
    elsif ($aa eq 'G') { return 7; }
    elsif ($aa eq 'H') { return 8; } 
    elsif ($aa eq 'I') { return 9; } 
    elsif ($aa eq 'L') { return 10;}
    elsif ($aa eq 'K') { return 11;} 
    elsif ($aa eq 'M') { return 12;}
    elsif ($aa eq 'F') { return 13;} 
    elsif ($aa eq 'P') { return 14;}
    elsif ($aa eq 'S') { return 15;}
    elsif ($aa eq 'T') { return 16;} 
    elsif ($aa eq 'W') { return 17;} 
    elsif ($aa eq 'Y') { return 18;} 
    elsif ($aa eq 'V') { return 19;} 
    elsif ($aa eq 'X') { return 20;}
    elsif ($aa eq 'J') { return 20;}
    elsif ($aa eq 'O') { return 20;}
    elsif ($aa eq 'U') { return 4; }
    elsif ($aa eq 'B') { return 3; }
    elsif ($aa eq 'Z') { return 6; }
    elsif ($aa eq '-') { return 20;}
    elsif ($aa eq '_') { return 20;}
    elsif ($aa eq '.') { return 20;}

    else {return 20;}
}


sub printSimilarityMatrix {
    my $self = shift;

    for (my $i=0; $i<20; $i++) {
	    for (my $j=0; $j<20; $j++) {
	        my $me = sprintf("%.3f", $SIM[$i][$j]);
	        print "$me ";
	    }
	    print "\n";
    }
}


sub setupSubstitutionMatrix {
    my $self = shift;
    my $m = shift;
    my $matrix = defined($m) ? $m : 0;

    my @pb;
    my @P;
    my @R;

    my $b;
    if ($matrix == 0) {
	    for (my $a=0; $a<20; ++$a) {
	        for ($pb[$a]=0.0, $b=0; $b<20; ++$b) {
		        $P[$a][$b] = 0.000001 * $GONNET[$a*20 + $b];
	        }
	    }
	    for (my $a=0; $a<20; ++$a) { 
            $P[$a][20] = $P[20][$a] = 1.0 
        }
    } else {
	    &setupBlosumMatrix($matrix, \@P, \@pb);
    }
    
    my $sumab=0.0;
    for (my $a=0; $a<20; $a++) {
	    for (my $b=0; $b<20; $b++) {
	        $sumab += $P[$a][$b];
	    }
    }
    for (my $a=0; $a<20; $a++) {
	    for (my $b=0; $b<20; $b++) {
	        $P[$a][$b] /= $sumab;
	    }
    }
    for (my $a=0; $a<20; $a++) {
	    for ($pb[$a]=0.0, my $b=0; $b<20; $b++) {
	        $pb[$a] += $P[$a][$b];
	    }
    }

    for (my $a=0; $a<20; ++$a) {
	    for (my $b=0; $b<20; ++$b) {   
	        $R[$a][$b] = $P[$a][$b]/$pb[$b];
	    }
    }
  
    ## Precompute matrix R for amino acid pseudocounts:
    for (my $a=0; $a<20; ++$a) {
	    for (my $b=0; $b<20; ++$b) {
	        $SIM[$a][$b] = log($R[$a][$b]/$pb[$a])/log(2);
	    }
    }
}


sub setupBlosumMatrix {
    my $self = shift;
    my $matrix = shift;
    my $Pref = shift;
    my $pbref = shift;

    my @P = @{$Pref};
    my @pb = @{$pbref};
    my $n = 0;

    for (my $a=0; $a<20; $a++) {
	    for ($pb[$a]=0.0, my $b=0; $b<=$a; $b++, $n++) {
	        $P[$a][$b] = $BLOSUM62[$n];
	    }
    }
    for (my $a=0; $a<19; $a++) {
	    for (my $b=$a+1; $b<20; $b++) {
	        $P[$a][$b] = $P[$b][$a]; 
	    }
    }
    for (my $a=0; $a<20; $a++) {
	    $P[$a][20] = 1.0;
	    $P[20][$a] = 1.0;
    }  
}

1;
     
