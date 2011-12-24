#!/usr/bin/env perl
# 
# Create a profile (.prf) from a given HMMER/HMMER3 file

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode
my $factor = 1;        # mix orig sequence and HMMER profile with this factor (1 = 50/50, 2 = 33/66, 0.5 = 66/33 (orig/HMMER))
my $scale = 1000;
my $log2 = log(2);

my $help="
Create a profile (.prf) from a given HMMER/HMMER3 file

Usage: perl create_profile_from_hmmer.pl -i <infile> [-o <outfile>]

Options:
  -i <infile>   Input file in HMMER/HMMER3 format
  -o <outfile>  Output file in prf-format (default: infile.prf)

  -v [0-5]      verbose mode (default: $v)
\n";

# Variable declarations
my $line;
my $infile;
my $outfile;
my $i;
my $a;

my @counts;     # count profile (normalised to 1)
my $name;       # name of HMMER profile
my $len;        # length of HMMER profile
my @hmmer_prof; # HMMER profile

                 #   A  C  D  E   F  G  H  I   K   L   M  N   P  Q  R   S   T   V   W   Y   
my @hmmeraa2csaa = ( 0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18);
my @aminoacids   = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');
my %aa2i;
for ($a = 0; $a < 20; $a++) {
    $aa2i{$aminoacids[$a]} = $a;
}

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

if ($options=~s/ -i\s+(\S+) //) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}

if ($options=~s/ -factor\s+(\S+) //) {$factor=$1;}

if ($options=~s/ -v\s+(\S+) //) {$v=$1;}

if (!$infile) {print($help); print "ERROR! No input file!\n"; exit(1);}

if (!$outfile) {
    $infile =~ /^(\S+)\.\S+?$/;
    $outfile = "$1.prf";
}

##############################################################################################
# Main part
##############################################################################################

######################################
# Read HMMER sequence and profile
######################################
open (IN, $infile);
$line = <IN>;

if ($line =~ /^HMMER3/) {
    while ($line = <IN>) {
	if ($line =~ /^NAME\s+(\S+)/) {
	    $name = $1;
	} elsif ($line =~ /^LENG\s+(\d+)/) {
	    $len = $1;
	} elsif ($line =~ /^HMM/) {
	    last;
	}
    }
    $line = <IN>; $line = <IN>; # Skip header lines
    if ($line =~ /^\s*COMPO/) {
	$line = <IN>; $line = <IN>; # Skip null model lines
    }
    # Read profiles and query seq
    $i = 0;
    while ($line = <IN>) {
	if ($line =~ /^\/\//) { last; }

	$line =~ s/^\s*\d+//; # sequence position
	for ($a = 0; $a < 20; $a++) {
	    $line =~ s/^\s*(\S+)\s/ /;
	    $hmmer_prof[$i][$hmmeraa2csaa[$a]] = exp(-1.0*$1);
	}

	# Read query char in count profile
	for ($a = 0; $a < 20; $a++) {
	    $counts[$i][$a] = 0;
	}
	$line =~ /^\s*\d+\s+(\S)/;
	$counts[$i][$aa2i{$1}] = 1;
	
	$line = <IN>; $line = <IN>;
	$i++;

    }
} elsif ($line =~ /^HMMER/) {
    my @pb;
    while ($line = <IN>) {
	if ($line =~ /^NAME\s+(\S+)/) {
	    $name = $1;
	} elsif ($line =~ /^LENG\s+(\d+)/) {
	    $len = $1;
	} elsif ($line =~ /^NULE/) {
	    $line =~ s/^NULE//; # sequence position
	    for ($a = 0; $a < 20; $a++) {
		$line =~ s/^\s*(\S+)\s/ /;
		#  pb[a] = (float) 0.05 * fpow2(float(ptr)/HMMSCALE);
		$pb[$a] = 0.05 * (2**($1/1000));
	    }
	} elsif ($line =~ /^HMM/) {
	    last;
	}
    }


    $line = <IN>; $line = <IN>; # Skip header lines
    # Read profiles and query seq
    $i = 0;
    while ($line = <IN>) {
	if ($line =~ /^\/\//) { last; }

	$line =~ s/^\s*\d+//; # sequence position
	for ($a = 0; $a < 20; $a++) {
	    $line =~ s/^\s*(\S+)\s/ /;
	    # prob = pb[a]*fpow2(float(ptr)/HMMSCALE);
	    $hmmer_prof[$i][$hmmeraa2csaa[$a]] = $pb[$a] * (2**($1/1000));
	}

	# Read query char in count profile
	$line = <IN>;
	for ($a = 0; $a < 20; $a++) {
	    $counts[$i][$a] = 0;
	}
	$line =~ /^\s*(\S)/;
	$counts[$i][$aa2i{$1}] = 1;
	
	$line = <IN>;
	$i++;
    }

} else {
    print($help); 
    print "ERROR! Unknown input format!\n";
    exit(1);
}
######################################
# build count_profile (mix orig sequence and HMMER-profile)
######################################
for ($i = 0; $i < $len; $i++) {
    my $sum = 0;
    for ($a = 0; $a < 20; $a++) {
	$counts[$i][$a] += $factor * $hmmer_prof[$i][$a];
	$sum += $counts[$i][$a];
    }
    # Normalize to one
    my $fac = 1 / $sum;
    for ($a = 0; $a < 20; $a++) {
	$counts[$i][$a] *= $fac;
    }
}

######################################
# write count_profile
######################################

open (OUT, ">$outfile");
# Write header
printf(OUT "CountProfile\n");
printf(OUT "NAME\t%s\n", $name);
printf(OUT "LENG\t%i\n", $len);
printf(OUT "ALPH\t20\n");              # 20 amino acid alphabet
printf(OUT "COUNTS");
for ($a = 0; $a < 20; $a++) {
    printf(OUT "\t%s", $aminoacids[$a]);
}
printf(OUT "\tNEFF\n");

# Write profile
for ($i = 0; $i < $len; $i++) {
    printf(OUT "%i", $i+1);
    for ($a = 0; $a < 20; $a++) {
	if ($counts[$i][$a] == 0) {
	    printf(OUT "\t*");
	} else {
	    printf(OUT "\t%i", int(-(log2($counts[$i][$a]) * $scale)+0.5));
	}
    }
    printf(OUT "\t%i\n", 1 * $scale); # set Neff to 1
}

printf(OUT "//\n");
close OUT;

exit;

sub log2()
{
    my $n = shift; 
    return (log($n)/$log2);
}
