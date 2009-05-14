#!/usr/bin/perl

use strict;

my $usage="
Calculate correlation for given diversities and lengths of query and template protein
Usage: calc_correlation.pl -qneff <float> -tneff <float> -qlen <int> -tlen <int>
Option:
-qneff : query diversity
-tneff : template diversity
-qlen  : query length
-tlen  : template length
\n";

our $qneff=0;
our $tneff=0;
our $qlen=0;
our $tlen=0;

our $inputs = 4;
our $hidden = 4;

my $ARGC=scalar(@ARGV);
if ($ARGC<8) {die ($usage);}

my $options="";
for (my $i=0; $i<$ARGC; $i++) {$options.=" $ARGV[$i] ";}

# General options
if ($options=~s/ -qneff\s*(\S+) / /) {$qneff=($1/10);}
if ($options=~s/ -tneff\s*(\S+) / /) {$tneff=($1/10);}
if ($options=~s/ -qlen\s*(\S+) / /) {$qlen=(log($1)/log(1000));}
if ($options=~s/ -tlen\s*(\S+) / /) {$tlen=(log($1)/log(1000));}

while ($options=~s/\s+(\S+)//) {print("WARNING: unknown option $1\n");}

if ($qneff == 0) { print "ERROR! query diversity not given!\n"; die ($usage);}
if ($tneff == 0) { print "ERROR! template diversity not given!\n"; die ($usage);}
if ($qlen == 0) { print "ERROR! query length not given!\n"; die ($usage);}
if ($tlen == 0) { print "ERROR! template length not given!\n"; die ($usage);}

my $alpha = &alpha_NN();
my $beta = &beta_NN();

print "Correlation factor:   alpha=$alpha   beta=$beta\n";

exit;

##############################################
# sub functions
##############################################

# Calculate output of hidden neural network units
sub calc_hidden_output()
{
    my $weights_ref = $_[0];
    my $bias = $_[1];
    my $res;
    # Calculate activation of hidden unit = sum of all inputs * weights + bias
    $res = $qlen*${$weights_ref}[0] + $tlen*${$weights_ref}[1] + $qneff*${$weights_ref}[2] + $tneff*${$weights_ref}[3] + $bias;
    $res = 1.0 / (1.0 + exp(-($res)));
    return $res;
}

# Neural network regressions of alpha for Evalue correction factor
sub alpha_NN()
{
    my @biases = (7.89636,3.68944,2.05448,3.69149);  # bias for all hidden units
    my $alpha_bias = 1.33439;
    # Weights for the neural networks (column = start unit, row = end unit)
    my @weights = (
	-6.72336, -4.73393, -2.15446, -4.75140,
	-14.54957, 4.05462, 0.57951, 3.55780,
	2.08289, -1.81976, -1.19936, -17.35097,
	1.53268, -8.13514, -2.50677, 1.51106,
	6.37397, -0.36254, 0.16279, -1.32174
	);
    my $alpha=0.0;
    for (my $h = 0; $h<$hidden; $h+=1) {
	my @tmp = @weights[($h*$inputs)..($h*$inputs+$hidden-1)];
	#print "alpha:  $alpha    tmp: @tmp\n";
	$alpha += &calc_hidden_output(\@tmp,$biases[$h]) * $weights[$hidden*$inputs+$h];
    }
    $alpha = 1.0 / (1.0 + exp(-($alpha + $alpha_bias)));
    return $alpha;
}

# Neural network regressions of beta for Evalue correction factor
sub beta_NN()
{
    my @biases = (7.89636,3.68944,2.05448,3.69149);
    my $beta_bias = 5.43347;
    my @weights = (
	-6.72336, -4.73393, -2.15446, -4.75140,
	-14.54957, 4.05462, 0.57951, 3.55780,
	2.08289, -1.81976, -1.19936, -17.35097,
	1.53268, -8.13514, -2.50677, 1.51106,
	-2.27841, -7.79426, -9.53092, 3.65717
	);
    my $beta=0.0;
    for (my $h = 0; $h<$hidden; $h+=1) {
	my @tmp = @weights[($h*$inputs)..($h*$inputs+$hidden-1)];
	#print "beta:  $beta    tmp: @tmp\n";
	$beta += &calc_hidden_output(\@tmp,$biases[$h]) * $weights[$hidden*$inputs+$h];
    }
    $beta = 1.0 / (1.0 + exp(-($beta + $beta_bias)));
    return $beta;
}
