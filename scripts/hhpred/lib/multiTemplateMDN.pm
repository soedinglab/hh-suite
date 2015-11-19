#!/user/bin/perl -w


package multiTemplateMDN;

use strict;
use Math::Trig;  ## for tanh

#our @ISA = qw(Exporter);
#our @EXPORT = qw(readMDNWeights getMixtureParametersFromMDN);


sub activationFunction
{
    my $value = shift;
    my $fct = shift;

    return 1.0/(1.0 + exp(-$value)) if ($fct =~ /L/);
    return tanh($value) if ($fct =~ /H/);
    return $value if ($fct =~ /N/);
}




##################################################################################
## input: reference to vector and reference to matrix, which are to be multiplied
##        dimensions must be correct
## output: matrix * vector => vector
##################################################################################
sub vectorTimesMatrix {
    my ($vecPtr, $matPtr) = @_;

    my @vec = @$vecPtr;
    my @matrix = @$matPtr;

    ## test dimensions
    if (scalar(@matrix) == 0) {
	print "matrix is empty!\n";
	exit(1);
    }
    if (scalar(@vec != scalar(@{$matrix[0]}))) {
	print "Error: vector and matrix have non-corresponding dim!\n";
	exit(1);
    }

    my @erg;
    for (my $i=0; $i<scalar(@matrix); $i++) {
	my $scaProd = 0;
	for (my $j=0; $j<scalar(@{$matrix[$i]}); $j++) {
	    $scaProd += $vec[$j] * $matrix[$i][$j];
	}
	push(@erg, $scaProd);
    }

    return @erg;
}


########################################################
## weights has the following structure:
## w_11 w_21 w_31 bias
## w_12 w_22 w_32 bias
## ...
## w_15 2_25 w_35 bias
##
## one row represents the weights w_ij 
## between node i and j (i is one layer before j)
########################################################
## requires two arguments:
## file for weights of layer1
## file for weights of layer2
##
########################################################
sub new {
    my $class = shift;

    my $self = {
	numberOfInputs => undef,
	numberOfHiddenNodes => undef,
	numberOfOutputs => undef,
	dimTarget => 1,
	hiddenLayerActivationFunction => undef,
	outputLayerActivationFunction => undef,
	weightsLayer1 => [],
	weightsLayer2 => []
    };

    bless($self, $class);

    $self->readMDNWeights(shift, shift);

    return $self;
}


## input: weight file name for layer 1 and layer 2
## output: fills weightsLayer1 and weightsLayer2 and some other
## variables
##
sub readMDNWeights {
    my $self = shift;
    my $MDNLayer1WeightsFile = shift;
    my $MDNLayer2WeightsFile = shift;

    my @weightsLayer1;
    my @weightsLayer2;

    open (WH1, "< $MDNLayer1WeightsFile") or die "Cant open $MDNLayer1WeightsFile! $!\n";

    my $weightLineRead = 0;

    while(<WH1>) {
	## ignore comments
	next if (/\s*#/);
	
	if(/hiddenLayerActivationFunction\s*=\s*(\S+)/) {
	    $self->{hiddenLayerActivationFunction} = $1;
	}
	## weights (at least 2)
	elsif (/^\s*(\S+\s+)+\S+/) {
	    my @weights = split(/\s+/);

	    $weightsLayer1[$weightLineRead] = \@weights;

       	    $weightLineRead++;
	}
    }

    close (WH1);

    my $numberOfInputs = scalar(@{$weightsLayer1[0]}) - 1;
    $self->{numberOfInputs} = $numberOfInputs;

    print "numberOfInputs=$numberOfInputs\n";

    my $numberOfHiddenNodes = $weightLineRead;
    $self->{numberOfHiddenNodes} = $numberOfHiddenNodes;

    ## reset
    $weightLineRead = 0;

    open(WH2, "< $MDNLayer2WeightsFile") or die "Cant open $MDNLayer2WeightsFile! $!\n";

    while(<WH2>) {
	next if (/\s*#/);

	if (/outputLayerActivationFunction\s*=\s*(\S+)/) {
	    $self->{outputLayerActivationFunction} = $1;
	}
	elsif (/^\s*(\S+\s+)+\S+/) {
	    my @weights = split(/\s+/);
	    $weightsLayer2[$weightLineRead] = \@weights;

	    $weightLineRead++;
	}
    }

    my $numberOfOutputs = $weightLineRead;
    $self->{numberOfOutputs} = $numberOfOutputs;

    close(WH2);

    $self->{weightsLayer1} = \@weightsLayer1;
    $self->{weightsLayer2} = \@weightsLayer2;
}


##############################################################
## input: vector containing input for MDN
## output: MDN result, i.e. parameters of Gaussian components
##         e.g. (pi1, pi2, mu1, mu2, sigma1, sigma2)
##############################################################

sub getMixtureParametersFromMDN {
    my $self = shift;
    ## vector containing input for MDN, e.g. (distance, post, sim)
    my @input = @_;

    if (scalar(@input) != $self->{numberOfInputs}) {
	print "number of inputs doesnt match weights in first layer!\n";
	exit(1);
    }

    ## since the weight matrices have bias as last column, 
    ## input must get a "1" as last entry
    push (@input, 1);

    my @inputLayer1 = &vectorTimesMatrix(\@input, $self->{weightsLayer1});
    my @outputLayer1;

    for (my $i=0; $i<@inputLayer1; $i++) {
	$outputLayer1[$i] = &activationFunction($inputLayer1[$i], $self->{hiddenLayerActivationFunction});
    }

    push (@outputLayer1, 1);
    my @outputLayer2 = &vectorTimesMatrix(\@outputLayer1, $self->{weightsLayer2});

    for (my $i=0; $i<@outputLayer2; $i++) {
	$outputLayer2[$i] = &activationFunction($outputLayer2[$i], $self->{outputLayerActivationFunction});
    }


    ## each component is decribed by 3 parameters:
    ## pi, mu, sigma (each one-dimensional)
    my $numberOfComponents = $self->{numberOfOutputs} / 3;

    ## from "raw output", calculate now the "real" parameters
    ## pi
    my $pi_normalizer = 0;

    for (my $i=0; $i<$numberOfComponents; $i++) {
	$outputLayer2[$i] = exp($outputLayer2[$i]);
	$pi_normalizer += $outputLayer2[$i];
    }
    for (my $i=0; $i<$numberOfComponents; $i++) {
	$outputLayer2[$i] /= $pi_normalizer;
    }

    ## mu stays unchanged, i.e. network emits mu directly

    ## sigma
    for (my $i=0; $i<$numberOfComponents; $i++) {
	$outputLayer2[$#outputLayer2 - $i] = sqrt( exp($outputLayer2[$#outputLayer2 - $i]) );
    }

    return @outputLayer2;
}

1;

## &readWeights();
## my @input = (5, 0.9, 0.7);
## my @erg = &getMixtureParametersFromMDN(@input);
## print "@erg\n";
