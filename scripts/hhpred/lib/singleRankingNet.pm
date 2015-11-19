#!/usr/bin/perl -w

package singleRankingNet;

require Exporter;

use strict;
use PDBResolution;
use utilities;

{
    my $nhidden = 5;
    my @bias = (2.849843e+00,2.111072e+00,-1.083691e+01,-4.135554e+00,-4.120570e-01,2.0);
    my @weights = ( -3.129355e-04, -1.589077e-02, -2.161848e+00 , 2.226650e-01,
		      -4.868254e-01,    4.035841e-03,  8.064615e-02, -1.659731e+01,
		      8.745763e-01, -1.147158e+00,    4.432926e-03 , 4.485627e-02,
		      1.612531e+01, -6.115551e+00, -3.739336e+00 ,    8.667732e-02,
		      2.469719e-01,  3.679384e+00 ,-8.063131e-01, -4.584766e+00 ,
		      -8.657036e-01, -3.492344e-01 ,-1.328472e-01 ,-5.744038e-01 , 3.360329e-02,
		      -7.531151e+00, -9.724150e-01 ,-2.660009e+00,  7.692027e-01,
		      1.438417e+00);


    sub _get_nhidden() {
		$nhidden;
    }

    sub _get_bias {
		my ($self, $i) = @_;
		$bias[$i];
    }

    sub _get_weight {
		my ($self, $i) = @_;
		$weights[$i];
    }

    sub _init {
		my ($self) = @_;
    }
}

sub new {
    my ($caller) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    no strict "refs";
    my $self = bless {}, $class;
    
    return $self;
}


sub read_from_file {
    my ($self, $paramFile) = @_;
    ## TODO
}


sub predict {
    my $self = shift;
    my $template = shift;
    my $queryLength = shift;

    my $TMscore = 0;

    my $templateName = $template->get_Hit();

    print "templateName=$templateName\n";
    print "score=" . $template->get_Score() . "\n";

    for (my $unit = 0; $unit<$self->_get_nhidden(); $unit++) { 
	#calculate input of hidden unit (sum of all inputs * weights)
	my $input = ($template->get_Score()) * $self->_get_weight(0 + $unit * 3) + 
	    $template->get_SS()/$queryLength * $self->_get_weight(1 + $unit * 3) + 
	    $template->get_SumProbL() * $self->_get_weight(2 + $unit * 3);

	#calculate output of hidden unit
	my $output = 1.0 / (1.0 + exp(-($input + $self->_get_bias($unit))));
	$TMscore += $output * $self->_get_weight( $self->_get_nhidden()*3 + $unit);	
    }		
    $TMscore = sprintf("%.6f", $TMscore);			
    #$template->set_predTM($TMscore);

    print "predict->$TMscore\n";
    return $TMscore;
}



sub predict_TMscore {
    my $self = shift;
    my $template = shift;
    my $QNeff = shift;

    $self->_init();

    my $TMscore = 0;
    my $numFeat = 5;

    my $templateName = $template->get_Hit();

    print "templateName=$templateName\n";
    print "score=" . $template->get_Score() . "\n";
    print "ssScore=" . $template->get_SS() . "\n";
    print "sumProbsL=" . $template->get_SumProbL() . "\n";
    print "sim=" . $template->get_Sim() . "\n";
    print "col=" . $template->get_Cols() . "\n";
    print "QNeff=" . $QNeff . "\n";

    my $unit = 0;
    for ($unit = 0; $unit<$self->_get_nhidden(); $unit++) { 
	#calculate input of hidden unit (sum of all inputs * weights)
	my $input = ($template->get_Score()) * $self->_get_weight(0 + $unit * $numFeat) + 
	    $template->get_SS() * $self->_get_weight(1 + $unit * $numFeat) + 
	    $template->get_SumProbL() * $self->_get_weight(2 + $unit * $numFeat) +
	    $template->get_Score()/$template->get_Cols() * $self->_get_weight(3 + $unit * $numFeat) + 
	    $template->get_Sim() * $self->_get_weight(4 + $unit*$numFeat);

	#calculate output of hidden unit
	my $output = 1.0 / (1.0 + exp(-($input + $self->_get_bias($unit))));
	$TMscore += $output * $self->_get_weight( $self->_get_nhidden()*$numFeat + $unit);	
    }		
    $TMscore += $self->_get_bias($unit);
    $TMscore = &sigmoid($TMscore);
    $TMscore = sprintf("%.6f", $TMscore);			

    print "predict->$TMscore\n";
    return $TMscore;
}


1;
