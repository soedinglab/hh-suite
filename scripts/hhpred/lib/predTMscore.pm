#!/usr/bin/perl -w

package TMscoreNN;

require Exporter;

use strict;
use config;
use PDBResolution;
use utilities;


{
    my $nhidden = 3;
    ## these are the original parameters
    my @bias = (-0.1635528802871704, -0.3720591366291046, 5.6198554039001465);
    my @weights = (-0.007, 11.739, 3.298, -0.008, 
     		   -4.479, -1.674, 1.605, -0.041, 
     		   -0.364, -17.33, -4.629, 0.161, 
     		   0.867, -0.405, -0.265);

    ## test parameters
    # my @bias = (-5.77800,-0.11515,-0.88430);
    # my @weights = ( 0.06798,6.66194,8.12909,-0.34054,
    # 		    0.17152,-8.21288,-3.86566,0.34767,
    # 		    0.28003,0.88483,13.91985,0.12512,
    # 		    0.30894,-0.25261,0.63520);

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
}

sub new {
    my ($caller, %arg) = @_;
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


## first predict TM score of each template in templateList with a NN
## which has as inputs score, sumProb/L, SS and template resolution (in A).
## Then sort templates wrt predicted TM score
sub rank_templates {
    my ($self, $templateList, $queryLength, $config) = @_;

#    print "reading resolutions of templates...\n";
    my %resolutions = &extractAllPDBResolutions();
    

    for (my $i=0; $i<$templateList->size(); $i++) {       
	my $template = $templateList->get($i);
	# print "single NN features of " . $template->get_Hit() . ":\n";
# 	print "single NN features:\n";
# 	print "score=" . $template->get_Score() . "\n";
# 	print "ss=" . $template->get_SS() . "\n";
# 	print "qLen=" . $templateList->get_queryLength() . "\n";
# 	print "sumProbL=" . $template->get_SumProbL . "\n";
	my $TMscore = 0;
	my $templateName = $template->get_Hit();

#	print "resolution of " . $template->get_Hit() . "=" . $resolutions{$templateName} . "\n";

	for (my $unit = 0; $unit<$self->_get_nhidden(); $unit++) { 
	    #calculate input of hidden unit (sum of all inputs * weights)
	    my $input = ($template->get_Score()/50) * $self->_get_weight(0 + $unit * 4) + 
		$template->get_SS()/$queryLength * $self->_get_weight(1 + $unit * 4) + 
		$template->get_SumProbL() * $self->_get_weight(2 + $unit * 4) + 
		$resolutions{$templateName}/10 * $self->_get_weight(3 + $unit*4);

		#calculate output of hidden unit
	    my $output = 1.0 / (1.0 + exp(-($input + $self->_get_bias($unit))));
	    $TMscore += $output * $self->_get_weight( $self->_get_nhidden()*4 + $unit);	
	}		
	$TMscore = sprintf("%.6f", $TMscore);			
	$template->set_predTM($TMscore);
    }
    $templateList->sort_by_predTM();
}


sub predict {
    my $self = shift;
    my $template = shift;
    my $queryLength = shift;
    my $config = shift;

    my $TMscore = 0;
    my %resolutions = &extractAllPDBResolutions();

    my $templateName = $template->get_Hit();

    print "templateName=$templateName\n";
    print "score=" . $template->get_Score() . "\n";

    for (my $unit = 0; $unit<$self->_get_nhidden(); $unit++) { 
	#calculate input of hidden unit (sum of all inputs * weights)
	my $input = ($template->get_Score()/50) * $self->_get_weight(0 + $unit * 4) + 
	    $template->get_SS()/$queryLength * $self->_get_weight(1 + $unit * 4) + 
	    $template->get_SumProbL() * $self->_get_weight(2 + $unit * 4) + 
	    $resolutions{$templateName}/10 * $self->_get_weight(3 + $unit*4);

	#calculate output of hidden unit
	my $output = 1.0 / (1.0 + exp(-($input + $self->_get_bias($unit))));
	$TMscore += $output * $self->_get_weight( $self->_get_nhidden()*4 + $unit);	
    }		
    $TMscore = sprintf("%.6f", $TMscore);			
    #$template->set_predTM($TMscore);

    print "predict->$TMscore\n";
    return $TMscore;
}


1;
