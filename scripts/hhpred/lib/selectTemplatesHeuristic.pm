#!/user/bin/perl -w

package selectTemplatesHeuristic;



use strict;
use config;
use Template;
use TemplateList;
use utilities;


our @ISA = qw(Exporter);
our @EXPORT = qw(preselectTemplates presetAccPosteriors ChooseTemplatesScoringHeuristic SingleTemplateSelection UpdateMaxProb);


#########################################################################
## input: inbase, i.e. path/name of query
##        outbase, i.e. path/name of outputs for query
##        bestByProp, i.e. number of best entries taken, e.g. by sim
## output: array of "best" templates wrt similarity, sumProb/len, score
#########################################################################
sub preselectTemplates {
    my $numBestByProp = shift;
    my $templateList = shift;

    my $verbose = 1;

    if ($verbose >= 2) {
	$templateList->print();
    }

    ## heuristic works as follows:
    ## get an hhr file and select the best X templates wrt.
    ## - similarity
    ## - sumprobs
    ## - probability
    ## and for the rest (about 15-20?) templates use the scoring procedure

    my $selectedTemplates = TemplateList->new();

    ## select most interesting 30-50 templates using a heuristic
    ## best numBestByProp wrt similarity
    $templateList->sort_by_sim();
    for (my $i=0; $i<$numBestByProp; $i++) {
	my $templ = $templateList->get($i);
	$selectedTemplates->check_and_add($templ);
    }

    ## best numBestByProp wrt sumOfProbs
    $templateList->sort_by_sumProbL();
    for (my $i=0; $i<$numBestByProp; $i++) {
	$selectedTemplates->check_and_add($templateList->get($i));
    }

    ## best numBestByProp wrt probability
    $templateList->sort_by_prob();
    for (my $i=0; $i<$numBestByProp; $i++) {
	$selectedTemplates->check_and_add($templateList->get($i));
    }

    return ($selectedTemplates);
}


########################################################
## input: array of templates (from CreateArrayInfo)
##        inbase: path/queryname to hhm file
##        outbase: path/queryname, where path is the 
##           location to write new files; all these new
##           files start with queryname
##        templateList: array of preselected templates
##           as created by CreateInfoArray         
##        
## output: hash of references to arrays containing
##         posteriori probabilities for each template
########################################################
sub presetAccPosteriors {
    print "calling presetAccPosteriors\n";

    my $queryLength = shift;
    my $outbase = shift;
    my $tlist = shift;
    my $config = shift;

    my $verbose = 3;

    my %posteriors;
 
    ## create tab files for all templates
    for (my $tl=0; $tl<$tlist->size(); $tl++)
    {
	my $template = $tlist->get($tl);

	my $querystart = $template->get_Qstart();
	my $queryend = $template->get_Qend();
	my $tempstart = $template->get_Tstart();
	my $tempend = $template->get_Tend();

	## filename for new tab file
	my $tabFileForHit = "$outbase." . $template->get_Hit() . ".HIT" . $template->get_No() . ".tab";

	## create tab file(s) if necessary (for posterior probabilities)
	if (not (-e $tabFileForHit))
	{
	    if ($verbose >= 1) {
		print "Creating $tabFileForHit...\n";
	    }
	    my $success = &BuildSingleTabFile("$outbase." . $template->get_Filt() . ".tab", $template->get_No(), $outbase);
	    if ( ! $success) {
		print "Warning: presetAccPosteriors could not build tab file $tabFileForHit!\n";
		next;
	    }
	}

	open (ALIH, "< $tabFileForHit") or die ("Cant open $tabFileForHit! $!\n");

	## save posterior probabilities for template at positions given in aligned range
	my @templatePosteriors;
	for(my $i=0; $i<$queryLength; $i++)
	{
	    push(@templatePosteriors, 0);
	}

	## go through all aligned positions
	while(<ALIH>)
	{
	    next if (/^\s*i\s+j\s+score\s+SS\s+probab/);

	    ## position in query, position in template, score, SS, posterior-probability, dssp
	    if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
	    {
		my $posAli = $1 - 1;
		my $posteriorProbab = $5;
		my $dssp = $6;

		if ($dssp eq '-') {
		    $templatePosteriors[$posAli] = 0;
		}
		else {
		    $templatePosteriors[$posAli] = $posteriorProbab;
		}

		next;
	    }

	    ## position in query, position in template, score, SS, posterior-probability
	    ## if case is only entered, if no dssp value is available
	    if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)
	    {
		if ($verbose >= 3) {
		    print "warning: no dssp value for " . $template->get_Hit() . ": i=$1, j=$2\n";
		}
		my $posAli = $1 - 1;
		my $posteriorProbab = $5;

		$templatePosteriors[$posAli] = $posteriorProbab;

	    }
	}
	close(ALIH);

	my $templateName = $template->get_Hit();
	my $hitnr = $template->get_No();

	## each accepted template gets an unique key in its posterior-probability table
	my $accKey = $templateName . "##" . $hitnr;

	## works because of reference counting
	$posteriors{$accKey} = \@templatePosteriors;

	if ($verbose >= 1) {
	    print "(created $tabFileForHit and) saved posteriors\n";
	}
    }
    
    return %posteriors;
}


sub UpdateMaxProb {
    my $maxProbPtr = shift;
    my $posteriorsPtr = shift;
    my $probHom = shift;
    my $similarity = 0;

    print "probHom=$probHom\n";
    my $increased = 0;

    for (my $i=0; $i<@{$maxProbPtr}; $i++) {
	my $probAtPos = $probHom * $posteriorsPtr->[$i] * ($similarity + 1);

	if ($probAtPos > $maxProbPtr->[$i]) {
	    $maxProbPtr->[$i] = $probAtPos;
	    $increased = 1;
	}

	my $r = sprintf("%.4f", $maxProbPtr->[$i]);
	if ($increased) {
	    print "$r(+) ";
	} else {
	    print "$r ";
	}
	$increased = 0;
    }
    print "\n";
}


## calculate score for each possible template as follows:
## S(t) = sum_{i \in A_t} { exp( alpha*(Prob_t*P_t(i) - P_{max}(i)) ) - 1}
##
## take template with greatest score as long as it is > 0
##
## TODO: split this function into 2 parts?
sub ChooseTemplatesScoringHeuristic 
{
    my $queryName = shift;
    my $queryPath = shift;
    my $queryLength = shift;
    my $outbase = shift;
    ## maximum number of templates to rank with this heuristic
    my $maxNumTemplates = shift;
    my $preselect = shift;   ## if preselect = 1 => preselect wrt score, sumprobs, and prob, 
                             ## if preselect = 2 => preselect first element in predictions
                             ## otherwise do not preselect anything
    my $templateList = shift;
    my $config = shift;


    my $verbose = 3;

    ## array for finally chosen templates
    my $templatesChosen = TemplateList->new();

    ## containing posterior probabilities for all templates
    my %accPosteriors;
    ## containing max_{t T_acc} {Prob_t * P_t(i) for all positions i
    my @maxProb;
    
    my %templateAlreadyAccepted;
    
    ## heuristic parameter
    my $alpha = 1;
    my $yshift = 0.95;
    my $overlapThreshold = 0.5;

    print "selectTempaltesHeuristic:\n";
    print "alpha=$alpha\n";
    print "yshift=$yshift\n\n";

    ## values for current best template
    my $maxScore = -9999999;
    my $bestTemplateIdx;
    my @posteriorsMax;
    my $probHomMax;
    my $similarityMax;

    my $needMoreTemplates = 1;

    my $bestRankedPtr = 0;
    my $bestRanked = "nnn";
    my %accepted;
    
    ## calculate posterior probabilities for all templates in list
    %accPosteriors = &presetAccPosteriors($queryLength, $outbase, $templateList, $config);

    ## set maxProb to zero array
    for (my $i=0; $i<$queryLength; $i++) {
	push(@maxProb, 0);
    }

    ## preselect templates wrt sumprobs, score and prob. Fill up according to scoring heuristic
    ## TODO: check for maxNumTemplates exceeding
    if ($preselect == 1) {
	## take best 3 wrt sumprobs, score and prob
	print "calling preselectTemplates...\n";
	my $preselectedTemplates = &preselectTemplates(3, $templateList);

	print "preselection wrt sumprobs, score, prob:\n";
	$preselectedTemplates->print();

	## for each preselected template, calculate posterior probabilities 
	print "calling presetAccPosteriors...\n";
	my %accPosteriorsPreset = &presetAccPosteriors($queryLength, $outbase, $preselectedTemplates, $config);

	## to be filled with all not preselected templates
	my $tmpTemplateList = TemplateList->new();

	## remove preselected templates from list and update maxProb array
	for (my $i=0; $i<$templateList->size(); $i++)
	{
	    my $template = $templateList->get($i);
	    my $hitnr = $template->get_No();
	    my $templateName = $template->get_Hit();
	    my $key = $templateName . "##" . $hitnr;

	    ## if template has been preselected: save it as 'chosen' and update maxProb
	    if ( exists($accPosteriorsPreset{$key}) ) {
		$templatesChosen->add_template($template);
		my @posteriors = @{$accPosteriorsPreset{$key}};
		my $probHom = $template->get_Prob()/100;
		my $similarity = $template->get_Sim();

		&UpdateMaxProb(\@maxProb, \@posteriors, $probHom, $similarity);
	    }
	    else {
		$tmpTemplateList->add_template($template);
	    }
	}

	## reset templateList
	$templateList = $tmpTemplateList;
    }


    ## preselect first template in list, update maxProb array, then go on with scoring heuristic
    if ($preselect == 2) {
	## first template in list is accepted for sure:
	$templateAlreadyAccepted{0} = 1;
	push(@{$accepted{$templateList->get(0)->get_Hit()}}, 0);

	$templatesChosen->add_template($templateList->get(0));
	my $keyFirst = $templateList->get(0)->get_Hit() . "##" . $templateList->get(0)->get_No();
	my @posteriorsFirst = @{$accPosteriors{$keyFirst}};
	my $probFirst = $templateList->get(0)->get_Prob() / 100;
	my $similarityFirst = $templateList->get(0)->get_Sim();
	
	print "accepting template name=" . $templateList->get(0)->get_Hit() . "hitnr=" . $templateList->get(0)->get_No() . "; score=first in list\n";
	
	## update 'maxProb' array with values from first template
	&UpdateMaxProb(\@maxProb, \@posteriorsFirst, $probFirst, $similarityFirst);
    }

    ## test remaining templates, whether to choose one more or not
    while($needMoreTemplates) {

	## if new best ranked template found => test overlap
	##
	if ($bestRankedPtr != 0) {
	    ## template with this id already accepted - test overlap
	    if (exists($accepted{$bestRanked}))
	    {
		my $overlapTooBig = 0;

		## test overlap
		foreach my $lineNr (@{$accepted{$bestRanked}})
		{
		    my $minLast = min($templateList->get($lineNr)->get_Tend(), $bestRankedPtr->get_Tend());
		    my $minFirst = min($templateList->get($lineNr)->get_Tstart(), $bestRankedPtr->get_Tstart());
		    my $maxLast = max($templateList->get($lineNr)->get_Tend(), $bestRankedPtr->get_Tend());
		    my $maxFirst = max($templateList->get($lineNr)->get_Tstart(), $bestRankedPtr->get_Tstart());

		    if (($minLast - $maxFirst)/($maxLast - $minFirst) > $overlapThreshold)
		    {
			## do not accept
			$overlapTooBig = 1;

			if ($verbose >= 1) {
			    print "=> $bestRanked overlaps " . $templateList->get($lineNr)->get_Hit() . "\n";
			}
			last;
		    }
		}
		if ($overlapTooBig == 0)
		{
		    print "=> accepted{$bestRanked} += $bestTemplateIdx\n";
		    push(@{$accepted{$bestRanked}}, $bestTemplateIdx);

		    $templatesChosen->add_template($bestRankedPtr);

		    &UpdateMaxProb(\@maxProb, \@posteriorsMax, $probHomMax, $similarityMax);

		    print "accepting template name=" . $templateList->get($bestTemplateIdx)->get_Hit() ."; hitnr=" . $templateList->get($bestTemplateIdx)->get_No() . "; score=$maxScore\n";
		}

	    }
	    else
	    {
		## new template
		if ($bestRanked ne "nnn")
		{
		    $accepted{$bestRanked} = ();
		    print "=> accepted{$bestRanked} = $bestTemplateIdx\n";

		    push(@{$accepted{$bestRanked}}, $bestTemplateIdx);
		    $templatesChosen->add_template($bestRankedPtr);

		    &UpdateMaxProb(\@maxProb, \@posteriorsMax, $probHomMax, $similarityMax);

		    print "accepting template name=" . $templateList->get($bestTemplateIdx)->get_Hit() ."; hitnr=" . $templateList->get($bestTemplateIdx)->get_No() . "; score=$maxScore\n";
		}
	    }
	}

	$bestRankedPtr = 0;
	$bestTemplateIdx = -1;
	$maxScore = -999999;

	## start scoring all templates
	for (my $tl=0; $tl<$templateList->size(); $tl++) {
	    next if (defined $templateAlreadyAccepted{$tl});

	    my $tscore = 0;
	    my $template = $templateList->get($tl);

	    ## probability of template being homologous
	    my $probHom = $template->get_Prob()/100;
	    my $similarity = $template->get_Sim();
	    
	    my $templateName = $template->get_Hit();
	    my $hitnr = $template->get_No();
	    
	    ## key should exist!?! not always, when hhsearch -tab file does not have entry for each hit - dont know why this happens...: same as in presetAccPosteriors
	    my $key = "$templateName" . "##" . "$hitnr";
	    next if (not exists($accPosteriors{$key}));
	    my @posteriors = @{$accPosteriors{$key}};
	    
	    for (my $p=0; $p<@posteriors; $p++) {
		next if ($posteriors[$p] == 0);
		
		#my $probTemp = $posteriors[$p] * $probHom * ($similarity + 1);
		my $probTemp = $posteriors[$p] * $probHom;
		
		my $exponent = $alpha * ($probTemp - $maxProb[$p]);
		$tscore += exp($exponent) - $yshift;	
	    }
	    if ($tscore > $maxScore) {
		$maxScore = sprintf("%.3f", $tscore);
		$bestTemplateIdx = $tl;
		$probHomMax = $probHom;
		$similarityMax = $similarity;
		@posteriorsMax = @posteriors;
	    }
	} ## all templates in list are scored now. 


	## test how to proceed:
	##
	if ($bestTemplateIdx == -1) {
	    print "NO new best template could be found! stopping.\n";
	    $needMoreTemplates = 0;
	    next;
	}
	## found a "valid" best template for "final" template selection
	if ($preselect == 2 && $maxScore > 0) {
	    $bestRanked = $templateList->get($bestTemplateIdx)->get_Hit();
	    $bestRankedPtr = $templateList->get($bestTemplateIdx);

	    #  push(@templatesChosen, $templateList[$maxTemplateIdx]);
	    $templateAlreadyAccepted{$bestTemplateIdx} = 1;
	}
	## no more "valid" template available
	elsif ($preselect == 2 and $maxScore <= 0) {
	    print "No more template has a score > 0! Stop searching for templates!\n";
	    $needMoreTemplates = 0;
	    next;
	}
	## fill preselection
	elsif ($preselect == 1) {
	    if ($templatesChosen->size() >= $maxNumTemplates) {
		$needMoreTemplates = 0;
		print "maximum number of templates ($maxNumTemplates) reached. stop selecting templates\n";
		next;
	    }
	    $bestRanked = $templateList->get($bestTemplateIdx)->get_Hit();
	    $bestRankedPtr = $templateList->get($bestTemplateIdx);

	    $templateAlreadyAccepted{$bestTemplateIdx} = 1;
	}

	## exceeding maximal number of templates
	if ($preselect == 2 && $templatesChosen->size() >= $maxNumTemplates) {
	    $needMoreTemplates = 0;
	    print "more templates with Score > 0 than needed!\n";
	}

	## no more templates in template list
	if( scalar(keys(%templateAlreadyAccepted)) >= $templateList->size() ) {
	    print "no more templates available ... stop with " . $templateList->size() . " templates.\n";
	    $needMoreTemplates = 0;   
	}
    }

    return $templatesChosen;
}





sub SingleTemplateSelection {
    my $templateList = shift;

    my $queryLength = $templateList->get_queryLength();
    my $templatesChosen = TemplateList->new();

    my $maxOverlap = 20;
    my $minNewCoverage = 40;
    my @coverage;

    for (my $i=0; $i<$queryLength; $i++) {
	$coverage[$i] = 0;
    }

    for (my $i=0; $i<$templateList->size(); $i++) {
	my $template = $templateList->get($i);
	my $qStart = $template->get_Qstart() - 1;
	my $qEnd = $template->get_Qend() - 1;

	my $unaligned = 0;
	my $aligned = 0;
	for (my $i=$qStart; $i<=$qEnd; $i++) {
	    if ($coverage[$i] == 0) {
		$unaligned++; 
	    } else {
		$aligned++;
	    }
	}
	if ($aligned < $maxOverlap && $unaligned > $minNewCoverage) {
	    for (my $j=$qStart; $j<$qEnd; $j++) { $coverage[$j] = 1; }
	    $templatesChosen->add_template($template);
	} 	    
    }
    return $templatesChosen;
}


sub checkCoverage {
    my $start = shift;
    my $end = shift;
    my $coveragePtr = shift;

    my $overlap = 0;
    my $uncovered = 0;

    for (my $i=$start; $i<$end; $i++) {
	$overlap++ if ($coveragePtr->[$i] == 1);
	$uncovered++ if ($coveragePtr->[$i] == 0);
    }
    return ($uncovered, $overlap);
}


sub updateCoverage {
    my $start = shift;
    my $end = shift;
    my $coveragePtr = shift;

    my $uncovered = 0;
    my $overlap = 0;

    for (my $i=$start; $i<$end; $i++) {
	if ($coveragePtr->[$i] == 0) {
	    $uncovered++;
	} elsif ($coveragePtr->[$i] == 1) {
	    $overlap++;
	}
	$coveragePtr->[$i] = 1;
    }
    return ($uncovered, $overlap);
}


1;
