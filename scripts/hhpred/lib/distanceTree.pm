#!/user/bin/perl -w

use strict;


package distanceTree;

use predTMscore;
use singleRankingNet;
use TemplateList;
use simpleTree;
use aminoAcid;
use utilities;
use config;

my $config = HHpredConfig->instance();


our @ISA = qw(Exporter);
our @EXPORT = qw(upgma readDistances printDistances calculateSequenceWeights calculatePairwiseDistances iterativelyWeightTemplates getProbSeqAndProfSimilarityBetween neighborJoin);


## pairwise distances between query and sequences in templateList
## 
## do this by:
## - calculating pairwise similarities
## - transforming these similarities into distances
sub calculatePairwiseDistances {
    my $queryName = shift;
    my $templateList = shift;
    my $hhmDir = shift;
    my $workingDir = shift;

    my $v = 1;

    my @simSeq;
    my @simProf;
    my @simPredTMscore;
    my @prob;
    my @dist;
    my @seqNames;

    ## get sequence names, ordering is important, see upgma!
    push (@seqNames, $queryName);
    for (my $i=0; $i<$templateList->size(); $i++) {
		push (@seqNames, $templateList->get($i)->get_Hit());
    }

    print "\nSimilarities:\n";
    print "-------------\n";
    
    for (my $i=0; $i<@seqNames; $i++) {
		for (my $j=0; $j<@seqNames; $j++) {
			my $qhhmFile = "$hhmDir/$seqNames[$i].hhm";
			my $thhmFile = "$hhmDir/$seqNames[$j].hhm";			
			$simPredTMscore[$i][$j] = &getPredTMscoreSimilarityBetween($i, $j, $qhhmFile, $thhmFile, $workingDir);	    
		}
    }
   
    ## calculate pairwise distances between sequences: 
    for (my $i=0; $i<@seqNames; $i++) {
		for (my $j=0; $j<@seqNames; $j++) {

			## NN TMscore based distance
			$dist[$i][$j] = -log( &min(0.99, &max($simPredTMscore[$i][$j], 0.01)) ) / log(2);
		   
			if ($i != $j and $dist[$i][$j] <= 0.001) {
				$dist[$i][$j] = 0.01;
			}

			## minimum distance between i,j and j,i (for templates only, not wrt query)
			if ($i > $j && $i != 0 && $j != 0) {
				if ($dist[$i][$j] < $dist[$j][$i]) {
					$dist[$j][$i] = $dist[$i][$j];
				} else {
					$dist[$i][$j] = $dist[$j][$i];
				}
			}
		}
    }

    if ($v >= 1) {
	for (my $i=0; $i<@seqNames; $i++) {
	    for (my $j=0; $j<@seqNames; $j++) {
		print sprintf("%.3f ", $dist[$i][$j]);
	    }
	    print "\n";
	}
	print "--------------------------------------------------------\n";
    }

    return (@dist);
}


## calculate sequence based similarity between to sequences
## (which are extracted from a pir alignment)
sub calculateSeqSimilarityBetween {
    my $seq1 = shift;
    my $seq2 = shift;
    my $aaO = shift; ## amino-acid object

    
    if (length($seq1) != length($seq2)) {
	print "Error: calculateSimilarityBetween: sequences differ in length!\n";
    }

    my $matchCols = 1;
    my $sim = 0;

    for (my $i=0; $i<length($seq1); $i++) {
	my $AA1 = substr($seq1, $i, 1);
	my $AA2 = substr($seq2, $i, 1);

	my $aaI1 = $aaO->aa2i($AA1);
	my $aaI2 = $aaO->aa2i($AA2);

	## if match-match
	if ($aaI1 >= 0 and $aaI1 < 20 and $aaI2 >= 0 and $aaI2 < 20) {
	    $sim += $SIM[$aaI1][$aaI2];
	    $matchCols++;
	}
    }

    return $sim/$matchCols;
}


## calculate profile based similarity between sequences
## (by hhalign)
sub getProbSeqAndProfSimilarityBetween {
    my $qhhmFile = shift;
    my $thhmFile = shift;
    my $workingDir = shift;
    
    my $hhalign = $config->get_hhalign();
    
    ## use full length sequences for profile similarity?? or only ranges which are aligned to query??
    ## unfortunately, $output = `cmd -o -` does not work on queue
    ## TODO: read hhalign parameters from config
    my $cmd = "$hhalign -i $qhhmFile -t $thhmFile -mact 0.05 -o $workingDir/hhalignPairwiseOutput.hhr";
    system("$cmd");
    print "calling hhalign(" .  &get_basename($qhhmFile) . "," . &get_basename($thhmFile) . ")\n";
    my $tlist = TemplateList->new();
    $tlist->hhr_to_TemplateList("$workingDir/hhalignPairwiseOutput.hhr");
    system("rm $workingDir/hhalignPairwiseOutput.hhr");
    
    ## sequence based similarity
    my $seqSim = $tlist->get(0)->get_Sim();
    
    ## profile based similarity
    my $score = $tlist->get(0)->get_Score();
    my $aliCols = $tlist->get(0)->get_Cols();
    my $profileSim = $score / $aliCols;
    my $prob = $tlist->get(0)->get_Prob();

    return ($prob, $seqSim, $profileSim);
}


sub getPredTMscoreSimilarityBetween {
    my $i = shift;
    my $j = shift;
    my $qhhmFile = shift;
    my $thhmFile = shift;
    my $workingDir = shift;
    
    my $hhalign = $config->get_hhalign();;
    my $params = " -mact 0.05 -ssm 4 ";
    if ($i == 0 || $j == 0) { $params = " -mact 0.05 "; }
    
    ## use full length sequences for profile similarity?? or only ranges which are aligned to query??
    ## unfortunately, $output = `cmd -o -` does not work on queue
    ## TODO: read hhalign parameters from config
    my $cmd = "$hhalign -i $qhhmFile -t $thhmFile $params -o $workingDir/hhalignPairwiseOutput.hhr";
    system("$cmd");
    print "calling hhalign(" .  &get_basename($qhhmFile) . "," . &get_basename($thhmFile) . ")\n";
    my $tlist = TemplateList->new();
    $tlist->hhr_to_TemplateList("$workingDir/hhalignPairwiseOutput.hhr");
    system("rm $workingDir/hhalignPairwiseOutput.hhr");

    my $TMscoreNet = singleRankingNet->new();
    my $QLen = $tlist->get_queryLength();
    my $QNeff = $tlist->get_neff();

    print "getPredTMscoreSimiliarityBetween: TMscoreNet->predict_n11)\n";
    my $predTMscore = $TMscoreNet->predict_TMscore($tlist->get(0), $QNeff);


#    my $Tstart = $tlist->get(0)->get_Tstart();
#    my $Tend = $tlist->get(0)->get_Tend();
 #   my $TLen = $Tend - $Tstart + 1;
    #my $TLen = &get_HMM_len($thhmFile);

    my $s = $predTMscore; # * $QLen / &min($QLen, $TLen);

    return $s;
}


## go through pariwise distance matrix (upper right triangle of matrix)
## and find pair with minimal distance
sub searchOpt {
    my $hasParentPtr = shift;
    my $distPtr = shift;

    my $minI = -1;
    my $minJ = -1;
    my $minDist = 999999999;

    for (my $i=0; $i<@{$distPtr}; $i++) {
	## nodes which already have a parent are skipped
	next if ($hasParentPtr->{$i});

	for (my $j=$i+1; $j<@{$distPtr->[$i]}; $j++) {
	    next if ($hasParentPtr->{$j});
	    if ($distPtr->[$i][$j] < $minDist) {
		$minDist = $distPtr->[$i][$j];
		$minI = $i;
		$minJ = $j;
	    }
	}
    }

    return ($minI, $minJ, $minDist);
}


## add a new node (parent of node i and j) to similarity matrix
sub addNodeAndConnect {
    my $hasParentPtr = shift;
    my $distPtr = shift;
    my $clusterSizePtr = shift;
    my $nodeI = shift;
    my $nodeJ = shift;

    ## set new parent k of nodeI and nodeJ
    my $parent = scalar(@{$distPtr->[0]});
    my $k = $parent;

    $hasParentPtr->{$nodeI} = $parent;
    $hasParentPtr->{$nodeJ} = $parent;

    $clusterSizePtr->{$parent} = $clusterSizePtr->{$nodeI} + $clusterSizePtr->{$nodeJ};

    ## calculate distances from all available nodes 
    ## (i.e. the ones without parent) to new node k
    for (my $el=0; $el<@{$distPtr}; $el++) {
	if ($hasParentPtr->{$el} or $el==$nodeI or $el==$nodeJ) {
	    $distPtr->[$el][$k] = 0;
	    $distPtr->[$k][$el] = 0;
	    next;
	}

	# print "distPtr->[$nodeI][$el]=" . $distPtr->[$nodeI][$el] . "\n";
# 	print "clusterSizePtr->{$nodeI}=" . $clusterSizePtr->{$nodeI} . "\n";
# 	print "distPtr->[$nodeJ][$el]=" . $distPtr->[$nodeJ][$el] . "\n";
# 	print "clusterSizePtr->{$nodeJ}=" . $clusterSizePtr->{$nodeJ} . "\n";

	$distPtr->[$el][$k] = ($distPtr->[$el][$nodeI]*$clusterSizePtr->{$nodeI} + $distPtr->[$el][$nodeJ]*$clusterSizePtr->{$nodeJ}) / ($clusterSizePtr->{$nodeI} + $clusterSizePtr->{$nodeJ});

	$distPtr->[$k][$el] = $distPtr->[$el][$k];
    }

    return $k;
}


sub connectNodes {
    my $hasParentPtr = shift;
    my $clusterSizePtr = shift;
    my $nodeI = shift;
    my $nodeJ = shift;
    my $parent = shift;


    $hasParentPtr->{$nodeI} = $parent;
    $hasParentPtr->{$nodeJ} = $parent;

#    print "clusterSize->{$nodeI} = " . $clusterSizePtr->{$nodeI} . "\n";
#    print "clusterSize->{$nodeJ} = " . $clusterSizePtr->{$nodeJ} . "\n";

    $clusterSizePtr->{$parent} = $clusterSizePtr->{$nodeI} + $clusterSizePtr->{$nodeJ};
}


sub printDistances {
    my $distPtr = shift;
    my @dist = @$distPtr;

    print "Distances\n";
    print "---------\n";
    for (my $i=0; $i<@dist; $i++) {
	my @row = @{$dist[$i]};
	for (my $j=0; $j<@row; $j++) {
	    print sprintf("%.3f ", $row[$j]);
	}
	print "\n";
    }
    print "\n";
}


sub readDistances {
    my $file = shift;

    my @sim;
    open(SH, "< $file") or die "Cant open $file: $!\n";

    while(my $line = <SH>) {
	my @toks = split(/\s+/, $line);
	push(@sim, \@toks);
    }
    close(SH);

    return @sim;
}


## build upgma tree from distance matrix
## the distance matrix is extended in every step (for every new node)
## leaves are numbered as rows (columns) in distances, ie.
## ordering is important, see calculatePairwiseDistances
sub upgma {
    my @distances = @_;

    my $tree = simpleTree->new();
    my $v = 2;

    my $numLeafs = scalar(@distances);

    my %clusterSize;
    for (my $i=0; $i<@distances; $i++) {
	$clusterSize{$i} = 1;
	$tree->addLeaf($i);
    }

    my %hasParent;

    for (my $step=0; $step<$numLeafs-1; $step++) {

	&printDistances(\@distances);

	## search for two nodes with minimal distance
	my ($nodeI, $nodeJ, $maxDist) = &searchOpt(\%hasParent, \@distances);
	if ($v >= 2) {
	    print "searchOpt=(nodeI=$nodeI,nodeJ=$nodeJ,maxDist=$maxDist)\n";
	}

	## add a new node connecting nodeI and nodeJ
	my $parent = &addNodeAndConnect(\%hasParent, \@distances, \%clusterSize, $nodeI, $nodeJ);
	if ($v >= 2) { print "parent=$parent\n"; }

	$tree->addNode($parent);
 	$tree->connectAdaptWeightToHeight($nodeI, $nodeJ, $maxDist/2, $maxDist/2);

    }

    return $tree;
}

sub calculateMatrixQ {
    my $distPtr = shift;
    my $LPtr = shift;

    ## calculate r_i
    my @idx = sort {$a <=> $b} keys(%{$LPtr});
    for (my $i=0; $i<@idx; $i++) {
	my $sum = 0;
	for (my $k=0; $k<@idx; $k++) {
	    $sum += $distPtr->[$idx[$i]][$idx[$k]];
	}
	$LPtr->{$idx[$i]} = $sum / (scalar(@idx) - 2);
    }
    ## calculate D_ij and save pair i,j which gives minimum
    my $minI;
    my $minJ;
    my $minD = 999999999;
    for (my $i=0; $i<@idx; $i++) {
	for (my $j=0; $j<@idx; $j++) {
	    next if ($idx[$i] == $idx[$j]);
	    my $D = $distPtr->[$idx[$i]][$idx[$j]] - ($LPtr->{$idx[$i]} + $LPtr->{$idx[$j]});
	    if ($D < $minD) {
		$minD = $D;
		$minI = $idx[$i];
		$minJ = $idx[$j];
	    }
	}
    }
    return ($minI, $minJ);
}

sub newNodeWithDistances {
    my $distPtr = shift;
    my $LPtr = shift;
    my $minI = shift;
    my $minJ = shift;

    my $k = scalar(@{$distPtr});
    my @idx = sort {$a <=> $b} (keys %{$LPtr});
    for (my $m=0; $m<@idx; $m++) {
	$distPtr->[$k][$idx[$m]] = 0.5*($distPtr->[$minI][$idx[$m]] + $distPtr->[$minJ][$idx[$m]] - $distPtr->[$minI][$minJ]);
	$distPtr->[$idx[$m]][$k] = $distPtr->[$k][$idx[$m]];
	
    }
    $distPtr->[$k][$k] = 0;
    return $k;
}

sub neighborJoin {
    my @distances = @_;
    &printMatrix(\@distances, "distances in neighborjoin");

    my $v = 2;
    my $tree = simpleTree->new();
    my %L;

    ## T is set of leaf nodes, L=T
    for (my $i=0; $i<@distances; $i++) {
	$tree->addLeaf($i);
	$L{$i} = 1;
	print "add leaf $i\n" if ($v >= 2);
    }

    
    while (scalar(keys(%L)) > 2) {
	## calculate matrix Q
	## and find nodes i,j to join
	my ($minI, $minJ) = &calculateMatrixQ(\@distances, \%L);
	printMatrix(\@distances, "MatrixQ") if ($v >= 2);

	print "minI=$minI, minJ=$minJ\n" if ($v >= 2);
	
	## define new node k and calculate its distances
	my $k = &newNodeWithDistances(\@distances, \%L, $minI, $minJ);
	print "new node k=$k\n" if ($v >= 2);
	
	## calculate edge lengths ik, jk
	my $distIK = 0.5*($distances[$minI][$minJ] + $L{$minI} - $L{$minJ});
	my $distJK = $distances[$minI][$minJ] - $distIK;
	print "distIK=$distIK, distJK=$distJK\n" if ($v >= 2);

	## add node k as parent of i,j  to tree
	$tree->connect($minI, $minJ, $distIK, $distJK, $k);
	print "connect($minI,$minJ,$distIK,$distJK,$k)\n" if ($v >= 2);

	## remove i and j from L, add k
	delete $L{$minI};
	delete $L{$minJ};
	$L{$k} = 1;
    }
    if (scalar(keys(%L)) == 2) {
	## only two leaves remain: add last "edge" (in fact two edges, each having half the weight since tree is rooted)
	my @final = keys(%L);
	my $i = $final[0];
	my $j = $final[1];

	print "final i=$i\n" if ($v>=2);
	print "final j=$j\n" if ($v>=2);

	my $finalWeight = $distances[$i][$j];
	print "finalWeight=$finalWeight\n" if ($v>=2);

	$tree->connect($i, $j, $finalWeight/2, $finalWeight/2);
    }
    return $tree;
}


sub iterativelyWeightTemplates {
    my $tree = shift;

    my $treeModify = $tree;  ## copy by value
    my $root = $treeModify->getRoot();

    my @leaves = $treeModify->getAllLeavesBelow($root);
    my %weights = map {$_ => 1} @leaves; ## initial weights
    
    ## preprocessing of tree to ensure that getNextAllBelowVisitedNode works correctly
    $treeModify->visitAllLeaves();
    my $v = 2;
    my $round = 0;
    my $curNode = -1;

    ## iteratively start with "lowest" nodes and go up till root is reached
    while ($curNode != $root) {
	if ($v>= 2) {
	    print "round=$round\n";
	    print "---------\n";
	    $round++;
	}
	## find node with all its children have already been visited
	$curNode = $treeModify->getNextAllBelowVisitedNode();      
	my $parent = $treeModify->{parent}->{$curNode};
	if ($v >= 2) {
	    print "curNode=$curNode\n";
	    print "parent=$parent\n";
	}
	
	## get its leaves 
	my @curLeaves = $treeModify->getAllLeavesBelow($curNode);
	if ($v >= 2) {
	    print "curLeaves=@curLeaves\n";
	}

	## update weights
	my $t0 = $treeModify->{weight}->{$curNode};
	if ($t0 < 0.00001) {
	    $t0 = 0.0001;
	}
	if ($v >= 2) {
	    print "t0=$t0\n";
	}
	my $normalizer = 1.0/$t0;
	foreach my $leaf (@curLeaves) {
	    $normalizer += $weights{$leaf} / $tree->{weight}->{$leaf};
	}
	foreach my $leaf (@curLeaves) {
	    my $distn = $treeModify->{weight}->{$leaf};
	    if ($distn < 0.00001) { $distn = 0.0001; }
	    $weights{$leaf} *= (1.0/$t0 + 1.0/$distn) / $normalizer;
	    if ($v >= 2) {
		print "weight{$leaf}=$weights{$leaf}\n";
	    }
	}

	## update distances
	foreach my $leaf (@curLeaves) {
	    $treeModify->{weight}->{$leaf} += $t0;
	}
    }

    foreach my $leaf (sort {$a <=> $b} keys(%weights)) {
	print "final weight for $leaf: $weights{$leaf}\n";
    }

    return %weights;
}


1;
