#!/user/bin/perl -w


## simple undirected binary tree

use strict;

package simpleTree;

sub new { 
    my $caller = shift;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    no strict "refs";

    my $self = bless({}, $class);

    $self->{nodes} = {};
    $self->{parent} = {};
    $self->{weight} = {};
    $self->{numNodes} = 0;
    $self->{left} = {};
    $self->{right} = {};
    $self->{visited} = {};

    if ($caller_is_obj) {
	$self->{nodes} = $caller->{nodes};
	$self->{parent} = $caller->{parent};
	$self->{weight} = $caller->{weight};
	$self->{numNodes} = $caller->{numNodes};
	$self->{left} = $caller->{left};
	$self->{right} = $caller->{right};
    }

    return $self;
}


sub addLeaf {
    my $self = shift;
    my $id = shift;

    $self->{nodes}{$id} = $self->{numNodes};

    $self->{parent}{$id} = -1;
    $self->{left}{$id} = -1;
    $self->{right}{$id} = -1;

    $self->{weight}{$id} = 0;
    $self->{numNodes}++;
}


sub addNode {
    my $self = shift;
    my $id =  shift;
    my $parent = shift;
    my $left = shift;
    my $right = shift;
    my $weight = shift;

    $self->{nodes}->{$id} = $id;
    $self->{parent}->{$id} = $parent;
    $self->{left}->{$id} = $left;
    $self->{right}->{$id} = $right;
    $self->{weight}->{$id} = $weight;
}


sub isConnected {
    my $self = shift;
    my $nodeNum1 = shift;
    my $nodeNum2 = shift;

    my $num = "$nodeNum1&$nodeNum2";
    my $numr = "$nodeNum2&$nodeNum1";

    return 1 if (exists($self->{nodes}->{$num}));
    return 1 if (exists($self->{nodes}->{$numr}));
    return 0;
}


sub connectAdaptWeightToHeight {
    my $self = shift;
    my $nodeNum1 = shift;
    my $nodeNum2 = shift;
    my $weight1 = shift;
    my $weight2 = shift;
    my $parentNum = shift;

    $parentNum = $self->{numNodes} if (not defined($parentNum));

    return if ($self->isConnected($nodeNum1, $nodeNum2));

    $self->addLeaf($parentNum);
    $self->{left}->{$parentNum} = $nodeNum1;
    $self->{right}->{$parentNum} = $nodeNum2;

    $self->{parent}->{$nodeNum1} = $parentNum;
    $self->{parent}->{$nodeNum2} = $parentNum;

    $self->{weight}->{$nodeNum1} = $weight1 - $self->getHeight($nodeNum1);
    $self->{weight}->{$nodeNum2} = $weight2 - $self->getHeight($nodeNum2);
}

sub connect {
    my $self = shift;
    my $nodeNum1 = shift;
    my $nodeNum2 = shift;
    my $weight1 = shift;
    my $weight2 = shift;
    my $parentNum = shift;

    $parentNum = $self->{numNodes} if (not defined($parentNum));

    return if ($self->isConnected($nodeNum1, $nodeNum2));

    $self->addLeaf($parentNum);
    $self->{left}->{$parentNum} = $nodeNum1;
    $self->{right}->{$parentNum} = $nodeNum2;

    $self->{parent}->{$nodeNum1} = $parentNum;
    $self->{parent}->{$nodeNum2} = $parentNum;

    $self->{weight}->{$nodeNum1} = $weight1;
    $self->{weight}->{$nodeNum2} = $weight2;
}


sub getHeight {
    my $self = shift;
    my $node = shift;

    my $height = 0;

    while ($self->{left}->{$node} != -1) {
	$height += $self->{weight}->{$self->{left}->{$node}};
	$node = $self->{left}->{$node};
    }

    return $height;
}



sub copyBranch {
    my $self = shift;
    my $newTree = shift;
    my $root = shift;

    my @stack = ();
    my %visited;
    push (@stack, $root);

    while(@stack > 0) {
	my $node = $stack[scalar(@stack)-1];

	## 
	if ($self->{left}->{$node} != -1 and not exists($visited{$self->{left}->{$node}})) {
	    push(@stack, $self->{left}->{$node});
	}
	else {
	    if ($self->{right}->{$node} != -1 and not exists($visited{$self->{right}->{$node}})) {
		push(@stack, $self->{right}->{$node});
	    }
	    else {
		
		$visited{$node} = 1;

		$newTree->{nodes}->{$node} = $self->{nodes}->{$node};
		$newTree->{parent}->{$node} = $self->{parent}->{$node};
		$newTree->{left}->{$node} = $self->{left}->{$node};
		$newTree->{right}->{$node} = $self->{right}->{$node};
		
		$newTree->{weight}->{$node} = $self->{weight}->{$node};

		pop(@stack);
	    }
	}
    }
}


sub copyNode {
    my $self = shift;
    my $newTree = shift;
    my $node = shift;

    return if ( exists($newTree->{nodes}->{$node}) );

    $newTree->{nodes}->{$node} = $self->{nodes}->{$node};
    $newTree->{parent}->{$node} = $self->{parent}->{$node};
    $newTree->{left}->{$node} = $self->{left}->{$node};
    $newTree->{right}->{$node} = $self->{right}->{$node};
    
    $newTree->{weight}->{$node} = $self->{weight}->{$node};

    $newTree->{numNodes}++;
}



## root the tree at one of its leafs "L" , i.e. rearrange the whole tree
## so that finally the parent of "L" is new root
sub reRoot {
    my $self = shift;
    my $newRoot = shift;

    if (not $self->isLeaf($newRoot)) {
	print "Error: $newRoot is not a leaf!\n";
	return -1;
    }

    my $newTree = simpleTree->new();

    my $prevNode = $newRoot;
    my $prevParent = -1;
    
    my $newParent = -1;
    my $parent = $self->{parent}->{$newRoot}; 
    my $weight = $self->{weight}->{$newRoot};

    ## save new "root node" for update at the end
    my $newRootWeight = $self->{weight}->{$newRoot};
    my $saveParent = $parent;

    #$newTree->addNode($newRoot, -1, -1, -1, $self->{weight}->{$newRoot});

    ## if parent of "new root" is already root: just set "new root" as root and update weight
    if ($parent == $self->getRoot()) {
	if ($self->{left}->{$parent} == $newRoot) {
	    $self->copyBranch($newTree, $self->{right}->{$parent});
	    $newTree->{parent}->{$self->{right}->{$parent}} = -1;
	    ## update weight of new root
	    $newTree->{weight}->{$self->{right}->{$parent}} += $weight;
	}
	else {
	    $self->copyBranch($newTree, $self->{left}->{$parent});
	    $newTree->{parent}->{$self->{left}->{$parent}} = -1;
	    ## update weight of new root
	    $newTree->{weight}->{$self->{left}->{$parent}} += $weight;
	}

	return $newTree;
    }

    ## if node which is to become root isnt directly below root in tree:
    ## set "new root" as root and go up the tree; at each level, copy branch not
    ## containing Q and reverse roles of current parent and child
    while ($parent != $self->getRoot()) {
	## insert current parent into new tree - set parent explicitly!
	$self->copyNode($newTree, $parent);
	$newTree->{parent}->{$parent} = $prevParent;


	## coming from left node, then copy right branch
	if ($self->{left}->{$parent} == $prevNode) {
	    $self->copyBranch($newTree, $self->{right}->{$parent});
	    $newParent = $self->{parent}->{$parent};

	    ## root of tree is reached: copy branch not containing Q
	    if ($newParent == $self->getRoot()) {
		## if Q is in right branch
		if ($self->{right}->{$newParent} == $parent) {
		    $newTree->{left}->{$parent} = $self->{left}->{$newParent};
		    $self->copyBranch($newTree, $self->{left}->{$newParent});
		    $newTree->{parent}->{$self->{left}->{$newParent}} = $parent;
		    ## update weight
		    $newTree->{weight}->{$self->{left}->{$newParent}} += $self->{weight}->{$parent};
		}
		## if Q is in left branch
		else {
		    $newTree->{left}->{$parent} = $self->{right}->{$newParent};
		    $self->copyBranch($newTree, $self->{right}->{$newParent});
		    $newTree->{parent}->{$self->{right}->{$newParent}} = $parent; 
		    ## update weight
		    $newTree->{weight}->{$self->{right}->{$newParent}} += $self->{weight}->{$parent};
		}
	    }
	    ## change role of parent and child
	    else {
		$newTree->{left}->{$parent} = $newParent;
		$newTree->{parent}->{$newParent} = $parent;
	    }

	    ## go up one level
	    $prevNode = $parent;
	    $prevParent = $parent;
	    $parent = $newParent;
	}
	## coming from right node, then copy left branch
	elsif ($self->{right}->{$parent} == $prevNode) {
	    $self->copyBranch($newTree, $self->{left}->{$parent});
	    $newParent = $self->{parent}->{$parent};

	    if ($newParent == $self->getRoot()) {
		if ($self->{right}->{$newParent} == $parent) {
		    $newTree->{right}->{$parent} = $self->{left}->{$newParent};
		    $self->copyBranch($newTree, $self->{left}->{$newParent});
		    $newTree->{parent}->{$self->{left}->{$newParent}} = $parent;
		    ## update weight
		    $newTree->{weight}->{$self->{left}->{$newParent}} += $self->{weight}->{$parent};
		}
		else {
		    $newTree->{right}->{$parent} = $self->{right}->{$newParent};
		    $self->copyBranch($newTree, $self->{right}->{$newParent});
		    $newTree->{parent}->{$self->{right}->{$newParent}} = $parent; 
		    ## update weight
		    $newTree->{weight}->{$self->{right}->{$newParent}} += $self->{weight}->{$parent};
		}
	    }
	    else {
		$newTree->{right}->{$parent} = $newParent;
		$newTree->{parent}->{$newParent} = $parent;
	    }

	    $prevNode = $parent;
	    $prevParent = $parent;
	    $parent = $newParent;
	}
    }

    ## finally update weight of new root (i.e. parent of "newRoot")
    $newTree->{weight}->{$saveParent} = $newRootWeight;

    return $newTree;
}



sub getRoot {
    my $self = shift;
    
    foreach my $node (keys %{$self->{nodes}}) {
	if ($self->{parent}->{$node} == -1) {
	    return $node;
	}
    }

    return -1;
}


## weight each leaf in tree with formula (5.9), Durbin
## i.e. do postorder traversal
sub weightScheme {
    my $self = shift;

    my %alreadyVisited;
    my %leafWeight;

    my @stack;


    my $root = $self->getRoot();

    ## initialize leaf weights with edge lengths of leafs
    my @leafs = $self->getAllLeavesBelow($root);

    foreach my $leaf (@leafs) {
	$leafWeight{$leaf} = $self->{weight}->{$leaf};
    }


    ## reweight edge length according to formula (5.9):
    ## do postorder traversal and for each inner node (these are visited in appropriate manner so that
    ## all nodes below have already been visited) share its weight among its leafs 
    push(@stack, $root);

    while (scalar(@stack) != 0) {
	my $currentNode = $stack[$#stack];

	my $left = $self->{left}->{$currentNode};
	my $right = $self->{right}->{$currentNode};

	## if unvisited inner node
	if ($left != -1 and not exists($alreadyVisited{$left})) {
	    push(@stack, $left);
	}
	else {
	    if ($right != -1 and not exists($alreadyVisited{$right})) {
		push(@stack, $right);
	    }
	    else {
		#print "$currentNode\n";

		if ($self->isLeaf($currentNode)) {
		    ## skip
		}
		## postorder traversal: all nodes below currentNode have been visited
		else {
		    ## share weight at parent among all leafs
		    my $parentWeight = $self->{weight}->{$currentNode};
		    my @myLeafs = $self->getAllLeavesBelow($currentNode);
		    
		    my $normalizer = 0;
		    foreach my $leaf (@myLeafs) {
			$normalizer += $leafWeight{$leaf};
		    }

		    foreach my $leaf (@myLeafs) {
			$leafWeight{$leaf} = $leafWeight{$leaf} + $parentWeight * $leafWeight{$leaf}/$normalizer;
			print "new leaf weight for $leaf: $leafWeight{$leaf} ($parentWeight, $normalizer)\n";
		    }
		    
		}

		$alreadyVisited{$currentNode} = 1;
		pop(@stack);
	    }
	}
    }

    return %leafWeight;
}


## calculate mean branch length from all leafs to root
sub meanBranchLength {
    my $self = shift;

    my @leafs = $self->getAllLeavesBelow($self->getRoot());
    my $totalLen = 0;

    for (my $i=0; $i<@leafs; $i++) {
	my $parent = $self->{parent}->{$leafs[$i]};

	my $branchLen = $self->{weight}->{$leafs[$i]};	

	while ($parent != -1) {
	    $branchLen += $self->{weight}->{$parent};
	    $parent = $self->{parent}->{$parent};
	}

	$totalLen += $branchLen;
    }

    my $meanLen = $totalLen / scalar(@leafs);

    return ($meanLen);
}


sub isLeaf {
    my $self = shift;
    my $node = shift;

    return 1 if $self->{left}->{$node} == -1 and $self->{right}->{$node} == -1;
    return 0;
}


## return all leaves below node "node"
sub getAllLeavesBelow {
    my $self = shift;
    my $node = shift;

    my @leaves;

    my @stack;

    push(@stack, $node);

    while(scalar(@stack) != 0) {
	my $currentNode = pop(@stack);

	## new leaf found
	if ($self->{left}->{$currentNode} == -1 and $self->{right}->{$currentNode} == -1) {
	    push(@leaves, $currentNode);
	}
	if ($self->{left}->{$currentNode} != -1) {
	    push(@stack, $self->{left}->{$currentNode});
	}
	if ($self->{right}->{$currentNode} != -1) {
	    push(@stack, $self->{right}->{$currentNode});
	}

    }

    return @leaves;
}


sub visitAllLeaves {
    my $self = shift;

    my $root = $self->getRoot();
    my @leaves = $self->getAllLeavesBelow($root);

    foreach my $node (keys %{$self->{nodes}}) {
	$self->{visited}->{$node} = 0;
    }

    foreach my $leaf (@leaves) {
	$self->{visited}->{$leaf} = 1;
    }
}


## get next node which has only "leaves" attached
sub getNextAllBelowVisitedNode {
    my $self = shift;
    
    my $root = $self->getRoot();

    while ( $self->{visited}->{$root} != 1) {
	if ( ! $self->{visited}->{$self->{left}->{$root}}) {
	    $root = $self->{left}->{$root}; 
	} else {
	    if ( ! $self->{visited}->{$self->{right}->{$root}}) {
		$root = $self->{right}->{$root};
	    } else {
		$self->{visited}->{$root} = 1;
	    }
	}
    }
    return $root;
}


sub lca {
    my $self = shift;
   
}


sub print {
    my $self = shift;

    foreach my $node (keys %{$self->{nodes}}) {
	#print "node $node\n";
	print ("node: $node has parent " . $self->{parent}->{$node} . " with weight " . $self->{weight}->{$node} . " left=" . $self->{left}->{$node} . " right=" . $self->{right}->{$node} ."\n");
    }
}

1;
