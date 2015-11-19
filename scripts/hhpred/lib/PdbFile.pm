#!/usr/bin/perl -w


package PdbFile;

use strict;
use utilities;

sub new {
    my ($caller, $filename) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    no strict "refs";

    my $self = bless {}, $class;
    $self->{residues} = {}; ## for a residue number => atom lines concatenated
    $self->{residOrder} = {}; ## i-th residue in "seq" is residOrder-th residue in "residues"
    $self->{atm2resid} = {}; ## atom i is in residue j
    $self->{atoms} = {};
    $self->{comment} = "";
    $self->{seq} = "";
    
    if (defined($filename)) {
	$self->read($filename);
    }

    return $self;
}


sub read {
    my $self = shift;
    my $filename = shift;

    $self->clear();

    my $seq = "";
    my $nresPrev = -1e4;
    my $internResidCount = 0;

    open(PDB, $filename) or die "Cant open $filename: $!\n";
    while (my $line = <PDB>) {
	my $residue;
	my $nres;
	# ATOM           1        N         PRO          1              -29.477        -9.021        66.175  1.00 75.79
	if ($line =~ /ATOM.{2}\s*(\d+)\s.([\S|\s]{3}).(\S{3})\s.\s*(\d+)[\s|\D]\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*/) {
	    my $atm = $1;
	    $residue = $3;
	    $nres = $4;

	    $self->{atm2resid}->{$atm} = $nres;
	    if (not exists($self->{residues}->{$nres})) { $self->{residues}->{$nres} = ""; }
	    $self->{residues}->{$nres} .= $line;
	    $self->{atoms}->{$atm} = $line;

	    if ($nresPrev != $nres) {
		$seq .= &Three2OneLetter($residue);
		$nresPrev = $nres;
		$self->{residOrder}->{$internResidCount} = $nres;
		$internResidCount++;
	    }
	}
	elsif ($line =~ /^COMMENT/) {
	    $self->{comment} .= $line;
	}
    }

    $self->{seq} = $seq;

    close(PDB);
}


sub clear {
    my $self = shift;
    %{$self} = ();
}


sub residue_to_string {
    my $self = shift;
    my $nres = shift;

    return $self->{residues}->{$nres} if exists($self->{residues}->{$nres});
    return "";
}


sub print_residue {
    my $self = shift;
    my $nres = shift;

    print $self->residue_to_string($nres) . "\n";
}


sub print_seq {
    my $self = shift;
    print $self->{seq} . "\n";
}


sub get_residue_for_atom {
    my $self = shift;
    my $atm = shift;
    return $self->{atm2resid}->{$atm};
}


sub get_atom_type {
    my $self = shift;
    my $atm = shift;
    
    my $atmLine = $self->{atoms}->{$atm};
    $atmLine =~ /^ATOM.{2}\s*\d+\s.([\S|\s]{3})/;
    my $type = &trim($1);
    return $type;
} 


sub get_coordinates_of {
    my $self = shift;
    my $resid = shift;
    my $atmType = shift;
    
    my $v = 0;
    
    my @coords;

#    print "resid=$resid\n";

    if (not exists($self->{residues}->{$resid})) {
	return ((-999999));
    }

    $resid = $self->{residues}->{$resid};
    my @atoms = split(/\n/, $resid);
    
    for (my $i=0; $i<@atoms; $i++) {
	if ($atoms[$i] =~ /ATOM.{2}\s*\d+\s.([\S|\s]{3}).\S{3}\s.\s*\d+[\s|\D]\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*(-?\d+.\d+)/ && &trim($1) eq $atmType)  {
	    push(@coords, $2, $3, $4);
	    last;
	}
    }

    if (scalar(@coords) != 3) {
	if ($v >= 1) {
	    print "atom coordinates dont exist for residue=$resid, atomType=$atmType\n";
	}
	return ((-999999));
    }
    return @coords;
}


sub get_num_of_residues {
    my $self = shift;
    
    my $numOfResidues = scalar(keys(%{$self->{residues}}));
    return $numOfResidues;
}


sub get_startIdx_of_seq {
    my $self = shift;
    my $seq = shift;

    my $startIdx = &KMP($self->{seq}, $seq);
    return $startIdx;
}


## having a sequence from a pir file (aseq), build a new pdb file
## which contains only those residues in "aseq" which have "pattern" 1
## aseq must be a substring of 
##
## residue numbers in pdb-file and intern
## pdb    3 4 5 8 9 ...
## intern 1 2 3 4 5 ...
##
## intern numbers the residues in self->{seq}
## when searching aseq in self->seq with KMP one gets an intern residue number as result
## 
## residOrder{intern} = pdb, e.g. residOrder{2} = 4
sub rebuild_pdb_file {
    my $self = shift;
    my $aseq = shift;
    my $apattern = shift;
    my $outfile = shift;
    
    my @seq = split(//, $aseq);
    my @pattern = split(//, $apattern);

    if (scalar(@seq) != scalar(@pattern)) {
	print "ERROR in rebuild_pdb_file: seq and pattern have different length!\n";
	return;
    }

    ## search for start index of seq in pdbSeq => intern idx
    my $startIdx = &KMP($self->{seq}, $aseq);

    if ($startIdx == -1) {
	print "ERROR in rebuild_pdb_file: KMP could not find $aseq in $self->{seq}!\n";
	return;
    }
    
    my $newPdb = "";
    $newPdb .= $self->{comment};
    $newPdb .= "COMMENT\n";
    $newPdb .= "COMMENT artificial pdb file for multitemplates\n";

    for (my $i=0; $i<@seq; $i++) {
	## build new pdb file only with residues in pattern (i.e. pattern=1)
	next if ($pattern[$i] == 0);

	## get pdb residue number from intern residue number
	my $internResid = $startIdx + $i;
	my $res = $self->{residOrder}->{$internResid};

	$newPdb .= $self->{residues}->{$res};
    }
    $newPdb .= "END";

    ## renumber atom indices (by removing residues, their order is not any longer valid)
    my @pdbLines = split(/\n/, $newPdb);
    my $atomCounter = 1;
    for (my $i=0; $i<@pdbLines; $i++) {
	if ($pdbLines[$i] =~ /^(ATOM.{2})([\s|\d]{5})(.*)/) {
	    my $atmIdx = sprintf("%5d", $atomCounter);
	    $pdbLines[$i] = "$1$atmIdx$3";
	    $atomCounter++;
	}
    }

    $newPdb = join("\n", @pdbLines);
    open(OH, "> $outfile") or die ("Cant write $outfile: $!\n");
    print (OH $newPdb);
    close(OH);
}


## excise residues "start" till "end" from current pdb-file (in self)
## and write new pdb into outfile
sub excise_pdb_file {
    my $self = shift;
    my $start = shift;
    my $end = shift;
    my $outfile = shift;

    my $newPdb = "";
    $newPdb .= "COMMENT excised residues $start-$end\n";

    for (my $i=$start; $i<=$end; $i++) {
	if (exists($self->{residues}->{$i})) {
	    $newPdb .= $self->{residues}->{$i};
	}
    }
    $newPdb .= "END";

    open(OH, "> $outfile") or die ("Cant write $outfile: $!\n");
    print (OH $newPdb);
    close(OH);
}


## calculate distance between CA-atom of residue1 and residue2
sub distance_between {
    my $self = shift;
    my $resid1 = shift;
    my $resid2 = shift;

    if (not exists($self->{residues}->{$resid1}) || not exists($self->{residues}->{$resid2})) {
	print "ERROR: PdbFile::distance_between cant find residue ($resid1, $resid2)!\n";
	return -1;
    }
    my @coord1 = $self->get_CA_coordinates($resid1);
    my @coord2 = $self->get_CA_coordinates($resid2);

    if ($coord1[0] == -999999 || $coord2[0] == -999999) {
	return -999999;
    }

    if ($#coord1 != $#coord2) {
	print "ERROR: distance_between: coord1 and coord2 differ in length (resid1=$resid1,resid2=$resid2)!\n";
    }

    my $dist = &euklid_dist(\@coord1, \@coord2);
    return $dist;
}


sub get_CA_coordinates {
    my $self = shift;
    my $resid = shift;
    my $v = shift || 0;

    if (not exists($self->{residues}->{$resid})) {
	if ($v >= 1) {
	    print "CA coordinates dont exist for residue\n$resid\n";
	}
	return ((-999999));
    }
    $resid = $self->{residues}->{$resid};
    my @atoms = split(/\n/, $resid);

    my @coords;
    
    for (my $i=0; $i<@atoms; $i++) {
	if ($atoms[$i] =~ /ATOM.{2}\s*\d+\s.CA\s.\S{3}\s.\s*\d+[\s|\D]\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*(-?\d+.\d+)/) {
	    push(@coords, $1, $2, $3);
	    last;
	}
    }

    if (scalar(@coords) != 3) {
	if ($v >= 1) {
	    print "CA coordinates dont exist for residue\n$resid\n";
	}
	return ((-999999));
    }

    return @coords;
}


sub radius_of_gyration {
    my $self = shift;

    my @coords;
    my @center = (0,0,0);
    my $numValidResidues = 0;
    my $gyrationRadius2 = 0;

    ## calculate center:
    foreach my $residue (keys(%{$self->{residues}})) {
	@coords = $self->get_CA_coordinates($residue);
	next if ($coords[0] == -999999);
	
	for (my $j=0; $j<3; $j++) {
	    $center[$j] += $coords[$j];
	}
	$numValidResidues++;
    }
    
    for (my $j=0; $j<3; $j++) {
	$center[$j] /= $numValidResidues;
    }

    ## calculate distances to center:
    foreach my $residue (keys(%{$self->{residues}})) {
	@coords = $self->get_CA_coordinates($residue);
	next if ($coords[0] == -999999);

	my $dist = &euklid_dist(\@coords, \@center);
	$gyrationRadius2 += $dist*$dist;
    }
    $gyrationRadius2 /= $numValidResidues;
    return sqrt($gyrationRadius2);
}


## calculate distances between residues
## if a pdb file is not complete, ie. there are missing residues
## then two "subsequent" (solved) residues are not subsequent in sequence
## eg. solved residues: 4 5 8 9
## => with inbetween==0: 4-5 5-8 8-9
## if only sequence-subsequent residues shall be considered use seqsub=1
sub pairwise_distances {
    my $self = shift;
    ## how many residues between residue x and y 
    my $inbetween = shift || 0;
    my $seqsub = shift || 0;

    my %distances;

    my @residues = sort {$a <=> $b} keys(%{$self->{residues}});
    for (my $i=0; $i<scalar(@residues)-$inbetween-1; $i++) {
	my $xResidue = $residues[$i];
	my $yResidue = $residues[$i+$inbetween+1];

	next if ($seqsub==1 && $yResidue - $xResidue != $inbetween+1);
	my $dist = $self->distance_between($xResidue, $yResidue);
	next if ($dist == -999999);
	$distances{"$xResidue-$yResidue"} = $dist;
    }
    return %distances;
}

1;
