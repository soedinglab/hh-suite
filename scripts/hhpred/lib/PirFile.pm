package PirFile;

use FastaFile;
use utilities;
use PdbFile;
use strict;

sub new {
    my ($caller, $filename) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    no strict "refs";
    my $self = bless [], $class;
    
    if (defined($filename)) {
	$self->read_from_file("$filename");
    }
    return $self;
}


sub read_from_file {
    my ($self, $filename) = @_;

    my @lines;
    open(FH, "< $filename") or die ("Cant open $filename: $!\n");
    @lines = <FH>;
    close(FH);

    my $readSeq = 0;
    my $seq = "";
    my $id;

    my $idxEntry = 0;

    for (my $i=0; $i<@lines; $i++) {
	my $curLine = $lines[$i];       

	if ($curLine =~ /^\s*>\S+;(\S+)/) {
	    if ($readSeq == 1) {
		$self->[$idxEntry]->{seq} = $seq;
		$idxEntry++;
	    }

	    $id = $1;

	    ## structure or sequence description line
	    $i++;
	    chomp($lines[$i]);
	    $curLine = $lines[$i];
	    $self->[$idxEntry] = {id=>$id, descr=>$curLine};

	    $readSeq = 1;
	    $seq = "";
	    next;
	}

	if ($readSeq == 1) {
	    chomp($curLine);
	    $seq .= $curLine;
	}
    }
    $self->[$idxEntry]->{seq} = $seq;
}


sub size {
    my $self = shift;
    scalar(@{$self});
}


sub get_seq {
    my ($self, $idx) = @_;
    $self->[$idx]->{seq};
}


sub get_id {
    my ($self, $idx) = @_;
    $self->[$idx]->{id};
}


sub get_descr {
    my ($self, $idx) = @_;
    $self->[$idx]->{descr};
}


sub add_entry {
    my ($self, $id, $descr, $seq) = @_;
    
    my $idx = $self->size();
    $self->[$idx]->{id} = $id;
    $self->[$idx]->{descr} = $descr;
    $self->[$idx]->{seq} = $seq;
}

sub to_string {
    my $self = shift;

    my $res = "";

    for (my $i=0; $i<$self->size(); $i++) {
	## skip gap lines (my be introduced in reduce_to_cores)
	next if ($self->[$i]->{seq} =~ /^-+\*/);
#	my $hit = &get_basename($self->[$i]->{id});
#	my $chain = &get_PDB_chain($hit);
#	my $len = &get_seq_len($self->[$i]->{seq});
	$res .= ">P1;" . $self->[$i]->{id} . "\n";
	$res .= $self->[$i]->{descr} . "\n" if ($self->[$i]->{descr} ne "");
#	$res .= "sequence:$hit: 1: : $len: : : : : \n" if ($self->[$i]->{descr} =~ /\S:\@-\@/ && $i==0);
#	$res .= "structureX:$hit: : $chain: : $chain: dummy : dummy : : \n" if ($self->[$i]->{descr} eq "" && $i>0);
	$res .= $self->[$i]->{seq} . "\n";
    }

    return $res;
}


sub print {
    my $self = shift;

    print $self->to_string();
}


sub write_to_file {
    my ($self, $outfile) = @_;

    open (OH, "> $outfile") or die "Cant write to $outfile: $!\n";
    print (OH $self->to_string());
    close(OH);
}


sub to_fasta {
    my ($self) = @_;
    
    my $fasta = FastaFile->new();

    for (my $i=0; $i<$self->size(); $i++) {
	my $seq = $self->[$i]->{seq};
	## remove final '*'
	$seq =~ s/\*$//;
	$fasta->add_entry($self->[$i]->{id}, $self->[$i]->{descr}, $seq);
    }
    
    return $fasta;
}



## transform current pir file into a new one (written to outdir/queryId.pir)
## in this new one, only alignment positions covered by the first template (the cores)
## are used in the remaining templates, i.e. residues uncovered by the first template
## but by the 2nd, 3rd,... will be deleted, i.e. replaced by gaps.
## In order to provide Modeller with apropriate pdb files for these "new" templates,
## new pdb-files are generated and written to outdir
sub reduce_to_cores {
    my $self = shift;
    my $outdir = shift;
    my $pdbdir = shift;
    my $firstTemplEntry = shift; # which entry in pir file is first (NN-ranked) template 
                                 # this is needed because cores are wrt first template

    my $reducedPirFile = $self;

    my @pattern;
    if (@{$self} < 2) {
	print "ERROR in reduce_to_cores: no template in alignment: $!\n";
	return;
    }

    ## get first template
    print "first entry id=" . $self->[$firstTemplEntry]->{seq} . "\n";
    my $firstTemplate = $self->[$firstTemplEntry]->{seq};
    $firstTemplate =~ s/\*$//;
    my @firstTemp = split(//, $firstTemplate);

    ## create pattern of first template, i.e. 1 if residue, 0 if gap
    for (my $i=0; $i<@firstTemp; $i++) {
	my $pat = 1;
	$pat = 0 if ($firstTemp[$i] eq "-");
	push(@pattern, $pat);
    }

    my %ids;

    ## all other templates are not allowed to cover new residues, i.e.
    ## have pat=1 where pattern==0 (but vice versa is allowed)
    for (my $i=1; $i<@{$self}; $i++) {

	my $tempSeq = $self->[$i]->{seq};
	my $idOrig = $self->[$i]->{id};
	my $idAct = $idOrig;

	if ($i==$firstTemplEntry) {
	    ## handle case if firstTemplate is coming up several times, see below
	    next;
	}
	
	## a template can be used several times (ie. occur several times in pir file)
	## thus the name of the new pdb file has to change by adding an index
	## because otherwise the file will be overwritten for each same template
	if (exists($ids{$idOrig})) {
	    my $occurrence = $ids{$idOrig};
	    $idAct = $idOrig . "" . $occurrence;
	    $reducedPirFile->[$i]->{id} = $idAct;
	    if ($reducedPirFile->[$i]->{descr} =~ /(.+:)$idOrig(:.+)/) {
	       $reducedPirFile->[$i]->{descr} = "$1$idAct$2";
	    } else {
		print "No match for: " . $reducedPirFile->[$i]->{descr} . "\n";
	}
	    $ids{$idOrig}++;
	} else {
	    $ids{$idOrig} = 0;
	}

	$tempSeq =~ s/\*$//;
	my @template = split(//, $tempSeq);

	my $tempPattern = "";

	## test each residue in current template
	for (my $j=0; $j<@template; $j++) {
	    if ($template[$j] ne "-" && $pattern[$j]==0) {
		$template[$j] = "-";
		$tempPattern .= "0";
	    } elsif ($template[$j] ne "-") {
		$tempPattern .= "1";
	    }
	}

	$reducedPirFile->[$i]->{seq} = join("", @template) . "*";

	## read in original pdb file
	my $templatePdb = PdbFile->new();
	my $pdbFileName = "$pdbdir/$idOrig.pdb";
	$templatePdb->read($pdbFileName);

	## create new adapted pdb file
	my $newTempPdbFile = "$outdir/$idAct.pdb";
	$tempSeq =~ s/-//g;
	$templatePdb->rebuild_pdb_file($tempSeq, $tempPattern, $newTempPdbFile);
    }
    my $queryId = $reducedPirFile->[0]->{id};

    ## find gaps-only columns (may have been created during core generation)
    my @gapPos = &match_all_positions("-+", $reducedPirFile->[0]->{seq});
    ## check if gap is in first template => gap column, otherwise => insert
    my @gapColPos;
    for (my $i=0; $i<@gapPos; $i++) {
	my $start = $gapPos[$i]->[0];
	my $end = $gapPos[$i]->[1];
	my $rangeInTemp1 = substr($reducedPirFile->[$firstTemplEntry]->{seq}, $start, $end-$start);
	if ($rangeInTemp1 =~ /^-+$/) {
	    push(@gapColPos, $gapPos[$i]);
	}
    }
    ## remove gap columns
    for (my $i=0; $i<$reducedPirFile->size(); $i++) {
	$reducedPirFile->[$i]->{seq} = &remove_ranges($reducedPirFile->[$i]->{seq}, @gapColPos);
    }

    $reducedPirFile->write_to_file("$outdir/$queryId.pir");
}

## this is a slight modification of the upper reduce_to_cores function:
## the former one reduces all templates (but the first) to the cores defined by the first template
## this one reduces only the first template to cores by removing residues given in hash %remmoveResidues
sub remove_blocks_from_first {
    my $self = shift;
    my $outdir = shift;
    my $pdbdir = shift;
    my $firstTemplEntry = shift; # which entry in pir file is first (NN-ranked) template 
                                 # this is needed because cores are wrt first template

    my $removeResiduesPtr = shift;

    my %removeResidues = %{$removeResiduesPtr};

    my $reducedPirFile = $self;

    my @pattern;
    if (@{$self} < 2) {
	print "ERROR in reduce_to_cores: no template in alignment: $!\n";
	return;
    }

    ## get first template
    print "first entry id=" . $self->[$firstTemplEntry]->{seq} . "\n";
    my $firstTemplate = $self->[$firstTemplEntry]->{seq};
    $firstTemplate =~ s/\*$//;
    my @firstTemp = split(//, $firstTemplate);

    ## create pattern of first template, i.e. 1 if residue, 0 if gap
    my $residueNr = 1;
    my $pattern = "";    ## pattern for building new pdb file (1 if residue, 0 if gap)
    my $newFirstTemplSeq = ""; ## modified first template sequence in pir file (the additional gaps

    for (my $i=0; $i<@firstTemp; $i++) {
	if ($firstTemp[$i] eq "-" ) { $newFirstTemplSeq .= "-"; }
	else {	    
	    if (exists($removeResidues{$residueNr})) {
		$pattern .= "0";		
		$newFirstTemplSeq .= "-";
	    }
	    else {
		$pattern .= "1";
		$newFirstTemplSeq .= $firstTemp[$i];
	    }

	    $residueNr++;
	}
    }

    my $oldFirstSeq = $reducedPirFile->[$firstTemplEntry]->{seq};
    $reducedPirFile->[$firstTemplEntry]->{seq} =  "$newFirstTemplSeq*";


    ## create new pdb file where residues are removed
    ## MODELLER needs this new pdb file

    ## read in original pdb file
    my $firstTemplName = $self->[$firstTemplEntry]->{id};
    my $templatePdb = PdbFile->new();
    my $pdbFileName = "$pdbdir/$firstTemplName.pdb";
    $templatePdb->read($pdbFileName);
    
    my $newTempPdbFile = "$outdir/$firstTemplName.pdb";
    my $tempSeq = $oldFirstSeq;
    $tempSeq =~ s/-//g;
    $tempSeq =~ s/\*$//;

    print "pattern=$pattern\n";
    print "        $tempSeq\n";


    $templatePdb->rebuild_pdb_file($tempSeq, $pattern, $newTempPdbFile);

    ## copy remaining pdb files to outdir (needed for MODELLER)
    for (my $i=1; $i<$self->size(); $i++) {
	next if ($i==$firstTemplEntry);
	my $templName = $self->[$i]->{id};
	my $pdbFileName = "$pdbdir/$templName.pdb";
	system("cp $pdbFileName $outdir");
    }

    my $queryId = $reducedPirFile->[0]->{id};

    ## find gaps-only columns (may have been created during core generation)
    my @gapPos = &match_all_positions("-+", $reducedPirFile->[0]->{seq});
    ## check if gap is in all sequences
    my @gapColPos;
    for (my $i=0; $i<@gapPos; $i++) {
	my $start = $gapPos[$i]->[0];
	my $end = $gapPos[$i]->[1];
	my $success = 1;
	for (my $j=0; $j<$reducedPirFile->size(); $j++) {
	    my $rangeInTemp1 = substr($reducedPirFile->[$j]->{seq}, $start, $end-$start);
	    if ($rangeInTemp1 !~ /^-+$/) {
		$success = 0;
		last;
	    }
	}
	if ($success) {
	    push(@gapColPos, $gapPos[$i]);
	}
    }
    ## remove gap columns
    for (my $i=0; $i<$reducedPirFile->size(); $i++) {
	$reducedPirFile->[$i]->{seq} = &remove_ranges($reducedPirFile->[$i]->{seq}, @gapColPos);
    }

    $reducedPirFile->write_to_file("$outdir/$queryId.pir");

    ## now, all necessary pdb files and the final pir file should be placed in outdir so that MODELLER can be 
    ## started
}


## find displaced gaps in the query, i.e.
## find situations where in the query there are two adjacent residues
## aligned to two template positions j1 and j2 (j2 > j1+1) and template
## distance is > 7A
sub num_of_displaced_gaps {
    my $self = shift;
    my $pdbdir = shift; 
    my $v = 1; 

    my @gapPos;
    my $numDisplacedGaps = 0;

    my $qseq = $self->[0]->{seq};
    ## find gaps in query
    @gapPos = &match_all_positions("-+", $qseq);
    if ($v >= 1) {
	print "gaps in query:\n";
	for (my $i=0; $i<@gapPos; $i++) {
	    print "$gapPos[$i]->[0] - $gapPos[$i]->[1]\n";
	}
    }

    ## skip front gaps (should not happen)
    if (scalar(@gapPos) > 0 && $gapPos[0]->[0] == 0) {
	shift @gapPos;
    }

    ## check for each template whether j2 > j1+1 and whether distance is > 7A
    for (my $i=1; $i<$self->size(); $i++) {
	my $tseq = $self->[$i]->{seq};
	my $tpdb = PdbFile->new();
	my $tpdbFileName = "$pdbdir/" . $self->[$i]->{id} . ".pdb";
	if ($v >= 1) {
	    print "checking template $tpdbFileName\n";
	}
	$tpdb->read($tpdbFileName);
	
	for (my $j=0; $j<@gapPos; $j++) {
	    my $start = $gapPos[$j]->[0];
	    my $end = $gapPos[$j]->[1];
	    ## get subsequence of template where there is a gap in query
	    my $excerpt = substr($tseq, $start, $end-$start);
	    ## test whether it is a gap => continue
	    next if ($excerpt =~ /^-+$/);
	    next if (substr($tseq, $start-1, 1) eq "-" || substr($tseq, $end, 1) eq "-");
	    ## calculate distance between aligned residues
	    # get residues indices in pdb-file of residues right before and after the gap (in query)
	    my $tseqWoGaps = $tseq;
	    $tseqWoGaps =~ s/-//g;
	    $tseqWoGaps =~ s/\*$//;
	    my $startIdx = $tpdb->get_startIdx_of_seq($tseqWoGaps);
	    print "startIdx=$startIdx\n";
	    my $prefix = substr($tseq, 0, $start);
	    if ($v >= 1) { print "prefix=$prefix\n"; }
	    $prefix =~ s/-//g;
	    my $offset = length($prefix);
	    my $res1Idx = $startIdx + $offset;
	    my $res2Idx = $res1Idx + $end - $start + 1;
	    if ($v >= 1) {
		print "resid1=$res1Idx, resid2=$res2Idx\n";
	    }
	    my $dist = $tpdb->distance_between($res1Idx, $res2Idx);
	    if ($v >= 1) { print "dist=$dist\n"; }
	    if ($dist >= 7) {
		$numDisplacedGaps++;
	    }
	}
    }
    return $numDisplacedGaps;
}


## are characters at pos1 and pos2 aligned in seq1 and seq2, i.e.
## are they all residues or gaps
sub is_aligned {
    my $self = shift;
    my $pos1 = shift;
    my $pos2 = shift;
    my $seq1 = shift;
    my $seq2 = shift;

    my $s1p1 = substr($self->[$seq1]->{seq}, $pos1, 1);
    my $s1p2 = substr($self->[$seq1]->{seq}, $pos2, 1);
    my $s2p1 = substr($self->[$seq2]->{seq}, $pos1, 1);
    my $s2p2 = substr($self->[$seq2]->{seq}, $pos2, 1);

    if ($s1p1 eq '-' || $s1p2 eq '-' || $s2p1 eq '-' || $s2p2 eq '-') {
	return 0;
    }
    return 1;
}


sub num_aligned_contacts {
    my $self = shift;
    my $pdbdir = shift;
    my $minSeqDist = shift || 10;  # minimum number of "between residues" for a contact
    my $contactThreshold = shift || 8; # max distance between residues to be a contact

    my $aliLen = $self->get_ali_length();
    my $numTemplates = $self->size() -1;
    my $num_contacts = 0;

    my $v = 0;
    
    ## do contact check for pairs Q-T1, Q-T2, ..., T1-T2, T1-T3, ... 
    for (my $i=1; $i<$self->size(); $i++) {
	for (my $j=1; $j<$self->size(); $j++) {
	    ## read pdb-file
	    my $pdbFile = "$pdbdir/" . $self->[$j]->{id} . ".pdb";
	    if ($v >= 1) { print "i=$i, j=$j, pdb=" . $self->[$j]->{id} . ".pdb\n"; }
	    my $pdb = PdbFile->new();
	    $pdb->read($pdbFile);
	    my $seqWoGaps = $self->[$j]->{seq};
	    $seqWoGaps =~ s/\*$//;
	    $seqWoGaps =~ s/-//g;
	    
	    ## residue-nr of first occurence of seqWoGaps in pdb-file
	    my $startIdx = $pdb->get_startIdx_of_seq($seqWoGaps);
	    if ($v >= 1) { print "startIdx=$startIdx\n"; }
	    ## test residues in pairwise alignment, if they are aligned and have right distance (in sequence)
	    ## => calculate distance and see if they are in contact
	    for (my $k=0; $k<$aliLen-$minSeqDist; $k++) {
		for (my $l=$k+$minSeqDist+1; $l<$aliLen; $l++) {
		    if ($v >= 2) { print "k=$k, l=$l, i-1=" . ($i-1) . ", j=$j\n"; }
		    next if (not $self->is_aligned($k, $l, $i-1, $j));
		    ## get prefixes of aligned residues
		    my $prefix1 = substr($self->[$j]->{seq}, 0, $k+1);
		    my $prefix2 = substr($self->[$j]->{seq}, 0, $l);
		    ## test whether there are gaps in between, i.e. they are distant enough in sequence
		    ## get substring between aligned residues => remove gaps => check length
		    my $between = $prefix2;
		    $between =~ s/.$//;
		    $between =~ s/^$prefix1//;
		    $between =~ s/-//g;
		    if ($v>=2) { print "between=$between\n"; }
		    next if (length($between) < $minSeqDist);
		    ## calculate residue numbers in pdb-file => remove gaps from prefixes
		    $prefix1 =~ s/-//g;
		    $prefix2 =~ s/-//g;
		    my $offset1 = length($prefix1);
		    my $offset2 = length($prefix2);
		    my $res1Idx = $startIdx + $offset1;
		    my $res2Idx = $startIdx + $offset2;
		    if ($v>=2) { print "res1Idx=$res1Idx, res2Idx=$res2Idx\n"; }
		    my $dist = $pdb->distance_between($res1Idx, $res2Idx);
		    if ($v>=2) { print "dist=$dist\n"; }
		    if ($dist < $contactThreshold) {
			$num_contacts++;
		    }
		}
	    }
	}
    }
    my $numOfPairs = $numTemplates*($numTemplates-1)/2;
    
    return $num_contacts/$numOfPairs;
}


## number of residues (in query + all templates)
sub num_residues {
    my $self = shift;

    my $numRes = 0;
    for (my $i=0; $i<$self->size(); $i++) {
	my $seq = $self->[$i]->{seq};
	$seq =~ s/\*$//;
	$seq =~ s/-//g;
	$numRes += length($seq);
    }
    return $numRes;
}


sub num_residues_woQuery {
    my $self = shift;
    my $query = $self->[0]->{seq};
    $query =~ s/\*$//;
    $query =~ s/-//g;
    my $resInQuery = length($query);
    
    my $totNumRes = $self->num_residues();
    return $totNumRes - $resInQuery;
}


sub get_ali_length {
    my $self = shift;
    return length($self->[0]->{seq})-1;
}


sub mean_sd_template_len {
    my $self = shift;
    
    my @lengths;
    for (my $i=1; $i<$self->size(); $i++) {
    }
}

sub score_gaps {
    my $self = shift;
    my $gop = shift || 10;  # gap open penalty
    my $gep = shift || 1;   # gap extension penalty

    my $aliLen = $self->get_ali_length();
    
    my @gaps;
    for (my $i=0; $i<$self->size(); $i++) {
	my @tmp = &match_all_positions("-+", $self->[$i]->{seq});
	push(@gaps, \@tmp);
    }

    my $gapScore = 0;
    for (my $i=0; $i<@gaps; $i++) {
	my @gapsInSeq = @{$gaps[$i]};
	for (my $j=0; $j<@gapsInSeq; $j++) {
	    my $start = $gapsInSeq[$j]->[0];
	    my $end = $gapsInSeq[$j]->[1];
	    next if ($start == 0 || $end == $aliLen); # skip terminal gaps
	    $gapScore += $gop + ($end-$start)*$gep;
	}
    }
    return $gapScore;
}


sub mean_sd_gap_len {
    my $self = shift;
    
    my @gapLens;
    my $mean = 0;
    my $sd = 0;
    my $aliLen = $self->get_ali_length();

    for (my $i=0; $i<$self->size(); $i++) {
	my @gaps = &match_all_positions("-+", $self->[$i]->{seq});
	
	for (my $j=0; $j<@gaps; $j++) {
	    my $start = $gaps[$j]->[0]; 
	    my $end = $gaps[$j]->[1];
	    next if ($start == 0 || $end == $aliLen);
	    push(@gapLens, $end-$start);
	}
    }
    if (scalar(@gapLens) == 0) { return (0,0); }
    $mean = &mean(@gapLens);
    $sd = &sd(@gapLens);

    return ($mean, $sd);
}


## calculate coverage of query by templates
## a residue is covered, if a template residue is aligned to it
## gaps in query are ignored (not counted)
## "coverage" is the coverage by all templates
## "coverageIncrease" is the increase in coverage by the last template in hhmakemodel list (-m option)
sub coverage {
    my $self = shift;
    my $lastTemplIdx = shift;                                                      
                          # if coverage of lastly added template is wanted:
                          # hhmakemodel option -m 4 3 creates a pir file with templates in order 3 4 (ascending), 
                          # but last added template is 3 => lastTemplIdx = 1 (i.e. first template in pir file (0 is query))

    my $coverageByLastTemplate = 0;
    my @coverage;

    ## initialization
    for (my $i=0; $i<$self->get_ali_length(); $i++) {
	$coverage[$i] = 0;
    }

    my $query = $self->[0]->{seq};
    $query =~ s/\*$//;
    my @qRes = split(//, $query);
    ## for all but the very last template
    for (my $i=1; $i<$self->size(); $i++) {
	next if ($i == $lastTemplIdx);
	my $tseq = $self->[$i]->{seq};
	$tseq =~ s/\*$//;
	my @tRes = split(//, $tseq);

	for (my $j=0; $j<@tRes; $j++) {
	    ## both query and template must not be a gap
	    if ($tRes[$j] ne '-' && $qRes[$j] ne '-') {
		$coverage[$j] = 1;
	    }
	}
    }
    my $qLen = $query;
    $qLen =~ s/-//g;
    $qLen = length($qLen);
    #print "qLen=$qLen\n";
        
    ## additional coverage of very last template
    my $newlyCovered = 0;
    my $lastTemplate = $self->[$lastTemplIdx]->{seq};
    $lastTemplate =~ s/\*$//;
    my @lastRes = split(//, $lastTemplate);
    
    for (my $i=0; $i<@lastRes; $i++) {
	if ($lastRes[$i] ne '-' && $coverage[$i] != 1 && $qRes[$i] ne '-') {
	    $newlyCovered++;
	    $coverage[$i] = 1;
	}
    }
    my $coverageIncrease = $newlyCovered/$qLen;
    my $coveredResidues = 0;
    
    ## total coverage
    for (my $i=0; $i<@coverage; $i++) {
	if ($qRes[$i] ne '-' && $coverage[$i] == 1) {
	    $coveredResidues++;
	}
    }
    my $coverage = $coveredResidues/$qLen;

    return ($coverage, $coverageIncrease);
}





## calculate number of newly covered residues for each hitNo
sub allCoverages {
    my $self = shift;
    my @hitNos = @_;  # hit numbers in the order given by template selection strategy                                                    
                          # if coverage of lastly added template is wanted:
                          # hhmakemodel option -m 4 3 creates a pir file with templates in order 3 4 (ascending), 
                          # but last added template is 3 => lastTemplIdx = 1 (i.e. first template in pir file (0 is query))

    my @coverage;
    my %newlyCovered;

    ## initialization
    for (my $i=0; $i<$self->get_ali_length(); $i++) {
	$coverage[$i] = 0;
    }

    
    my $query = $self->[0]->{seq};
    $query =~ s/\*$//;
    my @qRes = split(//, $query);

    ## for all templates
    for (my $i=0; $i<@hitNos; $i++) {
	my $hitNo = $hitNos[$i];
	## find index in pir file
	my $pirLineIdx = &rankForHitNo($hitNo, @hitNos);
	print "pirLineIdx($hitNo)=$pirLineIdx\n";
	$newlyCovered{$hitNo} = 0;

	## get sequence
	my $tseq = $self->[$pirLineIdx]->{seq};
	$tseq =~ s/\*$//;
	my @tRes = split(//, $tseq);

	## check its coverage
	for (my $j=0; $j<@tRes; $j++) {
	    ## both query and template must not be a gap
	    if ($tRes[$j] ne '-' && $qRes[$j] ne '-') {
		if ($coverage[$j] != 1) {
		    $newlyCovered{$hitNo}++;
		}
		$coverage[$j] = 1;
	    }
	}
    }
    

    return %newlyCovered;
}


sub rankForHitNo {
    my $hitNo = shift;
    my @allHitNos = @_;
    
    my $rank = 1;
    for (my $i=0; $i<@allHitNos; $i++) {
	if ($hitNo > $allHitNos[$i]) {
	    $rank++;
	}
    }
    return $rank;
}

1;
