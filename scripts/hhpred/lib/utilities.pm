#!/usr/bin/perl -w
package utilities;

use strict;
use config;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(System normalizeToOne sigmoid normalizeToOneHash Three2OneLetter One2ThreeLetter getTMscoreFrom BuildTabFiles BuildSingleTabFile max min KMP match_all_positions remove_ranges mean sd sample euklid_dist verbose get_dirname get_basename sumProbLen get_PDB_chain get_seq_len get_HMM_len trim ltrim rtrim TMalignBetween TMscoreBetween TMalignIDBetween hashToStr get_neff_from_hhm getSSPredFromHHM getCoverageApprox randBetween printMatrix printHash symmetricMinMatrix logistic scalarProduct PosteriorsFromTabFile getRandomString);

my $config = HHpredConfig->instance();

sub randBetween {
    my $lower = shift;
    my $upper = shift;

    my $rnum = int($lower + rand($upper-$lower+1));
    return $rnum;
}
    
sub isInRange {
    my $number = shift;
    my $lower = shift;
    my $upper = shift;

    if ($number >= $lower && $number <= $upper) { return 1; }
    return 0;
}


sub logistic {
    my $x = shift;
    return 1.0/(1.0 + exp(-$x));
}


sub scalarProduct {
    my ($xPtr, $yPtr) = @_;
    return -1 if (scalar(@$xPtr) != scalar(@$yPtr));
    
    my $result = 0;
    for (my $i=0; $i<scalar(@$xPtr); $i++) {
	$result += $xPtr->[$i] * $yPtr->[$i];
    }
    return $result;
}

## normalize entries in array to one
## assumes that all entries are non-negative
sub normalizeToOne {
    my @a = @_;

    my $sum = 0;

    for (my $i=0; $i<@a; $i++) {
	$sum += $a[$i];
    }

    if ($sum == 0) {
	print "normalizeToOne: Warning: normalizer equals zero!\n";
	return @a;
    }

    for (my $i=0; $i<@a; $i++) {
	$a[$i] /= $sum;
    }

    return @a;
}


sub printMatrix {
    my $matrixPtr = shift;
    my $matrixName = shift || "";

    my @matrix = @{$matrixPtr};

    if ($matrixName ne "") {
	print "$matrixName:\n";
    }

    for (my $i=0; $i<scalar(@matrix); $i++) {
	for (my $j=0; $j<scalar(@{$matrix[$i]}); $j++) {
	    print "$matrix[$i][$j] ";
	}
	print "\n";
    }
}


sub printHash {
    my $hashPtr = shift;
    my $inRow = 0;
    if (defined $_[1]) { $inRow = 1; }

    my %hash = %{$hashPtr};

    foreach my $key (sort keys %hash) {
	print "$key => $hash{$key}";
	($inRow == 0) ? print "\n" : print " ";
    }
}


## transform matrix into a symmetric one by always
## taking the minimum of two corresponding entries
sub symmetricMinMatrix {
    my $matrixPtr = shift;

    for (my $i=0; $i<scalar(@{$matrixPtr}); $i++) {
	for (my $j=$i+1; $j<@{$matrixPtr->[$i]}; $j++) {
#	    next if ($i == $j);
	    $matrixPtr->[$j][$i] = $matrixPtr->[$i][$j] if ($matrixPtr->[$i][$j] < $matrixPtr->[$j][$i]);
	    $matrixPtr->[$i][$j] = $matrixPtr->[$j][$i] if ($matrixPtr->[$j][$i] < $matrixPtr->[$i][$j]);
	}
    }
}

sub normalizeToOneHash {
    my $hashref = shift;

    my %hash = %$hashref;

    my $sum = 0;

    foreach my $key (keys %hash) {
	$sum += $hash{$key};
    }

    if ($sum == 0) {
	print "normalizeToOneHash: Warning: normalizer equals zero!\n";
	return %hash;
    }

    foreach my $key (keys %hash) {
	$hash{$key} /= $sum;
    }

    return %hash;
}


sub sigmoid {
    my $val = shift;
    
    if ($val < -15.0) {
	return 0.0;
    }
    elsif ($val > 15.0) {
	return 1.0;
    }
    else {
	return (1.0 / (1.0 + exp(-$val)));
    }
}



##################################################################################
# Convert three-letter amino acid code into one-letter code
##################################################################################
sub Three2OneLetter {
    my $res = uc($_[0]);

    if    ($res eq "GLY") {return "G";}
    elsif ($res eq "ALA") {return "A";}
    elsif ($res eq "VAL") {return "V";}
    elsif ($res eq "LEU") {return "L";}
    elsif ($res eq "ILE") {return "I";}
    elsif ($res eq "MET") {return "M";}
    elsif ($res eq "PHE") {return "F";}
    elsif ($res eq "TYR") {return "Y";}
    elsif ($res eq "TRP") {return "W";}
    elsif ($res eq "ASN") {return "N";}
    elsif ($res eq "ASP") {return "D";}
    elsif ($res eq "GLN") {return "Q";}
    elsif ($res eq "GLU") {return "E";}
    elsif ($res eq "CYS") {return "C";}
    elsif ($res eq "PRO") {return "P";}
    elsif ($res eq "SER") {return "S";}
    elsif ($res eq "THR") {return "T";}
    elsif ($res eq "LYS") {return "K";}
    elsif ($res eq "HIS") {return "H";}
    elsif ($res eq "ARG") {return "R";}
    elsif ($res eq "ASX") {return "D";}
    elsif ($res eq "GLX") {return "E";}
    elsif ($res eq "MSE") {return "M";} # SELENOMETHIONINE 
    elsif ($res eq "SEP") {return "S";} # PHOSPHOSERINE 
    elsif ($res eq "SEC") {return "C";} # SELENOCYSTEINE
    elsif ($res eq "TPO") {return "T";} # PHOSPHOTHREONINE 
    elsif ($res eq "TYS") {return "Y";} # SULFONATED TYROSINE 
    elsif ($res eq "KCX") {return "K";} # LYSINE NZ-CARBOXYLIC ACID
    else                  {return "X";}
}

    
sub One2ThreeLetter {
    my $res = uc($_[0]);

    if    ($res eq "G") {return "GLY";}
    elsif ($res eq "A") {return "ALA";}
    elsif ($res eq "V") {return "VAL";}
    elsif ($res eq "L") {return "LEU";}
    elsif ($res eq "I") {return "ILE";}
    elsif ($res eq "M") {return "MET";}
    elsif ($res eq "F") {return "PHE";}
    elsif ($res eq "Y") {return "TYR";}
    elsif ($res eq "W") {return "TRP";}
    elsif ($res eq "N") {return "ASN";}
    elsif ($res eq "D") {return "ASP";}
    elsif ($res eq "Q") {return "GLN";}
    elsif ($res eq "E") {return "GLU";}
    elsif ($res eq "C") {return "CYS";}
    elsif ($res eq "P") {return "PRO";}
    elsif ($res eq "S") {return "SER";}
    elsif ($res eq "T") {return "THR";}
    elsif ($res eq "K") {return "LYS";}
    elsif ($res eq "H") {return "HIS";}
    elsif ($res eq "R") {return "ARG";}
    elsif ($res eq "U") {return "SEC";}
    elsif ($res eq "B") {return "ASX";}
    elsif ($res eq "Z") {return "GLX";}
    else                {return "UNK";}
}


## get (approximated) coverage of query by template at index idx in tlist
sub getCoverageApprox {
    my $tlist = shift;
    my $idx = shift;

    if ($idx < 0 || $idx >= $tlist->get_queryLength()) {
	return -1;
    }
    my $range = $tlist->get($idx)->get_Qend() - $tlist->get($idx)->get_Qstart();
    ## gaps within this region are not considered => approximation
    my $coverageApprox = $range / $tlist->get_queryLength();

    return $coverageApprox;
}


## extract TMscore from TMscore/TMalign-output file
sub getTMscoreFrom {
    my $file = shift;

    my $tmscore = -1;

    open (TM, "< $file") or die "Cant open $file";
    while(my $line = <TM>) {
	if ($line =~ /TM-score\s*=\s*(\d+(\.\d+)?)/) {
	    $tmscore = $1;
	    last;
	}
    }
    close(TM);

    return $tmscore;
}
 

## TMalign between two given structures
sub TMalignBetween {
    my $strucFileOne = shift;
    my $strucFileTwo = shift;
    my $options = shift || " ";
    my $TMalign = shift || $config->get_TMalign();

    my $TMalignOutput = `$TMalign $strucFileOne $strucFileTwo $options`;
    my $TMscore = -1;
    if ($TMalignOutput =~ /TM-score\s*=\s*(\S+),/) {
	$TMscore = $1;
    }
    if ($TMscore == -1) {
	print "WARNING: TMalignBetween could not find TMscore in TMalign output!\n";
    }

    return $TMscore;
}


## TMalign (TMscore and sequence identity) between two given structures
sub TMalignIDBetween {
    my $strucFileOne = shift;
    my $strucFileTwo = shift;
    my $options = shift || " ";
    my $TMalign = shift || $config->get_TMalign();

    my $TMalignOutput = `$TMalign $strucFileOne $strucFileTwo $options`;
    my $TMID = -1;
    my $TMscore = -1;
    if ($TMalignOutput =~ /TM-score\s*=\s*(\S+),\s*ID\s*=\s*(\S+)/) {
	$TMscore = $1;
	$TMID = $2;
    }
    if ($TMscore == -1) {
	print "WARNING: TMalignBetween could not find TMscore in TMalign output!\n";
    }

    return ($TMscore, $TMID);
}


## TMalign (TMscore and sequence identity) between two given structures
sub TMscoreBetween {
    my $modelFile = shift;
    my $nativeFile = shift;
    my $options = shift || " ";
    my $TMscore = shift || $config->get_TMscore();

    my $TMscoreOutput = `$TMscore $modelFile $nativeFile $options`;
    my $score = -1;
    if ($TMscoreOutput =~ /TM-score\s*=\s*(\d+(\.\d+)?)/) {
	$score = $1;
    }
    if ($score == -1) {
	print "WARNING: TMscoreBetween could not find TMscore in TMscore output!\n";
    }

    return $score;
}

##################################################################################
## tabfile: tabfile containing all tabs-results (created by hhsearch options -atab)
## outbase: where to write the separate tab-files
## maxhits: how many hits to write (def: in fact all)
## 
## the input tabfile contains tab-entries (i.e. i j sim probab) for each template
## of the initial hhsearch. The order is as in initial hhr file.
##
## the created tab files (containing posteriori-prob-
## abilities) are saved in outbase.HITtemplateStartStop.tab,
## where start: first residue
##       stop: last residue
## it might be that a template is aligned more than one times at different
## positions. One can 
##################################################################################
sub BuildTabFiles {
    my $tabfile = shift;
    my $outbase = shift;
    my $maxHits = shift;

    $maxHits = defined($maxHits) ? $maxHits : 1000;


    open (TH, "< $tabfile") or die "Cant open $tabfile: $!\n";

    my $hitnr = 1;

    while (my $line = <TH>) {

	next if ($line =~ /^\s*i\s+j/);

	## new template
	if ($line =~ />(\S+)/) {
	    if ($hitnr > 1) { close(HH); }


	    ## write a new tabfile
	    my $singleTabFile = "$outbase.$1.tab";

	    if ($hitnr > $maxHits) { last; }
	    $hitnr++;

	    open (HH, "> $singleTabFile") or die "Cant open $singleTabFile: $!\n";
	    next;
	}

	## i j score ss probab [dssp]
	if ($line =~ /^\s*\S+\s+\S+\s+\S+\s+\S+\s+\S+/) {
	    print (HH $line);
	}
    }

    close(TH);
}


############################################################################
## see subroutine BuildTabFiles
## this one creates only a single tab-file - the hitnr-th in the tabfile
############################################################################
sub BuildSingleTabFile {
    my $tabfile = shift;
    my $hitnr = shift;
    my $outbase = shift;

    open (TH, "< $tabfile") or die "Cant open $tabfile: $!\n";

    my $hit = 1;
    my $template = "";
    my $found = 0;

    while (my $line = <TH>) {

	next if ($line =~ /^\s*i\s+j/);

	## begin of a new template
	if ($line =~ />(\S+)/) {

	    ## already found => stop
	    if ($found == 1) {
		close (HH);
		last;
	    }

	    ## found
	    if ($found == 0 and $hitnr == $hit) {
		$found = 1;
		$template = $1;
		open (HH, "> $outbase.$1.HIT$hitnr.tab") or die "Cant open $outbase.$1.HIT$hitnr.tab: $!\n";
		next;
	    }

	    $hit++;
	    
	}

	## not yet found
	if ($found == 0) {
	    next;
	}
	## found
	else {
	    ## i j score ss probab [dssp]
	    if ($line =~ /^\s*\S+\s+\S+\s+\S+\s+\S+\s+\S+/) {
		print (HH "$line");
	    }
	}
    }

    close(HH);
    close (TH);
    return $found;
}


############################################################################
## see subroutine BuildTabFiles
## returns hash for one hit with query-residue => pp entries
############################################################################
sub PosteriorsFromTabFile {
    my $tabfile = shift;
    my $hitnr = shift;

    open (TH, "< $tabfile") or die "Cant open $tabfile: $!\n";

    my $hit = 1;
    my $template = "";
    my $found = 0;
    my %QidxToPP;

    while (my $line = <TH>) {
	next if ($line =~ /^\s*i\s+j/);
	## begin of a new template
	if ($line =~ />(\S+)/) {
	    ## already found => stop
	    if ($found == 1) {
		last;
	    }
	    ## found
	    if ($found == 0 and $hitnr == $hit) {
		$found = 1;
		$template = $1;
		next;
	    }
	    $hit++;	   
	}
	## not yet found
	if ($found == 0) {
	    next;
	}
	## found
	else {
	    ## i j score ss probab [dssp]
	    if ($line =~ /^\s*(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)/) {
		$QidxToPP{$1} = $2;
	    }
	}
    }
    close (TH);
    
    return %QidxToPP;
}



sub max {
    my $max = shift;
    foreach (@_) {
	$max = $_ if ($_ > $max);
    }
    return $max;
}


sub min {
    my $min = shift;
    foreach (@_) {
	$min = $_ if ($_ < $min);
    }
    return $min;
}


sub mean {
    my @array = @_;
    my $sum = 0;
    
    for (my $i=0; $i<@array; $i++) {
	$sum += $array[$i];
    }
    return $sum/scalar(@array);
}


sub sd {
    my @array = @_;
    my $mean = &mean(@array);
    my $N = scalar(@array);

    my $var = 0;
    for (my $i=0; $i<@array; $i++) {
	$var += ($array[$i] - $mean)*($array[$i] - $mean)
    }
    $var *= 1.0/$N;
    
    return sqrt($var);
}

## sample a bin where probability of each bin is given in "probs"
## sum(probs) must be 1
sub sample {
    my @probs = @_;
    
    my $rand = rand();

    my $sum = 0;
    for (my $bin=0; $bin<@probs; $bin++) {
	$sum += $probs[$bin];
	if ($rand <= $sum) {
	    return $bin;
	}
    }
    return 0;
}


sub hashToStr {
    my $hashPtr = shift;

    my %myhash = %{$hashPtr};
    my $str = "";

    foreach my $key (sort keys(%myhash)) {
	$str .= "$key=$myhash{$key}\n";
    }

    return $str;
}


## calculate number of residues on which calculation
## of sumProbs is based (the ones which are aligned and
## which have dssp)
sub sumProbLen {
    my $ss_dssp = shift;
    my $conf = shift;

    my $len = length($ss_dssp);
    if (length($ss_dssp) != length($conf)) {
	print "WARNING: sumProbLen length(ss_dssp) != length(conf)!\n";
	$len = &min(length($ss_dssp), length($conf));
    }

    my @ssDsspTok = split(//, $ss_dssp);
    my @confTok = split(//, $conf);

    my $sumProbLen = 0;

    for (my $i=0; $i<$len; $i++) {
	if ($ssDsspTok[$i] ne '-' && $confTok[$i] ne " ") {
	    $sumProbLen++;
	}
    }
    return $sumProbLen;
}



sub euklid_dist {
    my $ref1 = shift;
    my $ref2 = shift;
    my $v = shift || 2;

    my @vec1 = @$ref1;
    my @vec2 = @$ref2;

    if ($v>=2) {
	if ($#vec1 != $#vec2) {
	    print "ERROR: euklid_dist: vec1 and vec2 differ in length!\n";
	}
    }

    my $sum = 0;

    for (my $i=0; $i<@vec1; $i++) {
	$sum += ($vec1[$i] - $vec2[$i]) * ($vec1[$i] - $vec2[$i]);
    }

    return (sqrt($sum));
}


## Knuth-Morris-Pratt algorithm 
## returns index of first occurence of ss in st
## or -1 otherwise
sub KMP {
    my $st = shift; ## text
    my $ss = shift; ## search string
    my $cs = shift || 0; ## case sensitivity

    if ($cs != 0) {
	$st = uc($st);
	$ss = uc($ss);
    }

    my @t = split(//, $st);
    my @s = split(//, $ss);

    my $n = scalar(@t);
    my $m = scalar(@s);

    ## compute borders
    my @borders;
    $borders[0] = -1;
    my $i = 0;
    $borders[1] = 0;

    for (my $j=2; $j<=$m; $j++) {
	while(($i>=0) && ($s[$i] ne $s[$j-1])) {
	    $i = $borders[$i];
	}
	$i++;
	$borders[$j] = $i
    }

    ## search routine
    $i = 0;
    my $j = 0;

    while ($i <= $n - $m) {
	while($t[$i+$j] eq $s[$j]) {
	    $j++;
	    if ($j == $m) {
		return $i;
	    }
	}
	$i = $i + ($j - $borders[$j]);
	$j = &max(0, $borders[$j]);
    }
    return -1;
}


## given a string and a regex,
## give back all positions (start and end) where regex matches string
sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [ $-[0], $+[0] ];
    }
    return @ret;
}


sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Left trim function to remove leading whitespace
sub ltrim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    return $string;
}

# Right trim function to remove trailing whitespace
sub rtrim($) {
    my $string = shift;
    $string =~ s/\s+$//;
    return $string;
}




## setsPtr is a ptr to an array as
## generated by match_all_positions, i.e. an array of arrays
## with 2 elements (start, end)
## and returns start if: start <= idx < end
## or -1 if no such array is found
sub set_of_idx {
    my $idx = shift;
    my $setsPtr = shift;

    my @sets = @$setsPtr;
    for (my $i=0; $i<@sets; $i++) {
	my $start = $sets[$i]->[0];
	my $end = $sets[$i]->[1];

	return ($start, $end) if ($start >= $idx && $end < $idx);
    }
    return (-1,-1);
}
	   


## remove ranges given in splits (generated by e.g. match_all_positions)
## from str and give back new 
sub remove_ranges {
    my $str = shift;
    my @splits = @_;

    my $result = "";
    
    my $start = 0;
    my $end = 0;
    for (my $i=0; $i<@splits; $i++) {
	my $gStart = $splits[$i]->[0];
	my $gEnd = $splits[$i]->[1];
	
	$end = $gStart;
	$result .= substr($str, $start, $end-$start);
	$start = $gEnd;
    }
    $result .= substr($str, $start);
    return $result;
}


sub get_basename {
    my $dirbasename = shift;
    $dirbasename =~ /^.*\/(\S+?)(\.\S+)?$/;
    my $basename = $1;
    return $basename;
}


sub get_dirname {
    my $dirbasename = shift;
    $dirbasename =~ /^(.*)\//;
    my $dirname = $1;
    return $dirname;
}


## very simple: assumes name to look like 1a7j_A
sub get_PDB_chain {
    my $name = shift;
    
    if ($name =~ /\S+\_(\S+)$/) {
	return $1;
    } else {
	print "WARNING utitlities.pm get_PDB_chain strange format!\n";
	return "";
    }
}



sub get_seq_len {
    my $seq = shift;
    chomp($seq);
    $seq =~ s/[\*-]//g;
    return length($seq);
}


sub get_neff_from_hhm {
    my $hhmFile = shift;
    my $neff = -1;

    open(HH, "< $hhmFile") or die ("Cant open $hhmFile: $!\n");
    while(my $line = <HH>) {
	if ($line =~ /^Neff\s+(\S+)/i) {
	    $neff = $1;
	    last;
	}	    
    }
    close(HH);

    return $neff;
}


sub get_HMM_len {
    my $hhmFile = shift;
    my $len = -1;

    open(HH, "< $hhmFile") or die ("Cant open $hhmFile: $!\n");
    while(my $line = <HH>) {
	if ($line =~ /^LENG\s+(\S+)/i) {
	    $len = $1;
	    last;
	}	    
    }
    close(HH);

    return $len;
}


sub verbose {
    my $level = shift;
    my $actLevel = shift;
    my $message = shift;

    if ($actLevel >= $level) { print "$message\n"; }
}


sub getSSPredFromHHM {
    my $hhmFile = shift;
    open(HHM, "< $hhmFile") or die "Cant open $hhmFile";
    my $ssFound = 0;
    my $sspred = "";

    while(my $line = <HHM>) {
	chomp($line);
	if ($line =~ /^>ss\_pred/) {
	    $ssFound = 1;
	    next;
	}
	if ($ssFound && $line =~ /^>/) {
	    last;
	}
	next if ($ssFound == 0);
	$sspred .= $line;
	    
    }
    close(HHM);
    return $sspred;
}


sub getRandomString {
	my $len = shift;
	
	my @chars = ("A".."Z", "a".."z");
	my $string = "";
	$string .= $chars[rand @chars] for 1..$len;
	return $string;
}


sub System {
    my $cmd = shift;
    print "$cmd\n";
    system("$cmd");
}

1;

