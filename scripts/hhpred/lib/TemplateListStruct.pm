#!/user/bin/perl -w
use strict;

package TemplateListStruct;

use config;
use utilities;
use TemplateStruct;

my $config = HHpredConfig->instance();

sub new {
    my ($caller, $filename) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    no strict "refs";
    my $self = bless {}, $class;

    $self->{templates} = [];
    $self->{queryLength} = -1;
    $self->{query} = "";
    $self->{neff} = -1;

    if ($caller_is_obj) {
	my $size = $caller->size();
	for (my $i=0; $i<$size; $i++) {
	    $self->{templates}->[$i] = $caller->{templates}->[$i];
	}
	$self->{queryLength} = $caller->{queryLength};
	$self->{query} = $caller->{caller};
	$self->{neff} = $caller->{neff};
    }

    if (defined($filename)) {
	$self->hhr_to_TemplateList($filename);
    }
    return $self;
}


sub add_template {
    my ($self, $template) = @_;    
    my $curSize = $self->size();
    $self->{templates}->[$curSize] = $template;
}


## before adding template, check whether it is already in list
sub check_and_add {
    my ($self, $template) = @_;

    for (my $i=0; $i<$self->size(); $i++) {
	if ($self->{templates}->[$i]->equals($template)) {
	    return;
	}
    }
    $self->add_template($template);
}
    

sub clear {
    my $self = shift;
    %{$self} = ();
    $self->{templates} = [];
    $self->{query} = "";
    $self->{queryLength} = -1;
    $self->{neff} = -1;
}


## delete hit with No "No"
sub delete_No {
    my $self = shift;
    my $No = shift;

    ## get idx for hit with No "No"
    my $deleteIdx = -1; 
    for (my $i=0; $i<$self->size(); $i++) {
	if ($self->{templates}->[$i]->get_No() == $No) {
	    $deleteIdx = $i;
	    last;
	}
    }
    print "deleting No=$No, idx=$deleteIdx\n";
    if ($deleteIdx != -1) {
	splice(@{$self->{templates}}, $deleteIdx, 1);
    }
}


sub size {
    my $self = shift;
    return scalar(@{$self->{templates}});
}


sub get {
    my ($self, $i) = @_;
    $self->{templates}->[$i];
}


sub get_last {
    my $self = shift;
    $self->{templates}->[$self->size()-1];
}



sub to_string {
    my $self = shift;
    my $spacer = shift;
    my $out = "";
    for (my $i=0; $i<$self->size(); $i++) {
	$out .= $self->{templates}->[$i]->to_string($spacer) . "\n";
    }
    return $out;
}


sub print {
    my $self = shift;
    my $out = $self->to_string();
    print $out;
}


sub to_TemplateList_helper {
    my $self = shift;
    my $hhrFile = shift;
    my @lines = @_;

    my $matchC;
    my $No;
    my $filtnr = "start";  ## filter step (start means no filtering)
    my $spaceLen = 12;

    if ($hhrFile =~ /\.(\d+)\.hhr/) {
	$filtnr = $1;
    }

    for (my $i=0; $i<@lines; $i++) {
	my $line = $lines[$i];

	if ($line =~ /^Match_columns\s*(\S+)/) {
	    $matchC = $1;
	    $self->_set_queryLength($matchC);
	}
	if ($line =~ /^Query\s+(\S+)/) {
	    my $query = $1;
	    $self->_set_query($query);
	}
	if ($line =~ /^Neff\s+(\S+)/) {
	    my $neff = $1;
	    $self->_set_neff($neff);
	}
	## No     Hit       Prob E-val  P-val  Score    SS      Cols  Query(start end) Template(start end) HMM
	elsif ($line=~/^\s*(\d+)\s+(\S+).+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)-(\d+)\s+(\d+)-(\d+)\s*\((\S+)\)$/) {
	    my $No       = $1;
	    my $Hit      = $2; 
	    my $Prob     = $3;	    			
	    my $Eval     = $4;
	    my $Pval     = $5;
	    my $Score    = $6;
	    my $SS       = $7;
	    my $Cols     = $8;
	    my $Qstart   = $9;
	    my $Qend     = $10;
	    my $Tstart   = $11;
	    my $Tend     = $12;
	    my $HMM      = $13;

	    my $SSL      = $SS/$matchC; 
	    $SSL = sprintf("%.4f", $SSL);

	    my $template = TemplateStruct->new(Filt => $filtnr,
					 No => $No,
					 Hit => $Hit,
					 Prob => $Prob,
					 Eval => $Eval,
					 Pval => $Pval,
					 Score => $Score,
					 SS => $SS,
					 Cols => $Cols,
					 Qstart => $Qstart,
					 Qend => $Qend,
					 Tstart => $Tstart,
					 Tend => $Tend,
					 HMM => $HMM);
	    $self->add_template($template);
	}
	elsif($line =~ /^No\s+(\d+)/) {
	    $No = $1;			
	    $line = $lines[++$i];

	    if ($line !~ /^>(\S+)\s/) {
		die("Error:: wrong format in \"$line\"\n");
	    }

	    my $hit = $1;
	    $line = $lines[++$i];

	    if ($line =~ /COMPACTNESS=(\S+)/i) {
		$self->get($No-1)->set_Compactness($1);
	    } else {
		#print "missing compactness in hhr-file!\n";
	    }

	    if ($line =~ /CSS=(\S+)/i) {
		$self->get($No-1)->set_Css($1);
	    } else {
		#print "missing css in hhr-file!\n";
	    }

	    if ($line =~ /CONTACT=(\S+)/i) {
		$self->get($No-1)->set_Contact($1);
	    } else {
		#print "missing contact in hhr-file!\n";
	    }

	    if ($line =~ /CONTACT_REALIGN=(\S+)/i) {
		$self->get($No-1)->set_ContactRealign($1);
	    } else {
		#print "missing contact realign in hhr-file!\n";
	    }

	    if ($line !~ /Similarity=(\S+)\s+Sum_probs=(\S+)\s*/) {
		die("Error: wrong format in \"$line\"\n");
	    }	    

	    my $Similarity = $1;
	    my $SumProbL = $2/$matchC;
	    $SumProbL = sprintf("%.4f" , $SumProbL);	 

	    if ($line =~ /Identities=(\S+)%\s/) {
		$self->get($No-1)->set_Ident($1);
	    }

	    $self->get($No-1)->set_Sim($Similarity);
	    $self->get($No-1)->set_SumProbL($SumProbL);
	}
	elsif ($line =~ /^T\s+ss_dssp(\s+)(\S+)/) {
	    $spaceLen = length($1)-1;
	    my $ss_dssp = $self->get($No-1)->get_ss_dssp();
	    $self->get($No-1)->set_ss_dssp("$ss_dssp" . $2);
	}
	## Confidence line may contain spaces => read number of spaces from ss_dssp line
	elsif ($line =~ /^Confidence\s{$spaceLen}(.*)\n/) {
	    my $conf = $self->get($No-1)->get_conf();
	    $self->get($No-1)->set_conf("$conf" . $1);
	}
    }	
}


sub str_to_TemplateList {
    my $self = shift;
    my $str = shift;

    my @lines;
    @lines = split(/\n/, $str);
    
    $self->to_TemplateList_helper("dummy", @lines);
}


sub hhr_to_TemplateList {
    my ($self, $hhrFile) = @_;

    my @lines;   
    open(HHR,"< $hhrFile") or die("Cant open $hhrFile: $!\n");
    @lines = <HHR>;
    close(HHR);

    $self->to_TemplateList_helper($hhrFile, @lines);
}


sub write_to_file {
    my ($self, $outfile) = @_;
    
    open (OH, "> $outfile") or die("Cant write to $outfile: $!\n");
    my $out = $self->to_string('===');
    print(OH $out);
    close(OH);
}


sub read_from_file {
    my ($self, $infile) = @_;
    my $append = 0;
    ## append template(s) to already existing ones
    $append = 1 if (scalar(@_) > 2 && $_[2] == 1);

    $self->clear() if (! $append);
    open(IH, "< $infile") or die("Cant open $infile: $!\n");
    while(<IH>) {
	chomp;
	if (/(\S+===)+/) {
	    my @entry = split(/===/, $_);
	    my $template = TemplateStruct->new(Filt => $entry[0],
					 No => $entry[1],
					 Hit => $entry[2],
					 Prob => $entry[3],
					 Eval => $entry[4],
					 Pval => $entry[5],
					 Score => $entry[6],
					 SS => $entry[7],
					 Cols => $entry[8],
					 Qstart => $entry[9],
					 Qend => $entry[10],
					 Tstart => $entry[11],
					 Tend => $entry[12],
					 HMM => $entry[13],
					 Ident => $entry[14],
					 Sim => $entry[15],
					 SumProbL => $entry[16],
					 predTM => $entry[17],
					 Compactness => $entry[18],
					 Css => $entry[19],
					 Contact => $entry[20],
					 ContactRealign => $entry[21]);
	    
	    $self->add_template($template);
	}
    }
    close(IH);
}


sub set_queryLength {
    my ($self, $len) = @_;
    $self->{queryLength} = $len;
}

sub get_queryLength {
    my $self = shift;
    $self->{queryLength};
}

sub set_query {
    my ($self, $query) = @_;
    $self->{query} = $query;
}

sub get_query {
    my $self = shift;
    $self->{query};
}

sub get_neff {
    my $self = shift;
    $self->{neff};
}

sub set_neff {
    my ($self, $neff) = @_;
    $self->{neff} = $neff;
}

## for backward compatibility ##
sub _set_queryLength {
    my ($self, $len) = @_;
    $self->{queryLength} = $len;
}

sub _get_queryLength {
    my $self = shift;
    $self->{queryLength};
}

sub _set_query {
    my ($self, $query) = @_;
    $self->{query} = $query;
}

sub _get_query {
    my $self = shift;
    $self->{query};
}

sub _get_neff {
    my $self = shift;
    $self->{neff};
}

sub _set_neff {
    my ($self, $neff) = @_;
    $self->{neff} = $neff;
}
######




sub sort_by_sim {
    my $self = shift;
    @{$self->{templates}} = sort {$b->get_Sim() <=> $a->get_Sim()} @{$self->{templates}};
}


sub sort_by_prob {
    my $self = shift;
    @{$self->{templates}} = sort {$b->get_Prob() <=> $a->get_Prob()} @{$self->{templates}};
}


sub sort_by_sumProbL {
    my $self = shift;
    @{$self->{templates}} = sort {$b->get_SumProbL() <=> $a->get_SumProbL()} @{$self->{templates}};
}


sub sort_by_predTM {
    my $self = shift;
    @{$self->{templates}} = sort {$b->get_predTM() <=> $a->get_predTM()} @{$self->{templates}};
}


sub templateList_to_hhr {
    my $self = shift;
    my $outbase = shift;

    my @hhrContent = ();
    my $hh = $config->get_hh();

    open(HHR, "> $outbase.hhr") or die ("Error in templateList_to_hhr: Cant write $outbase.hhr: $!\n");

    for (my $i=0; $i<$self->size(); $i++) {
	my $template = $self->get($i);

	## open apropriate hhr file (wrt filter step)
	my $infile = "$outbase." . $template->get_Filt() . ".hhr";
	open (IN, "< $infile") or die ("Error: cannot open $infile!\n");

	my $checkedHeader = 0;
	my $begin;
	my $e = 0;
	my $end;							
	my $line;
	my $hitnr = $i+1;

	while ($line = <IN>) {	    
	    ## copy first header lines:
	    if (($checkedHeader==0) && ($i==0) && ($line !~ /^\s*\d+\s+\S+.+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\d+-\d+\s+\S+\s*\(\S+\)$/)) {
		if ($line=~ /^Command/) {
		    $line=~ s/(^Command\s*)(.*)$/$1$hh\/hhsearch artificial hhr file/;
		}

		## replace P-value against TMscore
		if ($line=~ /\s+No\s+Hit\s+Prob\s+E-value\s+P-value\s+Score\s+SS\s+Cols\s+Query\s+HMM\s+Template\s+HMM\s*/) {
		    $line =~ s/(\s*No\s+Hit\s+Prob\s+E-value\s+)(P-value)(\s+Score\s+SS\s+Cols\s+Query\s+HMM\s+Template\s+HMM\s+)/$1TMScore$3/;
		}	
		print (HHR "$line");				
	    }
	    else {
		$checkedHeader = 1;
	    }
	    
	    ## get hit Info:
	    my $No = $template->get_No();
	    if ($line =~ /^\s*$No(\s+\S+.+\s+\S+\s+\S+)\s+\S+(\s+\S+\s+\S+\s+\S+\s+\d+-\d+\s+\S+\s*\(\S+\)$)/)	{		
		## replace P-value by TMScore in hit info
		$line = sprintf("%3s$1  %1.4f$2\n", $hitnr, $template->get_predTM());		
		print (HHR "$line");
		last;
	    }
	}

	## skip all lines up to alignment block
	## Find beginning of alignment and replace hit index by new one
	while ($line = <IN>){ 
	    my $No = $template->get_No();
	    if ($line =~ /^No\s+$No/) {
		last;
	    }
	} 
       
	$line =~ s/^No\s+\d+/No $hitnr/;
	push(@hhrContent, $line);

	## Push alignment block onto array				
	while ($line = <IN>) {
	    if(($line =~ /^No\s/)) {
		last;
	    }
	    if ($line =~ /Done!/) {}
	    else {
		push(@hhrContent, $line);
	    }
	}	
	close (IN);

	## create associated tab file
	&BuildSingleTabFile("$outbase." . $template->get_Filt() . ".tab", $template->get_No(), $outbase);
    }
    print(HHR "\n");
    print(HHR @hhrContent);
    print(HHR "Done!\n");
    close (HHR);
}


## starting from current hhr file, extract some features and save them into resultfile
## this is needed for benchmark set compilation
sub createBenchmarkInfoFile {
    my ($self, $resultFile, $pdbdir) = @_;

    my $TMalign = $config->get_TMalign();
    
    my $query = $self->_get_query();
    my $queryPDB = "$pdbdir/$query.pdb";

    my $res = "";
    $res .= "queryName"."\t"."TMID"."\t"."coverage"."\t"."queryLen"."\t"."templateName"."\t"."TMscore\n";

    ## extract information from max first 50 templates
    for (my $i=0; $i<50 && $i<$self->size(); $i++) {
	my $template = $self->get($i);

	my $TMscore = 0;
	my $TMid = 0;

	my $templatePDB = "$pdbdir/" . $template->get_Hit() . ".pdb";
	my $tmalignResult = `$TMalign $templatePDB $queryPDB`;
	if ($tmalignResult =~ /TM-score\s*=\s*(\S+),\s+ID\s*=\s*(\S+)/) {
	    $TMscore = $1;
	    $TMid= int(($2*100)+0.5);
	}

	my $queryLen = $self->_get_queryLength();
	my $coverage = int(($template->get_Cols()*100/$queryLen)+0.5);
	my $templateName = $template->get_Hit();

	$res .= "$query\t$TMid\t$coverage\t$queryLen\t$templateName\t$TMscore\n";
    }

    open(OH, "> $resultFile") or die "Cant write $resultFile: $!\n";
    print (OH $res);
    close(OH);
}

1;
