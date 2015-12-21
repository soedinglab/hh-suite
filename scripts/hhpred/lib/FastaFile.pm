package FastaFile;

use PirFile;

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
    my $comment = "";
    my $id;

    my $idxEntry = 0;

    for (my $i=0; $i<@lines; $i++) {
	my $curLine = $lines[$i];       

	if ($curLine =~ /^\s*>(\S+)/) {
	    if ($readSeq == 1) {
		$self->[$idxEntry] = {id=>$id, seq=>$seq, comment=>$comment};
		$idxEntry++;
	    }

	    $id = $1;
	    $comment = "";

	    if ($curLine =~ /^\s*>(\S+)\s+(.+)/) {
		$comment = $2;
	    }	    
		
	    $readSeq = 1;
	    $seq = "";
	    next;
	}

	if ($readSeq == 1) {
	    chomp($curLine);
	    $seq .= $curLine;
	}
    }
    $self->[$idxEntry] = {id=>$id, seq=>$seq, comment=>$comment};
}


sub size {
    my $self = shift;
    scalar(@{$self});
}


sub to_string {
    my ($self) = @_;

    my $res = "";

    for (my $i=0; $i<$self->size(); $i++) {
	$res .= ">" . $self->[$i]->{id} . " " . $self->[$i]->{comment} . "\n";
	$res .= $self->[$i]->{seq} . "\n";
    }	

    $res;
}




sub print {
    my $self = shift;

    print $self->to_string();
}


sub write_to_file {
    my ($self, $filename) = @_;

    open(OH, "> $filename") or die ("Cant write to $filename: $!\n");
    print OH $self->to_string();
    close(OH);
}


sub add_entry {
    my ($self, $id, $comment, $seq) = @_;

    my $idx = $self->size();

    $self->[$idx]->{id} = $id;
    $self->[$idx]->{seq} = $seq;
    $self->[$idx]->{comment} = $comment;
}


sub get_seq {
    my ($self, $idx) = @_;
   
    $self->[$idx]->{seq};
}


sub get_comment {
    my ($self, $idx) = @_;

    $self->[$idx]->{comment};
}


sub to_pir {
    my ($self) = @_;
    
    my $pir = PirFile->new();
    
    for (my $i=0; $i<$self->size(); $i++) {
	my $seq = $self->[$i]->{seq};
	$seq .= "*";
	$pir->add_entry($self->[$i]->{id}, $self->[$i]->{comment}, $seq);
    }

    return $pir;
}

1;
