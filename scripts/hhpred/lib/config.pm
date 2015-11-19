package HHpredConfig;

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;
use strict;
use vars qw($AUTOLOAD);
use Carp;
use base 'Class::Singleton';

{
    ## default values
    my %_attr_data = (
	_cpus						=> [4, 'read/write'],
	_hhsearch_mact				=> [0.05,	'read/write'],
	_hhblits_mact				=> [0.5,	'read/write'],
	_hhblits_rounds				=> [3, 'read/write'],
	_numberOfGeneratedModels	=> [3, 'read/write'],
	_maxNumOfTemplates			=> [8, 'read/write'],
	_doFiltering				=> [1, 'read/write'],
	_parallelFiltering			=> [0, 'read/write'],
	_replaceDistanceRestraints	=> [1, 'read/write'],
	_multiTemplate				=> [1, 'read/write'],
	_templateWeightStrategy		=> [1, 'read/write'],
	_preselectTemplates			=> [1, 'read/write'],
	_rankTemplates				=> [1, 'read/write'],
	_realignProbcons			=> [0, 'read/write'],
	_assessModel				=> [1, 'read/write'],
	_doParallelModeller			=> [0, 'read/write'],
	
	_hhlib			=> ["$hhlib",					'read'],
	_hhalign		=> ["$hhbin/hhalign",				'read'],
	_hhsearch		=> ["$hhbin/hhsearch",				'read'],
	_hhblits		=> ["$hhbin/hhblits",				'read'],
	_hhmake			=> ["$hhbin/hhmake",				'read'],
	_hhfilter		=> ["$hhbin/hhfilter",				'read'],
	_hhmakemodel	=> ["$hhlib/scripts/hhpred/dependencies/hhmakemodel.pl",	'read'],
	_addss			=> ["$hhscripts/addss.pl",		'read'],
	_multithread	=> ["$hhscripts/multithread.pl",	'read'],

	_TMscore	=> ["$hhlib/scripts/hhpred/bin/TMscore",	'read'],
	_TMalign 	=> ["$hhlib/scripts/hhpred/bin/TMalign",	'read'],

	_repairPDB			=> ["$hhlib/scripts/hhpred/bin/repair_pdb.pl",							'read'],
	_modellerParallel	=> ["$hhlib/scripts/hhpred/bin/modeller9.13/bin/modpy.sh python2.7 ",	'read'],
	_modeller			=> ["$hhlib/scripts/hhpred/bin/modeller9.13/bin/modpy.sh python2.7",	'read'],
	
	_pdbdir			=> ["XXXXX",								'read'],
	_uniprot20		=> ["XXXXX",	'read'],
	
	_MDNWeightsLayer1CACA	=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer1CACAminP.dat", 'read'],
	_MDNWeightsLayer2CACA	=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer2CACAminP.dat", 'read'],
	_MDNWeightsLayer1NO		=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer1NOminP.dat",	'read'],
	_MDNWeightsLayer2NO		=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer2NOminP.dat",	'read'],
	_MDNWeightsLayer1SCMC	=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer1SCMCminP.dat", 'read'],
	_MDNWeightsLayer2SCMC	=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer2SCMCminP.dat", 'read'],
	_MDNWeightsLayer1SCSC	=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer1SCSCminP.dat", 'read'],
	_MDNWeightsLayer2SCSC	=> ["$hhlib/scripts/hhpred/share/neural-net/MDNWeightsLayer2SCSCminP.dat", 'read'],
	);

    sub _accessible {
		my ($self, $attr, $mode) = @_;
		$_attr_data{$attr}[1] =~ /$mode/;
    }

    sub _default_for {
		my ($self, $attr) = @_;
		$_attr_data{$attr}[0];
    }

    sub _standard_keys {
		keys %_attr_data;
    }
}

my $matchKeyValue = qr/^\s*(\S+)\s*=\s*(.+?)\s*$/;

## configFile is contains entries of the form
## key = value
sub _new_instance  {
    my ($caller) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = bless {}, $class;
    
    my %argsInFile;
    my $configFile = $ENV{'HHPRED_CONFIG'};
    if (defined($configFile) && $configFile ne "") {
	open(CH, "< $configFile") or die("Cant open $configFile: $!\n");
	while(my $line = <CH>) {
	    next if ($line =~ /^\s*#/); ## skip comments
	    if ($line =~ /$matchKeyValue/) {
		my $key = $1;
		my $value = $2;
		$argsInFile{$key} = $value;
	    }
	}
	close(CH);
    }

    foreach my $attrname ( $self->_standard_keys() ) {
	my ($argname) = ($attrname =~ /^_(.*)/);
	## take value from file (if available)
	if (exists $argsInFile{$argname}) {
	    $self->{$attrname} = $argsInFile{$argname};	
	} elsif ($caller_is_obj) { ## copy values
	    $self->{$attrname} = $caller->{$attrname};
	} else { ## take default values
	    $self->{$attrname} = $self->_default_for($attrname) 
	}
    }
    return $self;
}


sub read_from_file {
    my ($self, $configFile) = @_;

    my %argsInFile;
    open(CH, "< $configFile") or die("Cant open $configFile: $!\n");
    while(my $line = <CH>) {
	next if ($line =~ /^\s*#/); ## skip comments
	if ($line =~ /$matchKeyValue/) {
	    my $key = $1;
	    my $value = $2;
	    $argsInFile{$key} = $value;
	}
    }
    close(CH);

    foreach my $attrname ( $self->_standard_keys() ) {
	my ($argname) = ($attrname =~ /^_(.*)/);

	if (exists $argsInFile{$argname}) {
	    $self->{$attrname} = $argsInFile{$argname};	
	}
    }
}


sub write_to_file {
    my ($self, $outFile) = @_;
    open(OH, "> $outFile") or die("Cant write to $outFile: $!\n");
    my $out = $self->to_string();
    print(OH $out);
    close(OH);
}


sub print {
    my $self = shift;
    my $out = $self->to_string();
    print $out;
}


sub to_string {
    my $self = shift;

    ## find longest key and value
    my $maxKeyLen = -99999;
    my $maxValLen = -99999;
    foreach my $attrname ( $self->_standard_keys() ) {
	$maxKeyLen = length($attrname) if (length($attrname) > $maxKeyLen);
	$maxValLen = length($self->{$attrname}) if (length($self->{$attrname}) > $maxValLen);
    }

    my $out = "";
    $out .= "-" x ($maxKeyLen + $maxValLen + 4) . "\n";
    $out .= "HHpred configuration parameters:\n";
    $out .= "-" x ($maxKeyLen + $maxValLen + 4) . "\n";

    foreach my $attrname ( sort $self->_standard_keys() ) {
	my ($argname) = ($attrname =~ /^_(.*)/);
	my $entry = sprintf("%-*s => %-s\n", $maxKeyLen, $argname, $self->{$attrname});
	$out .= $entry;
    }
    $out .= "-" x ($maxKeyLen + $maxValLen + 4) . "\n";

    return $out;
}


## automatically generated getters and setters:
## $AUTOLOAD contains full name of a routine and is checked for name/accessiblity
## then an anonymous routine (names e.g. get_name) is created and stored 
## in table for future use
sub AUTOLOAD {
    no strict "refs";
    my ($self, $newval) = @_;

    if ($AUTOLOAD =~ /.*::get(_\w+)/ && $self->_accessible($1, 'read')) {
	my $attr_name = $1;
	*{$AUTOLOAD} = sub { return $_[0]->{$attr_name} };
	return $self->{$attr_name}
    }
    if ($AUTOLOAD =~ /.*::set(_\w+)/ && $self->_accessible($1, 'write')) {
	my $attr_name = $1;
	*{$AUTOLOAD} = sub { $_[0]->{$attr_name} = $_[1]; return };
	$self->{$1} = $newval;
	return
    }
    ## mistaken?
    croak("No such method: $AUTOLOAD");
}


sub DESTROY {
}

1;
