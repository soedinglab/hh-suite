package QualityAssessModel;

## assesses a 3d model built from the templates in 'tList' based on hhsearch posterior probabilities
## for each query residue take the max. posterior probability over all templates as a quality measure

use strict;
use TemplateList;
use Template;
use utilities;

use vars qw($AUTOLOAD);
use Carp;

{
    my %_attr_data = #  DEFAULT  ACCESS
	( 
	  _tList    => ['###', 'read/write'],
	  _outbase  => ['###', 'read/write'],
	  _pdbFile  => ['###', 'read/write'],
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


## constructor
sub new {
    my ($caller, %arg) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = bless {}, $class;
    foreach my $attrname ($self->_standard_keys() ) {
	my ($argname) = ($attrname =~ /^_(.*)/);
	if (exists $arg{$argname}) {
	    $self->{$attrname} = $arg{$argname};
	} elsif ($caller_is_obj) {
	    $self->{$attrname} = $caller->{$attrname};
	} else {
	    $self->{$attrname} = $self->_default_for($attrname);
	}
    }
    return $self;
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


sub run {
    my $self = shift;
    my $pos = shift;

    my @templatePPs;
    ## generate all posteriors for all templates
    for (my $i=0; $i<$self->get_tList()->size(); $i++) {
	my $tabFile = $self->get_outbase() . "." . $self->get_tList()->get($i)->get_Filt() . ".tab";
	my %PP = &PosteriorsFromTabFile( $tabFile, $self->get_tList()->get($i)->get_No() );
	push (@templatePPs, \%PP);
    }

    ## find maximum pp for each query residue
    my @maxPP = (0) x $self->get_tList()->get_queryLength();
    for (my $i=0; $i<@maxPP; $i++) {
	for (my $j=0; $j<@templatePPs; $j++)  {
	    if (exists($templatePPs[$j]->{$i})) {
		$maxPP[$i] = &max($maxPP[$i], $templatePPs[$j]->{$i});
	    }
	}
    }

    ## write result in pdb file
    my $pdbFile = $self->get_pdbFile();
    open(PDB, "< $pdbFile");
    my @pdb = <PDB>;
    close(PDB);
    
    my $maxDisplacement = 100;
    my $minDisplacement = 5;

    # ATOM      1  N   SER A   2      31.507  77.672  34.913  1.00 39.58           N
    for (my $i=0; $i<@pdb; $i++) {
	if ($pdb[$i] =~ /^ATOM\s+\S+\s+\S+\s+\S+\s+(\d+)/) {
	    my $residue = $1;
	    my $temperatureFactor = $maxDisplacement - $maxPP[$residue-1]*$maxPP[$residue-1] * ($maxDisplacement - $minDisplacement);
	    $temperatureFactor = sprintf("%6.2f", $temperatureFactor);
	    $pdb[$i] =~ s/(.{60})(.{6})(.*)/$1$temperatureFactor$3/;
	}
    }
    
    open(PDB, "> $pdbFile");
    print (PDB @pdb);
    close(PDB);
}
