package Template;

use strict;
use utilities;

use vars qw($AUTOLOAD);
use Carp;

{
    my %_attr_data = #  DEFAULT  ACCESS
	( 
	  _Filt    => ['###',   'read/write'],
	  _No      => ['###',   'read'],
	  _Hit     => ['###',   'read'],
	  _Prob    => ['###',   'read'],
	  _Eval    => ['###',   'read'],
	  _Pval    => ['###',   'read'],
	  _Score   => ['###',   'read/write'],
	  _SS      => ['###',   'read'],
	  _Cols    => ['###',   'read'],
	  _Qstart  => ['###',   'read'],
	  _Qend    => ['###',   'read'],
	  _Tstart  => ['###',   'read'],
	  _Tend    => ['###',   'read'],
	  _HMM     => ['###',   'read'],
	  _Sim     => ['###',   'read/write'],
	  _Ident   => ['###',   'read/write'],
	  _SumProbL=> ['###',   'read/write'],
	  _ss_dssp => ['',      'read/write'],
	  _conf    => ['',      'read/write'],
	  _predTM  => ['###',   'read/write']
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


sub to_string {
    my $self = shift;
    my $spacer = shift;
    my $out = "";
    $out .= defined($spacer) ? "$self->{_Filt}$spacer" : "$self->{_Filt}\t";
    $out .= defined($spacer) ? "$self->{_No}$spacer" : "$self->{_No}\t";
    $out .= defined($spacer) ? "$self->{_Hit}$spacer" : "$self->{_Hit}\t";
    $out .= defined($spacer) ? "$self->{_Prob}$spacer" : "$self->{_Prob}\t";
    $out .= defined($spacer) ? "$self->{_Eval}$spacer" : "$self->{_Eval}\t";
    $out .= defined($spacer) ? "$self->{_Pval}$spacer" : "$self->{_Pval}\t";
    $out .= defined($spacer) ? "$self->{_Score}$spacer" : "$self->{_Score}\t";
    $out .= defined($spacer) ? "$self->{_SS}$spacer" : "$self->{_SS}\t";
    $out .= defined($spacer) ? "$self->{_Cols}$spacer" : "$self->{_Cols}\t";
    $out .= defined($spacer) ? "$self->{_Qstart}$spacer" : "$self->{_Qstart}\t";
    $out .= defined($spacer) ? "$self->{_Qend}$spacer" : "$self->{_Qend}\t";
    $out .= defined($spacer) ? "$self->{_Tstart}$spacer" : "$self->{_Tstart}\t";
    $out .= defined($spacer) ? "$self->{_Tend}$spacer" : "$self->{_Tend}\t";
    $out .= defined($spacer) ? "$self->{_HMM}$spacer" : "$self->{_HMM}\t";
    $out .= defined($spacer) ? "$self->{_Ident}$spacer" : "$self->{_Ident}\t";
    $out .= defined($spacer) ? "$self->{_Sim}$spacer" : "$self->{_Sim}\t";
    $out .= defined($spacer) ? "$self->{_SumProbL}$spacer" : "$self->{_SumProbL}\t";
    $out .= defined($spacer) ? "$self->{_predTM}$spacer" : "$self->{_predTM}\t";

    return $out;
}


## check whether two templates have same keys and values
sub equals {
    my ($self, $template) = @_;

    my %cmp = map { $_ => 1 } keys %{$self};
    for my $key (keys %{$template}) {
        last unless exists $cmp{$key};
        last unless $self->{$key} eq $template->{$key};
        delete $cmp{$key};
    }
    if (%cmp) {
	return 0;
    } else {
	return 1;
    }
}

1;
