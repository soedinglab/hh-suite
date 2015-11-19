#============================================================================
#
# Class::Singleton.pm
#
# Implementation of a "singleton" module which ensures that a class has
# only one instance and provides global access to it.  For a description 
# of the Singleton class, see "Design Patterns", Gamma et al, Addison-
# Wesley, 1995, ISBN 0-201-63361-2
#
# Written by Andy Wardley <abw@wardley.org>
#
# Copyright (C) 1998-2008 Andy Wardley.  All Rights Reserved.
# Copyright (C) 1998 Canon Research Centre Europe Ltd.
#
#============================================================================

package Class::Singleton;
require 5.004;
use strict;
use warnings;

our $VERSION = 1.4;


#========================================================================
#
# instance()
#
# Module constructor.  Creates an Class::Singleton (or derived) instance 
# if one doesn't already exist.  The instance reference is stored in the
# _instance variable of the $class package.  This means that classes 
# derived from Class::Singleton will have the variables defined in *THEIR*
# package, rather than the Class::Singleton package.  The impact of this is
# that you can create any number of classes derived from Class::Singleton
# and create a single instance of each one.  If the _instance variable
# was stored in the Class::Singleton package, you could only instantiate 
# *ONE* object of *ANY* class derived from Class::Singleton.  The first
# time the instance is created, the _new_instance() constructor is called 
# which simply returns a reference to a blessed hash.  This can be 
# overloaded for custom constructors.  Any addtional parameters passed to 
# instance() are forwarded to _new_instance().
#
# Returns a reference to the existing, or a newly created Class::Singleton
# object.  If the _new_instance() method returns an undefined value
# then the constructer is deemed to have failed.
#
#========================================================================

sub instance {
    my $class = shift;
    
    # already got an object
    return $class if ref $class;

    # we store the instance in the _instance variable in the $class package.
    no strict 'refs';
    my $instance = \${ "$class\::_instance" };
    defined $$instance
        ? $$instance
        : ($$instance = $class->_new_instance(@_));
}


#=======================================================================
# has_instance()
#
# Public method to return the current instance if it exists.
#=======================================================================

sub has_instance {
    my $class = shift;
    $class = ref $class || $class;
    no strict 'refs';
    return ${"$class\::_instance"};
}


#========================================================================
# _new_instance(...)
#
# Simple constructor which returns a hash reference blessed into the 
# current class.  May be overloaded to create non-hash objects or 
# handle any specific initialisation required.
#========================================================================

sub _new_instance {
    my $class = shift;
    my %args  = @_ && ref $_[0] eq 'HASH' ? %{ $_[0] } : @_;
    bless { %args }, $class;
}



1;

__END__

=head1 NAME

Class::Singleton - Implementation of a "Singleton" class 

=head1 SYNOPSIS

    use Class::Singleton;
    
    my $one = Class::Singleton->instance();   # returns a new instance
    my $two = Class::Singleton->instance();   # returns same instance

=head1 DESCRIPTION

This is the C<Class::Singleton> module.  A Singleton describes an object class
that can have only one instance in any system.  An example of a Singleton
might be a print spooler or system registry.  This module implements a
Singleton class from which other classes can be derived.  By itself, the
C<Class::Singleton> module does very little other than manage the instantiation
of a single object.  In deriving a class from C<Class::Singleton>, your module 
will inherit the Singleton instantiation method and can implement whatever
specific functionality is required.

For a description and discussion of the Singleton class, see 
"Design Patterns", Gamma et al, Addison-Wesley, 1995, ISBN 0-201-63361-2.

=head1 PREREQUISITES

C<Class::Singleton> requires Perl version 5.004 or later. If you have an older
version of Perl, please upgrade to latest version, available from your nearest
CPAN site (see L<INSTALLATION> below).

=head1 INSTALLATION

The C<Class::Singleton> module is available from CPAN. As the 'perlmod' man
page explains:

    CPAN stands for the Comprehensive Perl Archive Network.
    This is a globally replicated collection of all known Perl
    materials, including hundreds of unbunded modules.
    
    [...]
    
    For an up-to-date listing of CPAN sites, see
    http://www.perl.com/perl/ or ftp://ftp.perl.com/perl/ .

The module is available in the following directories:

    /modules/by-module/Class/Class-Singleton-<version>.tar.gz
    /authors/id/ABW/Class-Singleton-<version>.tar.gz

C<Class::Singleton> is distributed as a single gzipped tar archive file:

    Class-Singleton-<version>.tar.gz

Note that "<version>" represents the current version number, of the 
form "C<1.23>".  See L<VERSION> below to determine the current version 
number for C<Class::Singleton>.

Unpack the archive to create an installation directory:

    gunzip Class-Singleton-<version>.tar.gz
    tar xvf Class-Singleton-<version>.tar

'cd' into that directory, make, test and install the module:

    cd Class-Singleton-<version>
    perl Makefile.PL
    make
    make test
    make install

The 'C<make install>' will install the module on your system.  You may need 
root access to perform this task.  If you install the module in a local 
directory (for example, by executing "C<perl Makefile.PL LIB=~/lib>" in the 
above - see C<perldoc MakeMaker> for full details), you will need to ensure 
that the C<PERL5LIB> environment variable is set to include the location, or 
add a line to your scripts explicitly naming the library location:

    use lib '/local/path/to/lib';

=head1 USING THE CLASS::SINGLETON MODULE

To import and use the C<Class::Singleton> module the following line should 
appear in your Perl program:

    use Class::Singleton;

The L<instance()> method is used to create a new C<Class::Singleton> instance,
or return a reference to an existing instance. Using this method, it is only
possible to have a single instance of the class in any system.

    my $highlander = Class::Singleton->instance();

Assuming that no C<Class::Singleton> object currently exists, this first call
to L<instance()> will create a new C<Class::Singleton> and return a reference
to it. Future invocations of L<instance()> will return the same reference.

    my $macleod    = Class::Singleton->instance();

In the above example, both C<$highlander> and C<$macleod> contain the same
reference to a C<Class::Singleton> instance.  There can be only one.

=head1 DERIVING SINGLETON CLASSES

A module class may be derived from C<Class::Singleton> and will inherit the 
L<instance()> method that correctly instantiates only one object.

    package PrintSpooler;
    use base 'Class::Singleton';
    
    # derived class specific code
    sub submit_job {
        ...
    }
    
    sub cancel_job {
        ...
    }

The C<PrintSpooler> class defined above could be used as follows:

    use PrintSpooler;
    
    my $spooler = PrintSpooler->instance();
    
    $spooler->submit_job(...);

The L<instance()> method calls the L<_new_instance()> constructor method the
first and only time a new instance is created. All parameters passed to the
L<instance()> method are forwarded to L<_new_instance()>. In the base class
the L<_new_instance()> method returns a blessed reference to a hash array
containing any arguments passed as either a hash reference or list of named 
parameters. 

    package MyConfig;
    use base 'Class::Singleton';
    
    sub foo {
        shift->{ foo };
    }
    
    sub bar {
        shift->{ bar };
    }
    
    package main;
    
    # either: hash reference of named parameters
    my $config = MyConfig->instance({ foo => 10, bar => 20 });
    
    # or: list of named parameters
    my $config = MyConfig->instance( foo => 10, bar => 20 );
    
    print $config->foo();   # 10
    print $config->bar();   # 20

Derived classes may redefine the L<_new_instance()> method to provide more
specific object initialisation or change the underlying object type (to a list
reference, for example).

    package MyApp::Database;
    use base 'Class::Singleton';
    use DBI;
    
    # this only gets called the first time instance() is called
    sub _new_instance {
        my $class = shift;
        my $self  = bless { }, $class;
        my $db    = shift || "myappdb";    
        my $host  = shift || "localhost";
        
        $self->{ DB } = DBI->connect("DBI:mSQL:$db:$host")
            || die "Cannot connect to database: $DBI::errstr";
        
        # any other initialisation...
        
        return $self;
    }

The above example might be used as follows:

    use MyApp::Database;
    
    # first use - database gets initialised
    my $database = MyApp::Database->instance();

Some time later on in a module far, far away...

    package MyApp::FooBar
    use MyApp::Database;
    
    # this FooBar object needs access to the database; the Singleton
    # approach gives a nice wrapper around global variables.
    
    sub new {
        my $class = shift;
        bless {
            database => MyApp::Database->instance(),
        }, $class;
    }

The C<Class::Singleton> L<instance()> method uses a package variable to store
a reference to any existing instance of the object. This variable,
"C<_instance>", is coerced into the derived class package rather than the base
class package.

Thus, in the C<MyApp::Database> example above, the instance variable would
be:

    $MyApp::Database::_instance;

This allows different classes to be derived from C<Class::Singleton> that can
co-exist in the same system, while still allowing only one instance of any one
class to exists. For example, it would be possible to derive both
'C<PrintSpooler>' and 'C<MyApp::Database>' from C<Class::Singleton> and have a
single instance of I<each> in a system, rather than a single instance of
I<either>.

You can use the L<has_instance()> method to find out if a particular class 
already has an instance defined.  A reference to the instance is returned or
C<undef> if none is currently defined.

    my $instance = MyApp::Database->has_instance()
        || warn "No instance is defined yet";

=head1 METHODS

=head2 instance()

This method is called to return a current object instance or create a new
one by calling L<_new_instance()>.

=head2 has_instance()

This method returns a reference to any existing instance or C<undef> if none
is defined.

    my $testing = MySingleton1->has_instance()
        || warn "No instance defined for MySingleton1";

=head2 _new_instance()

This "private" method is called by L<instance()> to create a new object
instance if one doesn't already exist. It is not intended to be called
directly (although there's nothing to stop you from calling it if you're
really determined to do so).

It creates a blessed hash reference containing any arguments passed to the
method as either a hash reference or list of named parameters.

    # either: hash reference of named parameters
    my $example1 = MySingleton1->new({ pi => 3.14, e => 2.718 });

    # or: list of named parameters
    my $example2 = MySingleton2->new( pi => 3.14, e => 2.718 );

It is important to remember that the L<instance()> method will I<only> call
the I<_new_instance()> method once, so any arguments you pass may be silently
ignored if an instance already exists. You can use the L<has_instance()>
method to determine if an instance is already defined.

=head1 AUTHOR

Andy Wardley E<lt>abw@wardley.orgE<gt> L<http://wardley.org/>

Thanks to Andreas Koenig for providing some significant speedup patches and
other ideas.

=head1 VERSION

This is version 1.4, released September 2007

=head1 COPYRIGHT

Copyright Andy Wardley 1998-2007.  All Rights Reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut
