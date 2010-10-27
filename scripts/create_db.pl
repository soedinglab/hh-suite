#!/usr/bin/perl
# 
# Creates HHblits databases from HMMER-files, HMM-files or A3M-files

################################################################################################################################
# Update the following variables

$hh_dir = "/cluster/bioprogs/hhblits";         # path to needed tools (cstranslate, hhmake)

################################################################################################################################

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode
my $a3mext = "a3m";       # default A3M-file extension
my $hhmext = "hhm";       # default HHM-file extension
my $x = 0.3;
my $c = 4;
my $cs_lib = "$lib_dir/cs219.lib";
my $context_lib = "$lib_dir/context_data.lib";

my $append = 0;

my $help="
Creates HHblits databases from HMMER-files, HHM-files or A3M-files.

The recommended way to use this script is to start with a directory
of A3M-files (-a3mdir <DIR>) and let this script generates an A3M-
database (-oa3m <FILE>) and an HHM-database (-ohhm <FILE>).
If you already have HHM-models for your A3M-files, you can use them 
as additional input (-hhmdir <DIR>).
If you don't need the A3M-database, you can also start this script
with an directory of HHM-files (-hhmdir <DIR>) and as output only
the HHM-database (-ohhm <FILE>).

Usage: perl create_db.pl -i <dir> [options]

Options:
  -a3mdir <dir>  Input directory (directories) with A3M-files
  -hhmdir <dir>  Input directory (directories) with HHM- or HMMER-files 
                 (WARNING! Using HMMER databases could result in a decreased sensitivity!)

  -oa3m  <FILE>  Output filename for the A3M database 
                 (if not given, no A3M database will be build)
  -ohhm  <FILE>  Output filename for the HHM database 
                 (if not given, no HHM database will be build)

  -a3mext        Extension of A3M-files (default: $a3mext)
  -hhmext        Extension of HHM- or HMMER-files (default: $hhmext)

  -append        If the output file exists, append new files (default: overwrite)

  -v [0-5]       verbose mode (default: $v)

Examples:

   perl create_db.pl -a3mdir '/databases/scop_a3ms/*' -oa3m /databases/scop_a3m.db -ohhm /databases/scop_hhm.db

   perl create_db.pl -a3mdir /databases/scop_a3ms -hhmdir /databases/scop_hhms -oa3m /databases/scop_a3m.db -ohhm /databases/scop_hhm.db

   perl create_db.pl -hhmdir /databases/scop_hhms -ohhm /databases/scop_hhm.db
\n";

# Variable declarations
my $line;
my $command;
my $a3mdir  = "";
my $hhmdir  = "";
my $a3mfile = "";
my $hhmfile = "";
my @a3mfiles;
my @hhmfiles;

my $file;
my @files;
my $dir;
my @dirs;

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

if ($options=~s/ -a3mdir\s+(\S+) //) {$a3mdir=$1;}
if ($options=~s/ -hhmdir\s+(\S+) //) {$hhmdir=$1;}

if ($options=~s/ -oa3m\s+(\S+) //) {$a3mfile=$1;}
if ($options=~s/ -ohhm\s+(\S+) //) {$hhmfile=$1;}

if ($options=~s/ -a3mext\s+(\S+) //) {$a3mext=$1;}
if ($options=~s/ -hhmext\s+(\S+) //) {$hhmext=$1;}

if ($options=~s/ -append //) {$append=1;}

if ($options=~s/ -v\s+(\S+) //) {$v=$1;}


# Check input
if ($a3mdir eq "" && $hhmdir eq "") {
    print($help); print "ERROR! At least one input directory must be given!\n"; exit(1);
}

if ($a3mfile eq "" && $hhmfile eq "") {
    print($help); print "ERROR! At least one database (-ao3m or -ohhm) must be created!\n"; exit(1);
}

if ($a3mfile ne "") {
    if ($a3mdir eq "") {
	print($help); print "ERROR! Input directory with A3M-files needed for A3M database!\n"; exit(1);
    }
}

# Create tmp directory (plus path, if necessary)
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 

##############################################################################################
# Main part
##############################################################################################

if ($a3mfile ne "") {
    print "Creating A3M database $a3mfile ...\n";

    # check, if directory contains only a3m-files
    my $check = 1;
    @dirs = glob($a3mdir);
    foreach $dir (@dirs) {
	my @a3mfiles = glob("$dir/*.$a3mext");
	@files = glob("$dir/*");
	if (scalar(@files) != scalar(@a3mfiles)) {
	    $check = 0;
	    last;
	}
    }

    if ($check == 0) {   # write only a3mfiles in tmp-dir
	&System("mkdir $tmpdir/a3ms");
	@dirs = glob($a3mdir);
	foreach $dir (@dirs) {
	    &System("cp $dir/*.$a3mext $tmpdir/a3ms/");
	}
	$a3mdir = "$tmpdir/a3ms";
    }

    if ($append) {
	$command = "$hh_dir/ffindex_build -as $a3mfile $a3mfile.index $a3mdir";
    } else {
	$command = "$hh_dir/ffindex_build -s $a3mfile $a3mfile.index $a3mdir";
    }
    if (&System($command) != 0) {
	print "WARNING! Error with command $command!\n";
    }
} 
if ($hhmfile ne "") {

    if ($hhmdir eq "") {
	# Build HHMs from A3Ms
	print "Generate HHMs from A3Ms ...\n";
	$hhmdir = "$tmpdir/hhms";
	$hhmext = "hhm";
	&System("mkdir $hhmdir");
	
	@dirs = glob($a3mdir);
	foreach $dir (@dirs) {
	    @files = glob("$dir/*.$a3mext");
	    foreach $file (@files) {
		$file =~ /^\S+\/(\S+?)\.$a3mext$/;
		$command = "$hh_dir/hhmake -i $file -o $hhmdir/$1.hhm";
		if (&System($command) != 0) {
		    print "WARNING! Error with command $command!\n";
		}
	    }
	}
    } else {
	# check, if directory contains only hhm-files
	my $check = 1;
	@dirs = glob($hhmdir);
	foreach $dir (@dirs) {
	    my @hhmfiles = glob("$dir/*.$hhmext");
	    @files = glob("$dir/*");
	    if (scalar(@files) != scalar(@hhmfiles)) {
		$check = 0;
		last;
	    }
	}
	
	if ($check == 0) {   # write only a3mfiles in tmp-dir
	    &System("mkdir $tmpdir/hhms");
	    @dirs = glob($hhmdir);
	    foreach $dir (@dirs) {
		&System("cp $dir/*.$hhmext $tmpdir/hhms/");
	    }
	    $hhmdir = "$tmpdir/hhms";
	}
    }

    print "Creating HHM database $hhmfile ...\n";

    if ($append) {
	$command = "$hh_dir/ffindex_build -as $hhmfile $hhmfile.index $hhmdir";
    } else {
	$command = "$hh_dir/ffindex_build -s $hhmfile $hhmfile.index $hhmdir";
    }
    if (&System($command) != 0) {
	print "WARNING! Error with command $command!\n";
    }

}

if ($v < 4) {
    $command = "rm -rf $tmpdir";
    &System($command);
}

exit;


################################################################################################
### System command with return value parsed from output
################################################################################################
sub System()
{
    if ($v>2) {printf("\$ %s\n",$_[0]);} 
    return system($_[0])/256;
}
