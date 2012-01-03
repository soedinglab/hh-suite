#!/usr/bin/env perl
#
# create_hhblits.pl
# Creates HHblits database files from A3M and HHM/HMMER-formatted files 
# Usage: Usage: perl create_hhblitsdb.pl -o <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>] [more_options]
#
# (C) Johannes Soeding & Michael Remmert
# Available under GNU Public License Version 3
#
# Reference: 
# Remmert M., Biegert A., Hauser A., and Soding J.
# HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
# Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

use strict;

################################################################################################################################
# Check the following paths

my $script_dir=$ENV{"HHLIB"}."/scripts";# path to directory with scripts (create_profile_from_hhm.pl, create_profile_from_hmmer.pl)
my $hh_dir   = $ENV{"HHLIB"}."/bin";    # path to cstranslate, hhmake
my $data_dir = $ENV{"HHLIB"}."/data";   # path to context and cs files  (context_data.lib, cs219.lib)


################################################################################################################################

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=3;             # verbose mode
my $a_if_append = ""; # do not append by default
my $a3mext = "a3m";   # default A3M-file extension
my $hhmext = "hhm";   # default HHM-file extension
my $csext = "seq219";   # default HHM-file extension
my $cs_lib = "$data_dir/cs219.lib";
my $context_lib = "$data_dir/context_data.lib";
my $cpu = 4;

# Variable declarations
my $line;
my $command;
my $a3mdir  = "";     # name of input A3M directory
my $hhmdir  = "";     # name of input HHM/HMM directory
my $csdir  = "";      # name of input cs directory
my $a3mfile = "";     # name of packed ouput A3M file 
my $hhmfile = "";     # name of packed ouput HHM file
my $csfile = "";      # name of cs sequence db file
my $dbname = "";      # output db name
my $logfile = "/null/dev"; # log file 
my $file;
my $dir;
my $numcsfiles=0;
my $numa3mfiles=0;
my $numhhmfiles=0;
my $help="
Creates HHblits database files from MSA and HMM files 

Usage: perl create_hhblitsdb.pl -o <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>] [more_options]

Depending on the input directories, the script generates the following HHblits database files:
  <db_name>.cs219               column-state sequences, one for each MSA/HMM (used by HHblits for prefiltering)
  <db_name>.cs219.sizes         number of sequences and characters in <db_name>.cs219
  <db_name>_a3m_db              packed file containing A3M alignments read from <a3m_dir>
  <db_name>_a3m_db.index        index file for packed A3M file
  <db_name>_a3m.db.index.sizes  number of lines in <db_name>_a3m_db.index
  <db_name>_hhm_db              packed file containing HHM-formatted HMMs read from <hhm_dir>
  <db_name>_hhm_db.index        index file for packed HHM file
  <db_name>_hhm_db.index.sizes  number of lines in <db_name>_hhm_db.index

Options:
  -o <db_name>     name of database
  -ia3m <a3m_dir>  input directory (or glob of directories) with A3M-formatted files
  -ihhm <hhm_dir>  input directory (or glob of directories) with HHM- (or HMMER-) formatted files 
                   (WARNING! Using HMMER instead of HHM format will result in a decreased performance)
  -ics  <cs_dir>   input directory (or glob of directories) with column state sequences

  -log <logfile>   log file with output of cstranslate and hhmake commands
  -csext           extension of column state sequences (default: $csext)
  -a3mext          extension of A3M-formatted files (default: $a3mext)
  -hhmext          extension of HHM- or HMMER-formatted files (default: $hhmext)
  -append          if the packed db files exists, append input A3M/HHM files (default: overwrite)
  -v [1-3]         verbose mode (default: $v)
  -cpu <int>       numbers of threads to launch for generating cs219 and hhm files (default = $cpu)

 
Example 1: only -ia3m is given; cs sequences and hhm files are generated from a3m files
   perl create_hhblitsdb.pl -o hhblits_dbs/mydb -ia3m mydb/a3ms 

Example 2: only -ihhm is given; cs sequences are generated from hhm files, no a3m db file is generated
   perl create_hhblitsdb.pl -o hhblits_dbs/mydb -ihhm mydb/hhms 

Example 3: -ia3m and -ihhm are given; cs sequences are generated from a3m files
   perl create_hhblitsdb.pl -o hhblits_dbs/mydb -ia3m mydb/a3ms -ihhm mydb/hhms   

Example 4: -ics, -ia3m, and -ihhm are given; all db files are created 
   perl create_hhblitsdb.pl -o hhblits_dbs/mydb -ia3m mydb/a3ms -ihhm mydb/hhms -ics mydb/cs  

Example 5: using glob expression to specify several input databases
   perl create_hhblitsdb.pl -o hhblits_dbs/mydb -ihhm 'mydbs*/hhms'  
\n";




###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

for (my $i=0; $i<@ARGV; $i++) {

    if ($ARGV[$i] eq "-ics") {
	if (++$i<@ARGV) {
	    $csdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ics option!\n");
	}
    } elsif ($ARGV[$i] eq "-ia3m") {
	if (++$i<@ARGV) {
	    $a3mdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ia3m option!\n");
	}
    } elsif ($ARGV[$i] eq "-ihhm") {
	if (++$i<@ARGV) {
	    $hhmdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ihhm option!\n");
	}
    } elsif ($ARGV[$i] eq "-o") {
	if (++$i<@ARGV) {
	    $dbname=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } elsif ($ARGV[$i] eq "-log") {
	if (++$i<@ARGV) {
	    $logfile=$ARGV[$i];
	    unlink $logfile;
	} else {
	    die ("$help\n\nERROR! Missing filename after -log option!\n");
	}
    } elsif ($ARGV[$i] eq "-csext") {
	if (++$i<@ARGV) {
	    $csext=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing extension after -csext option!\n");
	}
    } elsif ($ARGV[$i] eq "-a3mext") {
	if (++$i<@ARGV) {
	    $a3mext=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing extension after -a3mext option!\n");
	}
    } elsif ($ARGV[$i] eq "-hhmext") {
	if (++$i<@ARGV) {
	    $hhmext=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing extension after -hhmext option!\n");
	}
    } elsif ($ARGV[$i] eq "-v") {
	if (++$i<@ARGV) {
	    $v=$ARGV[$i];
	} else {
	    $v = 2;
	}
    } elsif ($ARGV[$i] eq "-append") {
	$a_if_append="a";
    } else {
	print "WARNING! Unknown option $ARGV[$i]!\n";
    }
}

# Check input
if (!$dbname) {print($help); die("ERROR! Name of database is missing! Use -o <db_name>\n");}
if ($a3mdir) {$a3mdir=~s/\/$//; $a3mfile = $dbname."_a3m_db"; $hhmfile = $dbname."_hhm_db";}
if ($hhmdir) {$hhmdir=~s/\/$//;}
$csfile = $dbname.".cs219";

if ($a3mdir eq "" && $hhmdir eq "" && $csdir eq "") {
    print($help); print "ERROR! At least one input directory must be given!\n"; exit(1);
}

# Create tmp directory (plus path, if necessary)
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 

##############################################################################################
# Generate column-state sequences in $tmpdir if no -ics directory given
##############################################################################################

if (!$csdir) 
{
    my $x = 0.3;    # parameters for cstranslate
    my $c = 4;      # parameters for cstranslate
    
    if ($a3mdir) {
	my @dirs = glob($a3mdir);
	foreach $dir (@dirs) {
	    print("\nGenerating seq219 files in $tmpdir/ from a3m files in $dir/\n\n");
	    $command = "$hh_dir/cstranslate -i \$file -o $tmpdir/\$base.seq219 -D $context_lib -A $cs_lib -x $x -c $c 1>/dev/null 2>>$logfile";
	    &System("$script_dir/multithread.pl '".$dir."/*.".$a3mext."' '$command' -cpu $cpu");
	    $numa3mfiles += scalar(glob("$dir/*.a3m"));
	}
	
    } elsif ($hhmdir) {
	
	my @dirs = glob($hhmdir);
	foreach $dir (@dirs) {

	    if ($hhmext == "hmm") {
		print("\nGenerating prf profile files in $tmpdir/ from hmm files in $dir/\n\n");
		$command = "$script_dir/create_profile_from_hmmer.pl -i \$file -o $tmpdir/\$base.prf 1>>$logfile 2>&1";
		&System("$script_dir/multithread.pl '".$dir."/*.".$hhmext."' '$command' -cpu $cpu");
		
		print("\nGenerating seq219 files in $tmpdir/ from prf files in $tmpdir/\n\n");
		$command = "$hh_dir/cstranslate -i \$file -o \$name.seq219 -A $cs_lib 1>>$logfile 2>&1";
		&System("$script_dir/multithread.pl '".$tmpdir."/*.prf' '$command' -cpu $cpu");
		
	    } else { # $hhmext == "hhm"
		print("\nGenerating prf profile files in $tmpdir/ from hhm files in $dir/\n\n");
		$command = "$script_dir/create_profile_from_hmmer.pl -i \$file -o $tmpdir/\$base.prf 1>>$logfile 2>&1";
		&System("$script_dir/multithread.pl '".$dir."/*.".$hhmext."' '$command' -cpu $cpu");
		
		print("\nGenerating seq219 files in $tmpdir/ from prf files in $tmpdir/\n\n");
		$command = "$hh_dir/cstranslate -i \$file -o \$name.seq219 -A $cs_lib -D $context_lib -x $x -c $c 1>>$logfile 2>&1";
		&System("$script_dir/multithread.pl '".$tmpdir."/*.prf' '$command' -cpu $cpu");
	    }
	}
    }
}

##############################################################################################
# Generate column-state database file
##############################################################################################

# Concatenate columns state sequences into cs database file
if ($a_if_append == "" && -e $csfile) {unlink $csfile;}
&System("find $tmpdir -name '*.seq219' -exec cat '{}' + >> $csfile"); # 
#&System("cd $tmpdir; cat *.seq219 >> ".cwd()."/$csfile"); # might lead to command line buffer overflow
$numcsfiles = scalar(glob("$tmpdir/*.seq219"));

# Count number of sequences and characters in cs database file
$numcsfiles = 0;
my $num_chars = 0;
open (IN, "<$csfile");
while ( <IN>) {
    if (/^>/) {$numcsfiles++;}
    else {
	$num_chars += length;
    }
}
close(IN);

open (OUT, ">$csfile.sizes");
print OUT "$numcsfiles $num_chars\n";
close OUT;


##############################################################################################
# Generate hhm files with hhmake from a3m files if no -ihhm directory given
##############################################################################################

if (!$hhmdir) 
{
    if ($a3mdir) {
	my @dirs = glob($a3mdir);
	foreach $dir (@dirs) {
	    print("\nGenerating hhm files in $tmpdir/ from a3m files in $dir/\n\n");
	    $command = "$hh_dir/hhmake -i \$file -o $tmpdir/\$base.hhm  2>>$logfile 1>/dev/null";
	    &System("$script_dir/multithread.pl '".$dir."/*.".$a3mext."' '$command' -cpu $cpu");	
	}
	$hhmdir = $tmpdir;
	$numhhmfiles = scalar(glob("$tmpdir/*.hhm"));
    }
}


##############################################################################################
# Generate packed A3M and HMM files and index files
##############################################################################################

# Generate packed A3M file and index file?
if ($a3mfile ne "") {
    print "Creating packed A3M database file $a3mfile ...\n";

    open (OUT, ">$tmpdir/a3m.filelist");
    $numa3mfiles = 0;
    my @dirs = glob($a3mdir);
    foreach $dir (@dirs) {
	my @files = glob("$dir/*.$a3mext");
	$numa3mfiles += scalar(@files);
	foreach $file (@files) {
	    print OUT "$file\n";
	}
    }
    close OUT;
    
    $command = "$hh_dir/ffindex_build -".$a_if_append."s -f $tmpdir/a3m.filelist $a3mfile $a3mfile.index";
    &System($command);
 
    open (OUT, ">$a3mfile.index.sizes");
    print OUT "$numa3mfiles\n";
    close OUT;
} 

# Generate packed HHMM file and index file?
if ($hhmfile ne "") {
    print "Creating packed HHM database file $hhmfile ...\n";

    open (OUT, ">$tmpdir/hhm.filelist");
    $numhhmfiles = 0;
    my @dirs = glob($hhmdir);
    foreach $dir (@dirs) {
	my @files = glob("$dir/*.$hhmext");
	$numhhmfiles += scalar(@files);
	foreach $file (@files) {
	    print OUT "$file\n";
	}
    }
    close OUT;

    $command = "$hh_dir/ffindex_build -".$a_if_append."s -f $tmpdir/hhm.filelist $hhmfile $hhmfile.index";
    &System($command);
 
    open (OUT, ">$hhmfile.index.sizes");
    print OUT "$numhhmfiles\n";
    close OUT;
}

printf("Number of $a3mext files: %i\n",$numa3mfiles);
printf("Number of $hhmext files: %i\n",$numhhmfiles);
printf("Number of $csext  files: %i\n",$numcsfiles);

my $err=0;
if ($numa3mfiles && $numhhmfiles && $numa3mfiles != $numhhmfiles) {
    print("WARNING: Number of $a3mext files not equal to number of $hhmext files\n"); $err=1;
}
if ($numcsfiles && $numhhmfiles && $numcsfiles != $numhhmfiles) {
    print("WARNING: Number of $csext files not equal to number of $hhmext files\n"); $err=1;
}
if ($numcsfiles && $numa3mfiles && $numcsfiles != $numa3mfiles) {
    print("WARNING: Number of $csext files not equal to number of $a3mext files\n"); $err=1;
}

if ($err==1) {print("$tmpdir will not be removed to check for missing files\n");}
elsif ($v<3) {
    $command = "rm -rf $tmpdir";
    &System($command);
}
wait;
exit;


################################################################################################
### System command with return value parsed from output
################################################################################################
sub System()
{
    if ($v>=3) {printf("\$ %s\n",$_[0]);} 
    system($_[0]); # ==0 or print(STDERR "ERROR: command $_[0] failed with error code $?\n");
}
