#!/usr/bin/env perl
#
# hhbsuite.pl
# Creates HH-suite database files from A3M and HHM/HMMER-formatted files 
# Usage: Usage: perl hhbsuitedb.pl -o <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>] [more_options]
#
#     HHsuite version 3.0.0 (15-03-2015)
#
#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#     (C) Johannes Soeding, 2012

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use strict;
#use File::Glob 'bsd_glob'; # splits patterns delimited by spaces into multiple patterns and applies them using OR

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;             # verbose mode
my $a_if_append = ""; # do not append by default (default: create new db)
my $remove = 0;       # do not remove by default (default: create new db)
my $hhmext = "hhm";   # default HHM-file extension
my $csext = "cs219"; # default cs sequence file extension
my $cpu = 8;

# Variable declarations
my $line;
my $command;
my $a3mglob = "";      # name of input A3M directory or glob expression
my $hhmglob = "";      # name of input HHM/HMM directory or glob expression
my $csglob  = "";      # name of input cs directory or glob expression
my $dbname = "";      # output db name
my $logfile = "/dev/null"; # log file 
my $file;
my $numcsfiles= 0;
my $num_chars = 0;
my $numa3mfiles=0;
my $numhhmfiles=0;
my $remove_expr="this should never match";
my $help="
hhsuitedb.pl from HHsuite $VERSION  
Builds HH-suite database from a3m formatted MSAs and/or from HMMs (-o).
MSAs and HMMs can also be added (-a) to or removed (-r) from an existing database. 

Usage: hhblitsdb.pl -o|-a|-r <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>]... 

With option -o, the following HH-suite database files can be generated:
 <db_name>_cs219.ffdata    packed file containing column-state sequences read from <cs_dir>
 <db_name>_cs219.ffindex   index file for packed column-state sequence file
 <db_name>_cs219.ffsize    number of lines in <db_name>_cs219.ffindex
 <db_name>_a3m.ffdata      packed file containing A3M alignments read from <a3m_dir>
 <db_name>_a3m.ffindex     index file for packed A3M file
 <db_name>_a3m.ffsize      number of lines in <db_name>_a3m.ffindex
 <db_name>_hhm.ffdata      packed file containing HHM-formatted HMMs read from <hhm_dir>
 <db_name>_hhm.ffindex     index file for packed HHM file
 <db_name>_hhm.ffsize      number of lines in <db_name>_hhm.ffindex

Options:
 -o <db_name>    create database with this name
 -a <db_name>    append files to database with this name
 -r <db_name>    remove files from database with this name
 -ia3m <a3m_dir> input directory (or glob of directories) with A3M-formatted files
                 These files MUST have extension 'a3m'.
 -ihhm <hhm_dir> input directory (or glob of directories) with HHM (or HMMER) files 
                 These files MUST have extension 'hhm' (HHsuite) or 'hmm' (HMMER3). 
 -ics  <cs_dir>  input directory (or glob of directories) with column state sequences
 -csext <ext>    extension of column state sequences (default: $csext)
 -log <logfile>  log file recording stderr stream of cstranslate and hhmake commands
 -hmm            use HMMER-formatted files. These MUST have extension hmm
                 (WARNING! HMMER format results in decreased performance over HHM format)
 -v [1-3]        verbose mode (default: $v)
 -cpu <int>      number of threads to generate cs219 and hhm files (default = $cpu)
 -f 'file_glob'  string with list of glob expressions of files to remove
 
Example 1: only -ia3m given; cs sequences and hhm files are generated from a3m files
   perl hhblitsdb.pl -o databases/mydb -ia3m mydb/a3ms/ 

Example 2: only -ihhm given; cs sequences are generated from hhm files, but no a3m db file 
   perl hhblitsdb.pl -o databases/mydb -ihhm mydb/hhms/ 

Example 3: -ia3m and -ihhm given; cs sequences are generated from a3m files
   perl hhblitsdb.pl -o databases/mydb -ia3m mydb/a3ms/ -ihhm mydb/hhms/   

Example 4: -ics, -ia3m, and -ihhm given; all db files are created 
   perl hhblitsdb.pl -o databases/mydb -ia3m mydb/a3ms/ -ihhm mydb/hhms/ -ics mydb/cs/  

Example 5: using glob expression to specify files (note the single quotes)
   perl hhblitsdb.pl -o databases/mydb -ihhm 'mydbs*/hhms/*.hhm'  

Example 6: add files to database; cs sequences and hhm files are generated from a3m files
   perl hhblitsdb.pl -a databases/mydb -ia3m 'mydbs/a3ms/g1a*.a3m'  

Example 7: remove files from database
   perl hhblitsdb.pl -r databases/mydb -f 'mydbs/a3ms/g1a* mydbs2/*'  
\n";


###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

for (my $i=0; $i<@ARGV; $i++) {
    if ($ARGV[$i] eq "-ics") {
	if (++$i<@ARGV) {
	    $csglob=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ics option!\n");
	}
    } elsif ($ARGV[$i] eq "-ia3m") {
	if (++$i<@ARGV) {
	    $a3mglob=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ia3m option!\n");
	}
    } elsif ($ARGV[$i] eq "-ihhm") {
	if (++$i<@ARGV) {
	    $hhmglob=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ihhm option!\n");
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
    } elsif ($ARGV[$i] eq "-hmm") {
	$hhmext="hmm";
	print("\nWARNING! HMMER format results in decreased performance over HHM format. We recommend to generate hhm files directly from multiple sequence alignments using hmake.\n");
    } elsif ($ARGV[$i] eq "-v") {
	if (++$i<@ARGV) {
	    $v=$ARGV[$i];
	} else {
	    $v = 2;
	}
    } elsif ($ARGV[$i] eq "-cpu") {
	if (++$i<@ARGV) {
	    $cpu=$ARGV[$i];
	}
    } elsif ($ARGV[$i] eq "-f") {
	if (++$i<@ARGV) {
	    $remove_expr=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing expression after -f option!\n");
	}
   } elsif ($ARGV[$i] eq "-r") {
	if (++$i<@ARGV) {
	    if ($dbname!="") {die("$help\n\nERROR! options -o and -r not compatible!\n");}
	    $dbname=$ARGV[$i];
	    $remove=1;
	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } elsif ($ARGV[$i] eq "-a") {
	if (++$i<@ARGV) {
	    if ($remove==1) {die("$help\n\nERROR! options -r and -a not compatible!\n");}
	    if ($dbname!="") {die("$help\n\nERROR! options -o and -a not compatible!\n");}
	    $dbname=$ARGV[$i];
	    $a_if_append="a";

	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } elsif ($ARGV[$i] eq "-o") {
	if (++$i<@ARGV) {
	    if ($remove==1) {die("$help\n\nERROR! options -r and -o not compatible!\n");}
	    if ($a_if_append) {die("$help\n\nERROR! options -a and -o not compatible!\n");}
	    $dbname=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } else {
	if ($dbname=="") {
	    $dbname=$ARGV[$i];
	} else {
	    print "WARNING! Unknown option $ARGV[$i]!\n";
	}
    }
}

# Check input
if (!$dbname) {print($help); die("ERROR! Name of database is missing! Use -o <db_name>\n");}
$a3mglob=~s/\/$//; # remove trailing '/'
if ($hhmglob) {$hhmglob=~s/\/$//;} # remove trailing '/'
if ($csglob) {$csglob=~s/\/$//;} # remove trailing '/'
my $a3mdata  = $dbname."_a3m.ffdata"; 
my $a3mindex = $dbname."_a3m.ffindex"; 
my $a3msize  = $dbname."_a3m.ffsize"; 
my $hhmdata  = $dbname."_hhm.ffdata"; 
my $hhmindex = $dbname."_hhm.index"; 
my $hhmsize  = $dbname."_hhm.ffsize"; 
my $csdata   = $dbname."_cs.ffdata"; 
my $csindex  = $dbname."_cs.index"; 
my $cssize   = $dbname."_cs.ffsize"; 

if ($a_if_append eq "" && $remove==0) {
    unlink $a3mdata, $a3mindex, $a3msize, $hhmdata, $hhmindex, $hhmsize, $csdata, $csindex, $cssize
}

if ($a3mglob eq "" && $hhmglob eq "" && $csglob eq "" && $remove==0) {
    print($help); print "ERROR! At least one input directory must be given!\n"; exit(1);
}

# If $csglob is simple directory instead of glob expression, turn it into glob expression
$a3mglob = &TurnDirIntoGlob($a3mglob,"a3m");
$hhmglob = &TurnDirIntoGlob($hhmglob,$hhmext);
$csglob  = &TurnDirIntoGlob($csglob,$csext);


# Create tmp directory (plus path, if necessary)
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 
unlink glob("$tmpdir/*"); # clean up directory if it already exists
unlink $logfile;
    


##############################################################################################
# Remove files?
##############################################################################################
if ($remove==1) {

    # Transform glob into regex
    if ($remove_expr =~ /\\/ || $remove_expr =~ /\//) {
	die("Error: / and \\ symbols not allowed in -f <glob>. Give glob for file names **without path**\n");
    }
    $remove_expr =~ s/\./\\\./g;   # . => \.
    $remove_expr =~ s/\*/\.\*/g;   # * => \w*
    $remove_expr =~ s/\?/\./g;    # ? => .

    # Unlink files
    &UnlinkFilesFromFfindex($a3mindex,$a3msize,$remove_expr);
    &UnlinkFilesFromFfindex($hhmindex,$hhmsize,$remove_expr);
    &UnlinkFilesFromFfindex($csindex, $cssize, $remove_expr);
    printf("\n");

} else {

##############################################################################################
# Generate new db or append to old
##############################################################################################


##############################################################################################
# Generate column-state database file
##############################################################################################



    # Generate column-state sequences in $tmpdir if no -ics directory given
    if (!$csglob) 
    {
	my $x = 0.3;    # parameters for cstranslate
	my $c = 4;      # parameters for cstranslate
	
	if ($a3mglob) {
	    print("Generating cs219 files in $tmpdir/ from a3m files $a3mglob\n\n");
	    $command = "$hhbin/cstranslate -i \$file -o $tmpdir/\$base.cs219 -D $context_lib -A $cs_lib -x $x -c $c 1>>$logfile 2>>$logfile";
	    &HHPaths::System("$hhscripts/multithread.pl '".$a3mglob."' '$command' -cpu $cpu");
	    
	} elsif ($hhmglob) {
	    
	    if ($hhmext eq "hmm") {
		print("\nGenerating prf profile files in $tmpdir/ from hmm files $hhmglob/\n\n");
		$command = "$hhscripts/create_profile_from_hmmer.pl -i \$file -o $tmpdir/\$base.prf 1>/dev/null 2>>$logfile";
		&HHPaths::System("$hhscripts/multithread.pl '".$hhmglob."' '$command' -cpu $cpu");
	    } else { # $hhmext eq "hhm"
		print("\nGenerating prf profile files in $tmpdir/ from hhm files $hhmglob/\n\n");
		$command = "$hhscripts/create_profile_from_hhm.pl -i \$file -o $tmpdir/\$base.prf 1>/dev/null 2>>$logfile";
		&HHPaths::System("$hhscripts/multithread.pl '".$hhmglob."' '$command' -cpu $cpu");
	    }
	    
	    print("\nGenerating cs219 files in $tmpdir/ from prf files in $tmpdir/\n\n");
	    if ($hhmext eq "hmm") {
		$command = "$hhbin/cstranslate -i \$file -o \$name.cs219 -A $cs_lib 1>>$logfile 2>>$logfile";	    
	    } else { # $hhmext eq "hhm"
		$command = "$hhbin/cstranslate -i \$file -o \$name.cs219 -A $cs_lib -D $context_lib -x $x -c $c 1>>$logfile 2>>$logfile";
	    }
	    &HHPaths::System("$hhscripts/multithread.pl '".$tmpdir."/*.prf' '$command' -cpu $cpu");
	}
	print("\n");
	
	$csglob = $tmpdir."/*.$csext";
    }
        
    
##############################################################################################
# Generate hhm files with hhmake from a3m files if no -ihhm directory given
##############################################################################################
    
    if (!$hhmglob) 
    {
	if ($a3mglob) {
	    print("\nGenerating hhm files in $tmpdir/ from a3m files $a3mglob/\n\n");
	    $command = "hhmake -i \$file -o $tmpdir/\$base.hhm  1>/dev/null 2>>$logfile";
	    &HHPaths::System("$hhscripts/multithread.pl '".$a3mglob."' '$command' -cpu $cpu");	
	    $hhmglob = $tmpdir."/*.$hhmext";;
	    $numhhmfiles += scalar(glob("$hhmglob"));
	    print("\n");
	}
    }
    
    
##############################################################################################
# Generate packed ffdata files, ffindex and ffsize files for A3M, HMM, and cs files 
##############################################################################################
    
    if ($a3mglob ne "") {
	
	&GenerateFfindexDatabaseFiles($a3mdata,$a3mindex,$a3msize,$a3mglob);
    } 
    
    if ($hhmglob ne "") {
	&GenerateFfindexDatabaseFiles($hhmdata,$hhmindex,$hhmsize,$hhmglob);
    }

    if (1) {
	&GenerateFfindexDatabaseFiles($csdata,$csindex,$cssize,$csglob);
    }

} # end if $remove==0


##############################################################################################
# If $csglob is simple directory instead of glob expression, turn it into glob expression
##############################################################################################
sub TurnDirIntoGlob
{
    my $glob = $_[0];
    my $ext = $_[1];
    if ($glob) {
	if ($glob !~ /\*/ && $glob !~ /\?/ && $glob !~ / /) { 
	    $glob .= "/*.".$ext;
	}
    }
    return $glob;
}


##############################################################################################
# Unlink files in A3M, HHM, or cs index file 
##############################################################################################
sub UnlinkFilesFromFfindex
{
    my $ffindex = shift @_;
    my $ffsize = shift @_;
    my $remove_expr = shift @_;
    my $numfiles_to_unlink = 0;
    my $files_to_unlink = "";
    my $numfiles = 0;

    open (IN, "<$ffindex") || die("Error: can't open $ffindex: $!");
    while($line=<IN>) {
	if ($line =~ /^$remove_expr/) { 
	    $line =~ /(\S*)/; 
	    $files_to_unlink .= " ".$1;
	    $numfiles_to_unlink++;
	}
	else {$numfiles++;}
    }
    close IN;
    if ($v>=2) {	
#	printf("remove_expr = %s\n",$remove_expr);
	printf("Removing %i files from $ffindex => new number of files will be %i\n",$numfiles_to_unlink,$numfiles);}
    if ($numfiles_to_unlink > 0) {
	&HHPaths::System("ffindex_modify -su $ffindex ".$files_to_unlink); 
    }
    return $numfiles_to_unlink;
}

##############################################################################################
# Unlink files in A3M, HHM, or cs index file 
##############################################################################################
sub GenerateFfindexDatabaseFiles
{
    my $ffdata = shift @_;
    my $ffindex = shift @_;
    my $ffsize = shift @_;
    my $glob  = shift @_;
    my $numremoved = 0;
    print "Creating database files $ffdata, $ffindex, and $ffsize ...\n";

    my @files = glob("$glob");
    my $numfiles = scalar(@files);
    open (OUT, ">$tmpdir/filelist.txt");
    for (my $i=0; $i< scalar(@files); $i++) {
	printf(OUT "%s\n",$files[$i]);
	$files[$i] =~ s/.*\///;  # remove everything up to last /
	$files[$i] =~ s/.*\\//;  # remove everything up to last \
	$files[$i] =~ s/\./\./g; # . => \.
    }
    close OUT;
    
    # Build packed file (concatenated with '\0' as delimiters) and index file from files in file list
    # The ffindex binaries are contained in <install_dir>/lib/ffindex/bin/
    if ($a_if_append) {

	my $remove_expr = join("|",@files);
	$numremoved = &UnlinkFilesFromFfindex($ffindex,$ffsize,$remove_expr);
	# printf("numfiles=%i  numremoved=%i\n",$numfiles,$numremoved);

   }

    $command = "ffindex_build -".$a_if_append."s -f $tmpdir/filelist.txt $ffdata $ffindex";
    &HHPaths::System($command);

    if ($numremoved>0) { 
	printf("Warning: %i files in $ffindex were overwritten by appended files with identical names\n",$numremoved);
    }
 
    print("\n");
    return $numfiles;
}
