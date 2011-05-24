#!/usr/bin/perl
# 
# Create a HHblits CS-database from HMMER-files, HMM-files or A3M-files

################################################################################################################################
# Update the following variables

$script_dir = "/cluster/bioprogs/hhblits/scripts";  # path to directory with scripts (create_profile_from_hmmer.pl, create_profile_from_hhm.pl)
$lib_dir = "/cluster/bioprogs/hhblits";             # path to needed libraries (context_data.lib, cs219.lib)
$hh_dir = "/cluster/bioprogs/hhblits";              # path to needed HH-tools (cstranslate)

################################################################################################################################

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode
my $ext = "a3m";       # default file extension
my $x = 0.3;
my $c = 4;
my $cs_lib = "$lib_dir/cs219.lib";
my $context_lib = "$lib_dir/context_data.lib";
my $append = 0;

my $help="
Create a HHblits CS-database from HMMER-files, HMM-files or A3M-files

Usage: perl create_cs_db.pl -i <dir> [options]

Options:
  -i <dir>    Input directory with HMMER-, HMM- or A3M-files
  -o <file>   Output file for the CS-database (default: <indir>.cs_db)

  -ext <ext>  File extension, which identifies file-typ (default: $ext)
              - A3M-type  : a3m
              - HHM-type  : hhm
              - HMMER-type: hmmer or hmm

  -append     append to CS-database (default: overwrite file, if existing)

  -v [0-5]    verbose mode (default: $v)
\n";

# Variable declarations
my $line;
my $command;
my $indir;
my $outfile;
my $file;
my @files;

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

if ($options=~s/ -i\s+(\S+) //) {$indir=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}

if ($options=~s/ -ext\s+(\S+) //) {$ext=$1;}

if ($options=~s/ -append //) {$append=1;}

if ($options=~s/ -v\s+(\S+) //) {$v=$1;}

@files = glob("$indir/*.$ext");
if (scalar(@files) == 0) {print($help); print "ERROR! No files in $indir with extension $ext!\n"; exit(1);}

if ($append == 0 && -e $outfile) {
    system("rm $outfile");
}

if (!$outfile) {
    if ($indir =~ /^\S+\/(\S+?)$/) {
	$outfile = "$1.cs_db";
    } else {
	$outfile = "$indir.cs_db";
    }
    print("Create HHblits database with CS-database in file $outfile!\n");
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

print "Computing " . scalar(@files) . " files ...\n";
my $count = 0;

if ($ext eq "a3m") {

    foreach $file (@files) {

	$count++;
	if (($count % 1000) == 0) { print "$count\n"; }
	elsif (($count % 50) == 0) { print ". "; }

	# Create CS-sequence for CS-database
	$command = "$hh_dir/cstranslate -i $file -a $outfile -D $context_lib -A $cs_lib -x $x -c $c > /dev/null 2>&1";
	if (&System($command) != 0) {
	    print "WARNING! Error with command $command!\n";
	}
    }

} elsif ($ext eq "hhm") {

    foreach $file (@files) {

	$count++;
	if (($count % 1000) == 0) { print "$count\n"; }
	elsif (($count % 50) == 0) { print ". "; }

	# Extract profile from HHM-file
	$command = "$script_dir/create_profile_from_hhm.pl -i $file -o $tmpdir/file.prf > /dev/null 2>&1";
	if (&System($command) != 0) {
	    print "WARNING! Error with command $command!\n";
	}

	# Create CS-sequence for CS-database
	$command = "$hh_dir/cstranslate -i $tmpdir/file.prf -a $outfile -D $context_lib -A $cs_lib -x $x -c $c > /dev/null 2>&1";
	if (&System($command) != 0) {
	    print "WARNING! Error with command $command!\n";
	}
    }

} elsif ($ext eq "hmmer" || $ext eq "hmm") {

    foreach $file (@files) {

	$count++;
	if (($count % 1000) == 0) { print "$count\n"; }
	elsif (($count % 50) == 0) { print ". "; }

	# Extract profile from HMMER-file
	$command = "$script_dir/create_profile_from_hmmer.pl -i $file -o $tmpdir/file.prf > /dev/null 2>&1";
	if (&System($command) != 0) {
	    print "WARNING! Error with command $command!\n";
	}

	# Create CS-sequence for CS-database
	$command = "$hh_dir/cstranslate -i $tmpdir/file.prf -a $outfile -A $cs_lib > /dev/null 2>&1";
	if (&System($command) != 0) {
	    print "WARNING! Error with command $command!\n";
	}
    }

} else {
    print($help);
    print "ERROR! Unknown extension $ext!\n";
}

# Get number of sequences and characters from CS-database
$num_seqs = 0;
$num_chars = 0;

open (IN, "cat $outfile |grep \"^>\" |wc -l |");
$line = <IN>;
close IN;
$line =~ /^\s*(\d+)/;
$num_seqs = $1;

open (IN, "cat $outfile |grep -v \"^>\" |wc -c |");
$line = <IN>;
close IN;
$line =~ /^\s*(\d+)/;
$num_chars = $1;

open (OUT, ">$outfile.sizes");
print OUT "$num_seqs   $num_chars\n";
close OUT;

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
