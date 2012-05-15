#!/usr/bin/env perl
#
# multithread.pl: 
# Run a command with different file names as arguments on multiple threads in parallel 
# Usage:   multithread.pl <fileglob> '<command>' [-cpu <int>] 
#
#
#     HHsuite version 2.0.14 (May 2012)
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

#     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use strict;
use POSIX;

# Variables 
my $cpu=8;         # number of cpus to use
my $parent_pid=$$; # main process id
my $pid;           # process id of child
my %pid=();        # hash has all running PIDs as keys and the file name as data
my $children=0;    # number of child processes running
my $options="";
my $file;
my $ifile=0;
my $v=1;
my $numerr=0;

if (scalar(@ARGV)<2) {
    die("
multithread.pl from HHsuite $VERSION  
Run a command for many files in parallel using multiple threads
Usage: multithread.pl '<fileglob>' '<command>' [-cpu <int>] [-v {0,1,2}]

<command> can include symbol 
   \$file for the full filename, e.g. /tmp/hh/1c1g_A.a3m, 
   \$name the filename without extension, e.g. /tmp/hh/1c1g_A, and 
   \$base for the filename without extension and path, e.g. 1c1g_A.

 -cpu <int>  number of threads to launch (default = $cpu)
 -v {0,1,2}  verbose mode (default = $v)

Example: multithread.pl '*.a3m' 'hhmake -i \$file 1>\$name.log 2>>error.log' -cpu 16
\n"); 
}

$|=1; # autoflush on


my @files=glob($ARGV[0]); 
my $command=$ARGV[1]; 
$SIG{'CHLD'}='IGNORE';
$SIG{'USR1'}=\&ChildFinished;
$SIG{'INT'} =\&KillAllProcesses;

if (@ARGV>2) {
    $options.=join(" ",@ARGV[2..$#ARGV]);
}
# Set number of cpus to use
if ($options=~s/-cpu\s*(\d+)\s*//g) {$cpu=$1;}
if ($options=~s/-v\s*(\d+)\s*//g) {$v=$1;}

# Warn if unknown options found
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; print("WARNING: unknown options '$options'\n");}

if ($v>=1) {print (scalar(@files)." files read in ...\n");}

foreach $file (@files) {   
    $ifile++;

    # All cpus occupied? -> wait for a cpu to become free
    if ($children>=$cpu) {
	if ($v>=2) {print("\nParent $$ is sleeping (children=$children) ");}
	my $count=0;
	while ($children>=$cpu) {
	    if ($count++>=10) {
		$count=0;
		if ($v>=2) {print("\nProcesses running:");}
		$children=0;
		foreach $pid (keys(%pid)) { 
		    if (! kill(0,$pid)) {    # kill($pid,0) returns false if process is dead (finished)
			if ($v>=2) {printf("\nPID %5.5s: %s is removed from process table",$pid,$pid{$pid});}
			delete($pid{$pid});  # remove process from hash of PIDs
		    } else {
			if ($v>=2) {printf("\nPID %5.5s: %s",$pid,$pid{$pid});}
			$children++;   # In case a USR1 signal was caught twice (??)
		    }
		}
		if ($v>=2) {print("\n");}
	    } else {
		if ($v==1) {print(".");}
		select(undef, undef, undef, 0.1); # sleep 0.1 seconds
	    }
	} 
    }
    
    if ($pid=fork()) {
	# Main process
	$children++;
	$pid{$pid}="$file ($ifile)";

	# Print out running processes and remove defunct ones
	select(undef, undef, undef, 0.1); # sleep 0.1 seconds

    } elsif (defined $pid) {
	# Child process
	my $name;   # filename without extension
	my $base;   # basename without path
	if ($file =~/(.*)\..*?$/) {$name=$1;} else {$name=$file;}  
	if ($name =~/.*\/(.*?)$/) {$base=$1;} else {$base=$name;} 
	my $lcommand = $command; # need local variable for thread
	$lcommand=~s/\$file/$file/g;
	$lcommand=~s/\$name/$name/g;
	$lcommand=~s/\$base/$base/g;
	&System("$lcommand");
	if ($v>=2) {printf("\nProcess $$ for file %s (%i) finished.",$file,$ifile);}
	kill(USR1 => $parent_pid);
	$SIG{'CHLD'}='IGNORE';
	exit;
    } else {
	die("\nError: fork returned undefined PID: $!\n");
    }
}

# Wait for all children to finish
while (wait() != -1) {}

if ($v>=1) {print ("\nAll processes should be finished now\n");}
if ($numerr>0) {print(STDERR "WARNING: $numerr commands returned with error code.\n");}
exit(0);


sub ChildFinished() {
    $children--;
    $SIG{'USR1'}=\&ChildFinished;
    if ($v>=2) {printf("\nChildren counter reduced to children=$children",$file,$ifile);}
    return;
}

sub KillAllProcesses() 
{
    foreach $pid (keys(%pid)) { 
	if ($v>=2) {printf("\nKill process $pid: returned %i\n",kill(-9,$pid));}
    }
    die ("\nInterrupt: Killed main process $$\n");
}

################################################################################################
### System command
################################################################################################
sub System {
    if ($v>=2) {print("\n");} 
    if ($v>=1) {print("\n".$_[0]," ");}
    if (system($_[0])) {
# Why is always -1 returned???
#       print(STDERR "\nERROR: command '$command' returned error code $?\n");
#       $numerr++;
    };
}
