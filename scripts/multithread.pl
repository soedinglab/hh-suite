#! /usr/bin/perl -w
# multithread.pl: 
# Run a command with different file names as arguments on multiple threads in parallel 
# Usage:   multithread.pl <fileglob> '<command>' [-cpu <int>] 
# (C) Johannes Soeding
# Available under GNU Public License Version 3
#
# Reference: 
# Remmert M., Biegert A., Hauser A., and Soding J.
# HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
# Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

use strict;
use POSIX;
$|=1; # autoflush on

# Variables 
my $cpu=4;        # number of cpus to use
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
    print("
 Run a command with different file names as arguments on multiple threads in parallel 
 Usage: multithread.pl '<fileglob>' '<command>' [-cpu <int>] [-v {0,1,2}]

 <command> can include symbol 
    \$file for the full filename, e.g. /tmp/hh/1c1g_A.a3m, 
    \$name the filename without extension, e.g. /tmp/hh/1c1g_A, and 
    \$base for the filename without extension and path, e.g. 1c1g_A.

  -cpu <int>  number of threads to launch (default = $cpu)
  -v {0,1,2}  verbose mode (default = $v)

 Example: multithread.pl '*.a3m' 'hhmake -i \$file 1>\$name.log 2>>error.log' -cpu 16
\n"); 
    exit;
}

my @files=glob($ARGV[0]); 
my $command=$ARGV[1]; 
$SIG{'CHLD'}='IGNORE';
$SIG{'USR1'}=\&ChildFinished;
$SIG{'INT'} =\&KillAllProcesses;

if (@ARGV>2) {
    $options.=join(" ",@ARGV[2..$#ARGV]);
}
# Set number of cpus to use
if ($options=~s/-cpu\s*(\d)\s*//g) {$cpu=$1;}
if ($options=~s/-v\s*(\d)\s*//g) {$v=$1;}

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
wait;
if ($v>=1) {print ("\nAll processes should be finished now\n");}
if ($numerr>0) {print(STDERR "WARNING: $numerr commands returned with error code >0\n");}
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

sub System {
    if ($v>=2) {print("\n");} 
    if ($v>=1) {print("\n".$_[0]," ");}
    if (system($_[0])) {
# Why is always -1 returned??
#	print(STDERR "\nERROR: command '$command' returned error code $?\n");
    };
}

