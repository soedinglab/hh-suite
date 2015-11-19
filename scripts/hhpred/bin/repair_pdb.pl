#! /usr/bin/perl -w
# Fix errors pdb-file generate by modeller
# Usage:   repair_pdb.pl pdb-file

use strict;

my @res;

# check arguments
if ( @ARGV < 1 ) {
	error();
	exit(1);
}

open( IN, "$ARGV[0]" ) or die("Cannot open!");
@res = <IN>;
close IN;


for (my $a = 0; $a < scalar(@res); $a++) {
    if ($res[$a] =~ /^(ATOM.{50})(.{6})/) {
	my $begin = $1;
	if ($2 ne "  1.00") {
	    $res[$a] =~ s/$begin.{6}/$begin  1.00/;
	}
    }
}


my $i = scalar(@res) - 2;

my $line = $res[$i];

$line =~ s/^ATOM(\s+\d+)\s+\S+\s+(\S+\s+\d+).*\s+(\S+)\s*$/TER $1      $2                                               $3/;

$line =~ /^\S+\s+(\d+)/;
my $tmp = $1;
my $tmp1 = $tmp + 1;
my $tmp2 = $tmp + 2;
$line =~ s/$tmp1/$tmp2/;
$line =~ s/$tmp/$tmp1/;

$i++;
$res[$i] = "$line\n";
$i++;
$res[$i] = "END";

open( OUT, ">$ARGV[0]" ) or die("Cannot open!");
print OUT @res;
close OUT;

print "FINISHED!!!!\n";

exit(0);

#####################################################
#### sub functions

sub error {
	print("Usage: repair_pdb.pl pdb-file \n");
	print("\n");
}
