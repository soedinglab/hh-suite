#!/usr/bin/env perl

package FixOccupancy;

use strict;
use warnings;

use File::Copy;

sub FixOccupancy {
	my $file = shift;

	open my $input, '<', $file or die "Could not open $file: $!";
	open my $output, '>', $file . '.occupancy_fix' or die "Could not open $file: $!";
	while ( my $line = <$input> ) {
		if ( $line !~ /^ATOM/ ) {
			print $output $line;
			next;
		}
		my $occupancy = substr($line, 54, 5);
		if ( $occupancy > 0.0 ) {
			print $output $line;
			next;
		}

		print $output substr($line, 0, 55) . " 0.01" . substr($line, 60, 21);
		
	}
	close $input;
	close $output;

	move($file . '.occupancy_fix', $file);
}

1;
