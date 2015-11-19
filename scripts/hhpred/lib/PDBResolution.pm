#!/user/bin/perl -w
use strict;

package PDBResolution;
require Exporter;
require File::Temp;

use config;
my $config = HHpredConfig->instance();

our @ISA = qw(Exporter);
our @EXPORT = qw(extractAllPDBResolutions getTMscore getResolutionFor);




## extract resolutions for all proteins in .fas file "fasFile"
## and save them in hash %results
sub extractAllPDBResolutions {
    my $pdbdb = $config->get_pdbdir() . "/db/pdb.hhm";

    my %result;

    ## read pdb database file and extract resolutions
    open (FIN, "< $pdbdb") or die "Cant open $pdbdb: $!\n";
    while (my $line = <FIN>) {
    	my $resolution = 10;
    	my $ID = "dummy";
	    if ($line =~ /^NAME\s+(\S+).*;\s(\d+\.\d\d?\d?)A\s/) {
			$ID = $1;
			$resolution = $2;
		} elsif ($line =~ /^NAME\s+(\S+).*;\s(\d+)A\s/) {
			$ID = $1;
			$resolution = $2;
		} elsif ($line =~ /^NAME\s+(\S+).*;\sNMR\s/) {
			$ID = $1;
			$resolution = 10;			
		}
		$result{$ID} = $resolution;
    } 
    
    close(FIN);

    return %result;
}


## read resolution of protein with id "id" from
## resolutions file (each row has the format: <id> <resolution>)
## and return it
sub getResolutionFor {
    my $id = shift;
    my $resolFile = shift;

    my $resol = 10;

    open (RES, "< $resolFile") or die "Cant open $resolFile: $!\n";
    while(my $line = <RES>) {
    if ($line =~ /^(\S+)\s+(\S+)/) {
        if ($1 eq $id) {
        $resol = $2;
        last;
        }
    }
    }
    close(RES);

    return $resol;
}



sub getTMscore {
    my $modelFile = shift;
    my $nativeFile = shift;

    my $tmscore = $config->get_TMscore();

    my $fh = File::Temp->new();
    my $outfile = $fh->filename;
    my $TM = 0;

    system("$tmscore $modelFile $nativeFile > $outfile");

    open (TM, "< $outfile") or die "Cant open $outfile: $!\n";
    while(my $line = <TM>) {
    if ($line =~ /^TM-score\s+=\s+(\S+)/) {
        $TM = $1;
        last;
    }
    }
    close(TM);
    return $TM;
}

1;
