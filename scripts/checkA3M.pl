#!/usr/bin/perl

use strict;
use warnings;

if(scalar @ARGV eq 0 or $ARGV[0] eq "-h" or $ARGV[0] eq "--help") {
  print "Checks the sanity of an a3m file!\n";
  print "USAGE: ./checkA3M.pl [a3mfile]\n";
  exit(0);
}

my $infile = $ARGV[0];

my $line_nr = 0;

my $id = "";
my $nr = 0;
my $seq = "";
my $header_line = 0;

my $first_nr_matchstates = -1;

my $first_mismatch_id = 0;

open IN, "<$infile";
while(my $line = <IN>) {
  $line_nr++;

  if($line =~ /^>(\S*)\s*/) {
    #process old data
    if($seq ne "") {
      my $matchstates = &countMatchStates($seq);
      if($first_nr_matchstates eq -1) {
        $first_nr_matchstates = $matchstates;
      }

      if($first_nr_matchstates ne $matchstates) {
        print "Error: Mismatching Number of Match States!";
        print "\tFirst Sequence has $first_nr_matchstates matchstates!\n";
        print "\tSeqence Nr. $nr ($id, line $header_line) has $matchstates matchstates!\n";
      }

      my @invalid_chars = @{&checkValidAlphabet($seq)};
      @invalid_chars = @{&uniq2(@invalid_chars)};
      @invalid_chars = @{addQuotes(@invalid_chars)};

      if(scalar @invalid_chars ne 0) {
        print "Error: Invalid characters in Seqence!\n";
        print "\tSequence Nr. $nr ($id, line $header_line) !\n";
        print "\tFound invalid characters ".join(",", @invalid_chars)."!\n";
      }
    }

    $header_line = $line_nr;
    $id = $1;
    $nr++;
    $seq = "";
  }
  elsif($line =~ /^#/) {
    #do nothing
  }
  else {
    chomp($line);
    $seq .= $line;
  }
}

if($seq ne "") {
  my $matchstates = &countMatchStates($seq);
  if($first_nr_matchstates eq -1) {
    $first_nr_matchstates = $matchstates;
  }

  if($first_nr_matchstates ne $matchstates) {
    print "Error: Mismatching Number of Match States!";
    print "\tFirst Sequence has $first_nr_matchstates matchstates!\n";
    print "\tSeqence Nr. $nr ($id, line $header_line) has $matchstates matchstates!\n";
  }

  my @invalid_chars = @{&checkValidAlphabet($seq)};
  @invalid_chars = @{&uniq2(@invalid_chars)};
  @invalid_chars = @{addQuotes(@invalid_chars)};

  if(scalar @invalid_chars ne 0) {
    print "Error: Invalid characters in Seqence!\n";
    print "\tSequence Nr. $nr ($id, line $header_line) !\n";
    print "\tFound invalid characters ".join(",", @invalid_chars)."!\n";
  }
}

close IN;



if($line_nr eq 0) {
  print STDERR "Error: $infile is empty!";
}

if($nr eq 0) {
  print STDERR "Error: $infile contains no headers/sequences!";
}

sub countMatchStates{
  my $seq = $_[0];
  my $matchstates = 0;
  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if($char eq uc($char) or $char eq "-") {
      $matchstates++;
    }
  }

  return $matchstates;
}

sub checkValidAlphabet{
  my $seq = $_[0];
  my @invalid_characters = ();

  my %valid_characters = (
   '-', 0, 
   '.', 0, 
   'A', 0, 
   'B', 0, 
   'C', 0, 
   'D', 0, 
   'E', 0, 
   'F', 0,  
   'G', 0, 
   'H', 0, 
   'I', 0, 
   'J', 0, 
   'K', 0, 
   'L', 0, 
   'M', 0, 
   'N', 0, 
   'O', 0, 
   'P', 0, 
   'Q', 0, 
   'R', 0, 
   'S', 0, 
   'T', 0, 
   'U', 0, 
   'V', 0, 
   'W', 0, 
   'X', 0, 
   'Y', 0, 
   'Z', 0, 
   'a', 0, 
   'b', 0, 
   'c', 0, 
   'd', 0, 
   'e', 0, 
   'f', 0, 
   'g', 0, 
   'h', 0, 
   'i', 0, 
   'j', 0, 
   'k', 0, 
   'l', 0, 
   'm', 0, 
   'n', 0, 
   'o', 0, 
   'p', 0, 
   'q', 0, 
   'r', 0, 
   's', 0, 
   't', 0, 
   'u', 0, 
   'v', 0, 
   'w', 0, 
   'x', 0, 
   'y', 0, 
   'z', 0
  );


  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if(! defined $valid_characters{$char} ) {
      push(@invalid_characters, $char);
    }
  }

  return \@invalid_characters;
}

#code from stackoverflow
sub uniq2 {
  my %seen = ();
  my @r = ();
  foreach my $a (@_) {
    unless ($seen{$a}) {
      push @r, $a;
      $seen{$a} = 1;
    }
  }
  return \@r;
}

sub addQuotes {
  my @r = ();
  foreach my $a (@_) {
    push @r, "'".$a."'";
  }

  return \@r;
}

