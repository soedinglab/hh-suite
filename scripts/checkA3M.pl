#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $infile = "";
my $silent = 0;
my $help = 0;

Getopt::Long::Configure(qw{no_auto_abbrev no_ignore_case_always});

GetOptions (
  "i=s" => \$infile,
  "silent|V" => \$silent,
  "help|h" => \$help);

if($help or $infile eq "") {
  print("USAGE: checkA3m.pl -i [a3m_file|stdin] [-V|--silent] [-h|--help]\n");
  exit(0);
}

my $line_nr = 0;

my $EXIT_VALUE = 0;
my $a3m_id = "";
my $id = "";
my $nr = 0;
my $nr_consensus = 0;
my $seq = "";
my $header_line = 0;

my $first_nr_matchstates = -1;

my $first_mismatch_id = 0;

my @lines;
if($infile eq "stdin") {
  @lines = <STDIN>;
}
else {
  if(! -e $infile) {
    print STDERR "Input file ($infile) does not exist!\n";
    exit(2);
  }

  open IN, "<$infile";
  @lines = <IN>;
  close IN;
}

foreach my $line(@lines) {
  $line_nr++;

  if($line =~ /^>(\S*)\s*/) {
    #process old data
    if($seq ne "") {
      if($id eq "ss_pred") {
        #get length
        my $matchstates = &countSSPredStates($seq);
        if($first_nr_matchstates eq -1) {
          $first_nr_matchstates = $matchstates;
        }

        if($matchstates eq 0) {
          if(not $silent) {
            print STDERR "Error: Empty ss_pred in $infile!\n";
            print STDERR "\t\n"
          }
        }

        if($first_nr_matchstates ne $matchstates) {
          if(not $silent) {
            print STDERR "Error: Mismatching Number of Match States!\n";
            print STDERR "\tFirst Sequence has $first_nr_matchstates matchstates!\n";
            print STDERR "\tSeqence Nr. $nr ($a3m_id, $id, line $header_line) has $matchstates matchstates!\n";
          }
          $EXIT_VALUE = 1;
        }

        #check alphabet
        my @invalid_chars = @{&checkValidSSPredAlphabet($seq)};
        @invalid_chars = @{&uniq2(@invalid_chars)};
        @invalid_chars = @{addQuotes(@invalid_chars)};

        if(scalar @invalid_chars ne 0) {
          if(not $silent) {
            print STDERR "Error: Invalid characters in Seqence!\n";
            print STDERR "\tSequence Nr. $nr ($a3m_id, $id, line $header_line) !\n";
            print STDERR "\tFound invalid characters ".join(",", @invalid_chars)."!\n";
          }
          $EXIT_VALUE = 1;
        }
      }
      elsif($id eq "ss_conf") {
        #get length
        my $matchstates = &countSSConfStates($seq);
        if($first_nr_matchstates eq -1) {
          $first_nr_matchstates = $matchstates;
        }

        if($matchstates eq 0) {
          if(not $silent) {
            print STDERR "Error: Empty ss_conf in $infile!\n";
          }
        }

        if($first_nr_matchstates ne $matchstates) {
          if(not $silent) {
            print STDERR "Error: Mismatching Number of Match States!\n";
            print STDERR "\tFirst Sequence has $first_nr_matchstates matchstates!\n";
            print STDERR "\tSeqence Nr. $nr ($a3m_id, $id, line $header_line) has $matchstates matchstates!\n";
          }
          $EXIT_VALUE = 1;
        }

        #check alphabet
        my @invalid_chars = @{&checkValidSSConfAlphabet($seq)};
        @invalid_chars = @{&uniq2(@invalid_chars)};
        @invalid_chars = @{addQuotes(@invalid_chars)};

        if(scalar @invalid_chars ne 0) {
          if(not $silent) {
            print STDERR "Error: Invalid characters in Seqence!\n";
            print STDERR "\tSequence Nr. $nr ($a3m_id, $id, line $header_line) !\n";
            print STDERR "\tFound invalid characters ".join(",", @invalid_chars)."!\n";
          }
          $EXIT_VALUE = 1;
        }
      }
      else {
        my $matchstates = &countMatchStates($seq);
        if($first_nr_matchstates eq -1) {
          $first_nr_matchstates = $matchstates;
        }

        if($first_nr_matchstates ne $matchstates) {
          if(not $silent) {
            print STDERR "Error: Mismatching Number of Match States!\n";
            print STDERR "\tFirst Sequence has $first_nr_matchstates matchstates!\n";
            print STDERR "\tSeqence Nr. $nr ($a3m_id, $id, line $header_line) has $matchstates matchstates!\n";
          }
          $EXIT_VALUE = 1;
        }

        my @invalid_chars = @{&checkValidAlphabet($seq)};
        @invalid_chars = @{&uniq2(@invalid_chars)};
        @invalid_chars = @{addQuotes(@invalid_chars)};

        if(scalar @invalid_chars ne 0) {
          if(not $silent) {
            print STDERR "Error: Invalid characters in Seqence!\n";
            print STDERR "\tSequence Nr. $nr ($a3m_id, $id, line $header_line) !\n";
            print STDERR "\tFound invalid characters ".join(",", @invalid_chars)."!\n";
          }
          $EXIT_VALUE = 1;
        }
      }
    }

    $header_line = $line_nr;
    $id = $1;
    if(length($a3m_id) == 0) {
      $a3m_id = $id;
      $a3m_id =~ s/_consensus//;
    }

    if($id =~ s/_consensus//) {
      $nr_consensus++;
    }

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
    if(not $silent) {
      print STDERR "Error: Mismatching Number of Match States!\n";
      print STDERR "\tFirst Sequence has $first_nr_matchstates matchstates!\n";
      print STDERR "\tSeqence Nr. $nr ($a3m_id, $id, line $header_line) has $matchstates matchstates!\n";
    }
    $EXIT_VALUE = 1;
  }

  my @invalid_chars = @{&checkValidAlphabet($seq)};
  @invalid_chars = @{&uniq2(@invalid_chars)};
  @invalid_chars = @{addQuotes(@invalid_chars)};

  if(scalar @invalid_chars ne 0) {
    if(not $silent) {
      print STDERR "Error: Invalid characters in Seqence!\n";
      print STDERR "\tSequence Nr. $nr ($a3m_id, $id, line $header_line) !\n";
      print STDERR "\tFound invalid characters ".join(",", @invalid_chars)."!\n";
    }
    $EXIT_VALUE = 1;
  }
}

close IN;



if($line_nr eq 0) {
  if(not $silent) {
    print STDERR "Error: $infile is empty!\n";
  }
  $EXIT_VALUE = 1;
}

if(($nr - $nr_consensus) eq 0) {
  if(not $silent) {
    print STDERR "Error: $infile contains no headers/sequences!\n";
  }
  $EXIT_VALUE = 1;
}

if($nr_consensus > 2) {
  if(not $silent) {
    print STDERR "Error: $infile contains several headers!\n";
  }
  $EXIT_VALUE = 1;
}

exit($EXIT_VALUE);




sub countMatchStates{
  my $seq = $_[0];
  my $matchstates = 0;
  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if(($char eq uc($char) or $char eq "-") and $char ne "\0") {
      $matchstates++;
    }
  }

  return $matchstates;
}

sub countSSConfStates{
  my $seq = $_[0];
  my $matchstates = 0;
  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if($char =~ /^\d$/) {
      $matchstates++;
    }
  }

  return $matchstates;
}

sub countSSPredStates{
  my $seq = $_[0];
  my $matchstates = 0;
  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if($char eq 'E' or $char eq 'C' or $char eq 'H') {
      $matchstates++;
    }
  }

  return $matchstates;
}

sub checkValidSSConfAlphabet {
  my $seq = $_[0];
  my @invalid_characters = ();

  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if(not $char =~ /^\d$/) {
      push(@invalid_characters, $char);
    }
  }

  return \@invalid_characters;
}

sub checkValidSSPredAlphabet {
  my $seq = $_[0];
  my @invalid_characters = ();

  my %valid_characters = (
   'E', 0,
   'C', 0,
   'H', 0
  );

  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if(! defined $valid_characters{$char} ) {
      push(@invalid_characters, $char);
    }
  }

  return \@invalid_characters;
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
   'z', 0,
  );


  for(my $i = 0; $i < length($seq); $i++) {
    my $char = substr($seq, $i, 1);
    if(! defined $valid_characters{$char} and $char ne "\0") {
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


