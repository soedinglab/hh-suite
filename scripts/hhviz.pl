#!/usr/bin/env perl
#
# @author Andreas Biegert
# @time 2005-05-20
#
my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";

use MyPaths;
use warnings;
use GD;

my $resfile;
my $imgfile;
my $htmlfile;
my @hits;
my %hhcluster;
my $id;
my $line;
my $matchcols;
my $xscale;
my $im;
my $ypos;
my $x1;
my $y1;
my $x2;
my $y2;
my $white;
my $black;
my $red;
my $blue;
my $col;
my $leftopen;
my $rightopen;
my @rows;
my $row;
my $basedir;
my $rgbval;
my $probcutoff=40;
my $imgwidth=800;
my $imgheight;
my $barheight=10;
my $barsep=8;
my $border=10;
my $headsep=25;
my $headheight=5;
my $edgerad=3;
my $scraplen=5;
my $marklen=5;
my $imgdir;

# check arguments
if((scalar @ARGV)!=3){
    print "Usage: hhviz.pl ID BASEDIR IMGDIR\n";
    exit(1);
}

$id = $ARGV[0];
$basedir = $ARGV[1];
$imgdir = $ARGV[2];
$resfile=$id.".hhr";
$htmlfile = $id.".html";	
$imgfile = $id.".png";


# Read hit list
open (RESFILE,"<$basedir/$resfile") or die "unable to open hhrfile $resfile for reading\n";
while($line=<RESFILE>){
	chomp $line;
	if ($line =~ /Match_columns\s+(\d+)/) {
	    # From header of results file: length of query, needed further down
	    $matchcols=$1;
	} 
	elsif ($line=~/^\s*(\d+)\s+(\S.*?)\s+(\d+\.\d+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)-(\d+)\s+(\d+)-(\d+)\s*\((\d+)\)\S*$/)
	{
	    if ($3 >= $probcutoff) {
		my %tmp;
		my $m=$1;
		my $header = $2;
		my @header_words = split(' ', $header);		
		my $template="";
		my $keyword="";
		$template = $header_words[0];		
		if (scalar(@header_words) > 1) {
			$keyword = $header_words[1];
		}		
		my $prob=$3;
		my $Eval=$4;
		my $qbeg=$5;
		my $qend=$6;
		my $tbeg=$7;
		my $tend=$8;
		my $tms=$9;
		$tmp{'m'}=$1;
		if ($template =~ /^[defgh](\d[a-z0-9]{3})[a-z0-9_.][a-z0-9_]$/) {
		    # SCOP id? show family code
		    $tmp{'tag'}=$keyword;
		} elsif ($template =~ /^PF\d{4}/ || $template=~/pfam\d{4}/) {
		    # Pfam id? show short description
		    $tmp{'tag'}=$keyword;
		} elsif ($template =~ /^SM\d{4}/ || $template=~/smart\d{4}/) {
		    # SMART id? show short description
		    $tmp{'tag'}=$keyword;
		} elsif ($template =~ /^COG\d{4}/) {
		    # COG id? show short description
		    $tmp{'tag'}=$keyword;
		} else {
		    # All others: show id 
		    $tmp{'tag'}=$template;
		}
		$tmp{'title'}="Prob=$prob\% E=$Eval  NaMeLiNe";
		$tmp{'template'}=$template;
		$tmp{'prob'}=$prob;
		$tmp{'qbeg'}=$qbeg;
		$tmp{'qend'}=$qend;
		$tmp{'tbeg'}=$tbeg;
		$tmp{'tend'}=$tend;
		$tmp{'tms'}=$tms;
		$row=&calcRow($qbeg,$qend,\@rows);
		&putBar($row,$qbeg,$qend,\@rows);
		$tmp{'row'}=$row;
		push(@hits, \%tmp);
#		print "m=$m descr=".$hits[$m-1]{'title'}."\n";
	    }
	}
	
	# Start of next alignments found
	elsif ($line=~/^No (\d+)/) {
	    my $m=$1;
	    if (defined $hits[$m-1]) {
		$line=<RESFILE>;
		$line=~/^>(.{1,80})/;
		my $nameline=$1;
		$nameline=~tr/<>//d;
		${$hits[$m-1]}{'title'}=~s/NaMeLiNe/$nameline/;
	    }
	} 
}
close RESFILE;

$xscale=($imgwidth-2*$border)/$matchcols;
$imgheight=2*$border+$headheight+$headsep+scalar(@rows)*($barsep+$barheight)-$barheight; 
#print hits
foreach my $i (@hits) {
    my %hit = %$i;
    print $hit{'m'}."\t".$hit{'template'}."\t".$hit{'prob'}."\t".$hit{'qbeg'}."\t".$hit{'qend'}."\t".$hit{'tbeg'}."\t".$hit{'tend'}."\t".$hit{'tms'}."\t".$hit{'row'}."\n"; 
}

# visualize hits
$im = new GD::Image($imgwidth,$imgheight);
$white = $im->colorAllocate(255,255,255);
$black = $im->colorAllocate(0,0,0);       
$red = $im->colorAllocate(255,0,0);      
$blue = $im->colorAllocate(0,0,255);
$im->transparent($white);

$x1=$border+$xscale;
$y1=$border;
$x2=$border+$matchcols*$xscale;
$y2=$border+$headheight;
$im->filledRectangle($x1,$y1,$x2,$y2,$black);

for(my $i=10;$i<=$matchcols;$i+=10) {
    if ($i%100==0) {
	my $xoffset=length($i)*(5/2);
	$im->line($border+$i*$xscale,$border,$border+$i*$xscale,$border+$headheight+2*$marklen,$black);
	$im->string(gdSmallFont,$border+$i*$xscale-$xoffset,$border+$headheight+2*$marklen,$i,$black);
    } else {
	$im->line($border+$i*$xscale,$border,$border+$i*$xscale,$border+$headheight+$marklen,$black);
    }
}

open(HTMLFILE,">$basedir/$htmlfile") or die "unable to open htmlfile $htmlfile for writing\n";
print HTMLFILE "<map name=\"hhpredmap\">\n";

for(my $i=0; $i<scalar(@hits); $i++) {
    my %hit = %{$hits[$i]};
    print "Processing hit No ".$hit{'m'}." with template '".$hit{'template'}."' ...\n";
    
    $ypos = $border+$headheight+$headsep+$hit{'row'}*($barheight+$barsep);
    $x1=$border+$xscale*$hit{'qbeg'};
    $y1=$ypos;
    $x2=$border+$xscale*$hit{'qend'};
    $y2=$ypos+$barheight;
    
    &printMapEntry(HTMLFILE,$hit{'m'},$hit{'title'},$x1,$y1,$x2,$y2);
    
    if ($hit{'tbeg'}>1) {
	$leftopen=1;
    } else {
	$leftopen=0;
    }
    if ($hit{'tend'}-$hit{'tbeg'}<$hit{'tms'}-$hit{'tbeg'}) {
	$rightopen=1;
    } else {
	$rightopen=0;
    }
    print $hit{'tag'}."\n";
#    print "( $x1 / $y1 )  ( $x2 / $y2 )\n";
#    $rgbval=$basergb+($hit{'prob'}-$probcutoff)*$rgbscale;
#    $col = $im->colorAllocate(255-$rgbval,255-$rgbval,$rgbval);
#    &drawBar($x1,$y1,$x2,$y2,$col,$leftopen,$rightopen);

    my $minsat=0.8;           # minimum saturation
    my $signif = ($hit{'prob'}-$probcutoff)/(100-$probcutoff);
    my $satur  = 1.0-(1.0-$minsat)*(1.0-$signif);  
    my ($red,$grn,$blu);
    my $col=4*$signif;
    if ($col>3) {
	# red (1,0,0)-> dark yellow (1,0.7,0) transition
	$col-=3.0;
	$red = 1;
	$grn = 0.7*(1-$col);
	$blu = 0;
    } elsif ($col>2) {
	# dark yellow (1,0.7,0) -> green (0,1,0) transition
	$col-=2.0;
	$red = $col;
	$grn = 0.7+0.3*(1-$col);
	$blu = 0;
    } elsif ($col>1) {
	# green (0,1,0) -> cyan (0,0.7,1) transition
	$col-=1.0;
	$red = 0;
	$grn = 1-0.3*(1-$col);
	$blu = 1-$col;
    } else {
	# cyan (0,0.7,1) -> blue (0,0,1) transition
	$red = 0;
	$grn = 0.7*$col;
	$blu = 1;
    }
    
    $col = $im->colorAllocate(255*$red,255*$grn,255*$blu);
    &drawBar($x1,$y1,$x2,$y2,$col,$leftopen,$rightopen,$hit{'tag'}); # print tag (not template name)
}
print HTMLFILE "</map>\n";
print HTMLFILE "<p><img src=\"$imgdir/$imgfile\" border=\"0\" alt=\"hhpredhits\" usemap=\"#hhpredmap\"></p>";
close(HTMLFILE);

# write image to file
open (OUTFILE,">$basedir/$imgfile") or die "unable to open imgfile $imgfile for writing\n";
print OUTFILE $im->png;
close OUTFILE;

print "\nDone! Results written to '$htmlfile' and '$imgfile'.\n";

exit(0);

#############################################################################
# functions
sub drawBar() {
    my $x1 = $_[0];
    my $y1 = $_[1];
    my $x2 = $_[2];
    my $y2 = $_[3];
    my $col = $_[4];
    my $leftopen = $_[5];
    my $rightopen = $_[6];
    my $tag = $_[7];

    $im->filledRectangle($x1+$edgerad,$y1,$x2-$edgerad,$y2,$col);

    if ($leftopen) {
	# Ragged border: local alignment
	$poly = new GD::Polygon;
        $poly->addPt($x1+$edgerad,$y1);
        $poly->addPt($x1+$edgerad,$y2);
        $poly->addPt($x1+$edgerad-$scraplen,$y1+0.75*($y2-$y1));
	$poly->addPt($x1+$edgerad,$y1+0.5*($y2-$y1));
	$poly->addPt($x1+$edgerad-$scraplen,$y1+0.25*($y2-$y1));
	$im->filledPolygon($poly,$col);
    } else {
	# Rounded border: alignment up to end
	$im->filledArc($x1+$edgerad,$y1+$edgerad,2*$edgerad,2*$edgerad,0,360,$col);
	$im->filledArc($x1+$edgerad,$y2-$edgerad,2*$edgerad,2*$edgerad,0,360,$col);
	$im->filledRectangle($x1,$y1+$edgerad,$x1+$edgerad,$y2-$edgerad,$col);
    }
    if ($rightopen) {
	# Ragged border: local alignment
        $poly = new GD::Polygon;
        $poly->addPt($x2-$edgerad,$y1);
        $poly->addPt($x2-$edgerad,$y2);
        $poly->addPt($x2-$edgerad+$scraplen,$y1+0.75*($y2-$y1));
	$poly->addPt($x2-$edgerad,$y1+0.5*($y2-$y1));
	$poly->addPt($x2-$edgerad+$scraplen,$y1+0.25*($y2-$y1));
	$im->filledPolygon($poly,$col);
    } else {
	# Rounded border: alignment up to end
	$im->filledArc($x2-$edgerad,$y1+$edgerad,2*$edgerad,2*$edgerad,0,360,$col);
	$im->filledArc($x2-$edgerad,$y2-$edgerad,2*$edgerad,2*$edgerad,0,360,$col);
	$im->filledRectangle($x2-$edgerad,$y1+$edgerad,$x2,$y2-$edgerad,$col);
    }
    # add label
    if (length($tag)*6>($x2-$x1-2*$edgerad)) {
	$tag = substr($tag,0,($x2-$x1-2*$edgerad)/6);
    }
    my $xoffset=length($tag)*(6/2);
    $im->string(gdMediumBoldFont,(($x2+$x1)/2)-$xoffset,$y1-2,$tag,$white);
}

sub printMapEntry() {
    my $fh = $_[0];
    my $m = $_[1];
    my $mapstring = $_[2];
    my $x1 = $_[3];
    my $y1 = $_[4];
    my $x2 = $_[5];
    my $y2 = $_[6];

    if ($x1 =~ /(\d+)\.\d+/) {$x1 = $1}; 
    if ($y1 =~ /(\d+)\.\d+/) {$y1 = $1};
    if ($x2 =~ /(\d+)\.\d+/) {$x2 = $1};
    if ($y2 =~ /(\d+)\.\d+/) {$y2 = $1};
    print $fh "<area shape=\"rect\" coords=\"$x1,$y1,$x2,$y2\" href=\"#$m\" title=\"$mapstring\" />\n";
}

# Data structure:
# Each hit is recorded with its residue range in an array @bar=($beg,$end)
# The addresses of several @bar's may be contained in a @row. Each @row contains at least one @bar.
# The adresses of all occupied @row's are contained in array @rows.

# Determine index of @rows where to place new hit @bar=($beg,$end)
sub calcRow() {
    my $beg=$_[0];      # first residue of residue range
    my $end=$_[1];      # last residue 
    my $rowsref=$_[2];
    my $b;
    my $e;
    my $fits;
    my $bar;
    my @bars;

    for(my $r=0; $r<scalar(@$rowsref); $r++) {
	@bars=@{${$rowsref}[$r]};
	$fits=1;
	foreach $bar (@bars) {
	    ($b,$e)=@$bar;
	    if (!($beg>$e || $end<$b)) {
		$fits=0;
		last;
	    } 
	}
	if ($fits) {
	    return $r;
	}
    }
    return scalar(@$rowsref);
}

# Insert @bar=($beg,$end) for current hit into @rows
sub putBar() {
    my $pos=$_[0];      # row in which to insert @bar=($beg,$end)
    my $beg=$_[1];      # first residue of residue range
    my $end=$_[2];      # last residue
    my $rowsref=$_[3];  # data structure @rows
    my @bar=($beg, $end);
    
    # Do we place @bar into new row of @rows?
    if ($pos >= scalar(@$rowsref)) {
	my @newrow=(\@bar);
	${$rowsref}[$pos]=\@newrow;
    } else {
	my $rowref=$rows[$pos]; # get address of @row with index $pos 
	push(@$rowref, \@bar);  # add \@bar to @row
    }
    return;
}


