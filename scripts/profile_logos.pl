#!/usr/bin/perl -w
# Builds histograms for HHpred
#
# @author Andreas Biegert
# @time 2005-06-08
# Modified by Johannes 15.7.05, 1.9.2005

my $rootdir;
BEGIN {
   if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib"; 

use MyPaths;
use warnings;
use GD;
use POSIX;
use strict;

my $v=3;     # verbose mode
my $resfile;
my @hits;
my $q_hmm_ref;
my @t_hmms;
my $id;        # id for files $id.hhr, $id.hhm, $id.tar.gz
my $idpng;     # id for output files $idpng_1.png, $idpng_1.tt, $idpng_1.map etc.
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID

# Create directory path
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 

my $line;
my $matchcols;
my $xscale;
my $im;         # main image with profile logos
my $im_tt;      # tooltip image with additional information (such as sequence)
my $ypos;
my $x1;
my $y1;
my $x2;
my $y2;
my $col;
my $rows;
my $basedir;
my $tarfile;
my $zipfile;
my $hhpred_dir="$database_dir/hhpred/new_dbs";   # geaendert (Johannes 15.7.05)
my $hhpred_genomes_dir="$database_dir/hhpred/genomes";   # geaendert (Michael  23.10.07)
my $hhcluster_dir="$database_dir/hhcluster/";    # geaendert (Johannes 31.8.05)
my $rgbval;
my $comphash_ref;

my $imgwidth;
my $imgheight;
my $imgdir;
my $probcutoff=0.1;
my $alnwidth=80;                 # how many columns are shown in one line
my $barwidth=9;                 # bar width in pixels 
my $hbarsep=0;                   # horizontal space between adjacent bars in pixels
my $tborder=10;                  # top border in pixels
my $rborder=40;                  # right border in pixels
my $bborder=10;                  # bottom border in pixels
my $lborder=50;                  # left border in pixels
my $barheight=100;               # maximal bar height in pixels
my $linesep=5;                   # vertical space between lines in pixels
my $qtsep=15;                    # vertical space between query and template graph
my $stdfontwidth=5;              # bar font width in pixel
my $stdfontheight=8;             # bar font height in pixel
my $stdfont=gdSmallFont;         # font used for bar labels
my $barfontwidth=6;              # bar font width in pixel
my $barfontheight=8;             # bar font height in pixel
my $barfont=gdMediumBoldFont;    # font used for bar labels
my $ttreach=11;

# colors
my $white;
my $black; 
my $darkgreen;
my $lightgreen;
my $lightred;
my $darkred;
my $orange;
my $yellow;
my $yellowgreen;
my $grey;
my $lightgrey;
my $blue;
my $lightpink;
my $darkviolett;
my $violett;
my $lightviolett;
my %aa_colors;
my %ss_colors;

my %aa_order = ( 'A' => 4,
		 'C' => 3,
		 'D' => 6, 
		 'E' => 6,
		 'F' => 2,
		 'G' => 5,
		 'H' => 8,
		 'I' => 1,
		 'K' => 9,
		 'L' => 1,
		 'M' => 1,
		 'N' => 7,
		 'P' => 10,
		 'Q' => 7,
		 'R' => 9,
		 'S' => 4,
		 'T' => 4,
		 'V' => 1,
		 'W' => 2,
		 'Y' => 2,
		 );

# check arguments
if((scalar @ARGV)<3){
    print "Usage: profile_logos.pl IDENT BASEDIR IMGDIR\n";
    exit(1);
}else{
    $id = $ARGV[0];
    $basedir = $ARGV[1];
    $imgdir = $ARGV[2];
    $idpng = $id;        
    $idpng =~s/\.\d+$//; # remove all digits preceded by a dot, e.g. 63347.30 => 63347
    $resfile=$id.".hhr";
    $tarfile=$id.".tar";
    $zipfile=$tarfile.".gz";
}

# Read hit list
open (RESFILE,"<$basedir/$resfile") or die "unable to open hhrfile $basedir/$resfile for reading\n";
my $m=0;
my $qsspred='';   # PSIPRED SS for query 
my $tsspred='';   # PSIPRED SS for template 
my $qssconf='';   # PSIPRED confidens values for query
my $tssconf='';   # PSIPRED confidens values for template
my $qssdssp='';   # query consensus sequence
my $tssdssp='';   # template consensus sequence
my $qcons='';     # query consensus sequence
my $tcons='';     # template consensus sequence
my $qseq='';      # query seed sequence
my $tseq='';      # template seed sequence
my $match='';     # match quality symbols
my $program=0;    # 1: hhalign, 0: hhsearch; important to know in which databases to search for HMM

my $hhblits_db = "";

my $createmodel = 0; # 0: normal search, 1: create model

# Move up to Match columns line
while ($line=<RESFILE>) { if ($line=~/^Match.columns/) {last;} }
$line=~/^Match.columns\s+(\d+)/;
$matchcols=$1;

# Read 'Command' line
while ($line=<RESFILE>) { if ($line=~/^Command/ || $line=~/^\s*$/) {last;} }
if ($line=~/-aliw\s+(\d+)/) {$alnwidth=$1;} else {$alnwidth=80;}
if ($line=~/hhalign/) {$program=1;} 
elsif ($line=~/hhblits/) {
    $program=2;
    if ($line =~ /-d (\S+) /) {
	$hhblits_db = $1;
    }
}
else {$program=0;}
if ($line=~/$id\.db\.hhm/) {$createmodel=1;} else {$createmodel=0;}

# Move up to first empty line before summary hit list
while ($line=<RESFILE>) { if ($line=~/^\s*$/  || $line=~/^\s*No Hit/) {last;} }


########################################################################################################
# Read summary hit list

while($line=<RESFILE>){

    chomp $line;
# No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
#  1 gi|401465|sp|P31489|YDA1_YEREN 100.0       0       0  343.9  24.5  411    9-439    23-453 (455)
#  2 gi|401465|sp|P31489|YDA1_YEREN 100.0       0       0  291.2  16.6  319    5-333    33-367 (455)
#                     No     Hit      Prob        Eval     Pval Scor  SS    Cols    beg - end     beg-end
    if ($line=~/^\s*(\d+)\s+(.*\S)\s+(\d+\.\d+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+)-(\d+)\s+(\d+)-(\d+)\s*\((\d+)\)/)
    {
	my %tmp;
	my $m=$1;
	my $template=$2;
	my $prob=$3;
	my $Eval=$4;
	my $qbeg=$5;
	my $qend=$6;
	my $tbeg=$7;
	my $tend=$8;
	my $tms=$9;
	my $keyword=$template;
	my $desc="";
	if ($template=~/^\S+\s+(\S+)/) {
	    $keyword=$1;
	    if ($template=~/^\S+\s+\S+\s+(\S+)/) {
		$desc=$1;
	    }
	}
	$template=~s/^(\S+).*/$1/; # take only first word from name in hit list
	
	my $db;        # HMM database from which match originates
	my @dirs;
	$tmp{'m'}=$m;	

	print "Template: $template!\n";

	if ($program==1) { 
	    # hhalign: look in job dir for HMM
	    $db = $basedir;
	    $template = $id;

	    # HMM could not be found?
	    if (!$db || ! -e "$db/$template.hhm") {
		print("Error: '$template.hhm' could not be found in $db");
		next;
	    }
	} elsif ($createmodel==1) {
		$db = $basedir;
		# HMM could not be found?
		if (!$db || ! -e "$db/$template.hhm") {
			print("Error: '$template.hhm' could not be found in $db");
			next;
		}
	} else {
	    # hhsearch: look in hhpred databases for HMM

	    if ($template=~/pfam\d{4}/ || $template=~/smart\d{4}/ || $template =~ /^cd\d{5}/ || $template =~ /^LOAD\d/) {
		# Domain from CDD? show short description
		$tmp{'tag'}=$keyword;
		@dirs=glob "$hhpred_dir/cdd_*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^PF\d{4}/) {
		# Pfam id? show short description
		$tmp{'tag'}=$keyword;
		@dirs=glob "$hhpred_dir/pfamA*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^SM\d{4}/) {
		# SMART id? show short description
		$tmp{'tag'}=$keyword;
		@dirs=glob "$hhpred_dir/smart*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^COG\d{4}/) {
		# COG id? show short description
		$tmp{'tag'}=$keyword;
		@dirs=glob "$hhpred_dir/COG*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^KOG\d{4}/) {
		# KOG id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/KOG*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^PD\d{6}/) {
		# ProDom id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/pfamB*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^PTHR\d{5}/) {
		# Panther id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/panther*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^TIGR\d{5}/) {
		# TIGRFAM id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/tigrfam*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^PIRSF\d{6}/) {
		# PIRSF id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/pirsf*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^SUPFAM\d+/) {
		#SUPERFAMILY id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/supfam*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^\d[a-z0-9]{3}[A-Z]\d/) {
		#CATH id? show id
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/CATH*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^[defgh](\d[a-z0-9]{3})[a-z0-9_.][a-z0-9_]$/) {
		# SCOP id? show family code
		$tmp{'tag'}=$keyword;
		@dirs=glob "$hhpred_dir/scop*";
		$db=pop(@dirs);
	    } elsif ($template =~ /^\d[a-z0-9]{3}(_[A-Za-z0-9])?$/) {
		# PDB id? show id 
		$tmp{'tag'}=$template;
		@dirs=glob "$hhpred_dir/pdb*";
		$db=pop(@dirs);
	    } elsif ($template =~ /cl\|(\S+?)\|/) {
		# HHblits identifier
		$template = $1;
		$db = "$tmpdir";
		if ($hhblits_db ne "") {
		    system("$rootdir/bioprogs/hhblits/ffindex_get $hhblits_db"."_hhm_db $hhblits_db"."_hhm_db.index $template.hhm > $tmpdir/$template.hhm");
		}
		open (IN, "$tmpdir/$template.hhm");
		my $tmp = <IN>;
		close IN;
		if ($tmp !~ /^HH/ && $hhblits_db ne "") {
		    system("$rootdir/bioprogs/hhblits/ffindex_get $hhblits_db"."_a3m_db $hhblits_db"."_a3m_db.index $template.a3m > $tmpdir/$template.a3m");
		    system("$hh/hhmake -i $tmpdir/$template.a3m -o $tmpdir/$template.hhm");
		}
		$tmp{'tag'}=$template;
	    } 

	    if (!$db || (! -e "$db/$template.hhm" && ! -e "$db/$template.db")) {
		# Some other database?
		$tmp{'tag'}=$template;
		$template=~/(\w*)\d*/; # remove numbers
		@dirs=glob "$hhpred_dir/".lc($1)."*";
		$db=pop(@dirs);
		if (!$db || ! -e "$db/$template.hhm") {
		    # Look through all genome-databases if they contain the HMM
		    $template =~ s/\|/_/g;
		    $template =~ s/\./_/g;
		    my $template1="";
		    if ($template =~ /(gi_\d+)/) {
			$template1 = $1;
		    }
		    my @alldirs=((glob "$hhpred_genomes_dir/bacteria/*"),(glob "$hhpred_genomes_dir/archaea/*"),(glob "$hhpred_genomes_dir/eucarya/*"));
		    foreach my $db1 (@alldirs) {
			if (-e "$db1/$template".".hhm") {$db=$db1; last; print "DB: $db\n";}
			if (-e "$db1/$template1".".hhm") {$template=$template1; $db=$db1; last; print "DB: $db\n";}
		    }
		
		    # HMM could still not be found? Look through all databases if they contain the HMM
		    if (!$db || ! -e "$db/$template.hhm") {
			$template =~ /^(\w+)\d*/;
			$tmp{'tag'}=$template;
			my @alldirs=((glob "$hhpred_dir/*"),(glob "$hhcluster_dir/*"),$basedir);
			foreach my $db1 (@alldirs) {
			    if (-e "$db1/$template".".hhm") {$db=$db1; last;}
			}

			# HMM could not be found?
			if (!$db || ! -e "$db/$template.hhm") {
			    print("Error: '$template.hhm' could not be found in any of the databases in $hhpred_dir or $hhcluster_dir\n");
			    next;
			}
		    }
		}
	    }

	} # end $program==1 (hhsearch)

	$tmp{'hhmfile'}="$db/$template.hhm";
	$tmp{'title'}="Prob=$prob\% E=$Eval  NaMeLiNe";
	$tmp{'template'}=$template;
	$tmp{'prob'}=$prob;
	$tmp{'qbeg'}=$qbeg;
	$tmp{'qend'}=$qend;
	$tmp{'tbeg'}=$tbeg;
	$tmp{'tend'}=$tend;
	$tmp{'tms'}=$tms;
	push(@hits, \%tmp);
    }

    elsif ($line=~/^\s*$/) {last;}   # empty line found? => read alignments

}
########################################################################################################



# Search for beginning  of first alignment
while($line=<RESFILE>){
    if ($line=~/^No (\d+)/) {last;}
}
$m=0;

########################################################################################################
# Read alignments one by one

while (defined $line) {

    # Start of next alignment found
    if ($m>0) {
	${$hits[$m-1]}{'qsspred'}=$qsspred;
	${$hits[$m-1]}{'tsspred'}=$tsspred;
	${$hits[$m-1]}{'qssconf'}=$qssconf;
	${$hits[$m-1]}{'tssconf'}=$tssconf;
	${$hits[$m-1]}{'qssdssp'}=$qssdssp;
	${$hits[$m-1]}{'tssdssp'}=$tssdssp;
	${$hits[$m-1]}{'tcons'}=$tcons;
	${$hits[$m-1]}{'qcons'}=$qcons;
	${$hits[$m-1]}{'tcons'}=$tcons;
	$qseq=~tr/\.//d;
	$tseq=~tr/\.//d;
	${$hits[$m-1]}{'qseq'}=$qseq;
	${$hits[$m-1]}{'tseq'}=$tseq;
	${$hits[$m-1]}{'match'}=$match;
	if ($v>=3) {
	    print("\n$m\n");
	    print("ssdssp  $qssdssp\n");
	    print("sspred  $qsspred\n");
	    print("qseq    $qseq\n");
	    print("qcons   $qcons\n");
	    print("        $match\n");
	    print("tcons   $tcons\n");
	    print("tseq    $tseq\n");
	    print("ssdssp  $tssdssp\n");
	    print("sspred  $tsspred\n");
	}
    }
    
    if ($line!~/^No (\d+)/) {last;}

    $m=$1;

    $qsspred='';
    $tsspred='';
    $qssconf='';
    $tssconf='';
    $qssdssp='';
    $tssdssp='';
    $qcons='';
    $tcons='';
    $qseq='';
    $tseq='';
    $match='';
    
    $line=<RESFILE>;
    $line=~/^>(.{1,80})/;
    my $nameline=$1;
    $nameline=~tr/<>//d;
    ${$hits[$m-1]}{'title'}=~s/NaMeLiNe/$nameline/;

    # Read all blocks of this alignment
    while (1) {
	
	# Move up to first line of Query block
	while ($line && $line!~/^Q\s+/) {$line=<RESFILE>;}
	
	# Read lines in query block
	my $seed=0;
	while ($line) {
	    if    ($line=~/^Q\s+ss_pred\s+(\S+)/) {$qsspred.=$1;}
	    elsif ($line=~/^Q\s+ss_conf\s+(\S+)/) {$qssconf.=$1;}
	    elsif ($line=~/^Q\s+ss_dssp\s+(\S+)/) {$qssdssp.=$1;}
	    elsif ($line=~/^Q\s+Consens..\s+\d+\s+(\S+)/) {$qcons.=$1;}
	    elsif (!$seed && $line=~/^Q\s+\S+\s+\d+\s+(\S+)/) {$qseq.=$1; $seed=1;}
	    elsif ($line!~/^Q\s+/) {last;}
	    $line=<RESFILE>;
	}

	# Read match quality symbol line
	$line=~/^\s{22}(.+)\s$/;
	$match.=$1;
	$line=<RESFILE>;
	
	# Read lines in template block
	$seed=0;
	while ($line) {
	    if    ($line=~/^T\s+ss_pred\s+(\S+)/) {$tsspred.=$1;}
	    elsif ($line=~/^T\s+ss_conf\s+(\S+)/) {$tssconf.=$1;}
	    elsif ($line=~/^T\s+ss_dssp\s+(\S+)/) {$tssdssp.=$1;}
	    elsif ($line=~/^T\s+Consens..\s+\d+\s+(\S+)/) {$tcons.=$1;}
	    elsif (!$seed && $line=~/^T\s+\S+\s+\d+\s+(\S+)/) {$tseq.=$1; $seed=1;}
	    elsif ($line!~/^T\s+/) {last;}
	    $line=<RESFILE>;
	}

	# Move up to next non-empty line
	while ($line && ($line=~/^\s*$/ || $line=~/^Confidence/)) {$line=<RESFILE>;}

	if (!$line || $line=~/^No (\d+)/ || $line=~/^Done/) {last;}
    }
}
########################################################################################################

close RESFILE;

# read query hmm
$q_hmm_ref=&extractModel($basedir.'/'.$id.'.hhm');

# read template hmms
foreach my $i (@hits) {
    my %hit = %$i;
    push(@t_hmms, &extractModel($hit{'hhmfile'}));
}

# calc img width
$imgwidth=$lborder+$alnwidth*$barwidth+($alnwidth-1)*$hbarsep+$rborder; 

# for each hit: draw alignment graph
for(my $i=0; $i<scalar(@hits); $i++) {
    my $imgfile=$idpng.'_'.($i+1).'.png';
    my $imgfile_tt=$idpng.'_'.($i+1).'_tt.png';
    my $mapfile=$idpng.'_'.($i+1).'.map';
    my $ref=$hits[$i];
    my $len=length(${$ref}{'match'});
    my $nogaplen=0;
    my $mapname=$idpng."_".($i+1);
    my $imgname="hitimage".($i+1);
    
    # determine ungapped sequence length
    my @qcons = unpack("C*",${$ref}{'qcons'});
    my @tcons = unpack("C*",${$ref}{'tcons'});
    for (my $j=0; $j<$len; $j++) {
	if (( $qcons[$j]!=45 && $qcons[$j]!=46 ) || ($tcons[$j]!=45 && $tcons[$j]!=46) ) {
	    $nogaplen++;
	}
    }
    my $lines=ceil($nogaplen/$alnwidth);

    # calc img height
    $imgheight=$tborder+$bborder+$lines*(2*$barheight+$qtsep)+($lines-1)*$linesep;
    
    # create new image
    $im = new GD::Image($imgwidth,$imgheight);
    $im_tt = new GD::Image($imgwidth,$imgheight);
    
    # allocate colors
    &allocateColors();
    
    # draw profile logos
    my $mapfh;
    open ($mapfh,">$tmpdir/$mapfile") or die "unable to open mapfile $mapfile for writing\n";
    print $mapfh "<map name=\"".$mapname."\">\n";
    &drawProfile($ref,$t_hmms[$i],$mapfh,$imgfile,$imgfile_tt,$i+1);
    print $mapfh "</map>\n";
    print $mapfh "<img name=\"$imgname\" src=\"$imgdir/$imgfile\" border=\"0\" alt=\"hhpred profile ".($i+1)."\" usemap=\"#$mapname\" />\n";
    close $mapfh;
    
    # write main image to file
    open (OUTFILE,">$tmpdir/$imgfile") or die "unable to open imgfile $tmpdir/$imgfile for writing\n";
    print OUTFILE $im->png;
    close OUTFILE;
    # write tooltip image to file
    open (OUTFILE,">$tmpdir/$imgfile_tt") or die "unable to open imgfile $tmpdir/$imgfile_tt for writing\n";
    print OUTFILE $im_tt->png;
    close OUTFILE;

    if ($i) {
	# add files to existing tar archive
	&System("cd $tmpdir; tar -rf $tarfile $imgfile $imgfile_tt $mapfile");
    } else {
	# create new tar archive
    	&System("cd $tmpdir; tar -cf $tarfile $imgfile $imgfile_tt $mapfile;");
    }
    unlink ("$tmpdir/$mapfile", "$tmpdir/$imgfile", "$tmpdir/$imgfile_tt" );
}
# zip tarfile
system("gzip -c $tmpdir/$tarfile > $basedir/$zipfile");
system("rm -rf $tmpdir");

exit;

sub drawProfile {
    my $hit_ref=$_[0]; 
    my $hmm_ref=$_[1];
    my $mapfh=$_[2];
    my $imgfile=$_[3];
    my $imgfile_tt=$_[4];
    my $hitnum=$_[5];
    my $qpos=${$hit_ref}{'qbeg'}-1;
    my $tpos=${$hit_ref}{'tbeg'}-1;
    my $len=length(${$hit_ref}{'match'});

    if (scalar($hmm_ref)==0) { return; } # do not draw profile if HMM is empty! (e.g. hhm-file of template could not be found)

    if (${$hit_ref}{'qssdssp'}=~/^-*$/) {${$hit_ref}{'qssdssp'} = " " x $len;}
    if (${$hit_ref}{'qsspred'} eq "")   {${$hit_ref}{'qsspred'} = " " x $len;}
    if (${$hit_ref}{'qssconf'} eq "")   {${$hit_ref}{'qssconf'} = " " x $len;}
    if (${$hit_ref}{'qcons'} eq "")     {${$hit_ref}{'qcons'} = " " x $len;}
    if (${$hit_ref}{'tssdssp'}=~/^-*$/) {${$hit_ref}{'tssdssp'} = " " x $len;}
    if (${$hit_ref}{'tsspred'} eq "")   {${$hit_ref}{'tsspred'} = " " x $len;}
    if (${$hit_ref}{'tssconf'} eq "")   {${$hit_ref}{'tssconf'} = " " x $len;}
    if (${$hit_ref}{'tcons'} eq "")     {${$hit_ref}{'tcons'} = " " x $len;}
    my $i;                              # index of printed characters 
    my $j;                              # index of processed characters 
    my $gapcols=0;

    for ($i=0,$j=0; $j<$len; $j++,$i=$j-$gapcols) {
	my $line_i=$i%$alnwidth;
	my $xoff=$lborder+$line_i*($hbarsep+$barwidth);
	my $yoff=$tborder+floor($i/$alnwidth)*(2*$barheight+$qtsep)+(floor($i/$alnwidth))*$linesep+$barheight;
	my $x1=$xoff;
	my $y1;
	my $x2=$xoff+$barwidth;
	my $y2;
	my $prob;
	my $aa;
	my $qlab;
	my $tlab;
	my $match_chr=substr(${$hit_ref}{'match'},$j,1);
	my $qcons_chr=substr(${$hit_ref}{'qcons'},$j,1);
	my $tcons_chr=substr(${$hit_ref}{'tcons'},$j,1);
	my $qsspred_chr=substr(${$hit_ref}{'qsspred'},$j,1);
	my $tsspred_chr=substr(${$hit_ref}{'tsspred'},$j,1);
	my $qssconf_chr=substr(${$hit_ref}{'qssconf'},$j,1);
	my $tssconf_chr=substr(${$hit_ref}{'tssconf'},$j,1);
	my $qssdssp_chr=substr(${$hit_ref}{'qssdssp'},$j,1);
	my $tssdssp_chr=substr(${$hit_ref}{'tssdssp'},$j,1);
	my $qgap=0;
	my $tgap=0;

	if (substr(${$hit_ref}{'qcons'},$j,1) eq "-" || 
	    substr(${$hit_ref}{'qcons'},$j,1) eq ".") {
	    $qgap=1;
	}
	if (substr(${$hit_ref}{'tcons'},$j,1) eq "-" || 
	    substr(${$hit_ref}{'tcons'},$j,1) eq ".") {
	    $tgap=1;
	}
	
	# skip column entirely?
	if ($qgap && $tgap) {
	    $gapcols++;
	    next;
	}

	# print maplink
	if ($line_i==0) {
	    &printMapLink($hitnum,$mapfh,"$imgdir/$imgfile_tt","$imgdir/$imgfile",$xoff,$yoff-$qtsep,$xoff+$alnwidth*($barwidth+$hbarsep)-$hbarsep,$yoff+2*$qtsep);
	}

	# print sequence character to tooltip image
	my $qchar=substr(${$hit_ref}{'qseq'},$i,1);
	my $tchar=substr(${$hit_ref}{'tseq'},$i,1);
#	if ($qchar ne '-') {
#	    $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
#	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yoff-1.5*$qtsep-$stdfontheight*0.8,$qchar,$aa_colors{$qchar});
#	}
#	if ($tchar ne '-') {
#	    $im_tt->filledRectangle($xoff,$yoff+$qtsep,$xoff+$barwidth,$yoff+2*$qtsep,$black);
#	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yoff+2.5*$qtsep-$stdfontheight*0.8,$tchar,$aa_colors{$tchar});
#	}

	my $yqoff=$yoff;
	my $ytoff=$yoff;
	
	# draw conservation character
	if ($match_chr ne " ") {
	    # in tooltip image
	    $im->filledRectangle($xoff,$yoff,$xoff+$barwidth,$yoff+$qtsep,$black);
	    $im->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yoff+0.5*$qtsep-$stdfontheight*0.8,$match_chr,$white);
	    # in main image
#	    $im_tt->filledRectangle($xoff,$yoff,$xoff+$barwidth,$yoff+$qtsep,$black);
	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yoff+0.5*$qtsep-$stdfontheight*0.8,$match_chr,$black);
	} else {
	    # ...in main image
	    $im->filledRectangle($xoff+1,$yoff,$xoff+$barwidth,$yoff+$qtsep,$white);
	    # ...in tooltip image
	    $im_tt->filledRectangle($xoff+1,$yoff,$xoff+$barwidth,$yoff+$qtsep,$white);
	}

	# Draw qcons?
	if (0 && $qcons_chr ne " ") {
	    $yqoff-=$qtsep;
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yqoff+0.5*$qtsep-$stdfontheight*0.8,$qcons_chr,$black);
	}
	# Draw qchar?
	if ($qchar ne " ") {
	    $yqoff-=$qtsep;
#   	    $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yqoff+0.5*$qtsep-$stdfontheight*0.8,$qchar,$aa_colors{uc($qchar)});
	}
	# Draw qssconf?
	if ($qssconf_chr ne " ") {
	    $yqoff-=$qtsep;
#	    if ($qssconf_chr ne "-") {
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
		$im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yqoff+0.5*$qtsep-$stdfontheight*0.8,$qssconf_chr,$black);
#	    }
	}
	# Draw qsspred?
	if ($qsspred_chr ne " ") {
	    $yqoff-=$qtsep;
#	    if ($qsspred_chr ne "-") {
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
		$im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yqoff+0.5*$qtsep-$stdfontheight*0.8,$qsspred_chr,$ss_colors{$qsspred_chr});
#	    }
	}
	# Draw qssdssp?
	if ($qssdssp_chr ne " ") {
	    $yqoff-=$qtsep;
#	    if ($qssdssp_chr ne "-") {
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
		$im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$yqoff+0.5*$qtsep-$stdfontheight*0.8,$qssdssp_chr,$ss_colors{$qssdssp_chr});
#	    }
	}
	
	# Draw tcons?
	if (0 && $tcons_chr ne " ") {
	    $ytoff+=$qtsep;
#           $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$ytoff+0.5*$qtsep-$stdfontheight*0.8,$tcons_chr,$black);
	}
	# Draw tchar?
	if ($tchar ne " ") {
	    $ytoff+=$qtsep;
#   	    $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
	    $im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$ytoff+0.5*$qtsep-$stdfontheight*0.8,$tchar,$aa_colors{uc($tchar)});

	}
	# Draw tssdssp?
	if ($tssdssp_chr ne " ") {
	    $ytoff+=$qtsep;
#	    if ($tssdssp_chr ne "-") {
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
		$im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$ytoff+0.5*$qtsep-$stdfontheight*0.8,$tssdssp_chr,$ss_colors{$tssdssp_chr});
#	    }
	}
	# Draw tsspred?
	if ($tsspred_chr ne " ") {
	    $ytoff+=$qtsep;
#	    if ($tsspred_chr ne "-") {
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
		$im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$ytoff+0.5*$qtsep-$stdfontheight*0.8,$tsspred_chr,$ss_colors{$tsspred_chr});
#	    }
	}
	# Draw tssconf?
	if ($tssconf_chr ne " ") {
	    $ytoff+=$qtsep;
#	    if ($tssconf_chr ne "-") {
#   	        $im_tt->filledRectangle($xoff,$yoff-$qtsep,$xoff+$barwidth,$yoff,$black);
		$im_tt->string($stdfont,1+$x1+($x2-$x1)/2-$stdfontwidth/2,$ytoff+0.5*$qtsep-$stdfontheight*0.8,$tssconf_chr,$black);
#	    }
	}
	# draw profile bars for query sequence
	if (!$qgap) {
	    my $currheight=0.0;
	    $qpos++;
	    $comphash_ref=${$q_hmm_ref}[$qpos-1];
	    for $aa (sort byprob keys %$comphash_ref) {
		$prob=${$q_hmm_ref}[$qpos-1]{$aa};
		if ($prob>=$probcutoff) {
		    $y1=$yoff-$currheight-$prob*$barheight;
		    $y2=$yoff-$currheight;
		    $currheight+=$prob*$barheight;
		    &drawBar($x1,$y1,$x2,$y2,$aa);
		}
	    }
	    print "\n";
	}
	# draw profile bars for template sequence
	if (!$tgap) {
	    my $currheight=0.0;
	    $tpos++;
	    $comphash_ref=${$hmm_ref}[$tpos-1];
	    for $aa (sort byprob keys %$comphash_ref) {
		$prob=${$hmm_ref}[$tpos-1]{$aa};
		if ($prob>=$probcutoff) {
		    $y1=$yoff+$qtsep+$currheight;
		    $y2=$yoff+$qtsep+$currheight+$prob*$barheight;
		    $currheight+=$prob*$barheight;
		    &drawBar($x1,$y1,$x2,$y2,$aa);
		}
	    }
	}
	# draw state numbers?
	if ($line_i==0) {
	    if ($i>0 &&
		(substr(${$hit_ref}{'qcons'},$j-1,1) eq "-" || 
		substr(${$hit_ref}{'qcons'},$j-1,1) eq ".")) {
		$qlab='Q '.($qpos+1);
	    } else {
		$qlab='Q '.$qpos;
	    }
	    if ($i>0 &&
		(substr(${$hit_ref}{'tcons'},$j-1,1) eq "-" || 
		substr(${$hit_ref}{'tcons'},$j-1,1) eq ".")) {
		$tlab='T '.($tpos+1);
	    } else {
		$tlab='T '.$tpos;
	    }
	    $im->string($stdfont,$xoff-8*$stdfontwidth,$yoff-$stdfontheight*1.5,$qlab,$black);
	    $im->string($stdfont,$xoff-8*$stdfontwidth,$yoff+$qtsep+$stdfontheight*0.2,$tlab,$black);
	    $im_tt->string($stdfont,$xoff-8*$stdfontwidth,$yoff-$stdfontheight*1.5,$qlab,$black);
	    $im_tt->string($stdfont,$xoff-8*$stdfontwidth,$yoff+$qtsep+$stdfontheight*0.2,$tlab,$black);
	} elsif ($line_i+1==$alnwidth || $j+1==$len) {
	    $im->string($stdfont,$xoff+$barwidth+2*$stdfontwidth,$yoff-$stdfontheight*1.5,$qpos,$black);
	    $im->string($stdfont,$xoff+$barwidth+2*$stdfontwidth,$yoff+$qtsep+$stdfontheight*0.2,$tpos,$black);
	    $im_tt->string($stdfont,$xoff+$barwidth+2*$stdfontwidth,$yoff-$stdfontheight*1.5,$qpos,$black);
	    $im_tt->string($stdfont,$xoff+$barwidth+2*$stdfontwidth,$yoff+$qtsep+$stdfontheight*0.2,$tpos,$black);
	}

	
	# process tooltip
	#my $ttbeg=&max($i-$ttreach,0);
	#my $ttend=&min($i+$ttreach,$len);
	#my $lrange=$i-$ttreach<0 ? $i : $ttreach;
	#my $rrange=$i+$ttreach>$len ? $len-$i : $ttreach;
	#my $tooltiplen=$ttend-$ttbeg;
	#my $qtooltip='';
	#my $ttooltip='';
	#my $qres=' ';
	#my $tres=' ';
	
	#$qtooltip=$qtooltip.substr(${$hit_ref}{'qseq'},$ttbeg,$lrange);	
	#$qtooltip=$qtooltip."<span style=\\'background-color: #FFCC33;\\'>".substr(${$hit_ref}{'qseq'},$ttbeg+$lrange,1)."</span>";	
	#$qtooltip=$qtooltip.substr(${$hit_ref}{'qseq'},$ttbeg+$lrange+1,$rrange);	
	#$ttooltip=$ttooltip.substr(${$hit_ref}{'tseq'},$ttbeg,$lrange);	
	#$ttooltip=$ttooltip."<span style=\\'background-color: #FFCC33;\\'>".substr(${$hit_ref}{'tseq'},$ttbeg+$lrange,1)."</span>";	
	#$ttooltip=$ttooltip.substr(${$hit_ref}{'tseq'},$ttbeg+$lrange+1,$rrange);

	#if (!$qgap) {
	#    $qres=('-' x ($lrange-ceil(length($qpos)/2)+4)).$qpos;	    
	#} 
	#if (!$tgap) {
	#    $tres=('-' x ($lrange-ceil(length($tpos)/2)+4)).$tpos;
	#}
	#my $text=$qres."<br/>Q: ".$qtooltip."<br/>T: ".$ttooltip."<br/>".$tres;
	#&printTooltip($mapfh,$text,$xoff,$yoff-$barheight,$xoff+$barwidth,$yoff+$barheight+$qtsep);
    }
}

sub drawBar {
    my $x1 = $_[0];
    my $y1 = $_[1];
    my $x2 = $_[2];
    my $y2 = $_[3];
    my $aa = $_[4];

    $im->filledRectangle($x1,$y1,$x2,$y2,$aa_colors{$aa});
    $im->rectangle($x1,$y1,$x2,$y2,$black);
    $im->string($barfont,$x1+($x2-$x1)/2-$barfontwidth/2,$y1+($y2-$y1)/2-$barfontheight*0.8,$aa,$black);
}

sub printMapLink() {
    my $hitnum=$_[0];
    my $mapfh=$_[1];
    my $imgOn=$_[2];
    my $imgOff=$_[3];
    my $x1 = $_[4];
    my $y1 = $_[5];
    my $x2 = $_[6];
    my $y2 = $_[7];
    
    my $imgName="hitimage".$hitnum;
    my $onMouseOver="imgOn('$imgName','$imgOn');";
    my $onMouseOut="imgOff('$imgName','$imgOff');";

    if ($x1 =~ /(\d+)\.\d+/) {$x1 = $1}; 
    if ($y1 =~ /(\d+)\.\d+/) {$y1 = $1};
    if ($x2 =~ /(\d+)\.\d+/) {$x2 = $1};
    if ($y2 =~ /(\d+)\.\d+/) {$y2 = $1};

    print $mapfh "<area shape=\"rect\" coords=\"$x1,$y1,$x2,$y2\" href=\"#\" onmouseover=\"$onMouseOver\" onmouseout=\"$onMouseOut\" >\n";
}

sub printTooltip() {
    my $mapfh=$_[0];
    my $text=$_[1];
    my $x1 = $_[2];
    my $y1 = $_[3];
    my $x2 = $_[4];
    my $y2 = $_[5];

    if ($x1 =~ /(\d+)\.\d+/) {$x1 = $1}; 
    if ($y1 =~ /(\d+)\.\d+/) {$y1 = $1};
    if ($x2 =~ /(\d+)\.\d+/) {$x2 = $1};
    if ($y2 =~ /(\d+)\.\d+/) {$y2 = $1};
#    print $mapfh "<area shape=\"rect\" coords=\"$x1,$y1,$x2,$y2\" href=\"javascript:void(0);\" onmouseover=\"return overlib('$text', CSSSTYLE,TEXTSIZE,12,FGCOLOR,'#ffffff');\" onmouseout=\"return nd();\" >\n";
    print $mapfh "<area shape=\"rect\" coords=\"$x1,$y1,$x2,$y2\" href=\"#\" onmouseover=\"return ptt('$text');\" onmouseout=\"return nd();\" >\n";
}

sub extractModel {
    my $hmmfile = $_[0];
    my @hmm;
    my $line;
    
    if (-r $hmmfile) {
	open (HMMFILE,"<$hmmfile") or die "unable to open hhmfile $hmmfile for reading\n";
	while($line=<HMMFILE>) {
	    if ($line=~/^\w\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\d+\s*$/) {
		my %tmp;
		$tmp{'A'}=&calcProb($2);
		$tmp{'C'}=&calcProb($3);
		$tmp{'D'}=&calcProb($4);
		$tmp{'E'}=&calcProb($5);
		$tmp{'F'}=&calcProb($6);
		$tmp{'G'}=&calcProb($7);
		$tmp{'H'}=&calcProb($8);
		$tmp{'I'}=&calcProb($9);
		$tmp{'K'}=&calcProb($10);
		$tmp{'L'}=&calcProb($11);
		$tmp{'M'}=&calcProb($12);
		$tmp{'N'}=&calcProb($13);
		$tmp{'P'}=&calcProb($14);
		$tmp{'Q'}=&calcProb($15);
		$tmp{'R'}=&calcProb($16);
		$tmp{'S'}=&calcProb($17);
		$tmp{'T'}=&calcProb($18);
		$tmp{'V'}=&calcProb($19);
		$tmp{'W'}=&calcProb($20);
		$tmp{'Y'}=&calcProb($21);
		$hmm[$1-1]=\%tmp;
	    }
	}
	close HMMFILE;
    }
    
    return \@hmm;
}

sub calcProb {
    my $val = $_[0];
    
    if ($val eq '*') {
	return 0.0;
    } elsif($val == 0) {
	return 1.0;
    } else {
	return (2**((-1)*$val/1000))
    }
}

sub allocateColors
{
    # allocate colors for main image
    $white=$im->colorAllocate(255,255,255);
    $black=$im->colorAllocate(0,0,0);  
    $grey=$im->colorAllocate(144,144,144);
    $lightgrey=$im->colorAllocate(240,240,240);
    $darkgreen=$im->colorAllocate(0,150,0);
    $lightgreen=$im->colorAllocate(32,255,32);
    $lightred=$im->colorAllocate(255,0,0);
    $darkred=$im->colorAllocate(192,0,0);
    $yellow=$im->colorAllocate(255,255,0);
    $yellowgreen=$im->colorAllocate(240,255,0);
    $orange=$im->colorAllocate(255,128,0);
    $lightpink=$im->colorAllocate(255,180,220);
    $darkviolett=$im->colorAllocate(192,128,255);
    $violett=$im->colorAllocate(224,128,255);
    $lightviolett=$im->colorAllocate(255,192,255);
    $blue=$im->colorAllocate(92,128,255);
    $im->transparent($white);

    # allocate colors for tooltip image
    $white=$im_tt->colorAllocate(255,255,255);
    $black=$im_tt->colorAllocate(0,0,0);  
    $grey=$im_tt->colorAllocate(144,144,144);
    $lightgrey=$im_tt->colorAllocate(40,40,40);
    $darkgreen=$im_tt->colorAllocate(0,196,0);
    $lightgreen=$im_tt->colorAllocate(0,255,0);
    $lightred=$im_tt->colorAllocate(255,0,0);
    $darkred=$im_tt->colorAllocate(192,0,0);
    $yellow=$im_tt->colorAllocate(255,128,0);
    $yellowgreen=$im_tt->colorAllocate(240,255,0);
    $orange=$im_tt->colorAllocate(255,64,0);
    $lightpink=$im_tt->colorAllocate(255,180,220);
    $darkviolett=$im_tt->colorAllocate(192,128,255);
    $violett=$im_tt->colorAllocate(224,128,255);
    $lightviolett=$im_tt->colorAllocate(255,192,255);
    $blue=$im_tt->colorAllocate(0,0,255);
    $im_tt->transparent($white);

    %aa_colors = ( 'A' => $lightgrey,
		   'C' => $yellow,
		   'D' => $blue, 
		   'E' => $blue,
		   'F' => $darkgreen,
		   'G' => $lightgrey,
		   'H' => $orange,
		   'I' => $lightgreen,
		   'K' => $lightred,
		   'L' => $lightgreen,
		   'M' => $lightgreen,
		   'N' => $violett,
		   'P' => $grey,
		   'Q' => $violett,
		   'R' => $lightred,
		   'S' => $lightgrey,
		   'T' => $lightgrey,
		   'V' => $lightgreen,
		   'W' => $darkgreen,
		   'Y' => $darkgreen,
		   'B' => $blue,
		   'Z' => $blue,
		   'X' => $lightgrey,
		   '-' => $black,
		   '~' => $grey,
		   '.' => $grey,
		   ' ' => $black
		   );

    %ss_colors = ( 'H' => $lightred,
		   'h' => $lightred,
		   'E' => $blue,
		   'e' => $blue,
		   'C' => $black, 
		   'c' => $black, 
		   'G' => $black,
		   'I' => $black,
		   'S' => $black,
		   'T' => $black,
		   'B' => $black,
		   ' ' => $black,
		   '~' => $black,
		   '-' => $black
		   );
}

sub min {
    if ($_[0]<=$_[1]) {
	return $_[0];
    } else {
	return $_[1];
    }
}

sub max {
    if ($_[0]>=$_[1]) {
	return $_[0];
    } else {
	return $_[1];
    }
}

sub byprob {
    my $order = $aa_order{$a} <=> $aa_order{$b};
    if ($order==0) {
	return ${$comphash_ref}{$b} <=> ${$comphash_ref}{$a};
    } 
    return $order;
}

sub System()
{
    if ($v>=2) {printf("\$ %s\n",$_[0]);} 
    return system($_[0])/256;
}
