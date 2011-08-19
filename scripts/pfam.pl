#! /usr/bin/env perl
# Read the Pfam database file, reformat it into seperate a3m files,
# and generate hhm files.
# Usage:   pfam.pl Pfam-file
# Example: pfam.pl Pfam-A.full
#
# All sequences with pdb ids will be written at the top of the alignment
# The first sequence (=query) will be the (pdb-)sequence with the least
# number of deletions plus insertions
my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;                         # config file with path variables for nr, blast, psipred, pdb, dssp etc.

my $v=2;
my $format="a3m";
my $update=0;
my $Neff_min=0;
my $usage="
Read the Pfam database file and reformat it into seperate a3m or a2m files
Usage:   pfam.pl [options] file
Options: 
 -u            update: do not overwrite existing alignment files
 -v int        verbose mode
 -a2m          Extract alignments in a2m format instead of a3m
 -Neff         alignments with Neff < this threshold or Nseq < 2* this threshold will be skipped.
Example: pfam.pl -Neff 2.0 Pfam-A.full
\n";

# Set directory paths
my $dbname=$nr90;                 # The name of the BLAST data bank
my @iprodirs=glob("$newdbs/interpro_*");
my $iprodir= $iprodirs[$#iprodirs];

my $line;
my $infile;
my $ID="";           # Pfam id
my $AC="";           # Pfam accesion code
my $DE="";           # Pfam-A description
my $DR="";           # Pfam-A description
my $DC="";           # Pfam-B remark: This family was generated automatically from an alignment ...
my $basename;        # basename alignment file, psiblast file, and sequence file (PF00244.7 -> PF00244)
my @PDBs=();         # PDB codes in alignment
my $dssp=0;          # flag: is '#=GC SS_cons' included in alignment?
my $Ndom=0;          # number of domains read
my $name;            # name of sequence in alignment 
my $nameline;        # name of sequence + pdb codes
my $seq;             # residues of sequence in alignment
my $qname;           # first sequence of alignment, defined as first sequence with minimum number of indels
my $qnameline;       # nameline of first sequence of alignment
my $qseq;            # query sequence of alignment where '-' were substituted with 'X'
my $qmatch;          # number of deletions plus insertions in candidate query sequence
my $qpdb=0;          # 0:query is not pdb sequence  1:query is pdb sequence
my @names;           # names of sequences read in from alignment
my @namelines;       # namelines of sequences read in from alignment
my @seqs;            # sequences read in from alignment
my @inserts;         # insertions and deletions in first sequence of alignment 
my $i;               # counts columns
my @res;             # $res[$i] = amino acid in i'th column
my @match;           # $match[$i]=1 if $i is match column, else $match[$i]=0
my $dummydb="dummy"; # blast database consisting of just one sequence
my $titleline;       # first line in alignment: '#domain-AC domain-ID domain-DC; PDB-codes' 
my $nseq;            # number of sequences read from alignment
my %pdbids=();
my %id2ipr=();       # E.g. $id2ipr{"PF03445"}="IPR005105"
my %ipr2abstract=(); # E.g. $ipr2abstract{"IPR005105"}="This domain adopts a Rossman NAD binding fold. The C-termi..."
my %ipr2GO;          # E.g. $ipr2GO{"IPR001144"} = "GO: 0009405 catalytic activity, 0009405 pathogenesis, 0005576 extracellular region"

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

if ($options=~s/ -v\s+(\d) / /) {$v=$1;}
if ($options=~s/ -v / /) {$v=2;}
if ($options=~s/ -u / /) {$update=1;}
if ($options=~s/ -a2m / /g) {$format="a2m";}
if ($options=~s/ -Neff\s+(\S+) / /) {$Neff_min=$1;}
if ($options=~s/ -i\s+(\S+) / /) {$infile=$1;}

# Warn if unknown options found
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; print("WARNING: unknown options '$options'\n");}

if (!-e "$iprodir/interpro.xml") {
    print($usage); 
    die("Error in $0: $iprodir/interpro.xml could not be found.\nRun hhpred_update_ipro.pl to download file and check directory path in pfam.pl\n");
}
if (@ARGV<1) {die("$usage");}

if ($infile=~/Pfam.?B/i) {$format="a2m";}  # extract PfamB in A2M format to be able to redefine match states



&ReadInterProFiles();

open (INFILE,"<$infile") || die("Error: can't open $infile: $!\n");

while ($line=<INFILE>) {

#    print ("line=$line");
    if ($line=~/^\/\//) { 

	# End of alignment read

	if (-e "$basename.$format" && $update) {next;}
	if (@seqs<2*$Neff_min) {printf("Number of sequences=%i Skipping\n",scalar(@seqs)); next;} # not enogh sequence in alignment

	# Write into psipred query file
	$seq=$qseq;
	$seq=~tr/-a-z./X/d; # remove all inserts, transform gaps into 'X'
	open(SEQFILE,">$basename.seq") || die("Error: can't open $basename.seq: $!\n");
	printf(SEQFILE ">$qname ('-'=>'X')\n$seq\n");
	close(SEQFILE);
	    
	# Write into psipred jumpstart file
	open(PSIFILE,">$basename.psi") || die("Error: can't open $basename.psi: $!\n");
	printf(PSIFILE "%-27.27s %s\n",$qname,$seq);
	for (my $i=0; $i<@seqs && $i<30000; $i++) {
	    $seq=$seqs[$i];
	    $seq=~tr/a-z.//d;
	    if ($names[$i] ne $qname) {
		printf(PSIFILE "%-27.27s %s\n",$names[$i],$seq);
	    }
	}
	close(PSIFILE);
	
	# Do SS prediction with Psipred
	printf("Predicting secondary structure for $basename.seq \n");
	&RunPsipred("$basename.seq");
	
	# Read Psipred file
	open (PSIPREDFILE, "<$basename.horiz") || die("Error: could not open $basename.horiz: $!\n");  
	my $aa_pred="";
	my $ss_pred="";
	my $ss_conf="";
	while ($line=<PSIPREDFILE>) {
	    if    ($line=~/^Conf:\s+(\S+)/) {$ss_conf.=$1;}
	    elsif ($line=~/^Pred:\s+(\S+)/) {$ss_pred.=$1;}
	    elsif ($line=~/^  AA:\s+(\S+)/) {$aa_pred.=$1;}
	}
	close(PSIPREDFILE);

	$ss_conf =~ s/[^0-9]/0/g;
	
	if ($format ne "a3m") {
	    
	    # Change capitalization of psipred sequences to match that of the first sequence in alignment ($qseq)
	    my $i;   #position in $ss_pred sequence
	    my @ss_pred=split(//,$ss_pred);
	    my @ss_conf=split(//,$ss_conf);
	    my @aa_pred=split(//,$aa_pred);
	    push(@ss_pred,"");
	    push(@ss_conf,"");
	    push(@aa_pred,"");
	    $ss_pred="";
	    $ss_conf="";
	    $aa_pred="";
	    for ($i=0; $i<@inserts; $i++) {
		$inserts[$i]=~tr/a-z/........................../d; # substitute a-z by '.' 
		$ss_pred .= $inserts[$i].$ss_pred[$i];
		$ss_conf .= $inserts[$i].$ss_conf[$i];
		$aa_pred .= $inserts[$i].$aa_pred[$i];
	    }		
	    for (; $i<@ss_pred; $i++) {
		$ss_pred .= $ss_pred[$i];
		$ss_conf .= $ss_conf[$i];
		$aa_pred .= $aa_pred[$i];
	    }		
	    
	}		
	    

	# Write alignment file
	open(ALIFILE,">$basename.$format") || die("Error: can't open $basename.$format: $!\n");
	print(ALIFILE $titleline);
	if ($dssp ne "") {print(ALIFILE ">ss_dssp\n$dssp\n");}
	print(ALIFILE ">ss_pred\n$ss_pred\n");
	print(ALIFILE ">ss_conf\n$ss_conf\n");
#	print(ALIFILE ">aa_pred\n$aa_pred\n");
	print(ALIFILE ">".$qnameline."\n".$qseq."\n");
	for (my $i=0; $i<@seqs; $i++) {
	    if ($names[$i] ne $qname) {
		print(ALIFILE ">".$namelines[$i]."\n".$seqs[$i]."\n");
	    }
	}
	close(ALIFILE);

	# Generate HHM file
	if ($infile!~/Pfam.?B/i) {
	    &System("$hh/hhmake -M a3m -id 90 -diff 100 -v 2 -cons -i $basename.$format");   
	} else {  
	    # In Pfam-B, long inserts do not normally appear in one HSP. 
	    # Instead, two HSPs from the same sequence will be contained in the alignment.
	    # Therefore it makes sense to take the sequence with the most residues as the representative of the alignment
	    &System("$hh/hhfilter -M first -id 90 -v 2 -i $basename.$format -o $basename.a3m");   
	    print("$hh/hhmake -M a3m -id 100 -v 2 -cons -i $basename.a3m\n");
	    open(HHMAKE,"$hh/hhmake -M a3m -diff 100 -v 2 -i $basename.a3m |") or die("Error: can not open $hh/hhmake -M a3m -diff 100 -v 2 -i $basename.a3m |: $!\n");
	    while ($line=<HHMAKE>) {if ($line=~/^Effective number.*= (\S+)/) {last;}}
	    close(HHMAKE);
	    $line=~/= (\S+)/;
	    print("$line");
	    if (defined $1 && $1<$Neff_min) {&System("rm $basename.*");}
	}

	# Add consensus sequence from hhm file as first sequence of alignment
	if (-e "$basename.a3m") {
	    my $cname="";
	    my $cons="";
	    my @lines;
	    open(HHRFILE,"<$basename.hhm") or die("Error: can not open $basename.hhm: $!\n"); 
	    while($line=<HHRFILE>) {if ($line=~/^NAME\s+(\S+)/) {$cname=$1; last;}}
	    while($line=<HHRFILE>) {if ($line=~/^SEQ/) {last;}}
	    while($line=<HHRFILE>) {if ($line=~/^>($cname.cons\S*)/) {$cname=$1; last;}    }
	    while($line=<HHRFILE>) {if ($line!~/^>/) {chomp($line); $cons.=$line;} else {last;} }
	    close (HHRFILE);
	    open(ALIFILE,"<$basename.$format") or die("Error: can not open $basename.$format: $!\n"); 
	    @lines=<ALIFILE>;
	    close (ALIFILE);
	    open(ALIFILE,">$basename.$format") or die("Error: can not open $basename.$format: $!\n"); 
	    foreach $line (@lines) {
		if (defined $cname && $line=~/^>/ && $line!~/^>aa_/ && $line!~/^>ss_/) {
		    printf(ALIFILE ">%s\n",$cname);
		    printf(ALIFILE "%s\n",$cons);
		    undef($cname);
		}
		printf(ALIFILE "%s",$line);
	    }
	    close (ALIFILE);
	}
	unlink("$basename.log");
	unlink("$basename.seq");
	unlink("$basename.horiz");
	unlink("$basename.psi");

	# Initialize for reading next alignment
	@PDBs=(); 
	%pdbids=();
	$DC=""; $DE="";

    } elsif ($line=~/^\#=GF ID\s+(\S.*)/) {
	$ID=$1;
    } elsif ($line=~/^\#=GF AC\s+(\S.*)/) {
	$basename=$1;
	$basename=~s/\..*//;
	$AC=$basename;
    } elsif ($line=~/^\#=GF DR\s+PRODOM;\s+(\S+);/) {
	$DR=$1;
	if ($infile=~/Pfam.?B/i) {$basename=$DR;}
    } elsif ($line=~/^\#=GF DE\s+(\S.*)/) {
	$DE.=" $1";
    } elsif ($line=~/^\#=GF DC\s+(\S.*)/) {
	$DC.=" $1";
    } elsif ($line=~/^\#=GF DR\s+PFAMA;\s+(\S+);/) {
	$DC.=" $1";
    } elsif ($line=~/^\#=GS (\S+)\s+DR\s+PDB\;\s*(\S+)\s+(\S*)\;(.*)/) {
	$name=$1;
	my $pdbid=$2;
	my $chain=$3;
	$line=$3;
	if ($chain ne "") {$pdbid="$pdbid"."_$chain";}
	if ($line=~/^\s*(\S*)\;\s*(\S*)\;/) {
	    $pdbid="$pdbid($2-$1)";
	}
	push(@PDBs,$pdbid);
	$pdbids{$name}.=" $pdbid";
    } elsif ($line=~/^\#=GC SS_cons\s+(\S+)/) {
	# Read secondary structure consensus
	$dssp=$1;

	# Rectify capitalization of $dssp
	@res=();
	@res=split(//,$dssp); # unpack doesn't work because it generates strings instead of characters
	# Check length
	if (scalar(@res) != scalar(@match)) {
	    printf(STDERR "WARNING: in HMM $AC sequence $name has wrong length (%i instead of %i):\n%s\n%s\n",scalar(@res),scalar(@match),$qseq,$dssp);
	    next;
	}
	for ($i=0; $i<@match; $i++) {
	    if ($match[$i]) {
		$res[$i]=~tr/a-z./A-Z-/;
	    } else {
		$res[$i]=~tr/A-Z-/.........................../;
	    }
	}
	$dssp = join("",@res);
	$dssp=~tr/xX/.-/; # transform X (no SS) into - (unknown SS)
	if ($format eq "a3m") {$dssp=~tr/.//d;}

    } elsif ($line!~/^\#/) {

	# Read sequence?
	if ($line!~/^(\S+)\s+(\S*)/) {
	    chomp($line); print("Error in pfam-format. Expecting sequence in '$line'\n");
	} elsif ($ID ne "") {
	    # Read first sequence in a new alignment

	    $name=$1;
	    $seq=$2;
	    
	    if ($infile=~/Pfam.?B/i) {
		# Make everything a match state when reading PfamB
		$seq=~tr/./-/;
		@res=(); @match=();
		@res=split(//,$seq); # unpack doesn't work because it generates strings instead of characters
		for ($i=0; $i<@res; $i++) {
		    $match[$i]=1;
		}		
	    } else {
		# Store length and which columns are match states
		@res=(); @match=();
		@res=split(//,$seq); # unpack doesn't work because it generates strings instead of characters
		for ($i=0; $i<@res; $i++) {
		    $match[$i]= ($res[$i]=~tr/A-Z-/A-Z-/);
		}
	    }

	    @inserts=();
	    @inserts=split(/[A-Z-]/,$seq);

	    if ($format eq "a3m") {$seq=~tr/.//d;}
	    @namelines=();
	    @names=();
	    @seqs=();
	    if ($infile=~/Pfam.?B/i) {
		$titleline="#$DR $ID: $DE$DC";
	    } else {
		$titleline="#$AC $ID: $DE$DC";
	    }
	    $titleline =~s /\s+The following Pfam-B.*?\.//;
	    # Add InterPro reference and annotation?
	    if (defined $id2ipr{$AC}) {  
		$titleline .= ";  InterPro: ".$id2ipr{$AC}." ".$ipr2abstract{$id2ipr{$AC}};
		if (defined $ipr2GO{$id2ipr{$AC}} && $ipr2GO{$id2ipr{$AC}} ne "") {
		    $titleline .= "; ".$ipr2GO{$id2ipr{$AC}};
		}
	    }
	    # Add pdb ids?
	    if (@PDBs) {
		my %PDBs=();
		$i=0;
		$titleline.="; PDB:";
		foreach my $PDB (@PDBs) {
		    if ($i>=10) {$titleline.=" ..."; last;}  # if there are more than 10 pdb ids print ...
		    $PDB=~/^(....)/;
		    if (defined $PDBs{$1}) {next;}
		    $i++;
		    $PDBs{$1} = 1;
		    $titleline.=" ".$PDB;
		}
		$titleline.=".";
	    }
	    $titleline.="\n";
	    $ID=""; $dssp=""; $nseq=1;

	    if (defined $pdbids{$name}) {
		$nameline=$name." PDB: ".$pdbids{$name};
		unshift(@namelines,$nameline);
		unshift(@names,$name);
		unshift(@seqs,$seq);
		$qpdb=1;
	    } else {
		$nameline=$name;
		push(@namelines,$nameline);
		push(@names,$name);
		push(@seqs,$seq);
		$qpdb=0;
	    }
#	    print(">$AC ($1) $ID: $DE$DC$PDB\n$2\n");

	    # First sequence is candidate for query/seed sequence
	    $qseq=$seq;
	    $qname=$name;
	    $qnameline=$nameline;
	    if ($infile=~/pfam.?B/i) {
		$qmatch=($qseq=~tr/.-/.-/);     # number of gaps in candidate query sequence!!
	    } else {
		$qmatch=($qseq=~tr/a-z-/a-z-/); # number of deletions plus inserts in candidate query sequence!!
	    }
#	    printf("%-12.12s %-4i      %s\n",$name,$qmatch,$seq);


	    printf("%3i: $basename\n",++$Ndom);

	} else {
	    # Read next sequence of alignment
	    
	    $name=$1;
	    $seq=$2;
	    $nseq++;

	    # Rectify capitalization of $seq
	    @res=();
	    @res=split(//,$seq); # unpack doesn't work because it generates strings instead of characters
	    # Check length
	    if (scalar(@res) != scalar(@match)) {
		printf(STDERR "WARNING: in HMM $AC sequence $name has wrong length (%i instead of %i):\n%s\n%s\n",scalar(@res),scalar(@match),$qseq,$seq);
		next;
	    }
	    for ($i=0; $i<@match; $i++) {
		if ($match[$i]) {
		    $res[$i]=~tr/a-z./A-Z-/;
		} else {
		    $res[$i]=~tr/A-Z-/a-z./;
		}
	    }
	    $seq = join("",@res);
	    
	    if ($format eq "a3m") {$seq=~tr/.//d;}
	    
	    if (defined $pdbids{$name}) {
		$nameline=$name." PDB:".$pdbids{$name};
		unshift(@namelines,$nameline);
		unshift(@names,$name);
		unshift(@seqs,$seq);
	    } else {
		$nameline=$name;
		push(@namelines,$nameline);
		push(@names,$name);
		push(@seqs,$seq);
	    }
#	    print(">$name\n$seq\n");


	    # Does current sequence have fewer deletions than $qseq?
	    my $seqpdb;
	    if (defined $pdbids{$name}) {$seqpdb=1;} else {$seqpdb=0;}
	    my $tmatch;
	    if ($infile=~/pfam.?B/i) {
		$tmatch=($seq=~tr/.-/.-/);     # number of gaps in candidate query sequence!!
	    } else {
		$tmatch=($seq=~tr/a-z-/a-z-/); # number of deletions plus inserts in candidate query sequence!!
	    }
	    if ( (!$qpdb && $seqpdb) || ($qpdb==$seqpdb && $tmatch<$qmatch) ) {
		$qseq=$seq; 
		$qname=$name; 
		$qnameline=$nameline; 
		$qmatch=$tmatch;
		if (defined $pdbids{$name}) {$qpdb=1;}
	    }
#	    printf("%-12.12s %-3i %-3i %s\n",$qname,($qseq=~tr/A-Z/A-Z/),$qmatch,$qseq);
#	    printf("%-12.12s %-3i %-3i %s\n",$name,($seq=~tr/A-Z/A-Z/),$tmatch,$seq);

	}
    }
}
close (INFILE);






##############################################################################################
# Run SS prediction for $file starting from alignment in $file.psi
##############################################################################################
sub RunPsipred() {
    # This is a simple script which will carry out all of the basic steps
    # required to make a PSIPRED V2 prediction. Note that it assumes that the
    # following programs are in the appropriate directories:
    # blastpgp - PSIBLAST executable (from NCBI toolkit)
    # makemat - IMPALA utility (from NCBI toolkit)
    # psipred - PSIPRED V2 program
    # psipass2 - PSIPRED V2 program

    my $file=$_[0];
    my $basename;  #file name without extension
    my $rootname;  #basename without directory path
    
    if ($file =~/(.*)\..*/)    {$basename=$1;} else {$basename=$file;}
    if ($basename=~/.*\/(.*)/) {$rootname=$1;} else {$rootname=$basename;}

    system("cp -f $file $basename.fasta");
    
    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile

    # Create dummy database?
    if (!-e "$dummydb.phr") {&System("cp $file $dummydb"); &System("$ncbidir/formatdb -i $dummydb");}
    

#    &System("$blastpgp -b 0 -j 2 -h 0.001 -d $dbname -i $file -B $basename.psi -C $basename.chk >& $basename.log");
    &System("$blastpgp -a 2 -b 0 -j 1 -h 0.001 -d $dummydb -i $file -B $basename.psi -C $basename.chk >& $basename.log");
    
    system("echo $rootname.chk > $basename.pn\n");
    system("echo $rootname.fasta > $basename.sn\n");
    system("$ncbidir/makemat -P $basename");
    
    if (-e "$datadir/weights.dat4") { # Psipred version < 3.0
	&System("$execdir/psipred $basename.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $basename.ss");
    } else {
	&System("$execdir/psipred $basename.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $basename.ss");
    }

    &System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $basename.ss2 $basename.ss > $basename.horiz");
    
    # Remove temporary files
    system ("rm -f $basename.pn $basename.sn $basename.mn $basename.chk $basename.fasta $basename.mtx $basename.aux $basename.ss $basename.ss2");
    return;
}



################################################################################################
### Read InterPro file
################################################################################################
sub ReadInterProFiles()
{
    # Read tables IPR_ID -> abstracts and member_db_id -> IPR_ID 
    print("Read abstracts and InterPro ids from interpro.xml\n");
    my $ipr=""; 
    my $abstract;
    my $dbkey;
    open (INFILE, "<$iprodir/interpro.xml") || die("ERROR: cannot open: $!\n");
    while ($line=<INFILE>) {
	if    ($line=~/^\s*<interpro id=\"(\S+)\"/) {$ipr=$1;}
	elsif ($line=~/^\s*<abstract>/)  {
	    $abstract="";
	    while ($line=<INFILE>) {
		if ($line=~/^\s*<\/abstract>/) {last;}
		$abstract.=$line;
	    }
	}
	elsif ($line=~/^\s*<member_list>/)  {
	    while ($line=<INFILE>) {
		if ($line=~/dbkey=\"(\S+)\"/) {
		    $id2ipr{$1}=$ipr;
		} else {
		    last;
		}
	    }
	}
	elsif ($line=~/^\s*<\/interpro>/) {
	    chomp($abstract);
	    $abstract =~s/\s+/ /g;          # remove newlines and multiple spaces
	    $abstract =~s/<\/?taxon.*?>//g; # remove xml tags
	    $abstract =~s/<\/?cite.*?>//g;  # remove xml tags
	    $abstract =~s/<db_xref db=\"(\S+?)\" dbkey=\"(\S+)\"\/>/$2 from $1/g;  # remove xml tags
	    $abstract =~s/&lt;p&gt;/  /g;   # remove html tags
	    $abstract =~s/&lt;\/p&gt;//g;   # remove html tags
	    $abstract =~s/&lt;sub&gt;/_/g;  # remove html tags
	    $abstract =~s/&lt;sup&gt;/^/g;  # remove html tags
	    $abstract =~s/&lt;\/?\S{1,4}&gt;//g;  # remove html tags
	    $abstract =~s/<\/?\S{1,4}>//g;  # remove xml tags
	    $ipr2abstract{$ipr}=$abstract;
	}
    }
    close(INFILE);

    # Read interpro2go into hash
    print("Read interpro2go into hash\n");
    open (INFILE, "<$iprodir/interpro2go") || die("ERROR: cannot open: $!\n");
    while ($line=<INFILE>) {
	if ($line=~/^InterPro:(\S+)\s+.*?>\s+GO:(.*)\s+;\s+GO:(\d+)/) {
	    if ($3 ne "0005554") { # GO:0005554 "molecular function unknown"
		if (! defined $ipr2GO{$1}) {
		    $ipr2GO{$1}="GO: $3 $2";
		} else {
		    $ipr2GO{$1}.=", $3 $2";
		}
	    }
	}
    }
    close(INFILE);
}

################################################################################################
### System command
################################################################################################
sub System()
{
    my $command=$_[0];
    if ($v>=2) {print("$command\n");} 
    return system($command)/256; # && die("\nERROR: $!\n\n"); 
}

