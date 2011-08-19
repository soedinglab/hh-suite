#! /usr/bin/env perl

my $rootdir;
BEGIN {
    if (defined $ENV{TK_ROOT}) {$rootdir=$ENV{TK_ROOT};} else {$rootdir="/cluster";}
};
use lib "$rootdir/bioprogs/hhpred";
use lib "/cluster/soeding/perl";     # for cluster
use lib "/cluster/lib";              # for chimaera webserver: ConfigServer.pm

use strict;
use MyPaths;

# Default parameters
my $v=2;
my $help="
 Generate FASTA nonredundant sequence file from SEQRES records of globbed pdb files.
 For updating purposes, you can write only those sequences to pdb_new.fas 
 that are not already contained in an old file by giving as third argument the old pdb.fas file:

 Usage:   pdb2fasta.pl 'pdb-fileglob' pdb_newseqs.fas [options]
 Options:
  -u oldfile   update: write only those sequences to pdb_new.fas that are not contained in oldfile
  -dali dir    read FoldIndex.html and domain_definitions.txt in DALI directory and list DALI fold(s) in sequence name
  -scop file   read dir.cla.scop.txt_1.65 and list SCOP fold(s) in sequence name
  -v int       verbose mode
  -t MTH-YR    use only structures released until the given month and year, e.g. APR-04 or SEP-98
  -all         include all sequences instead of nonredundant set

 Example: pdb2fasta.pl '/raid/db/pdb/*.ent' /data/pdbfas/pdb_20Oct2004.fas
 Example: pdb2fasta.pl '/raid/db/pdb/*.ent' pdb_new.fas -u pdb.fas -dali /data/dali -scop /data/scop/dir.cla.scop.txt_1.65
\n";
my $TOTLEN=160;   # maximum length of name, description, and keywords
my $DESCLEN=80;   # maximum length of description


if (@ARGV<2) {die($help);}

my @pdbfiles;        
my $newseqfile;
my $oldfile="";
my $dalidir="";
my $scopfile="";
my $date="";
my %months=("JAN"=>1,"FEB"=>2,"MAR"=>3,"APR"=>4,"MAY"=>5,"JUN"=>6,"JUL"=>7,"AUG"=>8,"SEP"=>9,"OCT"=>10,"NOV"=>11,"DEC"=>12);
my %oldpdbids=();  # hash contains all pdbids in $oldfile
our $pdbfile;
my $pdbid;         # four-letter PDB identifier, e.g. 1hz4
my $resolution;    # experimental resolution in Angstrom
my $rvalue;        # R-value
my $free_rvalue;   # free R-value
my $molid=0;       # molecule id (for multichain structures)
my $length;        # number of residues in a chain
my @seqres=();     # three-letter code of chain currently read in 
my $seqres;        # 
my %descript;      # $descript{"A"} contains the description for chain A
my $descript;
my %organism;      # $organism{"A"} contains the organism for chain A
my $organism;    
my $organism_common;    
my @chain;         # $chain[$molid]
my $chain;         # either A for chain A or "" if no chain id 
my @chains;        # 
my @keywds;        # keywords for the structure
my $keywds;        # keywords for the structure
my $token;
my $synonym;       # read from COMPND SYNONYM records of pdb files
my @synonyms;      # read from COMPND SYNONYM records of pdb files
my $line;          # line read in from file
my @sequences=();  # contains all sequences of chains to be printed to outfile
my @resolution=(); # $resolution[$nc] contains resolution of $nc'th sequence, where $nc=$nchains{$seqres}  
my %nchains;       # $nchains{$seqres} is index in @sequences of sequence with these residues
my $nchains=0;     # number of chains written to $newseqfile
my $k=0;           # counts pdb files processed
my @equiv_pdbs;    # list of pdbids (including _chain) with identical residues (maximum one pdbid_chain per pdb file)
my %dalifamids=(); # $foldids{$pdbid} contains a list of (one or more) DALI or SCOP foldids
my %scopfamids=(); # $foldids{$pdbid} contains a list of (one or more) DALI or SCOP foldids
my $het;           # $het  contains list of hetero ligands with at least 10 atoms in current pdb file (e.g. "DAC,PTR") 
my %words;         # for debugging upper case -> lower case
my $nr=1;          # 1: create nonredundant set  0:do not eliminate redundant sequences
my %three2one=(
	       "ALA"=>"A","VAL"=>"V","PHE"=>"F","PRO"=>"P","MET"=>"M","ILE"=>"I","LEU"=>"L","ASP"=>"D","GLU"=>"E","LYS"=>"K",
	       "ARG"=>"R","SER"=>"S","THR"=>"T","TYR"=>"Y","HIS"=>"H","CYS"=>"C","ASN"=>"N","GLN"=>"Q","TRP"=>"W","GLY"=>"G",
	       "2AS"=>"D","3AH"=>"H","5HP"=>"E","ACL"=>"R","AIB"=>"A","ALM"=>"A","ALO"=>"T","ALY"=>"K","ARM"=>"R","ASA"=>"D",
	       "ASB"=>"D","ASK"=>"D","ASL"=>"D","ASQ"=>"D","AYA"=>"A","BCS"=>"C","BHD"=>"D","BMT"=>"T","BNN"=>"A","BUC"=>"C",
	       "BUG"=>"L","C5C"=>"C","C6C"=>"C","CCS"=>"C","CEA"=>"C","CHG"=>"A","CLE"=>"L","CME"=>"C","CSD"=>"A","CSO"=>"C",
	       "CSP"=>"C","CSS"=>"C","CSW"=>"C","CXM"=>"M","CY1"=>"C","CY3"=>"C","CYG"=>"C","CYM"=>"C","CYQ"=>"C","DAH"=>"F",
	       "DAL"=>"A","DAR"=>"R","DAS"=>"D","DCY"=>"C","DGL"=>"E","DGN"=>"Q","DHA"=>"A","DHI"=>"H","DIL"=>"I","DIV"=>"V",
	       "DLE"=>"L","DLY"=>"K","DNP"=>"A","DPN"=>"F","DPR"=>"P","DSN"=>"S","DSP"=>"D","DTH"=>"T","DTR"=>"W","DTY"=>"Y",
	       "DVA"=>"V","EFC"=>"C","FLA"=>"A","FME"=>"M","GGL"=>"E","GLZ"=>"G","GMA"=>"E","GSC"=>"G","HAC"=>"A","HAR"=>"R",
	       "HIC"=>"H","HIP"=>"H","HMR"=>"R","HPQ"=>"F","HTR"=>"W","HYP"=>"P","IIL"=>"I","IYR"=>"Y","KCX"=>"K","LLP"=>"K",
	       "LLY"=>"K","LTR"=>"W","LYM"=>"K","LYZ"=>"K","MAA"=>"A","MEN"=>"N","MHS"=>"H","MIS"=>"S","MLE"=>"L","MPQ"=>"G",
	       "MSA"=>"G","MSE"=>"M","MVA"=>"V","NEM"=>"H","NEP"=>"H","NLE"=>"L","NLN"=>"L","NLP"=>"L","NMC"=>"G","OAS"=>"S",
	       "OCS"=>"C","OMT"=>"M","PAQ"=>"Y","PCA"=>"E","PEC"=>"C","PHI"=>"F","PHL"=>"F","PR3"=>"C","PRR"=>"A","PTR"=>"Y",
	       "SAC"=>"S","SAR"=>"G","SCH"=>"C","SCS"=>"C","SCY"=>"C","SEL"=>"S","SEP"=>"S","SET"=>"S","SHC"=>"C","SHR"=>"K",
	       "SOC"=>"C","STY"=>"Y","SVA"=>"S","TIH"=>"A","TPL"=>"W","TPO"=>"T","TPQ"=>"A","TRG"=>"K","TRO"=>"W","TYB"=>"Y",
	       "TYQ"=>"Y","TYS"=>"Y","TYY"=>"Y","AGM"=>"R","GL3"=>"G","SMC"=>"C","ASX"=>"B","CGU"=>"E","CSX"=>"C","GLX"=>"Z",
	       "LED"=>"L"
	       );

# Read command line options
my $options="";
for (my $i=0; $i<=$#ARGV; $i++) {$options.=" $ARGV[$i]";}
if ($options=~s/ -u\s+(\S+)//) {$oldfile=$1;} 
if ($options=~s/ -dali\s+(\S+)//) {$dalidir=$1;} 
if ($options=~s/ -scop\s+(\S+)//) {$scopfile=$1;} 
if ($options=~s/ -v\s*(\d+)//) {$v=$1;} 
if ($options=~s/ -v//) {$v=2;} 
if ($options=~s/ -all//) {$nr=0;} 
if ($options=~s/ -t (\w\w\w)-(\d\d)//) {$date=($months{$1}-1)/12+$2+100*($2<50);} 
if ($options=~s/^\s*([^- ]\S+)\s*//) {
    if ($v>=2) {print("Globbing...")};
    @pdbfiles=glob($1);
    if ($v>=2) {print(" found ".scalar(@pdbfiles)." files\n")};
}
if ($options=~s/^\s*([^- ]\S+)\s*//) {$newseqfile=$1;} 

# Warn if unknown options found or no infile/newseqfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!@pdbfiles)   {print($help); print("Error: no input files given\n"); exit;}
if (!$newseqfile) {print($help); print("Error: no output file given\n"); exit;}

# Updating option?
if ($oldfile) {
    # Reading pdb codes from $oldfile
    if ($v>=3) {printf("Reading pdb codes from $oldfile ... \n");}
    open (OLDFILE,"<$oldfile") || die ("ERROR: cannot open $oldfile for writing: $!\n");
    while ($line=<OLDFILE>) {
	if ($line=~/^>(\S\S\S\S)/o) {$oldpdbids{$1}=1;}
    }
    close(OLDFILE);
}

# Add fold identifiers?
if ($dalidir) {&ReadDaliFiles();}  
if ($scopfile) {&ReadScopFile();}  

############################################################################################
# Read one pdb file after the other
foreach $pdbfile (@pdbfiles) {

    $k++;

    if ($pdbfile=~/^.*\/(.*?)$/) {$pdbid=$1;} else {$pdbid=$pdbfile;}  # remove path
    if ($pdbid=~/^(pdb)?(.*)\..*$/) {$pdbid=lc($2);} else {die("Error: globbed file $pdbfile has no extension\n");}
    if (exists $oldpdbids{$pdbid}) {next;}
    if ($v>=1) {
	print("."); 
	if (!($k%100)) {print("$k\n");}
    } elsif ($v>=2) {printf("Reading %4i %s\n",$k,$pdbfile);}

    open (PDBFILE, "<$pdbfile") || die ("Error: couldn't open $pdbfile: $!\n");
    if ($v>=4) {print("Reading $pdbfile...\n");}
    $line=<PDBFILE>;

    # Initialize before reading new pdb file
    $resolution=0;
    $rvalue=0;
    $free_rvalue=0;
    $molid=0;
    $token="MOLECULE";
    $descript="";
    $organism="";
    $organism_common="";
    $keywds="";    
    $chain="";
    %organism=();
    %descript=();
    @keywds=();
    @chain=();
    @seqres=();
    $synonym="";
    @synonyms=();
    $het="";   # will contain list of hetero groups (if found)

    # COMPND    ASPARTATE AMINOTRANSFERASE (E.C.2.6.1.1) WILD TYPE COMPLEXED  1ASA   3
    # COMPND   2 WITH PYRIDOXAL-5'-PHOSPHATE AND MALEATE                      1ASA   4
    #  or
    # COMPND    MOL_ID: 1;                                                            
    # COMPND   2 MOLECULE: ACTIN-LIKE PROTEIN 3;                                      
    # COMPND   3 CHAIN: A;                                                            
    # COMPND   4 SYNONYM: ARP3; ACTIN-RELATED PROTEIN 3; ACTIN-2;                     
    # COMPND   5 OTHER_DETAILS: PART OF THE ARP2/3 COMPLEX;                           
    # COMPND   6 MOL_ID: 2;                                                           
    # COMPND   7 MOLECULE: ACTIN-LIKE PROTEIN 2;                                      
    # COMPND   8 CHAIN: B;                                                            
    # COMPND   9 SYNONYM: ARP2; ACTIN-RELATED PROTEIN 2;                              
    # COMPND  10 OTHER_DETAILS: PART OF THE ARP2/3 COMPLEX;   
    #  or
    # COMPND    MOL_ID: 1;                                                            
    # COMPND   2 MOLECULE: PHOSPHATE SYSTEM POSITIVE REGULATORY PROTEIN               
    # COMPND   3 PHO4;                                                                
    # COMPND   4 CHAIN: A, B;                                                         
    # COMPND   5 FRAGMENT: DNA BINDING DOMAIN;                                        
    # COMPND   6 SYNONYM: BHLH;                                                       
    # COMPND   7 ENGINEERED: YES;                                                     
    # COMPND   8 BIOLOGICAL_UNIT: DIMER;                                              
    while ($line && $line!~/^COMPND /o && $line!~/^REMARK /o) {$line=<PDBFILE>;}			     
    if (!$line || $line!~/^COMPND /) 
    {
	if ($v>=2) {print("\n\nWarning: no COMPND line found in $pdbfile; skipping pdb file\n");}
	next;
    }
    while ($line && $line=~/^COMPND /o && $line!~/^REMARK /o) {
	$line=~s/^(.{70}).*/$1/;
	if ($line=~/^COMPND\s+\d*\s+MOL_ID:\s*(\d+)/) {
	    if ($molid>0) {
		&SetDescript();
	    } else {$molid=$1;}
	    $descript="";
	    $synonym="";
	}
	elsif ($line=~/^COMPND\s+\d*\s+MOLECULE:\s*(.*\S)/) {$descript=$1; $token="MOLECULE";}
	elsif ($line=~/^COMPND\s+\d*\s+CHAIN:\s*(.*\S)/)    {$chain=$1; $token="CHAIN";}
	elsif ($line=~/^COMPND\s+\d*\s+SYNONYM:\s*(.*\S)/)   {$synonym=$1;} 
	elsif ($line=~/^COMPND\s+\d*\s+FRAGMENT:/) {$token="";}	
	elsif ($line=~/^COMPND\s+\d*\s+EC:/) {$token="";}
	elsif ($line=~/^COMPND\s+\d*\s+ENGINEERED:/) {$token="";}
	elsif ($line=~/^COMPND\s+\d*\s+MUTATION:/) {$token="";}
	elsif ($line=~/^COMPND\s+\d*\s+BIOLOGICAL UNIT:/) {$token="";}
	elsif ($line=~/^COMPND\s+\d*\s+OTHER_DETAILS:/) {$token="";}
	else {
	    $line=~/^COMPND\s+\d*\s+(.*\S)/;
	    if ($token eq "MOLECULE") {$descript.=" ".$1;}
	    elsif ($token eq "SYNONYM") {$synonym.=" ".$1;}
	    elsif ($token eq "CHAIN") {$chain.=" ".$1;}
	}
	$line=<PDBFILE>;
    }
    &SetDescript();
    if (!$line) {if ($v>=2) {print("\nFormat error in $pdbfile. Skipping file ...\n");} next;} 


    # SOURCE    (ESCHERICHIA COLI)                                            1ASA   5
    #  or
    # SOURCE    MOL_ID: 1;                                                            
    # SOURCE   2 ORGANISM_SCIENTIFIC: BOS TAURUS;                                     
    # SOURCE   3 ORGANISM_COMMON: BOVINE;                                             
    # SOURCE   4 ORGAN: THYMUS;                                                       
    # SOURCE   5 MOL_ID: 2;                                                           
    # SOURCE   6 ORGANISM_SCIENTIFIC: BOS TAURUS;                                     
    # SOURCE   7 ORGANISM_COMMON: BOVINE;                                             
    # SOURCE   8 ORGAN: THYMUS;                                                       
    $molid=0;
    $token="ORGANISM";
 #   $organism="Synthetic?";
    while ($line && $line!~/^SOURCE /o && $line!~/^REMARK /o) {$line=<PDBFILE>;}			     
    if (!$line) {if ($v>=2) {print("\nFormat error in $pdbfile. Skipping file ...\n");} next;} 
    if ($v>=2 && $line!~/^SOURCE /) {print("\n\nWarning: no SOURCE line found in $pdbfile\n");}
    while ($line=~/^SOURCE /o && $line!~/^REMARK /o) {
	$line=~s/^(.{70}).*/$1/;
	if ($line=~/^SOURCE\s+\d*\s+MOL_ID:\s*(\d+)/) {
	    if ($molid>0) { 
		&SetOrganism();
	    } else {$molid=$1;}
	} 
	elsif ($line=~/^SOURCE\s+\d*\s+ORGANISM_SCIENTIFIC:\s*(.*\S)/) {$organism=$1; $token="ORGANISM";}
	elsif ($line=~/^SOURCE\s+\d*\s+SYNTHETIC/) {if ($organism eq "") {$organism="Synthetic"; $token="";}}
	elsif ($line=~/^SOURCE\s+\d*\s+FRAGMENT:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+ORGANISM_COMMON:\s+(.*\S)/) {$organism_common=$1; $token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+STRAIN:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+VARIANT:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+CELL_LINE:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+ATCC:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+ORGAN:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+TISSUE:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+CELL:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+ORGANELLE:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+SECRETION:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+CELLULAR_LOCATION:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+PLASMID:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+GENE:/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+EXPRESSION_/) {$token="";}
	elsif ($line=~/^SOURCE\s+\d*\s+OTHER_DETAILS:/) {$token="";}
	else {
	    $line=~/^SOURCE\s+\d*\s+(.*\S)/;
	    if ($token eq "ORGANISM") {$organism.=$1;}
	}
	$line=<PDBFILE>;
    }
    &SetOrganism();
    if (!$line) {if ($v>=2) {print("\nFormat error in $pdbfile. Skipping file ...\n");} next;} 


    # KEYWDS    KETOLISOMERASE, XYLOSE METABOLISM, GLUCOSE-FRUCTOSE
    # KEYWDS   2 INTERCONVERSION, HYDRIDE TRANSFER, ALPHA-BETA BARREL,
    # KEYWDS   3 METALLOENZYME, THERMOPHILE
    while ($line && $line!~/^KEYWDS /o && $line!~/^REMARK /o) {$line=<PDBFILE>;}			     
    if (!$line) {if ($v>=2) {print("\nFormat error in $pdbfile. Skipping file ...\n");} next;} 
    while ($line=~/^KEYWDS /o) {
	$line=~s/^(.{70}).*/$1/;
	if ($line=~/^KEYWDS\s+\d*\s+(.*\S)/) {$keywds.=" ".$1;}
	$line=<PDBFILE>;
    }
    $keywds=~s/CRYSTAL STRUCTURE,?\s*//i;
    $keywds=~s/X-RAY STRUCTURE,?\s*//i;
    $keywds=~s/THREE-DIMENSIONAL STR.CTURE,?\s*//i;
    $keywds=~s/NMR,?\s*//i;
    if ($keywds) {
	$keywds=~s/\s+/ /g;
	$keywds=~s/^\s+/ /;
	$keywds=~s/\s*;?$/;/;
	@keywds=split(/[,;]\s+/,$keywds);
    } else {@keywds=();}

    # Include keywords up to 50 chars taken
    my @these_keywds=@keywds;
    if (@these_keywds) {

	# Remove keywords that are substring of description or organism
	for (my $k=0; $k<@these_keywds; $k++) {
	    $keywds=$these_keywds[$k];
	    my $ddescript=$descript;
	    my $kkeywds=$keywds;
	    my $oorganism=$organism;
	    $ddescript=~tr/a-zA-Z//cd;
	    $kkeywds=~tr/a-zA-Z//cd;
	    $oorganism=~tr/a-zA-Z//cd;
	    if ($ddescript=~/$kkeywds/i) {splice(@these_keywds,$k,1); $k--;}
	    elsif ($kkeywds=~/$ddescript/i) {splice(@these_keywds,$k,1); $k--;}
	    elsif ($oorganism=~/$kkeywds/i) {splice(@these_keywds,$k,1); $k--;}
	   
	}

	# Add keywords until length limitation is exceeded ($TOTLEN chars)
	if (@these_keywds) {
	    $keywds=$these_keywds[0];
	    for (my $k=1; $k<@these_keywds && length($descript.$keywds.", ".$these_keywds[$k].$organism)<$TOTLEN; $k++) {
		$keywds.=", ".$these_keywds[$k];
	    }
	    $keywds=~s/^(\S)/ $1/;  # add space at first position
	} else {$keywds="";}
    } else {$keywds="";}
    foreach my $chain (keys(%descript)) {$descript{$chain}.=$keywds;}

    # Check date?
    if ($date) {
	while ($line && $line!~/^REVDAT   1/o && $line!~/^REMARK /o) {$line=<PDBFILE>;}			     
	if ($line=~/^REVDAT   1/) {
	    if ($line=~/^.{16}(\w\w\w)-(\d\d)/) {
		my $thisdate=($months{$1}-1)/12+$2+100*($2<50);
#	        print("This date: $thisdate  date=$date\n$line");
		if ($thisdate>$date) {next;}
	    } elsif ($v>=2) { 
		print("WARNING: no valid date in header: \n$line");
	    } 
	}
    }

    # REMARK   2 RESOLUTION. 2.00 ANGSTROMS.                                          
    while ($line && $line!~/^REMARK   2 /o && $line!~/^SEQRES /o) {$line=<PDBFILE>;}			     
    if (!$line) {if ($v>=2) {print("\nFormat error in $pdbfile. Skipping file ...\n");} next;} 
    if ($v>=2 && $line!~/^REMARK   2 /) {print("\n\nWarning: no REMARK   2 line found in $pdbfile\n");}
    while ($line=~/^REMARK   2 /o && $line!~/^SEQRES /o) {
	if ($line=~/^REMARK   2\s+RESOLUTION\.\s+(\d+\.?\d*)/) {$resolution=$1; last;}
	$line=<PDBFILE>;
    }


    # REMARK   3   R VALUE            (WORKING SET) : 0.216                           
    # REMARK   3   FREE R VALUE                     : 0.251                           
    while ($line && $line!~/^REMARK   3 /o && $line!~/^SEQRES /o) {$line=<PDBFILE>;}			     
    if (!$line) {if ($v>=2) {print("\nFormat error in $pdbfile. Skipping file ...\n");} next;} 
    if ($v>=2 && $line!~/^REMARK   3 /) {print("\n\nWarning: no REMARK   3 line found in $pdbfile\n");}
    while ($line=~/^REMARK   3 /o && $line!~/^SEQRES /o) {
	if ($line=~/^REMARK   3\s+R VALUE\s+\(WORKING SET\)\s*:\s*(\d+\.?\d*)/o) {$rvalue=$1;}
	if ($line=~/^REMARK   3\s+FREE R VALUE\s*:\s*(\d+\.?\d*)/o) {$free_rvalue=$1;}	
	$line=<PDBFILE>;
    }

    # Record current position in PDBFILE
    my $file_pos = tell(PDBFILE);
    
    # Search for hetero groups BEFORE adding seqeunces => read pdb file twice :(
    while ($line && $line!~/^HET /o && $line!~/^ATOM /o) {$line=<PDBFILE>;}			     
    if (defined $line && $line=~/^HET /o) {
	while ($line) {
#           ----+----1----+----2----+----3
#	    HET    DAC  A 172      18      
	    if ($line=~/^HET\s+(\S+) ..........\s*(\d+)/) {
		if ($2>=10 ) {
		    my $this_het=$1;
		    if ($het!~/$this_het/) { # don't list any hetgoup twice
			if ($het eq "") {$het=" HET: $1";} else {$het.=" $1";}
		    }
		}
	    } else {
		last;
	    }
	    $line=<PDBFILE>;
	}
    }
    if ($het ne "") {$het.=";";}

    # Rewind the current position in PDBFILE back before SEQRES records
    seek (PDBFILE,$file_pos,0);
 
    # SEQRES   1    396  MET PHE GLU ASN ILE THR ALA ALA PRO ALA ASP PRO ILE  1ASA  60
    # SEQRES   1 A  366  SER ARG MET PRO SER PRO PRO MET PRO VAL PRO PRO ALA          
    @seqres=();
    my $newchain="@"; # make sure that sequence is not printed before first chain has been read
    my $newlength;    # compare previous to current chain to find out when one chain is finished
    while ($line && $line!~/^SEQRES /o) {$line=<PDBFILE>;}			     
    if (!defined $line) {
	if ($v>=2) {print("\n\nWarning: no SEQRES line found in $pdbfile. Skipping file ...\n");}
	close (PDBFILE);
	next;
    }
    while ($line=~/^SEQRES /o) {
	if ($line=~/^SEQRES\s+\d+\s+(\S?)\s+(\d+)\s+(.{51})/o) {
	    $chain=$newchain;
	    $newchain=$1;
	    $length=$newlength;
	    $newlength=$2;
	    $seqres=$3;
	    $seqres=~s/\s*$//;
#	    print($line);
#	    printf("line $.: prevlen=%-4i counted=%-4i newchain=%1s  chain=%1s \n\n",$length,scalar(@seqres),$newchain,$chain);

	    # Compare previous to current chain to find out when one chain is finished
	    if ($chain ne $newchain && $chain ne "@") {
		if (scalar(@seqres)!=$length) {
		    if ($v>=2) {printf("\nWarning: in $pdbfile, line $., sequence length=$length, counted residues = %i\n",scalar(@seqres));}
		}
		&AddSequence();
		@seqres=();
	    } 
	    push(@seqres,split(/\s+/,$seqres));
	} else {
	    print("\nError: found invalid SEQRES record in $pdbfile, line $. : line=$line"); 
	    next;
	}
	$line=<PDBFILE>;
    }
    $chain=$newchain;
    $length=$newlength;
    &AddSequence();

    close (PDBFILE);

} # end foreach $pdbfile
############################################################################################

# Print all sequences
open (NEWSEQFILE,">$newseqfile") || die ("ERROR: cannot open $newseqfile for writing: $!\n");
for (my $nc=0; $nc<@sequences; $nc++) {
    if ($equiv_pdbs[$nc] ne "") {
	$sequences[$nc]=~s/^(.*)/$1 PDB:$equiv_pdbs[$nc]/; # Add list of equivalent pdb codes
    }
    printf(NEWSEQFILE "%s",$sequences[$nc]);
}
close(NEWSEQFILE);


foreach my $word (keys(%words)) {print("$word ");}
print("\n");
print("Written $nchains chains to $newseqfile\n");

exit;


##################################################################################
# Set description when a new MOL_ID line is found, or at the end of COMPND records
##################################################################################
sub SetDescript() {
    my $i;
#   print("chain=$chain\n");
    if ($chain eq "NULL;")  {$chain="";}
    $chain[$molid]=$chain;
    $molid=$1;

    if ($descript=~/^DNA[ ;]/) {$het=" HET: DNA";}

    # Cut description down to max. $DESCLEN letters
    $descript=~s/\s*;\s*$//;
    if (length($descript)>$DESCLEN) {
	$descript=~s/(.{$DESCLEN}\S*).*/$1.../; # remove everything after first comma
    }

    # Add synonyms with a maximum of 16 letters to description
    if ($synonym ne "") {
	@synonyms=split(/;\s+/,$synonym);
	for ($i=0; $i<scalar(@synonyms); $i++) {
	    if (length($synonyms[$i])>16) {
		splice(@synonyms,$i,1);
	    } else {
		$synonyms[$i]=~s/;\s*//;
		$token="SYNONYM";
	    }
	}
    }

    # Choose shortname (<=16 characters) from synonyms
    unshift(@synonyms,$descript);
    for ($i=0; $i<scalar(@synonyms); $i++) {
	if (length($synonyms[$i])<=16) {last;}
    }
    if ($i>=scalar(@synonyms)) {$i=0;}
    $descript=$synonyms[$i];
    $descript=~s/^\s*//;
    $descript=~s/\s*$//;
    splice(@synonyms,$i,1); 

    for $synonym (@synonyms) {$descript.=", ".$synonym;}

    if ($descript ne "") {
	$descript=~s/\s+/ /g;
	$descript=~s/^\s+//g;
	$descript=~s/\s*;*\s*$//; # remove ';'
    }
    $descript.=";";   # append a semicolon ';'

    if ($chain ne "") {
	$chain=~s/\s*;\s*$//;
	@chains=split(/[,; ]\s*/,$chain);
	foreach $chain (@chains) {
	    $descript{$chain}=$descript;
#	    printf("chain='$chain'   description='$descript'\n");
	} 
    } else {
	$descript{$chain}=$descript;
#	printf("chain='$chain'   description='$descript'\n");
    }
}


##################################################################################
# Set organism when a new MOL_ID line is found, or at the end of SOURCE records
##################################################################################
sub SetOrganism() {
    if ($organism eq "" && $organism_common ne "") {$organism=$organism_common;} 
    if (!exists $chain[$molid]) {$molid=1-$molid;}
    $chain=$chain[$molid];
    $molid=$1;
    $organism=~tr/$//d;
    $organism=~s/^\s*([^;:]*).*/$1/;
    if ($organism=~/^\S*\s*\(([\w ]*)/) {$organism=$1;} # bovine (Bos taurus)
    elsif ($organism=~/^([\w -]+)/)     {$organism=$1;} # BACTERIOPHAGE T4 (MUTANT GENE DERIVED ...)
    elsif ($organism=~/^(\S+\s+\S+)/)   {$organism=$1;} # maximum two words
    elsif ($organism=~/^(\S+)/)         {$organism=$1;} 
    $organism=~s/\s*$//g;
    if ($chain ne "") {
	$chain=~s/\s*;\s*$//;
	@chains=split(/[,; ]\s*/,$chain);
	foreach $chain (@chains) {
	    $organism{$chain}=$organism;
#	    printf("chain='%s'  organism{chain}='%s'\n",$chain,$organism{$chain});
	}
    } else {
	$organism{$chain}=$organism;
#	printf("chain='%s'  organism{chain}='%s'\n",$chain,$organism{$chain});
    }
}    



##################################################################################
# Print out sequence of last chain read in 
##################################################################################
sub AddSequence() {
    my $seqres="";
    my $pdbidchain;
    my $nc;   # $nc= either next chain number OR, if identical seq exists with better resolution, index of this seq

    foreach my $aa (@seqres) {$seqres.=&Three2OneLetter($aa);}
    if ($v>=3) {
	printf("CHAIN ='%s'\n",$chain);
	printf("DESCRP='%s'\n",$descript{$chain});
	printf("KEYWDS='$keywds'\n");
	printf("ORGANI='%s'\n",$organism{$chain});
	printf("SEQRES='%s'\n",$seqres);
    }
    if (length($seqres)<=20) {return 1;} # skip short protein/DNA chains (for DNA's ADGT &Three2OneLetter() returns "")

    if ($chain ne "") {$pdbidchain=$pdbid."_".$chain;} else {$pdbidchain=$pdbid;}

    # Check for nonredundancy
    if ($nr==1 && defined $nchains{$seqres} ) {
	$nc=$nchains{$seqres}; 
#	$sequences[$nchains{$seqres}]=~/^(\S+)/;
#	print("Sequence $pdbid"."_$chain redundant with $1. Res now: $resolution  Res before: $resolution[$nc]\n");
	if ($resolution==0 || $resolution[$nc]<=$resolution) {
	    # PDB identifier not yet contained in list of equivalent pdb ids? 
	    if ($equiv_pdbs[$nc]!~/$pdbid/ && $sequences[$nc]!~/^>$pdbid/) {
		$equiv_pdbs[$nc].=" $pdbidchain"; # Add new pdbid_chain to list $equiv_pdbs[$nc]
		if ($het ne "") {$equiv_pdbs[$nc].="*";}
	    }
	    return 1;
	} else {
	    # Sequence redundant
            # => Throw out earlier sequence and keep this one	    
	    # => Keep list $equiv_pdbs[$nc] from earlier sequence and append its pdbid
	    $sequences[$nc]=~/>(\S+)/;
	    $equiv_pdbs[$nc].=" $1";
	    if ($het ne "") {$equiv_pdbs[$nc].="*";}
	}

    } else {
	$nc=$nchains{$seqres}=$nchains;
	$nchains++;
	$equiv_pdbs[$nc]="";
    }
    $resolution[$nc]=$resolution;

    # If descript{chain} does not exist, it was not specified seperately for each chain
    if (exists $descript{$chain}) {$descript=$descript{$chain}}		    
    $descript=~s/;*$//;       # remove ; at the end
   
    # If organism{chain} does not exist, it was not specified seperately for each chain
    if (exists $organism{$chain}) {
	$organism=lc($organism{$chain});
    } else {
	if($organism=~/\((.*)\)/) {$organism=$1;}
	$organism=lc($organism);
    }
    
    if ($v>=3) {
	printf("Accept:\n",$chain);
	printf("CHAIN ='%s'\n",$chain);
	printf("DESCRP='%s'\n",$descript);
	printf("KEYWDS='$keywds'\n");
	printf("ORGANI='%s'\n",$organism);
	printf("SEQRES='%s'\n",$seqres);
    }
    
    # Correct upper/lower case
    $descript=" $descript ";
    $descript=~s/([a-zA-Z]{5,})/\L$1/g; # make everything longer than 5 word characters lower case
    $descript=~s/([\s]OF[\s])/\L$1/g;
    $descript=~s/([\s]OR[\s])/\L$1/g;
    $descript=~s/([\s]ON[\s])/\L$1/g;
    $descript=~s/([\s]NO[\s])/\L$1/g;
    $descript=~s/([\s]IN[\s])/\L$1/g;
    $descript=~s/([\s]IS[\s])/\L$1/g;
    $descript=~s/([\s]BY[\s])/\L$1/g;
    $descript=~s/([\s]AT[\s])/\L$1/g;
    $descript=~s/([\s]TO[\s])/\L$1/g;
    $descript=~s/([ -]ALL[ -])/\L$1/g;
    $descript=~s/([ -]AND[ -])/\L$1/g;
    $descript=~s/([ -]ARM[ -])/\L$1/g;
    $descript=~s/([ -]BOX[ -])/\L$1/g;
    $descript=~s/([ -]BOX[ -])/\L$1/g;
    $descript=~s/([ -]EGG[ -])/\L$1/g;
    $descript=~s/([ -]EYE[ -])/\L$1/g;
    $descript=~s/([ -]FOR[ -])/\L$1/g;
    $descript=~s/([ -]HAS[ -])/\L$1/g;
    $descript=~s/([ -]HEN[ -])/\L$1/g;
    $descript=~s/([ -]HOT[ -])/\L$1/g;
    $descript=~s/([ -]LOW[ -])/\L$1/g;
    $descript=~s/([ -]MOL[ -])/\L$1/g;
    $descript=~s/([ -]NON[ -])/\L$1/g;
    $descript=~s/([ -]ONE[ -])/\L$1/g;
    $descript=~s/([\s]OUT[\s])/\L$1/g;
    $descript=~s/([ -]SEX[ -])/\L$1/g;
    $descript=~s/([ -]SIX[ -])/\L$1/g;
    $descript=~s/([ -]TEN[ -])/\L$1/g;
    $descript=~s/([\s]THE[\s])/\L$1/g;
    $descript=~s/([ -]TWO[ -])/\L$1/g;
    $descript=~s/([ -]WAY[ -])/\L$1/g;
    $descript=~s/([\W]ACID[\W])/\L$1/g;
    $descript=~s/([\W]ACYL[\W])/\L$1/g;
    $descript=~s/([\W]ALDO[\W])/\L$1/g;
    $descript=~s/([\W]ANTI[\W])/\L$1/g;
    $descript=~s/([\W]AUTO[\W])/\L$1/g;
    $descript=~s/([\W]AXIN[\W])/\L$1/g;
    $descript=~s/([\W]BASE[\W])/\L$1/g;
    $descript=~s/([\W]BEAN[\W])/\L$1/g;
    $descript=~s/([\W]BETA[\W])/\L$1/g;
    $descript=~s/([\W]BILE[\W])/\L$1/g; 
    $descript=~s/([\W]BLUE[\W])/\L$1/g; 
    $descript=~s/([\W]BONE[\W])/\L$1/g;
    $descript=~s/([\W]BOND[\W])/\L$1/g;
    $descript=~s/([\W]CELL[\W])/\L$1/g;
    $descript=~s/([\W]COAT[\W])/\L$1/g;
    $descript=~s/([\W]COIL[\W])/\L$1/g;
    $descript=~s/([\W]COLD[\W])/\L$1/g;
    $descript=~s/([\W]COLI[\W])/\L$1/g;
    $descript=~s/([\W]CORE[\W])/\L$1/g;
    $descript=~s/([\W]CRYO[\W])/\L$1/g;
    $descript=~s/([\W]DRUG[\W])/\L$1/g;
    $descript=~s/([\W]DUAL[\W])/\L$1/g;
    $descript=~s/([\W]DUCK[\W])/\L$1/g;
    $descript=~s/([\W]ENDO[\W])/\L$1/g;
    $descript=~s/([\W]FAST[\W])/\L$1/g;
    $descript=~s/([\W]FIVE[\W])/\L$1/g;
    $descript=~s/([\W]FOLD[\W])/\L$1/g;
    $descript=~s/([\W]FOOT[\W])/\L$1/g;
    $descript=~s/([\W]FORM[\W])/\L$1/g;
    $descript=~s/([\W]FOUR[\W])/\L$1/g;
    $descript=~s/([\W]FROM[\W])/\L$1/g;
    $descript=~s/([\W]FLAP[\W])/\L$1/g;
    $descript=~s/([\W]FREE[\W])/\L$1/g;
    $descript=~s/([\W]GENE[\W])/\L$1/g;
    $descript=~s/([\W]GOOD[\W])/\L$1/g;
    $descript=~s/([\W]HALF[\W])/\L$1/g;
    $descript=~s/([\W]HAND[\W])/\L$1/g;
    $descript=~s/([\W]HAVE[\W])/\L$1/g;
    $descript=~s/([\W]HEAD[\W])/\L$1/g;
    $descript=~s/([\W]HEAT[\W])/\L$1/g;
    $descript=~s/([\W]HEME[\W])/\L$1/g;
    $descript=~s/([\W]HEXA[\W])/\L$1/g;
    $descript=~s/([\W]HIGH[\W])/\L$1/g;
    $descript=~s/([\W]HOLO[\W])/\L$1/g;
    $descript=~s/([\W]IRON[\W])/\L$1/g;
    $descript=~s/([\W]KETO[\W])/\L$1/g;
    $descript=~s/([\W]KNOT[\W])/\L$1/g;
    $descript=~s/([\W]LATE[\W])/\L$1/g;
    $descript=~s/([\W]LENS[\W])/\L$1/g;
    $descript=~s/([\W]LIKE[\W])/\L$1/g;
    $descript=~s/([\W]LONG[\W])/\L$1/g;
    $descript=~s/([\W]LOOP[\W])/\L$1/g;
    $descript=~s/([\W]MAIN[\W])/\L$1/g;
    $descript=~s/([\W]MEAN[\W])/\L$1/g;
    $descript=~s/([\W]MINI[\W])/\L$1/g;
    $descript=~s/([\W]MITE[\W])/\L$1/g;
    $descript=~s/([\W]MODE[\W])/\L$1/g;
    $descript=~s/([\W]MONO[\W])/\L$1/g;
    $descript=~s/([\W]MUCH[\W])/\L$1/g;
    $descript=~s/([\W]NINE[\W])/\L$1/g;
    $descript=~s/([\W]NULL[\W])/\L$1/g;
    $descript=~s/([\W]ONLY[\W])/\L$1/g; 
    $descript=~s/([\W]OPEN[\W])/\L$1/g; 
    $descript=~s/([\W]PARA[\W])/\L$1/g; 
    $descript=~s/([\W]PLUS[\W])/\L$1/g; 
    $descript=~s/([\W]POST[\W])/\L$1/g; 
    $descript=~s/([\W]POLY[\W])/\L$1/g; 
    $descript=~s/([\W]PORE[\W])/\L$1/g; 
    $descript=~s/([\W]PUMP[\W])/\L$1/g; 
    $descript=~s/([\W]RICH[\W])/\L$1/g; 
    $descript=~s/([\W]RING[\W])/\L$1/g; 
    $descript=~s/([\W]ROLE[\W])/\L$1/g; 
    $descript=~s/([\W]ROLL[\W])/\L$1/g; 
    $descript=~s/([\W]SALT[\W])/\L$1/g; 
    $descript=~s/([\W]SEMI[\W])/\L$1/g; 
    $descript=~s/([\W]SITE[\W])/\L$1/g; 
    $descript=~s/([\W]STEM[\W])/\L$1/g; 
    $descript=~s/([\W]TAIL[\W])/\L$1/g;
    $descript=~s/([\W]TATA[\W])/\L$1/g;
    $descript=~s/([\W]TRAP[\W])/\L$1/g; 
    $descript=~s/([\W]TUBE[\W])/\L$1/g;
    $descript=~s/([\W]TURN[\W])/\L$1/g;
    $descript=~s/([\W]TWIN[\W])/\L$1/g;
    $descript=~s/([\W]TYPE[\W])/\L$1/g;
    $descript=~s/([\W]WELL[\W])/\L$1/g;
    $descript=~s/([\W]WILD[\W])/\L$1/g;
    $descript=~s/([\W]WITH[\W])/\L$1/g;
    $descript=~s/([\W]WROM[\W])/\L$1/g;
    $descript=~s/([\W]ZETA[\W])/\L$1/g;
    $descript=~s/([\W]ZINC[\W])/\L$1/g;
    $descript=~s/DE NOVO/de novo/ig;
    $descript=~s/(\W)KDA(\W)/$1kDa$2/g;

    $descript=~s/(\S+[CBDFGJKLMNPQRTVWXZ]{4,}\S+)/\U$1/ig;   
    $descript=~s/(\W)(\S[CBDFGHJKLMNPQRSTVWXZ]{4,}\W)/$1\U$2/ig; # no vowels for at least 4 letters -> abbreviation -> upper case
    $descript=~s/(\W)([CBDFGHJKLMNPQRSTVWXZ]{4,}\S\W)/$1\U$2/ig; # no vowels for at least 4 letters -> abbreviation -> upper case
    $descript=~s/(\w+ii+\w+)/\U$1/ig;
    $descript=~s/([\W]rossman[\W])/\u$1/g;
    $descript=~s/([\W]nadph[\W])/\U$1/g;
    $descript=~s/([\W]gapdh[\W])/\U$1/g;
    $descript=~s/([\W]f[\W])/\u$1/g;
    $descript=~s/(\W)(\w{0,3})RNP(\W)/$1\L$2\URNP$3/ig;
    $descript=~s/(\W)(\w{0,3})RNA(\W)/$1\L$2\URNA$3/ig;
    $descript=~s/(\W)(\w{0,3})DNA(\W)/$1\L$2\UDNA$3/ig;
    $descript=~s/RNASE(\W)/RNAse$1/ig;
    $descript=~s/DNASE(\W)/DNAse$1/ig;
    $descript=~s/barnase/barnase/ig;
    $descript=~s/atpase(\W)/ATPase$1/ig;
    $descript=~s/gtpase(\W)/GTPase$1/ig;
        
    # Write amino acid three letter symbols with one capital letter
    foreach my $aa ("Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn","Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr","MSe") {
	$descript=~s/$aa([ -;:.+])/$aa$1/ig;
    }    
    # Write ions as chemical elements 
    foreach my $ion ("Zn","Mg","Na","Ka","Ca","Fe","Cu","Se","Al","Mn") {
	$descript=~s/([ -:;])$ion([ -;:.+])/$1$ion$2/ig;
	$descript=~s/([ -:;])$ion([ -;:.+])/$1$ion$2/ig;
    }    
    
    $descript=~s/^\s*//;
    $descript=~s/\s*$//;
    $descript="\u$descript";  # first letter of description upper case
    $organism="\u$organism";  # first letter of organism upper case
    $organism=~s/ ([\w\d]{0,2})$/ \U$1/;  # Influenza A virus
    $organism=~s/ (\w+\d+\w*)$/ \U$1/;  # Influenza A virus
    if ($v>=2 && $organism eq "") {print("\n\nWarning: no organism found for chain $chain in $pdbfile\n");}

    my $foldid="";
    if ($dalifamids{$pdbidchain}) {$foldid=" $dalifamids{$pdbidchain}".$foldid;}
    if ($scopfamids{$pdbidchain}) {$foldid=" $scopfamids{$pdbidchain}".$foldid;}

    my $res;
    if ($resolution>0) {$res=$resolution."A";} else {$res="NMR";}
    # Set sequence record
    $seqres=~s/(\S{1,100})/$1\n/g; # insert newlines after each 70 characters
    $sequences[$nc]=sprintf(">%-6.6s %s;%s %s {%s}%s\n",$pdbidchain,$descript,$het,$res,$organism,$foldid);
    $sequences[$nc].=sprintf("$seqres");
    if($v>=3) {print($sequences[$nc]);}
    return 0;
}

sub ReadScopFile() 
{
    # Read dir.cla.scop.txt_1.65
    if ($v>=2) {print("Reading $scopfile ... ");}
    open (SCOPFILE,"<$scopfile") || die ("ERROR: cannot open $scopfile: $!\n");
#d1dlwa_ 1dlw    A:      a.1.1.1 14982   cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=46461,px=14982

    my $n=0;
    my $scopfamid;
    my $chain;
    my $pdbidchain;
    while($line=<SCOPFILE>) {
	if ($line=~/^\S+\s+(\S\S\S\S)\s+(\S+)\s+(\w\.\d+\.\d+\.\d+)/) {
	    $pdbid=$1;
	    $chain=$2;
	    $scopfamid=$3;
	    $pdbidchain=$1;
#	    printf("chain=$chain\n");
	    if ($chain!~/([A-Za-z\d]):/) {
		    if ($scopfamids{$pdbidchain}) {
			$scopfamids{$pdbidchain}.=" ".$scopfamid;
#			printf("scopfamids($pdbidchain)=%s\n",$scopfamids{$pdbidchain});
		    } else {
			$scopfamids{$pdbidchain}="SCOP: ".$scopfamid;
#			printf("scopfamids($pdbidchain)=%s\n",$scopfamids{$pdbidchain});
		    }
		    $n++;
	    } else {
		my %chains=();
		while ($chain=~/([A-Za-z\d]):/) {
		    $chain=~s/([A-Za-z\d]):[^A-Za-z:]*//;
		    if ($chains{$1}) {next;}
		    if ($1 ne "") {$pdbidchain="$pdbid"."_$1";}
		    $chains{$1}=1;
		    if ($scopfamids{$pdbidchain}) {
			$scopfamids{$pdbidchain}.=" ".$scopfamid;
#			printf("scopfamids($pdbidchain)=%s\n",$scopfamids{$pdbidchain});
		    } else {
			$scopfamids{$pdbidchain}="SCOP: ".$scopfamid;
#			printf("scopfamids($pdbidchain)=%s\n",$scopfamids{$pdbidchain});
		    }
		    $n++;
		} 
	    }
	}
    }
    close(SCOPFILE);
    print(" found $n domains in SCOP file $scopfile\n");
}


sub ReadDaliFiles() 
{
    # Read FoldIndex.html
    if ($v>=2) {print("Reading $dalidir/FoldIndex.html ... ");}
    open (FOLDINEXFILE,"<$dalidir/FoldIndex.html") || die ("ERROR: cannot open $dalidir/FoldIndex.html: $!\n");
    # 1.1.1.1.1.1     1cunA_2 ...
    # 1.1.1.2.1.1     ___1lvfA_1 ...
    my $n=0;
    my $dalifamid;
    my %fold_for_repr;
    while($line=<FOLDINEXFILE>) {
	if ($line=~/(\d+\.\d+\.\d+\.\d+\.\d+\.\d+)\s+_*(\S+)(_\d+)/) {
	    $fold_for_repr{$2.$3}=$1;
	    if ($3 eq "_0") {
		$dalifamid=$1;
		$pdbid=$2;
		$pdbid=~s/^(\S\S\S\S)(\S)/$1_$2/;
		$dalifamids{$pdbid}="DALI: ".$dalifamid;	    
	    }
	    $n++;
	}
    }
    close(FOLDINEXFILE);
    if ($v>=1) {print(" found $n representative domains in DALI's FoldIndex.html\n");}

    # Read domain_definitions.txt
    my $domainfile="$dalidir/domain_definitions.txt";
    my $repr;
    if ($v>=2) {print("Reading $domainfile ... ");}
    open (DOMAINFILE,"<$domainfile") || die ("ERROR: cannot open $domainfile: $!\n");
    # 1cunA/1-106	1cunA_1	1	ALPHA SPECTRIN
    $n=0;
    while($line=<DOMAINFILE>) {
	if ($line=~/^(\S+)\/\S+\s+(\S+)/) {
	    if (!$fold_for_repr{$2}) {
		if ($v>=2) {print("WARNING: no fold for DALI representative $2 in $domainfile\n");}
		next;
	    }
	    $pdbid=$1;
	    $repr=$2;
	    $pdbid=~s/^(\S\S\S\S)(\S)/$1_$2/;
	    if ($dalifamids{$pdbid}) {
		$dalifamids{$pdbid}.=" ".$fold_for_repr{$repr};
	    } else {
		$dalifamids{$pdbid}=" DALI: ".$fold_for_repr{$repr};
	    }
	    $n++;
	}
    }
    close(DOMAINFILE);
    print(" found $n domains in DALI file $domainfile\n");
}


##################################################################################
# Convert three-letter amino acid code into one-letter code
##################################################################################
sub Three2OneLetter {
    my $res = $three2one{uc($_[0])};
    if (defined $res) {
	return $res;
    } elsif ($_[0] =~ /^\s+$/) {
	return "";
    } else {
	return "X";
    }
}
