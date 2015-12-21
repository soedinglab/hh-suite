
## read hhsearch results file (containing all filtered hits ranked by a predicted TM-Score)
## create pir alignment for MODELLER software
## if necessary create artificial pdb files (add 50 Angstrom)	

package checkTemplates;
require Exporter;


use strict;
use config;
use utilities;

our @ISA = qw(Exporter);
our @EXPORT = qw(CheckTemplates);

## THIS IS (not my own) UGLY CODE AND SHOULD BE REWRITTEN

## requires valid path to hhmakemodel

################################################################
## description: read selected hits in hhr-file
## check for overlap and repeat domains!
## input: hhrFile - the aritificial hhr file produced by 
##           CreateHHRFromPrediction
## a3mQuery: the a3m query alignment file-path
## pirfile: basename of pirfile; pir-file will be created now
##             by calling hhmakemodel, see below
################################################################
sub CheckTemplates {
    my $outbase = shift;
    my $workingDir = shift;
    my $config = shift;

    my $hhrFile = "$outbase.hhr";
    my $a3mQuery = "$outbase.a3m";
    my $pirFile = "$outbase.pir";

    my $pdbdir = $config->get_pdbdir();
    my $hhmakemodel = $config->get_hhmakemodel();

    my @hits_check;
    my $hit;
    my @hits;
    my $ID;
    my $qstart;
    my $qend;
    my $tstart;
    my $tend;

    my $verbose = 1;

    ## verbose mode for hh tools
    my $v = 1;

    open(HHR, "< $hhrFile") or die("Error checkTemplates: Cant open $hhrFile: $!\n");

    while (my $line = <HHR>) {
	if ($line =~ /^\s*(\d+)\s+(\S+).+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)-(\d+)\s+(\d+)-(\d+)\s*\((\S+)\)$/) {	
	    $hit = $1;
	    push (@hits, $hit);

	    $ID = $2;
	    $qstart = $9;
	    $qend = $10;
	    $tstart = $11;
	    $tend = $12;      

	    if ($ID =~ /^[defgh](\d[a-z0-9]{3})([a-z0-9_\.])[a-z0-9_]/) {
		$ID = $1;
		my $chain;

		if ($2 eq "_") {
		    $chain = "[A ]";
		} 
		else {
		    $chain = uc($2);
		}
		$ID .= "_$chain";
	    } 
	    push(@hits_check, [$hit, $ID, $qstart, $qend, $tstart, $tend]);      
	}  
	
    }
    close (HHR);


    if ($verbose >= 2) {
	foreach my $c (@hits_check){
	    foreach my $d (@{$c}){
		print "$d\t";
	    } 
	    print "\n";
	}
	print "\n";
    }

    #########################################
    ## check if same template is used several 
    ## times and in overlapping parts (e.g repeatproteins!!!) 
    ## -> because of problems with MODELLER restraints!
    #########################################
    my @aligned;
    my $donttake = "";
    my $added    = "";
    my $doubles  = 0;

    for (my $h=0; $h<@hits_check-1; $h++) {
	for (my $i=$h+1; $i<@hits_check; $i++) {
	    if (($donttake !~ /,$hits_check[$h][0],/) && ($donttake !~ /,$hits_check[$i][0],/)) {

		my $queryov=0;
		my $templov=0;			  
		#check: same templates?
		if ($hits_check[$h][1] eq $hits_check[$i][1]){
		    print "Hit: $hits_check[$h][0] & $hits_check[$i][0] same templates!\t";
		    #check: overlap?
		    #print "$hits_check[$h][2]&$hits_check[$i][2]\n$hits_check[$h][3]&$hits_check[$i][3]\n$hits_check[$h][4]&$hits_check[$i][4]\n$hits_check[$h][5]&$hits_check[$i][5]\n\n";
		    
		    #check: overlap in query!?!	
		    for (my $i=1; $i<=5000; $i++) {
			$aligned[$i]=0;
		    }

		    my $firstq=$hits_check[$h][2];
		    my $lastq=$hits_check[$h][3];						
		    for (my $i=$firstq; $i<=$lastq; $i++) {
			$aligned[$i]=1
		    }

		    $firstq=$hits_check[$i][2];
		    $lastq=$hits_check[$i][3];
		    my $queryov=0;
		    my $mlenq=0;

		    for (my $i=$firstq; $i<=$lastq; $i++) {
			if ($aligned[$i]==1) {
			    $queryov++;
			} 
			else {
			    $mlenq++;
			}
		    }
		    
		    #check: overlap in template!?!			
		    for (my $i=1; $i<=5000; $i++) { $aligned[$i]=0; }
		    my $firstt = $hits_check[$h][4];
		    my $lastt = $hits_check[$h][5];						

		    for (my $i=$firstt; $i<=$lastt; $i++) { $aligned[$i]=1; }
		    $firstt = $hits_check[$i][4];
		    $lastt = $hits_check[$i][5];

		    my $templov=0;
		    my $mlent=0;
		    for (my $i=$firstt; $i<=$lastt; $i++) {
			if ($aligned[$i]==1) {
			    $templov++;
			} 
			else {
			    $mlent++;
			}
		    }
		    
		    print "Overlap: in Query:  $queryov AA	in Template:  $templov AA\t";
		    
		    if ($queryov<=20 && $templov>20){			
			print "repeat-domain: $hits_check[$h][0] & $hits_check[$i][0]\t";		
			if (($added eq "") || ($added !~ /$hits_check[$h][1]/)){
			    $doubles++;
			    print " -> create $workingDir/$hits_check[$h][1]_$h.pdb\n";
			    system ("cp $pdbdir/$hits_check[$h][1].pdb	$workingDir/$hits_check[$h][1]_$h.pdb\n");
			    #add 50 Angstrom to x coordinate in pdb file:
			    open (PDB, "< $workingDir/$hits_check[$h][1]_$h.pdb") or die ("Error: cannot open $workingDir/$hits_check[$h][1]_$h.pdb\n");
			    my @lines = <PDB>;
			    close(PDB);    				 
			    for (my $i=0; $i<@lines; $i++) {
				if ($lines[$i]=~ /^(.{30})(\s*)(-?\d+.\d+)(\s*-?\d+.\d+\s*-?\d+.\d+\s*.*)/){
				    my $xcoord = $3+50;		   
				    $lines[$i] = sprintf("$1%8.3f$4\n", $xcoord);		
				}						
			    }						   
			    open (PDBnew, "> $workingDir/$hits_check[$h][1]_$h.pdb") or die ("Error: cannot open $workingDir/$hits_check[$h][1]_$h.pdb\n");
			    print(PDBnew @lines);	
			    close (PDBnew);	
			    $added .= " $hits_check[$h][1]_$h";	
			}
			elsif($added !~ /$hits_check[$h][1]_$h/){
			    $doubles++;
			    print " -> create $workingDir/$hits_check[$h][1]_$h.pdb\n";
			    my $copy;				
			    if ($added =~ /.*($hits_check[$h][1]_\d*).*?$/) { 
				$copy="$workingDir/$1.pdb";
			    }     
			    system ("cp $copy $workingDir/$hits_check[$h][1]_$h.pdb\n");
			    #add 50 Angstrom to x coordinate in pdb file:
			    open (PDB, "< $workingDir/$hits_check[$h][1]_$h.pdb") or die ("Error: cannot open $workingDir/$hits_check[$h][1]_$h.pdb\n");
			    my @lines = <PDB>;
			    close(PDB);    				 
			    for (my $i=0; $i<@lines; $i++) {
				if ($lines[$i]=~ /^(.{30})(\s*)(-?\d+.\d+)(\s*-?\d+.\d+\s*-?\d+.\d+\s*.*)/) {
				    my $xcoord = $3+50;		    
				    $lines[$i] = sprintf("$1%8.3f$4\n", $xcoord);   
				}						
			    }						   
			    open (PDBnew, "> $workingDir/$hits_check[$h][1]_$h.pdb") or die ("Error: cannot open $workingDir/$hits_check[$h][1]_$h.pdb\n");
			    print(PDBnew @lines);	
			    close (PDBnew);	
			    $added .= " $hits_check[$h][1]_$h";
			}
			else{print "\n";}																    
		    }
		    elsif ($queryov>20){
			print "ignore HIT $hits_check[$i][0] for modeling\n";						
			$donttake .= ",$hits_check[$i][0],";		      
		    }	
		    else{print "\n"}
		}
	    }
	}
    }
    my @hits_take;
    for (my $h=0; $h<@hits; $h++){
	if ($donttake !~ /,$hits[$h],/){
	    push(@hits_take,$hits[$h]);
	}
    }

    &System("$hhmakemodel $hhrFile -m @hits_take -q $a3mQuery -v $v -d $pdbdir -pir $outbase.pir"); 

    if ($doubles > 0) {
	my @pdbs = <$workingDir/*.pdb>;
	my $checkedPDBs = ""; 
	my $duplis = 1;	  

	open (PIR, "< $outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
	my @lines = <PIR>;
	close(PIR);

	for (my $i=0; $i<@lines; $i++) {    
	    if ($lines[$i] =~ /^structureX:/) {
		$lines[$i-1] =~ />P1;(\S+)/;
		my $pdbid = $1;

		for (my $j=0; $j<@pdbs; $j++) { 
		    if (($pdbs[$j] =~ /$workingDir\/$pdbid\_\d+\.pdb/) && ($checkedPDBs !~ /$pdbid/)) {
			$checkedPDBs .= "$pdbid ";
			last;
		    }
		    elsif (($pdbs[$j] =~ /$workingDir\/$pdbid\_\d+\.pdb/)) {
			print"\nReplace $pdbid by DUPLICATEPDB$duplis.$pdbid.pdb!\n";
			system "mv $pdbs[$j] $workingDir/DUPLICATEPDB$duplis.$pdbid.pdb";
			$lines[$i-1] =~ s/>P1;$pdbid/>P1;DUPLICATEPDB$duplis.$pdbid/;
			$lines[$i] =~ s/^structureX:$pdbid:(.*)/structureX:DUPLICATEPDB$duplis.$pdbid:$1/;
			$duplis++;
			@pdbs = grep !/$pdbs[$j]/, @pdbs;
			last;
		    }
		}
	    }	
	}

	open (OUT, "> $outbase.pir") or die ("Error: cannot write to pir-file: $!\n");
	print(OUT @lines);
	close(OUT);      
    }
}

