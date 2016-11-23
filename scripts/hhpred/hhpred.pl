#!/usr/bin/perl -w
use strict;
use warnings;

## general comments:
## this script is indirectly called by run_casp.pl, see comments therein
## it calls several functions and mainly functions as a "gluing" script
## its behavious is mainly controled by the config file read at the beginning
## in general it does the following:
## 1) build MSA from sequence via hhblits
## 2) build hhr-template-list via hhsearch
## 3) (pre-)select several templates via template selection strategy
## 4) filter preselected templates and query
## 5) rank templates with single template NN
## 6) select final templates 
## 7) generate final (artifical) hhr file with selected templates
## 8) optionally replace distance restraints
## 9) call MODELLER and build final model
## this final model in saved at outdir and will then be sent to 
## the CASP organizers via the run_casp.pl script


## more details:
## all temporary files are written to workingDir
## temporary files are eg. *.a3m, *.hhm, *.tab, *.hhr
## amd are named workingDir/queryName.*
## final results must be copied back to "dirbasename", e.g. outbase.pdb to dirbasename.pdb
##
## if filtering is on the filtered files are put in another temporary directory,
## see filterAlignments.pm
## the final pdb-file must be in outbase.pdb so that run_casp.pl can find this file and
## attach it to the reply-e-mail
## 

use File::Basename;
use File::Spec;


BEGIN {
    my $dirname = dirname(File::Spec->rel2abs(__FILE__));
    unshift @INC, $dirname . "/lib";
}

use File::Temp;

use config;
use utilities;
use Template;
use TemplateList;
use TemplateStruct;
use TemplateListStruct;
use selectTemplatesHeuristic;
use filterAlignments;
use predTMscore;
use checkTemplates;
use modeller;
use QualityAssessModel;

## called via run_casp.pl with these parameters
my $server      = "hhpred";
my $outFile     = "";
my $configFile  = undef;
my $queryEnding = "";
my $queryName   = "";
my $queryFile   = "";

#############################
## parse command line options

my $options     = "";

for (my $i=0; $i<@ARGV; $i++) {
    $options .= " $ARGV[$i] ";
}

if ($options =~ s/-i\s+(\S+)/ /) { 
    $queryFile = $1; 
    my $fpath;
    ($queryName, $fpath, $queryEnding) = fileparse("$queryFile", qr/\.[^.]*/);
}
if ($options =~ s/-o\s+(\S+)/ /) { $outFile = $1; }
if ($options =~ s/-config\s+(\S+)/ /) { $configFile = $1; }

##############################

print "\n";
print "==========================================================\n";
print "|                HHPRED structure predictor              |\n";
print "==========================================================\n";
print "\n";

## set default values

#my $queryFile = "$dirbasename" . "$queryEnding";

my $workingDir = "/tmp/$$"; 

## set configuration
my $config = HHpredConfig->instance($configFile);

my $preselectedTemplatesFile = "";
my $allRankedTemplates = TemplateList->new();

if ($queryFile eq "" or $outFile eq "" or $queryName eq "" or ! -e $queryFile) {
    print "usage: $0 -i <infile> -o <outfile> [-c <configfile>]!\n";
    print "   <infile>  .a3m or .seq file\n";
    print "   <outfile> resulting pdb file\n\n";
    exit 1;
}

## create working directory
if (! -e "$workingDir") {
    &System("mkdir -p $workingDir");
} else {
    my $dir = File::Temp->newdir("XXXXX", DIR => "$workingDir", CLEANUP => 0);
    $workingDir = $dir->dirname;
    &System("chmod 777 $workingDir");
}

if (! -e "$workingDir/$queryName") {
    &System("mkdir -p $workingDir/$queryName");
    $workingDir = "$workingDir/$queryName";
}

#############################
$config->print();
print "\n";

#############################
## set up starting files
## all files are written to workingDir
$queryName .= &getRandomString(7); # add random stuff so that queryName != template (needed for MODELLER)
my $outbase = "$workingDir/$queryName";

## either a3m or seq file
if ($queryEnding ne ".a3m") {
    &System("cp $queryFile $outbase.seq");
} else {
   &System("cp $queryFile $outbase.a3m");
}


#############################
## build query alignment if necessary
if (! -e "$outbase.a3m") {
    my $command = "HHLIB=" . $config->get_hhlib() 
                . " " 
                . $config->get_hhblits() 
                . " -i $outbase.seq"
                . " -oalis $outbase" 
                . " -ohhm $outbase.hhm"
                . " -n " . $config->get_hhblits_rounds() 
                . " -mact " . $config->get_hhblits_mact() 
                . " -oa3m $outbase.a3m"
                . " -d " . $config->get_uniprot20()
                . " -cpu " . $config->get_cpus();
    &System($command);
    sleep(1);

    &System($config->get_addss() . " $outbase.a3m");
    sleep(1);
}

&System($config->get_hhmake() . " -i $outbase.a3m -o $outbase.hhm");
sleep(1);

##############################
## search against database via hhsearch
my $pdbdir = $config->get_pdbdir();

$pdbdir =~ s/
    \/$ # trim trailing slash
    //x;
my $pdbdb = "$pdbdir/db/pdb.hhm";

&System("HHLIB=" . $config->get_hhlib() 
        . " "
        . $config->get_hhsearch()
        . " -i $outbase.hhm"
        . " -d $pdbdb"
        . " -mact " . $config->get_hhsearch_mact()
        . " -cpu " . $config->get_cpus()
        . " -atab $outbase.start.tab");  

while (!(-e "$outbase.hhr")) {
    sleep(1)
}

&System("cp $outbase.hhr $outbase.start.hhr");

##############################
## start handling of templates
my $tList = TemplateListStruct->new();
$tList->hhr_to_TemplateList("$outbase.hhr");
my $queryLength = $tList->get_queryLength();
$tList->print();

#############################
## preselect template (best according to similarity, sumprobs, probability
## and fill up rest with heuristic
if ($config->get_preselectTemplates()) {
    $tList = &ChooseTemplatesScoringHeuristic($workingDir, $queryLength, $outbase, 100, 1, $tList, $config);

    print "preselected templates:\n";
    $tList->print();
    print "====================================\n\n";
}

#############################
## filtering
if ($config->get_doFiltering()) {
    $tList =  &Filtering($queryName, $outbase, $tList, $server, $config);

    print "Filtered templates:\n";
    $tList->print();
    print "====================================\n\n";
}

#############################
## rank templates with NN
if ($config->get_rankTemplates()) {
    my $rankingNN = TMscoreNN->new();
    $rankingNN->rank_templates($tList, $queryLength, $config);
    
    $allRankedTemplates = $tList;

    print "TM score ranked templates\n";
    $tList->print();
    print "====================================\n\n";
}

#############################
## final template selection
if ($config->get_multiTemplate()) {
    my $maxNumTemplates = $config->get_maxNumOfTemplates();
    $tList = &ChooseTemplatesScoringHeuristic($workingDir, $queryLength, $outbase, $maxNumTemplates, 2, $tList, $config);
} else {
    $tList = &SingleTemplateSelection($tList);
}

## take same template(s) as in file (instead of the ones selected before);
## keep previous template(s) if preselected ones are not in template list
if ($preselectedTemplatesFile ne "") {
    my $preselectedTemplates = TemplateList->new();
    $preselectedTemplates->read_from_file($preselectedTemplatesFile);
    
    my $newTemplates = TemplateList->new();
    for (my $i=0; $i<$preselectedTemplates->size(); $i++) {
        for (my $j=0; $j<$allRankedTemplates->size(); $j++) {
            my $preTemplate = $preselectedTemplates->get($i);
            my $rankTemplate = $allRankedTemplates->get($j);

            if ($preTemplate->get_Hit() eq $rankTemplate->get_Hit()) {
                ## check overlap
                my $overlap = &min($preTemplate->get_Qend(), $rankTemplate->get_Qend()) - &max($preTemplate->get_Qstart(), $rankTemplate->get_Qstart());
                my $preLen = $preTemplate->get_Qend() - $preTemplate->get_Qstart();
                if ($overlap / $preLen > 0.5) {
                    $newTemplates->check_and_add($rankTemplate);
                    last;
                }
            }
        }
    }

    ## templates found
    if ($newTemplates->size() > 0) {
        print "preselected templates could be found\n";
        $tList = $newTemplates;
    } else {
        ## keep old templates
        print "preselected templates could NOT be found - keep old templates\n";
    }
}


print "final templates:\n";
$tList->print();
print "====================================\n\n";

$tList->write_to_file("$outbase.templates");

#############################
## create "artificial" hhr file for input, needed by hhmakemodel to build pir alignment
$tList->templateList_to_hhr($outbase);

#############################
## generate pir alignment
## the hhr file is generated/overwritten by the function call before
## the a3m is generated above using hhblast
## CheckTemplates possibly creates "corrected" pdb files and saves them
## in the "working" directory; then calls hhmakemodel to build outbase.pir
print "checking templates\n";
&CheckTemplates($outbase, $workingDir, $config);
print "\n====================================\n\n";

# if ($config->get_realignProbcons() && $tList->size() > 1) {
#     my $pirFile = PirFile->new("$outbase.pir");
#     my $fastaFile = $pirFile->to_fasta();
#     $fastaFile->write_to_file("$outbase.fas");
#     my $probcons = $config->get_probcons();

#     &System("$probcons $outbase.fas > $outbase.fas");
#     sleep(1);

#     my $fastaRealigned = FastaFile->new("$outbase.fas");
#     my $pirRealigned = $fastaRealigned->to_pir();
#     $pirRealigned->write_to_file("$outbase.pir");
# }

#############################
## Modeller initialization
print "Modeller initialization\n";
&ModellerRSR(
    queryName => $queryName,
    workingDir => $workingDir,
    outbase => $outbase,
    config => $config
    );

#############################
## replace distance restraints
## and call Modeller
if ($config->get_replaceDistanceRestraints()) {
    my $normalizer = &ChangeDistanceRestraints($outbase, $workingDir, $tList, $config);
    &ModellerNewDistance(
     rsrFile => "$outbase.new.rsr",
     queryName => $queryName,
     workingDir => $workingDir,
     outbase => $outbase,
     config => $config,
     normalizer => $normalizer);
} else {
     &ModellerRaw(
      rsrFile => "$outbase.rsr",
      queryName => $queryName,
      workingDir => $workingDir,
      outbase => $outbase,
      config => $config
     );
}

&WriteCASPpdbFile("$outbase.pdb", "$outbase.pdb", $queryName, $server, $tList);

if ($config->get_assessModel()) {
    $tList->set_queryLength($queryLength);
    my $assessor = QualityAssessModel->new(tList=>$tList, outbase=>$outbase, pdbFile=>"$outbase.pdb");
    $assessor->run();
}

#############################
## clean up
my @tabs = <$outbase*.tab>;
for (my $i=0; $i<@tabs; $i++) {
    system("rm $tabs[$i]");
}

#############################
## copy files of interest back to outdir , e.g. to the tempdir
&System("cp $outbase.pdb $outFile");

#&System("rm -r $workingDir");

## create a directory containing intermediate files
#my $resultDir = $dirbasename . "Results";
#&System("mkdir -p $resultDir");
#&System("chmod 777 $resultDir");
#&System("cp $outbase.* $resultDir");

print "END HHpred for $queryName\n";
print "=========================\n";
