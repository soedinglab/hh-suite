#!/usr/bin/perl -w

package filterAlignments;
require Exporter;

use strict;
use config;
use utilities;

use File::Temp qw/ tempfile tempdir /;


our @ISA = qw(Exporter);
our @EXPORT = qw(Filtering FilteringRanking RemoveSSFromA3m);

## TODO: what about scopdir??, see below


sub RemoveSSFromA3m {
    my $a3mFile = shift;

    open (A3M, "< $a3mFile") or die "Cant open $a3mFile: $!\n";
    my @a3m = <A3M>;
    close(A3M);

    ## assume sequence ID starts with a number
    ## remove all lines before ID
    while($a3m[0] !~ /^>\d/) {
    shift(@a3m);
    }

    open (A3M, "> $a3mFile") or die "Cant open $a3mFile: $!\n";
    print (A3M @a3m);
    close(A3M);
}

#######################################################
## input: inbase, i.e. basename of hhr, a3m and hhm
##        outbase, i.e. where to write output for query
##                 (filtered query a3m, hhm); also base
##                 for other output files (filtered a3m
##                 templates)                                                        
## output: list of templates with parameters
##         obtained by the following procedure
##
## description:
## first create a template info array
## take the first HITS templates and use them as
## "database" for the filter steps
## do filter-steps: for each filtering strength:
## create new filtered a3m files (templates and query)
## => create hhm files
## => do hhsearch => hhr file => CreateInfoArray
## append the templates from CreateInfoArray to
## list of templates, update best templates for
## each position in query
## do filtering until either last filter step is
## reached or filtering gives no more templates
## 
#######################################################
sub Filtering {
    my $outbase = shift;        ## basename of outputs wrt query during filtering (e.g. filt.a3m), remove these outputs later
    my $templateList = shift;
    my $config = shift;
    my $qscStep = shift || 0.1;

    my $outdir;                 ## directory where to write output during filtering, e.g. *.filt.a3m files

    ## setup
    my $verbose = 1;
    my $v = 1;

    if ($outbase =~ /^(.*?)\/([^\/]*)$/) {
    $outdir = $1; 
    } else {
    $outdir = "."; 
    }

    my $cpus = $config->get_cpus();
    my $options = "-mact " . $config->get_hhsearch_mact();

    my $qscmax = 0;     #-qsc (max similarity)

    ## TODO: remove limitation, all preselected hits are filtered
    my $HITS = 50;      #number of hits to be filtered
    my %bestscores;     #hash with raw scores for each hit (in order to update SCOREs correctly)    
    my $MAXOVLAP = 20;  #maximum allowable overlap to treat same templatehits independently

    my $hhmake   = $config->get_hhmake();
    my $hhfilter = $config->get_hhfilter();
    my $pdbdir   = $config->get_pdbdir();
          
    ## fill bestscores hash
    for (my $i=0; $i<$templateList->size(); $i++){
    my $template = $templateList->get($i);
    #new template:
    if (! $bestscores{$template->get_Hit()}){
        push (@{$bestscores{$template->get_Hit()}}, [$template->get_Qstart(),
                             $template->get_Qend(),
                             $template->get_Score(),
                             $template->get_SumProbL()]);
    }
    #template exists already: check if independent one      
    else {     
        if ($bestscores{$template->get_Hit()}[0][0] < $template->get_Qstart()) {
        #existing template "in front of" new one
        if (($bestscores{$template->get_Hit()}[0][1]-$template->get_Qstart())<=$MAXOVLAP) {
            #overlap <= 20 residues
            #push to key -> additional array
            push (@{$bestscores{$template->get_Hit()}}, [$template->get_Qstart(),
                                 $template->get_Qend(),
                                 $template->get_Score(), 
                                 $template->get_SumProbL()]); 
        } 
        else{
            if ($bestscores{$template->get_Hit()}[0][2]< $template->get_Score()){
            #if existing score is smaller than the score of 
            #another part of the query-template alignment: update it
            $bestscores{$template->get_Hit()}[0][2] = $template->get_Score();
            }
        }
        }       
        else {
        #existing template "behind" 
        if (($template->get_Qend() - $bestscores{$template->get_Hit()}[0][0]) <= $MAXOVLAP) {
            #overlap <= 20 residues
            #push to key -> additional array
            push (@{$bestscores{$template->get_Hit()}}, [$template->get_Qstart(),
                                 $template->get_Qend(),
                                 $template->get_Score(), 
                                 $template->get_SumProbL()]);   
        }
        else{
            if ($bestscores{$template->get_Hit()}[0][2]< $template->get_Score()){ 
            #if existing score is smaller than the score of another 
            #part of the query-template alignment: update it
            $bestscores{$template->get_Hit()}[0][2] = $template->get_Score();
            }
        } 
        }                                       
    }   
    }


    ####################
    #search for -qsc max    
    ####################

    for (my $i=0; $i<$templateList->size(); $i++) {
    my $tt = $templateList->get($i);
    if ($tt->get_Sim() > $qscmax) { $qscmax = $tt->get_Sim(); }             
    }
    
    my @hits  = (); #array with hits of prob>80%
    my $db    = ""; #contains hhms of templates (for filtering)
    my $a3msToFilter = ""; #contains a3ms of templates (for filtering)

#    my $pid   = getppid();

    ## copy a3m files into temporary directory
    my $tmpDir = tempdir( CLEANUP => 1 );    


    for (my $i=0; $i<$templateList->size(); $i++) {
    my $tt = $templateList->get($i);
    #if probability of hit is greater than 80%
    if ($tt->get_Prob() >= 80){
        my $hhm = "$tmpDir/" . $tt->get_Hit() . ".filt.hhm";
        #check for double entries!
        if ($db !~ /$hhm/){            
        $db .= "$hhm "; 
        $a3msToFilter .=  "$tmpDir/" . $tt->get_Hit() . ".filt.a3m ";

        &System("cp " . $config->get_pdbdir() . "/" . $tt->get_Hit() . ".a3m $tmpDir/" . $tt->get_Hit() . ".filt.a3m"); 
           
        push (@hits, $tt->get_Hit());
        }   
    }       
    }   

    #############
    #filtersteps:
    #############       

    &System("cp $outbase.a3m $tmpDir/query.filt.a3m");

    
    my $hhrnr = -1; 
    
    for (my $qsc=0; $qsc<$qscmax+0.1; $qsc+=$qscStep)
    {
    #if there are hits with more than 80% probability:
    if (@hits>0){           
        $hhrnr++;
        #filter all hits with prob>= 80% from previous search with "-qsc $qsc"
        
        my $filteredA3mFiles = "";

        if ($config->get_parallelFiltering()) {
        my $multithread = $config->get_multithread();
        my $cpus = $config->get_cpus();
        my $cmd = "$multithread '$a3msToFilter' '$hhfilter -i \$file -id 100 -diff 0 -qsc $qsc -o \$file.out &> /dev/null' -cpu $cpus";
        &System("$cmd");
        sleep(1);
        $cmd = "$multithread '$a3msToFilter' '$hhmake -i \$file.out -id 100 -diff 0 -qsc $qsc -o $tmpDir/\$base.hhm &> /dev/null' -cpu $cpus";
        &System("$cmd");
        sleep(1);
        } else {
        foreach my $fileToFilter (split(/\s+/, $a3msToFilter))  {
            &System("$hhfilter -i $fileToFilter -id 100 -diff 0 -qsc $qsc -o $fileToFilter.out -v $v &> /dev/null");

            sleep 1 while (! -e "$fileToFilter.out");
            $filteredA3mFiles .= "$fileToFilter.out ";
        }
        

        foreach my $fileToMake (split(/\s+/, $filteredA3mFiles)) {
            my $base = $fileToMake;
            $base =~ s/\.a3m\.out$//;
            
            &System("$hhmake -i $fileToMake -diff 100 -o $base.hhm -v $v &> /dev/null");
            sleep 1 while (! -e "$base.hhm");
        }
        }
        
        #filter query:
        &System("$hhfilter -i $tmpDir/query.filt.a3m -id 100 -diff 0 -qsc  $qsc -o $tmpDir/query.filt.a3m -v 1");      
        
        sleep 1 while ( !(-e "$tmpDir/query.filt.a3m") );      

        &System("$hhmake -i $tmpDir/query.filt.a3m -diff 100 -o $tmpDir/query.filt.hhm -v 1"); 
        sleep 1 while ( !(-e "$tmpDir/query.filt.hhm") );      
        
        #hhsearch:
        &System($config->get_hhsearch() . " -cpu " . $config->get_cpus() . " -i $tmpDir/query.filt.hhm -d \"$db\" -o $outbase.$hhrnr.hhr $options -Z $HITS -B $HITS -atab $outbase.$hhrnr.tab");       

        sleep 1 while ( !(-e "$outbase.$hhrnr.hhr") );      

        my $tListFiltered = TemplateList->new();
        $tListFiltered->hhr_to_TemplateList("$outbase.$hhrnr.hhr");

        ###########################################
        #update rawscores & SumProbs/Length & @hits 
        ###########################################
        @hits  = ();
        $db    = "";
        $a3msToFilter = "";

        for (my $i=0; $i<$tListFiltered->size(); $i++) {
        my $filtTemplate = $tListFiltered->get($i);
        
        my $difference = 1000;
        my $chooseScore;    
        my $templ = -1;
        my $choostempl;

        #to distinguish between "same" or "different" covered part of query 
        foreach my $e (@{$bestscores{$filtTemplate->get_Hit()}}) {
            $templ++;
            
            my $startdifference = abs(@{$e}[0] - $filtTemplate->get_Qstart());      
            #difference between start of observed and start of actual queryalignmentpositon 
            #score belongs to template with smallest difference between startpositions
            if ($startdifference < $difference){
            $chooseScore = @{$e}[2];
            $choostempl  = $templ;
            $difference  = $startdifference;
            }   
        }
        if ($chooseScore< $filtTemplate->get_Score()){
            #if existing score is smaller than the filtered score: update it
            $bestscores{$filtTemplate->get_Hit()}[$choostempl][2] = $filtTemplate->get_Score();         #update scores in  %bestscores                  
        }                               

        if ($filtTemplate->get_Prob() >= 80) {  
            my $hhm = "$tmpDir/" . $filtTemplate->get_Hit() . ".filt.hhm";
            #check for double entries!
            if ($db !~ /$hhm/) {        
            $db .= "$hhm "; 
            $a3msToFilter .=  "$tmpDir/" . $filtTemplate->get_Hit() . ".filt.a3m ";
            push (@hits, $filtTemplate->get_Hit());        
            }   
        }               
        }
    
        ## add filtered hits to final templates list
        for (my $j=0; $j<$tListFiltered->size(); $j++) {
        $templateList->add_template($tListFiltered->get($j));
        }
    }
    }
    
    #update rawscores in array with all possible templates
    for (my $i=0; $i<$templateList->size(); $i++) {
    my $template = $templateList->get($i);

    my $difference = 1000;  
    my $templ      = -1;
    my $choostempl; 
    
    foreach my $e (@{$bestscores{$template->get_Hit()}}){
        $templ++;
        #difference between start of observed and start of actual queryalignmentpositon 
        my $startdifference = abs(@{$e}[0] - $template->get_Qstart());      

        if ($startdifference < $difference) {
        $choostempl = $templ;
        $difference = $startdifference;
        }   #score belongs to template with smallest difference between startpositions
    }       
    #update scores in @prediction; omg this is ugly
    $template->set_Score( $bestscores{$template->get_Hit()}[$choostempl][2] );   
    }   
    
    &System("rm -R $tmpDir");

    return $templateList;
}


1;
