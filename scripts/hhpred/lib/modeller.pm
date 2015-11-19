##########################################################
## this package contains functions handling modeller, i.e.
## creating restraint file
## changing restraint file
## calling modeller program
##########################################################

package modeller;
require(Exporter);


use strict;
use config;
use utilities;
use TemplateList;
use multiTemplateMDN;
use distanceTree;
use FixOccupancy;

our @ISA = qw(Exporter);
our @EXPORT = qw(Modeller ModellerRSR ModellerRaw ModellerNewDistance ChangeDistanceRestraints WriteCASPpdbFile reweighModellerDistRestraints);


######################################################
## input: name of query
##        directory of query
##        inbase: basename of files read in (pir)
##        outbase: basename for outputs
##
## output: indirect: py file is written and modeller
##         is called generating some more files
##
## description: calls Modeller and creates restraints
##              file - needed when restraints file
##              is to be changed
######################################################
sub ModellerRSR{
    my %args = @_;
        
    my $queryName     = $args{queryName};
    my $workingDir    = $args{workingDir};
    my $outbase       = $args{outbase};
    my $config        = $args{config};
    my $pdbdir        = $args{pdbdir} || $config->get_pdbdir();

    my $tmpdir = $workingDir;

    my $modeller = $config->get_modeller();


    #Modeller constructs only restraints file
    my $offset;    
    #Change sequence name to id_temp and get templates -> $templates
    my $templates = "";
    my $knowns    = "";
    my %knowns;

    open (IN, "<$outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    my @lines = <IN>;
    close(IN);

    $lines[0] =~ s/^>P1;.*/>P1;$queryName/;

    for (my $i=0; $i<@lines; $i++) {
	if ($lines[$i] =~ /^sequence:.*?:\s*(\d+)\s*:/) {
	    $offset = $1-1;
	    $lines[$i] =~ s/^sequence:.*?:/sequence:$queryName:/; # replace sequence name with id
	}
	elsif ($lines[$i] =~ /^structureX:/) {
	    $lines[$i-1] =~ />P1;(\S+)/;
	    $knowns .= " '$1'," if (not exists($knowns{$1}));
	    $templates .= " $1";
	    $knowns{$1} = 1; # avoid double entry in knowns, Modeller cant handle
	} 
    }

    $templates =~ s/^ //;

    open (OUT, "> $outbase.pir") or die ("Error: cannt write  $outbase.pir: $!\n");
    print(OUT @lines);
    close(OUT);
    
    # Write the py-file
    open(PY,">$outbase.py") or die("Cannot open: $outbase.py");   
 
    print PY ("# Homology modeling by the automodel class\n");
    print PY ("from modeller import *		   # Load standard Modeller classes\n");
    print PY ("from modeller.automodel import *    # Load the automodel class\n");
    print PY ("log.verbose()\n");
    print PY ("env = environ()\n");			   
    print PY ("env.io.atom_files_directory = '$pdbdir:$workingDir'\n");    
    print PY ("a = automodel(env,\n");
    print PY ("	   	alnfile = '$outbase.pir',  # alignment filename\n");
    print PY ("		knowns  = ($knowns),        # codes of the templates\n");
    print PY ("		sequence = '$queryName')    # code of the target\n");
    print PY ("a.homcsr(0)\n");	    
    close (PY);
    
    #prepare for MODELLER:    
    system("chmod 777 $outbase.pir");
    system("chmod 777 $outbase.py");
    
    #run MODELLER: create *.ini & *.rsr 
    &System("cd $tmpdir; $modeller $outbase.py > $workingDir/$queryName.modeller.log");
    
    #TREAT FAILED MODELLER!
    return;	
}


sub ModellerRaw {
    my %args = @_;
    
    my $myrsr         = $args{rsrFile};
    my $queryName     = $args{queryName};
    my $workingDir    = $args{workingDir};
    my $outbase       = $args{outbase};
    my $config        = $args{config};

  #  print "myrsr=$myrsr\n";
  #  print "queryName=$queryName\n";
  #  print "workingDir=$workingDir\n";
  #  print "outbase=$outbase\n";

    my $NMODELS = $config->get_numberOfGeneratedModels();
    my $modeller = $config->get_modeller();
    my $pdbdir = $config->get_pdbdir();

    my $offset;    
    my $templates = "";
    my $knowns    = "";
    my %knowns;

    open (IN, "< $outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    my @lines = <IN>;
    close(IN);

    $lines[0] =~ s/^>P1;.*/>P1;$queryName/;

    for (my $i=0; $i<@lines; $i++) {
	if ($lines[$i] =~ /^sequence:.*?:\s*(\d+)\s*:/) {
	    $offset = $1-1;
	    $lines[$i] =~ s/^sequence:.*?:/sequence:$queryName:/; # replace sequence name with id
	}
	elsif ($lines[$i] =~ /^structureX:/) {
	    $lines[$i-1] =~ />P1;(\S+)/;
	    $knowns .= " '$1'," if (not exists($knowns{$1})); 
	    $templates .= " $1";
	    $knowns{$1} = 1; # avoid double entry in knowns, Modeller cant handle
	} 
    }

    $templates =~ s/^ //;

    open (OUT, ">$outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    print(OUT @lines);
    close(OUT);
    
    # Write the py-file
    open(PY, "> $outbase.py") or die("Cannot write $outbase.py");   
 
    print PY ("# Homology modeling by the automodel class\n");
    print PY ("from modeller import *			# Load standard Modeller classes\n");
    print PY ("from modeller.automodel import *   	# Load the automodel class\n");
    print PY ("log.verbose()\n");
    print PY ("env = environ()\n");
    print PY ("# directories for input atom files\n");
    print PY ("env.io.atom_files_directory = '$pdbdir:$workingDir'\n");    
    print PY ("a = automodel(env,\n");
    print PY ("	   	alnfile = '$outbase.pir',  # alignment filename\n");
    print PY ("		knowns  =($knowns),        # codes of the templates\n");
    print PY ("		sequence ='$queryName',    # code of the target\n");
    print PY ("         csrfile  = '$myrsr') \n");
    print PY ("a.starting_model= 1          	   # index of the first model\n");
    print PY ("a.ending_model = $NMODELS           # index of the last model\n");
    print PY ("a.make()				   # do the actual homology modeling\n");    
    close (PY);
    
    #prepare for MODELLER:    
    system("chmod 777 $outbase.pir");
    system("chmod 777 $outbase.py");
    
    #run MODELLER: create *.ini & *.rsr 
    &System("cd $workingDir; $modeller $outbase.py > $workingDir/$queryName.modeller.log");
    

    #TREAT FAILED MODELLER!
    #check model(s):
    my @modelsBuilt;

    for (my $i=1; $i<=$NMODELS; $i++) {
	my $modelID = "B" . (99990000 + $i);
	if (! -e "$outbase.$modelID.pdb") {
	    print("WARNING: MODELLER failed for $outbase.$modelID.pdb!\nSee log file $outbase.modeller.log\n");
	
	    #EXCEPTION does not work!!!!debug hhmakemodel.pl
	    #&System("perl $hh/hhmakemodel.pl $tmpname.hhr -m $singleTemp -q $tmpname.a3m -v -ts $tmpname.mod -d $pdbdir");
	}
	else {
	    push (@modelsBuilt, $i);
	}
    }

    # find model with minimal energy
    if (scalar(@modelsBuilt) > 0) {
	my $minEnergy = +1E8;
	my $imin = -1;
	my $line;

	my $bestModelID;

	## search for best model (lowest objective function value)
	## and save this model - remove all others
	foreach my $model (@modelsBuilt) {
	    my $modelID = "B" . (99990000 + $model);
	    open(IN, "<$outbase.$modelID.pdb") || die("Error: can't open $outbase.$modelID.pdb");
	    my $energy;
	    while ($line = <IN>) {
		if ($line =~ /MODELLER OBJECTIVE FUNCTION:\s*(\S+)/) {
		    $energy = $1; 
		    last;
		}	
	    }
	    if ($energy < $minEnergy) {
		$minEnergy = $energy; 
		$imin = $model;
	    }
	    close(IN);
	}
	if ($imin != -1) {
	    my $minModelID = "B" . (99990000 + $imin);
	    &System("cp $outbase.$minModelID.pdb $outbase.pdb");
	} else {
	    print "WARNING: MODELLER no minimal energy structure found!\n";
	}
    }
    
    &System("rm $outbase.sch $outbase.V999* $outbase.D000*");    
    return;	
}




sub Modeller {
    my %args = @_;
    
    my $queryName     = $args{queryName};
    my $workingDir    = $args{workingDir};
    my $outbase       = $args{outbase};
    my $config        = $args{config};
    my $physicalScale = 1;
    $physicalScale = $args{physicalScale} if (defined($args{physicalScale}));
    my $pdbdir        = $args{pdbdir} || $config->get_pdbdir();
    

    my $NMODELS = $config->get_numberOfGeneratedModels();
    my $modeller = $config->get_modeller();

    my $offset;    
    my $templates = "";
    my $knowns    = "";
    my %knowns;

    open (IN, "< $outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    my @lines = <IN>;
    close(IN);

    $lines[0] =~ s/^>P1;.*/>P1;$queryName/;

    for (my $i=0; $i<@lines; $i++) {
	if ($lines[$i] =~ /^sequence:.*?:\s*(\d+)\s*:/) {
	    $offset = $1-1;
	    $lines[$i] =~ s/^sequence:.*?:/sequence:$queryName:/; # replace sequence name with id
	}
	elsif ($lines[$i] =~ /^structureX:/) {
	    $lines[$i-1] =~ />P1;(\S+)/;
	    $knowns .= " '$1'," if (not exists($knowns{$1}));
	    $templates .= " $1";
	    $knowns{$1} = 1; # avoid double entry in knowns, Modeller cant handle
	} 
    }

    $templates =~ s/^ //;

    open (OUT, ">$outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    print(OUT @lines);
    close(OUT);
    
    # Write the py-file
    open(PY, "> $outbase.py") or die("Cannot write $outbase.py");   
 
    print PY ("# Homology modeling by the automodel class\n");
    print PY ("from modeller import *			# Load standard Modeller classes\n");
    print PY ("from modeller.automodel import *   	# Load the automodel class\n");
    print PY ("log.verbose()\n");
    print PY ("env = environ()\n");

    if ($physicalScale != 1) {
	print PY ("env.schedule_scale=physical.values(default=1.0, ca_distance=$physicalScale, n_o_distance=$physicalScale, sd_mn_distance=$physicalScale, sd_sd_distance=$physicalScale)\n");	
    }

    print PY ("# directories for input atom files\n");
    print PY ("env.io.atom_files_directory = '$pdbdir:$workingDir'\n");    
    print PY ("a = automodel(env,\n");
    print PY ("	   	alnfile = '$outbase.pir',  # alignment filename\n");
    print PY ("		knowns  = ($knowns),        # codes of the templates\n");
    print PY ("		sequence = '$queryName')    # code of the target\n");
    print PY ("a.starting_model= 1          	   # index of the first model\n");
    print PY ("a.ending_model = $NMODELS           # index of the last model\n");
 #   print PY ("a.deviation = 1.0\n");
 #   print PY ("a.spline_on_site = False            # do not use splines for approximation\n");
 #   print PY ("a.library_schedule = autosched.very_fast\n");
 #   print PY ("a.max_var_iterations = 300\n");
 #   print PY ("a.repeat_optimization = 2\n");
 #   print PY ("a.write_intermediates = True\n");
    print PY ("a.make()				   # do the actual homology modeling\n");    
    close (PY);
    
    #prepare for MODELLER:    
    system("chmod 777 $outbase.pir");
    system("chmod 777 $outbase.py");
    
    #run MODELLER: create *.ini & *.rsr 
    &System("cd $workingDir; $modeller $outbase.py > $workingDir/$queryName.modeller.log");
    

    my @modelsBuilt;

    #TREAT FAILED MODELLER!
    #check model(s):
    for (my $i=1; $i<=$NMODELS; $i++) {
	my $modelID = "B" . (99990000 + $i);
	if (! -e "$outbase.$modelID.pdb") {
	    print("WARNING: MODELLER failed for $outbase.$modelID.pdb!\nSee log file $outbase.modeller.log\n");
	
	    #EXCEPTION does not work!!!!debug hhmakemodel.pl
	    #&System("perl $hh/hhmakemodel.pl $tmpname.hhr -m $singleTemp -q $tmpname.a3m -v -ts $tmpname.mod -d $pdbdir");
	}
	else {
	    push (@modelsBuilt, $i);
	}
    }

    # find model with minimal energy
    if (scalar(@modelsBuilt) > 0) {
	my $minEnergy = +1E8;
	my $imin = -1;
	my $line;

	## search for best model (lowest objective function value)
	## and save this model - remove all others
	foreach my $model (@modelsBuilt) {
	    my $modelID = "B" . (99990000+$model);
	    open(IN, "<$outbase.$modelID.pdb") || die("Error: can't open $outbase.$modelID.pdb");
	    my $energy;
	    while ($line = <IN>) {
		if ($line =~ /MODELLER OBJECTIVE FUNCTION:\s*(\S+)/) {
		    $energy = $1; 
		    last;
		}	
	    }
	    if ($energy < $minEnergy) {
		$minEnergy = $energy; 
		$imin = $model;
	    }
	    close(IN);
	}
	if ($imin != -1) {
	    my $minModelID = "B" . (99990000+$imin);
	    &System("cp $outbase.$minModelID.pdb $outbase.pdb");
	} else {
	    print "WARNING: MODELLER no minimal energy structure found!\n";
	}
    }
    
    &System("rm $outbase.sch $outbase.V999* $outbase.D000* $outbase.ini ");    
    return;	
}



#############################################################
## input: myrsr:      rsr file
##        queryName:
##        workingDir: where save (temp) files
##        outbase:    outbase for pir files
##        config:     
##        
##        
## output: mod file containing the best modeller built model
############################################################# 
sub ModellerNewDistance {
    my %args = @_;
    
    my $myrsr         = $args{rsrFile};
    my $queryName     = $args{queryName};
    my $workingDir    = $args{workingDir};
    my $outbase       = $args{outbase};
    my $config        = $args{config};
    my $pdbdir        = $args{pdbdir} || $config->get_pdbdir();
    my $templWeightNormalizer = $args{normalizer};
    my $doTemplWeight = 1;
    $doTemplWeight = $args{doTemplWeightEqual} if (defined($args{doTemplWeightEqual}));
    my $physicalScale = 1;
    $physicalScale = $args{physicalScale} if (defined($args{physicalScale}));


    my $NMODELS = $config->get_numberOfGeneratedModels();
    my $modeller = $config->get_modeller();  

    my $numTemplates = 0;

    #Change sequence name to id_temp and get templates -> $templates
    my $templates = "";
    my $knowns = "";
    my %knowns;
    my @allTemplates; ## $knowns does not contain duplicates
    my %templateWeights;

    open (IN, "< $outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    my @lines = <IN>;
    close(IN);

    $lines[0] =~ s/^>P1;.*/>P1;$queryName/;

    for (my $i=0; $i<@lines; $i++) {
	if ($lines[$i] =~ /^sequence:.*?:\s*(\d+)\s*:/) {
	    $lines[$i] =~ s/^sequence:.*?:/sequence:$queryName:/; # replace sequence name with id
	}
	elsif ($lines[$i] =~ /^structureX:/) {
	    $lines[$i-1] =~ />P1;(\S+)/;
	    push (@allTemplates, $1);
	    $knowns .= " '$1'," if (not exists($knowns{$1}));
	    $templates .= " $1";
            $numTemplates++;
	    $knowns{$1} = 1;
	} 
    }

    $templates =~ s/^ //;

    ## weight for template based restraints
    ## default
    my $rsrWeight = 1.0/$numTemplates;
    $rsrWeight = sprintf("%.3f", $rsrWeight);
    my $templRsrWeight = 1.0/$templWeightNormalizer;
    ## TODO: adapt weight to new template-weighting strategy

    open (OUT, ">$outbase.pir") or die ("Error: cannot open $outbase.pir: $!\n");
    print(OUT @lines);
    close(OUT);
    
    # Write the py-file
    open(PY,"> $outbase.py") or die("Cannot open: $outbase.py");
    
    if ($config->get_templateWeightStrategy() == 1) {  ## WITH phylo weight
	print PY ("# Homology modeling by the automodel class\n");
	print PY ("from modeller import *			# Load standard Modeller classes\n");
	print PY ("from modeller.automodel import *   	        # Load the automodel class\n");

	if ($config->get_doParallelModeller()) {
	    print PY ("from modeller.parallel import *\n");
	}

	print PY ("from modeller.scripts import complete_pdb\n");
	
	print PY ("from modeller import _cuser_form54\n");
	print PY ("log.verbose()\n");
	print PY ("env = environ()\n");

	## downweight template based restraints
	if ($doTemplWeight) {
	    my $rsrWeight = $physicalScale;
	    print PY ("env.schedule_scale=physical.values(default=1.0, ca_distance=$rsrWeight, n_o_distance=$rsrWeight, sd_mn_distance=$rsrWeight, sd_sd_distance=$rsrWeight)\n"); #, chi1_dihedral=$templRsrWeight, chi2_dihedral=$templRsrWeight, chi3_dihedral=$templRsrWeight, chi4_dihedral=$templRsrWeight, phi_dihedral=$templRsrWeight, psi_dihedral=$templRsrWeight, omega_dihedral=$templRsrWeight, phi_psi_dihedral=$templRsrWeight)\n");	
	}

	print PY ("env.io.atom_files_directory = '$pdbdir:$workingDir'\n");
	print PY ("\n");
	
	print PY ("class MyGauss(forms.restraint_form):\n"); 
	print PY ("    \"\"\"Multiple Gaussian\"\"\"\n");
	print PY ("    _builtin_index = _cuser_form54.myform_create()\n");
	print PY ("    def __init__(self, group, feature, weights, means, stdevs, phyloWeight):\n");
	print PY ("        lv = -1\n");
	print PY ("        for var in (weights, means, stdevs):\n");
	print PY ("            if (lv >= 0 and lv != len(var)) \\\n");
	print PY ("               or not isinstance(var, (tuple, list)):\n");
	print PY ("                raise TypeError(\"weights, means and stdevs should all be \"\n");
	print PY ("                                + \"sequences of the same length\")\n");
	print PY ("            lv = len(var)\n");
	print PY ("        forms.restraint_form.__init__(self, group, feature, len(weights),\n");
	print PY ("                                tuple(weights) + tuple(means) + tuple(stdevs) + phyloWeight)\n"); 
	print PY ("a = automodel(env,\n");
	print PY ("	   	alnfile = '$outbase.pir',  # alignment filename\n");
	print PY ("		knowns  =($knowns),        # codes of the templates\n");
	print PY ("		sequence ='$queryName',    # code of the target\n");
	print PY ("		csrfile  = '$myrsr')       # use 'my' restraints file\n");
	print PY ("a.starting_model= 1          	   # index of the first model\n");
	print PY ("a.ending_model = $NMODELS           # index of the last model\n");
	#  print PY ("a.library_schedule = autosched.slow\n");
	#  print PY ("a.max_var_iterations = 300\n");
	#  print PY ("a.repeat_optimization = 2\n");
#	print PY ("a.deviation = 0.0\n");
#	print PY ("a.write_intermediates = True\n");

	if ($config->get_doParallelModeller()) {
	    ## all slaves have to have their own initialization with new restraints...
	    print PY ("j = job()\n");
	    print PY ("j.slave_startup_commands.append(\"from modeller import _cuser_form54\")\n"); 
	    print PY ("j.slave_startup_commands.append(\"class MyGauss(forms.restraint_form):\\n\" + \n");
            print PY ("                                \"      _builtin_index = _cuser_form54.myform_create()\")\n");
	    my $cpus = $config->get_cpus();
	    for (my $i=0; $i < &min($NMODELS, $cpus); $i++) {
		print PY ("j.append(local_slave())\n");
	    }
	    print PY ("a.use_parallel_job(j)\n");
	}

	print PY ("a.make()				# do the actual homology modeling\n");    
    } else {
	## NO PHYLO WEIGHT
	print PY ("# Homology modeling by the automodel class\n");
	print PY ("from modeller import *			# Load standard Modeller classes\n");
	print PY ("from modeller.automodel import *   	        # Load the automodel class\n");

	if ($config->get_doParallelModeller()) {
	    print PY ("from modeller.parallel import *\n");
	}

	print PY ("from modeller.scripts import complete_pdb\n");
	
	print PY ("from modeller import _cuser_form52\n");
	print PY ("log.verbose()\n");
	print PY ("env = environ()\n");

	## downweight template based restraints
	if ($doTemplWeight) {
	    ## !!! TEST !!!
	    my $rsrWeight = $physicalScale;
	    print PY ("env.schedule_scale=physical.values(default=1.0, ca_distance=$rsrWeight, n_o_distance=$rsrWeight, sd_mn_distance=$rsrWeight, sd_sd_distance=$rsrWeight)\n"); #, chi1_dihedral=$rsrWeight, chi2_dihedral=$rsrWeight, chi3_dihedral=$rsrWeight, chi4_dihedral=$rsrWeight, phi_dihedral=$rsrWeight, psi_dihedral=$rsrWeight, omega_dihedral=$rsrWeight, phi_psi_dihedral=$rsrWeight)\n"); 
	    ## !!! end TEST !!!
	}

	print PY ("env.io.atom_files_directory = '$pdbdir:$workingDir'\n");
	print PY ("\n");
	
	print PY ("class MyGauss(forms.restraint_form):\n"); 
	print PY ("    \"\"\"Multiple Gaussian\"\"\"\n");
	print PY ("    _builtin_index = _cuser_form52.myform_create()\n");
	print PY ("    def __init__(self, group, feature, weights, means, stdevs):\n");
	print PY ("        lv = -1\n");
	print PY ("        for var in (weights, means, stdevs):\n");
	print PY ("            if (lv >= 0 and lv != len(var)) \\\n");
	print PY ("               or not isinstance(var, (tuple, list)):\n");
	print PY ("                raise TypeError(\"weights, means and stdevs should all be \"\n");
	print PY ("                                + \"sequences of the same length\")\n");
	print PY ("            lv = len(var)\n");
	print PY ("        forms.restraint_form.__init__(self, group, feature, len(weights),\n");
	print PY ("                                tuple(weights) + tuple(means) + tuple(stdevs))\n"); 
	print PY ("a = automodel(env,\n");
	print PY ("	   	alnfile = '$outbase.pir',  # alignment filename\n");
	print PY ("		knowns  =($knowns),        # codes of the templates\n");
	print PY ("		sequence ='$queryName',    # code of the target\n");
	print PY ("		csrfile  = '$myrsr')       # use 'my' restraints file\n");
	print PY ("a.starting_model= 1          	   # index of the first model\n");
	print PY ("a.ending_model = $NMODELS           # index of the last model\n");
	#  print PY ("a.library_schedule = autosched.slow\n");
	#  print PY ("a.max_var_iterations = 300\n");
	#  print PY ("a.repeat_optimization = 2\n");
#	print PY ("a.deviation = 0.0\n");
#	print PY ("a.write_intermediates = True\n");

	if ($config->get_doParallelModeller()) {
	    ## all slaves have to have their own initialization with new restraints...
	    print PY ("j = job()\n");
	    print PY ("j.slave_startup_commands.append(\"from modeller import _cuser_form52\")\n"); 
	    print PY ("j.slave_startup_commands.append(\"class MyGauss(forms.restraint_form):\\n\" + \n");
            print PY ("                                \"      _builtin_index = _cuser_form52.myform_create()\")\n");
	    my $cpus = $config->get_cpus();
	    for (my $i=0; $i < &min($cpus,$NMODELS); $i++) {
		print PY ("j.append(local_slave())\n");
	    }
	    print PY ("a.use_parallel_job(j)\n");
	}

	print PY ("a.make()				# do the actual homology modeling\n");    
    }

    close (PY);
    
    system("chmod 777 $outbase.pir");
    system("chmod 777 $outbase.py");
    
    #run MODELLER: 
    if ($config->get_doParallelModeller()) {
	my $modellerParallel = $config->get_modellerParallel();
	&System("cd $workingDir; $modellerParallel $outbase.py > $workingDir/$queryName.modeller.log");
    } else {
	&System("cd $workingDir; $modeller $outbase.py > $workingDir/$queryName.modeller.log");
    }

    ## modelling results, e.g. pdb files
    my $resultbase = $outbase;

    my @modelsBuilt;

    #TREAT FAILED MODELLER!
    #check model(s):
    for (my $i=1; $i<=$NMODELS; $i++) {
	my $modelID = "B" . (99990000+$i);
	if (! -e "$outbase.$modelID.pdb") {
	    print("WARNING: MODELLER failed for $outbase.$modelID.pdb!\nSee log file $outbase.modeller.log\n");
	
	    #EXCEPTION does not work!!!!debug hhmakemodel.pl
	    #&System("perl $hh/hhmakemodel.pl $tmpname.hhr -m $singleTemp -q $tmpname.a3m -v -ts $tmpname.mod -d $pdbdir");
	}
	else {
	    push (@modelsBuilt, $i);
	}
    }

    # find model with minimal energy
    if (scalar(@modelsBuilt) > 0) {
	my $minEnergy = +1E8;
	my $imin = -1;
	my $line;

	## search for best model (lowest objective function value)
	## and save this model - remove all others
	foreach my $model (@modelsBuilt) {
	    my $modelID = "B" . (99990000+$model);
	    open(IN, "<$outbase.$modelID.pdb") || die("Error: can't open $outbase.$modelID.pdb");
	    my $energy;
	    while ($line = <IN>) {
		if ($line =~ /MODELLER OBJECTIVE FUNCTION:\s*(\S+)/) {
		    $energy = $1; 
		    last;
		}	
	    }
	    if ($energy < $minEnergy) {
		$minEnergy = $energy; 
		$imin = $model;
	    }
	    close(IN);
	}
	if ($imin != -1) {
	    my $minModelID = "B" . (99990000+$imin);
	    &System("cp $outbase.$minModelID.pdb $resultbase.pdb");
	} else {
	    print "WARNING: MODELLER no minimal energy structure found!\n";
	}
    }
    
    &System("rm $resultbase.sch $resultbase.V999* $resultbase.D000*");    
    
    #Add END line
    #my $repairpdb = $config->get_repairPDB();
    ## repairpdb is not necessary any more in newer Modeller version
   # &System("$repairpdb $resultbase.mod &> /dev/null");
    
    return;    
}


sub readTemplateWeightsFromFile {
    my $file = shift;

    my %weights;
    open(TW, "< $file") or die "Cant read $file: $!\n";
    while (my $line = <TW>) {
	if ($line =~ /(\S+)\s+(\S+)/) {
	    $weights{$1} = $2;
	}
    }
    close(TW);

    return %weights;
}


#############################################
## ugly but faster than passing by reference
## maybe avoid it by writing a class instead
my @PDBs;
my @atomsForResidue;


#################################################################
## input: rsrfile
##        inifile
##        inbase
##        array of templates
## output: side effect: write out a new restraint file:
##         rsrfile.new.rsr
#################################################################
sub ChangeDistanceRestraints {	
    my $outbase = shift;
    my $workingDir = shift;
    my $templateList = shift;
    my $config = shift;
    my $pdbdir = shift || $config->get_pdbdir();

    my $queryName = &get_basename($outbase);

    my $weightTemplates = $config->get_templateWeightStrategy();
    my %templateWeight;
    my $templWeightNormalizer = 0;

    @PDBs = ();
    @atomsForResidue = ();

    if ($weightTemplates == 1) {
	if ($templateList->size() == 1) {
	    $templateWeight{$templateList->get(0)->get_Hit()} = 1;	    
	    $templWeightNormalizer = 1;
	} else {
	    ## copy template hhm files to workingDir
	    for (my $i=0; $i<$templateList->size(); $i++) {
			my $tname = $templateList->get($i)->get_Hit();
			system("cp $pdbdir/$tname.hhm $workingDir");
	    }
	    ## build upgma tree, leaves are numbered according to row-number in templateList (0=query)
	    
	    my @distanceMatrix = &calculatePairwiseDistances($queryName, $templateList, $workingDir, $workingDir);
	    my $phyloTree = &upgma(@distanceMatrix);

	    ## root tree at query (= node 0)
	    print "PhyloTree:\n";
	    print "==========\n";
	    print $phyloTree->print() . "\n";
	    my $rerootedTree = $phyloTree->reRoot(0);
	    print "rerootedTree:\n";
	    print "=============\n";
	    $rerootedTree->print();
	    print "=============\n";

	    ## tweights containts leaf/row=>weight entries
	    my %tweights = &iterativelyWeightTemplates($rerootedTree);
	    
	    ## replace indices with template names
	    foreach my $leaf (keys %tweights) {
		my $tname = $templateList->get($leaf-1)->get_Hit();
		$templateWeight{$tname} = $tweights{$leaf};
	    }
	    print "ChangeDistanceRestraints template weights\n";
	    foreach my $template (keys %templateWeight) {
		print "$template=>$templateWeight{$template}\n";
		$templWeightNormalizer += $templateWeight{$template};
	    }
	    print "==============================\n";
	}
    } else {
	$templWeightNormalizer = $templateList->size();
    }

    ## these files must have been generated before
    my $rsrfile = "$outbase.rsr";
    my $inifile = "$outbase.ini";
    my $pirfile = "$outbase.pir";

    #get query length
    #my $querylength = $templateList->get_queryLength();

    ## set up mixture density networks
    ## these will be needed later in order to calculate parameters of new contraints
    my $mdnCACA = multiTemplateMDN->new($config->get_MDNWeightsLayer1CACA(), $config->get_MDNWeightsLayer2CACA());
    my $mdnNO = multiTemplateMDN->new($config->get_MDNWeightsLayer1NO(), $config->get_MDNWeightsLayer2NO());
    my $mdnSCMC = multiTemplateMDN->new($config->get_MDNWeightsLayer1SCMC(), $config->get_MDNWeightsLayer2SCMC());
    my $mdnSCSC = multiTemplateMDN->new($config->get_MDNWeightsLayer1SCSC(), $config->get_MDNWeightsLayer2SCSC());    

    
    #get similarity of each template and read TAB file(s)
    my @sim;	  #contains similarities of Templates
    my @templates;#contains PDB-Identifier of Templates	
    my @TABprobs; #contains probability of ij to according i from *.tab file			
    
    #use *.afh file (hit list with filter strength of each hit):
    #	
    #HHR   HIT  TEMPLATE   	PROB   SCORE   SS/Length   	SIM   	SumPr/L	QueryStart-End	COLS	TempStart-End	TMSCORE
    #1		1	1rwz_A		98.0	56.2	0.0833		0.164	0.5697		3	225		213			5	245		0.522695	
    #2		1	1rwz_A		97.3	56.2	0.0882		0.177	0.5412		2	225		211			9	245		0.518335	
    #0		1	1rwz_A		98.0	56.2	0.0794		0.173	0.5697		2	225		214			5	245		0.517708	
    #1		5	1plq		97.0	43.8	0.0877		0.094	0.5342		1	225		217			2	257		0.503807	
    #...
    
   	    
    ## for each template ...
    for(my $t=0; $t<$templateList->size(); $t++) {
	my $template = $templateList->get($t);

	$templates[$t] = $template->get_Hit();
	my $hitnr = $template->get_No();
	my $hhrnr = $template->get_Filt();

	$sim[$t] = $template->get_Sim();
	#my $templatelength=$4-$3+1;

	my $tab = "$outbase." . $templates[$t] . ".HIT$hitnr.tab";
	print "HIT: $t TAB:$tab SIM:$sim[$t]\n";
	

	## ... build TAB file, based on corresponding tab file
	if (! (-e "$tab")) {
	    &BuildSingleTabFile("$outbase.$hhrnr.tab", $hitnr, $outbase);
	}

	open(TAB,"< $tab") or die("Error reading $tab: $!\n");

	while (my $line = <TAB>) {									
	    #   i     j  score     SS  probab
	    #	4   400  -0.64  -0.28  0.0101
	    if ($line =~ /\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\S+)/) {
		my $i = $1;
		my $j = $2;
		my $probab = $3;
		if (!($TABprobs[$i])) {
		    $TABprobs[$1] = [$probab,$t,$j];
		}
		else { 
		    push(@{$TABprobs[$i]},$probab,$t,$j);
		}
	    }
	}
	close(TAB);

	my $seenValidLine = 0;

	## ... read template PDB
	open(PDBTemp, "< $pdbdir/$templates[$t].pdb") or die("Error reading $pdbdir/$templates[$t].pdb: $!\n");		
	while (my $line = <PDBTemp>){ 
	    #ATOM           1        N         PRO          1              -29.477        -9.021        66.175  1.00 75.79           N
	    #ATOM        3781       CB         MET    1    477              25.003       121.648        32.112  1.00 90.08           C
	    if ($line =~ /^ATOM.{2}\s*(\d+)\s.([\S|\s]{3}).(\S{3})\s.\s*(\d+)[\s|\D]\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*/){	    		  	    		
		$seenValidLine = 1;
		## save Atomname Residuename ResidueNr x y z	    	
		$PDBs[$t][$1]=[$2,$3,$4,$5,$6,$7]; 

		## speed up: save all atoms for a given residue (in a given template); will be needed later to compute the distance between "matching" atoms
	 	if (not defined(@{$atomsForResidue[$t][$4]})) {
 		    $atomsForResidue[$t][$4] = [$1];
 		}
 		else {
 		    push(@{$atomsForResidue[$t][$4]}, $1);
	        }
	    }
	}	

	if ($seenValidLine) {
	    for (my $r=0; $r < @{$atomsForResidue[$t]}; $r++) {
		my $rr = @{$atomsForResidue[$t]}[$r];
		next if (not defined($rr));
	    }
	}
#	$atomsForResiue[$t] = \%atomsForResidue;
	close(PDBTemp);			  
    } ## end for all templates


    #e.g.: print"$PDBs[4][1227][0]	$PDBs[4][1227][1]	$PDBs[4][1227][2]\n";	 
    #$TABprobs[x][0] x=atom 0=probabilityij
    #$TABprobs[x][1]        1=#of tabfile
    #$TABprobs[x][2]	    2=probabilityij(e.g. if more than one template is used!!!)
    #$TABprobs[x][3]        3=aligned residue number of template
    #$TABprobs[x][4]        4=#of tabfile
    #...		    5=probabilityij(e.g. if more than one template is used!!!)
    #...					...
    #################################
    #read INIPDB file of query: *.ini 
    #################################
    #INIPDB-FILE:
    #ATOM      1  N   PRO     1     -29.477  -9.021  66.175  1.00 75.79           N
    #ATOM      2  CA  PRO     1     -28.807  -9.573  65.023  1.00 75.45           C
    #...
    #1234567890123456789012345678901234567890...
    
    #find atom name in query ini-pdb:
    my @iniPDBresNr;
    open(iniPDB,"< $inifile") or die("Error reading file: $inifile: $!\n");

    while (my $line = <iniPDB>){
	if ($line =~ /^ATOM.{2}\s*(\d+)\s.([\S|\s]{3}).(\S{3})\s.\s*(\d+)[\s|\D]\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*/){
	    ## atomNr => residueNr, atomName, residueName
	    $iniPDBresNr[$1] = [$4,$2,$3];
	}
    }
    close(iniPDB);	

    #0:
    #1: 1    N       MET
    #2: 1    CA      MET
    #3: 1    CB      MET
    #4: 1    CG      MET
    #...									
    #####################
    #read restraints file
    #####################
    open(RSR,"< $rsrfile") or die("Error reading file: $rsrfile: $!\n");	

    my @rsrlines = <RSR>;					

    ## R    3   1   1   1   2   2   1     3     2       1.5380    0.0364 
    ## R    3   1   1  26   2   2   1  1394  1422       3.8982    6.2527 
    close(RSR);
    
    #Spline restraint
    #	Form 	Modality	Feature		Group		Numb_atoms	Numb_parameters		Numb_Feat	Atom_indices	Parameter(e.g. mean, standard deviation)
    #R   	10  	39   		1   		9   		2  			45   				1    		25   363       1.0000    0.0000   26.6000    0.7000   -0.6596...	
    
    ############################################
    #change old restraints in Modeller rsr-file:  
    ############################################
    ## comment: simply using Modellers mean for the template distance
    ##          works only for single template modelling (otherwise splines are used), therefore
    ##          one needs a subroutine to extract the distances in templates
    ##
    my @newrsrfile;

    ## if no distance found in template for a specific restraint, the old one can be used (remainRestr=1)
    ## or the restraint can be neglected (rmainRestr=0)
    my $remainRestr = 1;


    for (my $i=1; $i<@rsrlines; $i++) {		
	#	
	if ($i % 1000 == 0) {					                              
	    print "$i iterations in changeRestraints\n";
	}

	## 9|10|23|26
	if ($rsrlines[$i] =~ /^R\s+(3|10)\s+\d+\s+\d+\s+(9|10|23|26)\s+2\s+\d+\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)/) {   
	    #print "\nLINE:$i	$rsrlines[$i]";
	    #print "rsrType:$1	group:$2	atom1:$3	atom2:$4\n";
	    #search i, j, delta from *.ini 
	    my $Qi;
	    my $Qj;	
	    my $queryatomnameOne;			
	    my $queryatomnameTwo;			
	    my $queryresiduenameOne;
	    my $queryresiduenameTwo;	
	    
	    if ($iniPDBresNr[$3][0] && $iniPDBresNr[$4][0]){
		$Qi = $iniPDBresNr[$3][0];			#residue i of Query (number)
		$Qj = $iniPDBresNr[$4][0];			#residue j of Query (number)
		$queryatomnameOne = $iniPDBresNr[$3][1];	#atomname of residue i of Query
		$queryatomnameTwo = $iniPDBresNr[$4][1];	#atomname of residue j of Query	
		$queryresiduenameOne = $iniPDBresNr[$3][2];	#residuename of i of Query
		$queryresiduenameTwo = $iniPDBresNr[$4][2];	#residuename of j of Query
	    }
	    else { die("Error in file: $inifile"); }
	    
	    my $form = $1;
	    my $rsrType = $2;
	    my $atom1 = $3;
	    my $atom2 = $4;
	    my $checkdelta = $5;
	    ################################################################################################################
	    #search prob & sim for i-i' and j-j'from TAB file & calculate new values for new log-normal Gaussian distribution:			
	    if ($TABprobs[$Qi][0] && $TABprobs[$Qj][0]) {
		#e.g.:
		#      atom     	prob	Template	alternative probability from alternative Template
		#	55: 		0.1119      1
		#	56: 		0.0874      1
		#	57: 		0.0733      1
		#	58: 		0.3256      0       	0.0454  1
		#	59: 		0.3552      0       	0.0369  1
		#	60: 		0.4030      0       	0.0274  1
		#	61: 		0.4952      0       	0.0178  1
		#	62: 		0.5366      0
		#	63: 		0.5986      0
		#	64: 		0.6060      0
		
		#e.g. case2 (more than one template for i-i'):	atom 59: prob i-i' 0.3552 in Template 1 and 0.0369 in Template 2
		#						atom 63: prob j-j' 0.5986 in Template 1
		#     						so take prob i-i' 0.3552 from Template 1 (NOT from Template 2!!!)				

		# $TABprobs[$i][+0] = probab
		# $TABprobs[$i][+1] = t
		# $TABprobs[$i][+2] = j
		
		#if there is more than one template for i-i' AND j-j'
                #if there is more than one template for i-i' and template for j-j'
		# Q  |---------------------------------------------------------------|
		# T1              |--------------------------------|
		# T2         |-------------------------------|
		#                    j             i        
		if (@{$TABprobs[$Qi]}>3 && @{$TABprobs[$Qj]}>3){				
		    for (my $t=0;$t<@{$TABprobs[$Qi]}; $t=$t+3){

			my $probi = $TABprobs[$Qi][$t];						
			for (my $r=0;$r<@{$TABprobs[$Qj]}; $r=$r+3){
			    ## Qi and Qj aligned to residues in same template
			    if ($TABprobs[$Qi][$t+1] == $TABprobs[$Qj][$r+1]){
				my $probj =  $TABprobs[$Qj][$r];				#probability of j-j
				my $sim = $sim[$TABprobs[$Qj][$r+1]];				#similarity of according template
				my $dtemplate = &DistanceInTemplate($atom1,$atom2, $Qi,$Qj,$TABprobs[$Qi][$t+2],$TABprobs[$Qj][$r+2],$queryatomnameOne,$queryatomnameTwo,$queryresiduenameOne,$queryresiduenameTwo,$TABprobs[$Qj][$r+1], 1);
						
				if ( ($dtemplate ne "NULL")) {
				    my $posteriori = $probi*$probj;

				    if ($Qi - $TABprobs[$Qi][$t+2] == $Qj - $TABprobs[$Qj][$r+2]) {
					$posteriori = min($probi, $probj);
				    }
				    
				   # next if ($posteriori < 0.2);
				    ## skip too long distance-restraints
				    next if (($dtemplate > 28 and $rsrType==9) or ($dtemplate > 19 and $rsrType==10) or ($dtemplate > 11 and $rsrType==23) or ($dtemplate>10 and $rsrType==26));


				    my $newline = "";
				    my $hit = $templateList->get($TABprobs[$Qi][$t+1])->get_Hit();

				    if ($weightTemplates) {
					my $tweight = $templateWeight{$hit};
					$newline = &CalcRSRTweight($posteriori, $sim, $dtemplate, $rsrType, $atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC, $tweight);
				    }
				    else {
					$newline = &CalcRSR($posteriori, $sim, $dtemplate, $rsrType, $atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC);
				    }

				    $newline .= "  $hit $posteriori $sim $dtemplate $Qi $Qj $TABprobs[$Qi][$t+2] $TABprobs[$Qj][$r+2]\n";
				    push(@newrsrfile, $newline);									
				}
				elsif ($remainRestr) { push(@newrsrfile, $rsrlines[$i]); }
			    }
			}											
		    }
		}
		#if there is more than one template for i-i' and one template for j-j'
		# Q  |---------------------------------------------------------------|
		# T1                         |--------------------------------|
		# T2         |-------------------------------|
		#                    j             i        
		elsif (@{$TABprobs[$Qi]}>3 && @{$TABprobs[$Qj]}==3) {
		    #print "more than one template for i-i'		atom:	$atom1	$atom2\n";
		    my $probi;
		    my $probj = $TABprobs[$Qj][0];			#probability of j-j'
		    my $sim = $sim[$TABprobs[$Qj][1]];		#similarity of according template
		    my $dtemplate = "NULL";
		    my $posteriori = 0;

		    for (my $r=1; $r<@{$TABprobs[$Qi]}; $r=$r+3){
			if ($TABprobs[$Qj][1] == $TABprobs[$Qi][$r]){
			    $probi = $TABprobs[$Qi][$r-1];	#probability of i-i'
			    $dtemplate = &DistanceInTemplate($atom1,$atom2, $Qi,$Qj,$TABprobs[$Qi][$r+1],$TABprobs[$Qj][2],$queryatomnameOne,$queryatomnameTwo,$queryresiduenameOne,$queryresiduenameTwo,$TABprobs[$Qj][1], 2);

			    $posteriori = $probi*$probj;
			
			    if ($Qi - $TABprobs[$Qi][$r+1] == $Qj - $TABprobs[$Qj][2]) {
				$posteriori = min($probi, $probj);
			    }
			    
			    if ( ($dtemplate ne "NULL")) {
				#	next if ($posteriori < 0.2);
				## skip too long distance-restraints
				next if (($dtemplate > 28 and $rsrType==9) or ($dtemplate > 19 and $rsrType==10) or ($dtemplate > 11 and $rsrType==23) or ($dtemplate>10 and $rsrType==26));    

				my $newline = "";
				my $hit = $templateList->get($TABprobs[$Qj][1])->get_Hit();

				if ($weightTemplates) {
				    my $tweight = $templateWeight{$hit};
				    $newline = &CalcRSRTweight($posteriori, $sim, $dtemplate, $rsrType, $atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC, $tweight); 
				} 
				else {
				    $newline = &CalcRSR($posteriori, $sim, $dtemplate, $rsrType, $atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC); 
				}

				$newline .= "  $hit $posteriori $sim $dtemplate $Qi $Qj $TABprobs[$Qi][$r+1] $TABprobs[$Qj][2]\n";
				push(@newrsrfile,$newline);
			    }
			    elsif ($remainRestr) { push(@newrsrfile, $rsrlines[$i]); }
			}
		    }
		}				

		#if there is more than one template for j-j' and one template for i-i'
		# Q  |---------------------------------------------------------------|
		# T1      |------------------------------------|
		# T2                          |-----------------------------------|
		#                   i                    j
		elsif (@{$TABprobs[$Qi]}==3 && @{$TABprobs[$Qj]}>3){
		    #print "more than one template for j-j'		atom:	$atom1	$atom2\n";					
		    my $probi = $TABprobs[$Qi][0];			#probability of i-i'
		    my $probj;
		    my $sim = $sim[$TABprobs[$Qi][1]];		#similarity of according template					
		    my $dtemplate = "NULL";
		    my $posteriori = 0;


		    for (my $r=1;$r<@{$TABprobs[$Qj]}; $r=$r+3){
			if ($TABprobs[$Qi][1] == $TABprobs[$Qj][$r]){
			    $probj = $TABprobs[$Qj][$r-1];	#probability of j-j'
			    $dtemplate = &DistanceInTemplate($atom1,$atom2, $Qi,$Qj,$TABprobs[$Qi][2],$TABprobs[$Qj][$r+1],$queryatomnameOne,$queryatomnameTwo,$queryresiduenameOne,$queryresiduenameTwo,$TABprobs[$Qi][1], 3);

			    $posteriori = $probi*$probj;
			
			    if ($Qi - $TABprobs[$Qi][2] == $Qj - $TABprobs[$Qj][$r+1]) {
				$posteriori = min($probi, $probj);
			    }

			    if ( ($dtemplate ne "NULL")) {
				#	next if ($posteriori < 0.2);
				## skip too long distance-restraints
				next if (($dtemplate > 28 and $rsrType==9) or ($dtemplate > 19 and $rsrType==10) or ($dtemplate > 11 and $rsrType==23) or ($dtemplate>10 and $rsrType==26));

				my $newline = "";
				my $hit = $templateList->get($TABprobs[$Qi][1])->get_Hit();

				if ($weightTemplates) {
				    my $tweight = $templateWeight{$hit};
				    $newline = &CalcRSRTweight($posteriori, $sim, $dtemplate, $rsrType,$atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC, $tweight);					 
				}
				else {
				    $newline = &CalcRSR($posteriori, $sim, $dtemplate, $rsrType,$atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC);				    
				}

				$newline .= " $hit $posteriori $sim $dtemplate $Qi $Qj $TABprobs[$Qi][2] $TABprobs[$Qj][$r+1]\n";

				push(@newrsrfile, $newline);
			    }
			    elsif ($remainRestr) { push(@newrsrfile, $rsrlines[$i]); }			    
			}
		    }				    		 
		}		
		#only one template (per position; but it is possible that the two residues are (one-fold) covered by two 
		# different (single covering) templates) => test if same template
		# Q  |--------------------------------------------------------|
		# T1   |----------------------|
		# T2                                  |-------------------|
		#                    i                      j
		elsif (@{$TABprobs[$Qi]}==3 && @{$TABprobs[$Qj]}==3) {
		    my $probi = $TABprobs[$Qi][0];			#probability of i-i'
		    my $probj = $TABprobs[$Qj][0];			#probability of j-j'
		    my $sim = $sim[$TABprobs[$Qi][1]];		#similarity of according template

		    my $dtemplate = "NULL";

		    # same template???
		    if ($TABprobs[$Qi][1] == $TABprobs[$Qj][1]) {

			$dtemplate = &DistanceInTemplate($atom1,$atom2, $Qi,$Qj,$TABprobs[$Qi][2],$TABprobs[$Qj][2],$queryatomnameOne,$queryatomnameTwo,$queryresiduenameOne,$queryresiduenameTwo,$TABprobs[$Qi][1], 4);
			#if ($atom1==19) {
			#    print "DistanceInTemplate($atom1,$atom2, $Qi,$Qj,$TABprobs[$Qi][2],$TABprobs[$Qj][2],$queryatomnameOne,$queryatomnameTwo,$queryresiduenameOne,$queryresiduenameTwo,$TABprobs[$Qi][1], 4)=$dtemplate\n";
			#}
		    }		  	

		    if ( ($dtemplate ne "NULL")) {
			my $posteriori = $probi*$probj;
			
			if ($Qi - $TABprobs[$Qi][2] == $Qj - $TABprobs[$Qj][2]) {
			    $posteriori = min($probi, $probj);
			}
			
		#	next if ($posteriori < 0.2);
					       
			## skip too long distance-restraints
			next if (($dtemplate > 28 and $rsrType==9) or ($dtemplate > 19 and $rsrType==10) or ($dtemplate > 11 and $rsrType==23) or ($dtemplate>10 and $rsrType==26));

			my  $newline = "";
			my $hit = $templateList->get($TABprobs[$Qi][1])->get_Hit();

			if ($weightTemplates) {
			    my $tweight = $templateWeight{$hit};
			    $newline = &CalcRSRTweight($posteriori, $sim, $dtemplate, $rsrType,$atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC, $tweight);		
			}
			else {
			    $newline = &CalcRSR($posteriori, $sim, $dtemplate, $rsrType,$atom1,$atom2, $mdnCACA, $mdnNO, $mdnSCMC, $mdnSCSC); 
			}

			$newline .= "  $hit $posteriori $sim $dtemplate $Qi $Qj $TABprobs[$Qi][2] $TABprobs[$Qj][2]\n";
			push(@newrsrfile,$newline);						
		    }
		    elsif ($remainRestr) { push(@newrsrfile, $rsrlines[$i]); }
		    
		}			
		else{print "WARNING1: error in line: $rsrlines[$i] check PIR-Alignment: $outbase.pir and Tab!\n";}								
	    }
	    else{
		print "WARNING2: in line: $rsrlines[$i] check tabfile(s), residues Qi=$Qi, Qj=$Qj; keep old restraint!\n";
		if ($remainRestr) {push(@newrsrfile, $rsrlines[$i])};
	    }
	}
	else {push(@newrsrfile,$rsrlines[$i]);}
    }

    $rsrfile =~ s/\.rsr$//;
    open (OUT, ">$rsrfile.new.rsr") or die ("Error: cannot open $rsrfile.new.rsr: $!\n");    
    print(OUT @newrsrfile);
    close(OUT);  

    return $templWeightNormalizer;
}



########################################################
## input: self explanatory, see arguments
## output: distance between specified atoms in specified
##         residue
########################################################
sub DistanceInTemplate {
    my $queryatomOne = $_[0];		#209
    my $queryatomTwo = $_[1];		#222
    my $queryresidueOne = $_[2];		#27
    my $queryresidueTwo = $_[3];		#29
    my $tempresidueOne = $_[4];		#4
    my $tempresidueTwo = $_[5];		#6
    my $queryatomnameOne = $_[6];		#CA	
    my $queryatomnameTwo = $_[7];		#CA	
    my $queryresiduenameOne = $_[8];		#ILE
    my $queryresiduenameTwo = $_[9];		#LEU
    my $templatenr = $_[10];			#0
    my $case = $_[11];
    #@PDBs; global	contains: Atomname[0] Residuename[1] ResidueNr[2] x[3] y[4] z[5]
    my $dtemplate;
    
    #print"search Query:$queryatomnameOne,$queryresiduenameOne & $queryatomnameTwo,$queryresiduenameTwo	Temp:$tempresidueOne $tempresidueTwo\n";
    #print"OPEN: $templatename.pdb search:$tempresidueOne\n";		
    
    #find xyz in template pdb
    my @coordinatesOne = ();
    my @coordinatesTwo = ();

    my $a = 0; 
    my $b = 0; 
    my $c = 0; 
    my $d = 0; 
    my $e = 0; 
    my $f = 0;

    if (defined @{$atomsForResidue[$templatenr][$tempresidueOne]}) {
	foreach my $atom (@{$atomsForResidue[$templatenr][$tempresidueOne]}) {
	    push(@coordinatesOne, [$PDBs[$templatenr][$atom][0],$PDBs[$templatenr][$atom][1],$PDBs[$templatenr][$atom][3],$PDBs[$templatenr][$atom][4],$PDBs[$templatenr][$atom][5]]);
	}
    }
    if (defined @{$atomsForResidue[$templatenr][$tempresidueTwo]}) {
	foreach my $atom (@{$atomsForResidue[$templatenr][$tempresidueTwo]}) {
	    push(@coordinatesTwo, [$PDBs[$templatenr][$atom][0],$PDBs[$templatenr][$atom][1],$PDBs[$templatenr][$atom][3],$PDBs[$templatenr][$atom][4],$PDBs[$templatenr][$atom][5]]);
	}
    }

    ## find all atoms for one template belonging to the template residues given, e.g. 27 and 29
    ## TODO: avoid traversing whole array by recording all atoms for each residue for each template
    # for(my $i=0; $i<@{$PDBs[$templatenr]}; $i++) {
# 	if (($PDBs[$templatenr][$i][2]) && ($PDBs[$templatenr][$i][2] == $tempresidueOne)) {
# 	    #print"$PDBs[$templatenr][$i][0]	$PDBs[$templatenr][$i][1]	$PDBs[$templatenr][$i][2]	$PDBs[$templatenr][$i][3]	$PDBs[$templatenr][$i][4]	$PDBs[$templatenr][$i][5]\n";
# 	    push(@coordinatesOne,[$PDBs[$templatenr][$i][0],$PDBs[$templatenr][$i][1],$PDBs[$templatenr][$i][3],$PDBs[$templatenr][$i][4],$PDBs[$templatenr][$i][5]]);
# 	}
# 	if (($PDBs[$templatenr][$i][2])&&($PDBs[$templatenr][$i][2] == $tempresidueTwo)) {
# 	    #print"$PDBs[$templatenr][$i][0]	$PDBs[$templatenr][$i][1]	$PDBs[$templatenr][$i][2]	$PDBs[$templatenr][$i][3]	$PDBs[$templatenr][$i][4]	$PDBs[$templatenr][$i][5]\n";
# 	    push(@coordinatesTwo,[$PDBs[$templatenr][$i][0],$PDBs[$templatenr][$i][1],$PDBs[$templatenr][$i][3],$PDBs[$templatenr][$i][4],$PDBs[$templatenr][$i][5]]);
# 	}    	
#     }    

#    foreach my $atom {$atomsForResidue[$templatenr]}{$tempresidueOne}

    my $foundOne = 0;
    my $foundTwo = 0;

    for(my $k=0; $k<@coordinatesOne; $k++) {
	#print "$coordinatesOne[$k][0] eq $queryatomnameOne?\n";
	if ($coordinatesOne[$k][0] eq $queryatomnameOne) {
	    $a = $coordinatesOne[$k][2];
	    $b = $coordinatesOne[$k][3];
	    $c = $coordinatesOne[$k][4];
	    $foundOne = 1;    
	    last;
	}
    }    

    for(my $k=0; $k<@coordinatesTwo; $k++) {
	#print "$coordinatesTwo[$k][0] eq $queryatomnameTwo?\n";
	if ($coordinatesTwo[$k][0] eq $queryatomnameTwo) {
	    $d = $coordinatesTwo[$k][2];
	    $e = $coordinatesTwo[$k][3];
	    $f = $coordinatesTwo[$k][4];
	    $foundTwo = 1;    		
	    last;
	}
    }    


    if($foundOne==0 || $foundTwo==0) { 		
 	if (@coordinatesOne) {
	    if (@coordinatesTwo) {
		if($foundOne==0) {
		    my @alternatives = &searchAlternatives($queryresiduenameOne,$queryatomnameOne);
		    #print "QUERY1:$queryatomnameOne	in $queryresiduenameOne	TEMPLATE:\t";
		    #print"\nsearchAlternatives 1	zu $queryresiduenameOne $queryatomnameOne ($queryatomOne $queryatomTwo) found:\n";
		    for(my $k=0; $k<@coordinatesOne; $k++) {
			#print "$coordinatesOne[$k][0] in $coordinatesOne[$k][1]\nALTERNATIVE?\n";
			for (my $w=1; $w<@alternatives; $w++) {
			    #print "$alternatives[$w]\n";
			    if ($alternatives[$w] =~ /\s$coordinatesOne[$k][1]:$coordinatesOne[$k][0]\s/) {
				#print "HIT!!!: $coordinatesOne[$k][2]	$coordinatesOne[$k][3]	$coordinatesOne[$k][4]";
				$a += $coordinatesOne[$k][2];
				$b += $coordinatesOne[$k][3];
				$c += $coordinatesOne[$k][4];
				$foundOne++;								
			    }
			}
		    }
		}
		if($foundTwo==0) {
		    my @alternatives = &searchAlternatives($queryresiduenameTwo,$queryatomnameTwo);
		    #print "\nQUERY2:$queryatomnameTwo	in $queryresiduenameTwo	TEMPLATE:\t";  
		    #print"\nsearchAlternatives 2	zu $queryresiduenameTwo $queryatomnameTwo ($queryatomOne $queryatomTwo)	found:\n";
		    for(my $k=0; $k<@coordinatesTwo; $k++) {
			#print "$coordinatesTwo[$k][0] in $coordinatesTwo[$k][1]\nALTERNATIVE?\n";
			for (my $w=1; $w<@alternatives; $w++) {
			    #print "$alternatives[$w]\n";
			    if ($alternatives[$w] =~ /\s$coordinatesTwo[$k][1]:$coordinatesTwo[$k][0]\s/){
				#print "HIT!!!: $coordinatesTwo[$k][2]	$coordinatesTwo[$k][3]	$coordinatesTwo[$k][4]";
				$d += $coordinatesTwo[$k][2];
				$e += $coordinatesTwo[$k][3];
				$f += $coordinatesTwo[$k][4];
				$foundTwo++;								
			    }
			}
		    }
		} 
		if ($foundOne>0 && $foundTwo>0) {
		    $a /= $foundOne;
		    $b /= $foundOne;
		    $c /= $foundOne;
		    $d /= $foundTwo;
		    $e /= $foundTwo;
		    $f /= $foundTwo;					
		    $dtemplate= sqrt(($a-$d)*($a-$d) + ($b-$e)*($b-$e) + ($c-$f)*($c-$f));
		    print"FOUND ALTERNATIVE distance: $dtemplate\n"; 					  
		}
		else{$dtemplate = "NULL";}					    
	    }
	    else{$dtemplate = "NULL";}		  
	}
	else{$dtemplate = "NULL";}    	

    }
    else {
	$dtemplate= sqrt(($a-$d)*($a-$d) + ($b-$e)*($b-$e) + ($c-$f)*($c-$f));   
	if ($dtemplate < 1) { 
	    print "==> Distance < 1 : $dtemplate\n"; 
	    print "queryatomOne = $queryatomOne\n";
	    print "queryatomTwo = $queryatomTwo\n";
	    print "queryresidueOne = $queryresidueOne\n";
	    print "queryresidueTwo = $queryresidueTwo\n";
	    print "tempresidueOne = $tempresidueOne\n";
	    print "tempresidueTwo = $tempresidueTwo\n";
	    print "queryatomnameOne = $queryatomnameOne\n";
	    print "queryatomnameTwo = $queryatomnameTwo\n";
	    print "queryresidue = $queryresiduenameOne\n";
	    print "queryresiduenameTwo = $queryresiduenameTwo\n";
	    print "templatenr = $templatenr\n";
	    print "case = $case\n";
	}
    } 
    #print "DELTA: $dtemplate\n";
    return $dtemplate;
}


##################################################################
## input: self explanatory, see arguments
## output: constraint line as needed for modeller constraint file
##
## description:
## the parameters for the new constraint are calculated by a MDN
## parameters with index 2 are background parameters
## the new restraint is generated and returned
##################################################################
sub CalcRSR{	
    my $posteriori = shift;
    my $sim = shift;
    my $dtemplate = shift;
    my $rsrType	 = shift;
    my $atom1 = shift;
    my $atom2 = shift;
    my $mdnCACA = shift;
    my $mdnNO = shift;
    my $mdnSCMC = shift;
    my $mdnSCSC = shift;
        

    my @mdnOutput;
    my @mdnInput = ($dtemplate, $posteriori, $sim);

    if ($rsrType == 9) {
	@mdnOutput = $mdnCACA->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 10) {
	@mdnOutput = $mdnNO->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 23) {
	@mdnOutput = $mdnSCMC->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 26) {
	@mdnOutput = $mdnSCSC->getMixtureParametersFromMDN(@mdnInput);
    }

    my $p1 = $mdnOutput[0];
    my $p2 = $mdnOutput[1];
    my $mu1 = log($dtemplate) + $mdnOutput[2];
    my $mu2 = log($dtemplate) + $mdnOutput[3];
    my $sigma1 = $mdnOutput[4];
    my $sigma2 = $mdnOutput[5];

    #      Form	 Modality	Feature	       Group	 Numb_atoms	Numb_parameters		Numb_Feat	Atom_indices	Parameter(e.g. mean, standard deviation)	
    #R    	3    	1         	1      	  	1      	     2	   		2   		1     		3     2       	1.5000    0.0364 
    #R     50	1		1	    9|10|23|26	     2			6		1		..    .. 	...	  ...     ...    ...    ...   ...
    my $newline = sprintf("R 50 1 1 %3.0f 2 6 1 %5.0f %5.0f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f",$rsrType,$atom1,$atom2,$p1,$mu1,$sigma1,$p2,$mu2,$sigma2);			
    
   
    return $newline;	
}


##################################################################
## input: self explanatory, see arguments
## output: constraint line as needed for modeller constraint file
##
## description:
## the parameters for the new constraint, are calculated by a MDN
## parameters with index 2 are background parameters
##################################################################
sub CalcRSRTweight{	
    my $posteriori = shift;
    my $sim = shift;
    my $dtemplate = shift;
    my $rsrType	 = shift;
    my $atom1 = shift;
    my $atom2 = shift;
    my $mdnCACA = shift;
    my $mdnNO = shift;
    my $mdnSCMC = shift;
    my $mdnSCSC = shift;
    my $tweight = shift;
        

    my @mdnOutput;
    my @mdnInput = ($dtemplate, $posteriori, $sim);

    if ($rsrType == 9) {
	@mdnOutput = $mdnCACA->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 10) {
	@mdnOutput = $mdnNO->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 23) {
	@mdnOutput = $mdnSCMC->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 26) {
	@mdnOutput = $mdnSCSC->getMixtureParametersFromMDN(@mdnInput);
    }

    my $p1 = $mdnOutput[0];
    my $p2 = $mdnOutput[1];
    my $mu1 = log($dtemplate) + $mdnOutput[2];
    my $mu2 = log($dtemplate) + $mdnOutput[3];
    my $sigma1 = $mdnOutput[4];
    my $sigma2 = $mdnOutput[5];

    #      Form	 Modality	Feature	       Group	 Numb_atoms	Numb_parameters		Numb_Feat	Atom_indices	Parameter(e.g. mean, standard deviation)	
    #R     3       1         	1      	  	1      	     2	   		2   		1     		3     2       	1.5000    0.0364 
    #R     50	1		1	    9|10|23|26	     2			6		1		..    .. 	...	  ...     ...    ...    ...   ...
    my $newline = sprintf("R 50 1 1 %3.0f 2 7 1 %5.0f %5.0f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f",$rsrType,$atom1,$atom2,$p1,$mu1,$sigma1,$p2,$mu2,$sigma2, $tweight);
    
   
    return $newline;	
}


##################################################################
## input: self explanatory, see arguments
## output: constraint line as needed for modeller constraint file
##
## description:
## the parameters for the new constraint, are calculated by a MDN
## instead of writing the parameters mu, mu_BG, sigma, sigma_BG,
## w, w_BG into the restraints file, use precalculated parameters
## a,b,c in order to save computing steps when Modeller optimization
## runs
## parameters with index 2 are background parameters
##################################################################
sub CalcRSRefficient {	
    my $posteriori = shift;
    my $sim = shift;
    my $dtemplate = shift;
    my $rsrType	 = shift;
    my $atom1 = shift;
    my $atom2 = shift;
    my $mdnCACA = shift;
    my $mdnNO = shift;
    my $mdnSCMC = shift;
    my $mdnSCSC = shift;
        

    my @mdnOutput;
    my @mdnInput = ($dtemplate, $posteriori, $sim);

    if ($rsrType == 9) {
	@mdnOutput = $mdnCACA->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 10) {
	@mdnOutput = $mdnNO->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 23) {
	@mdnOutput = $mdnSCMC->getMixtureParametersFromMDN(@mdnInput);
    }
    elsif ($rsrType == 26) {
	@mdnOutput = $mdnSCSC->getMixtureParametersFromMDN(@mdnInput);
    }

    my $p1 = $mdnOutput[0];
    my $p2 = $mdnOutput[1];
    my $mu1 = log($dtemplate) + $mdnOutput[2];
    my $mu2 = log($dtemplate) + $mdnOutput[3];
    my $sigma1 = $mdnOutput[4];
    my $sigma2 = $mdnOutput[5];

    #      Form	 Modality	Feature	       Group	 Numb_atoms	Numb_parameters		Numb_Feat	Atom_indices	Parameter(e.g. mean, standard deviation)	
    #R    	3    	1         	1      	  	1      	     2	   		2   		1     		3     2       	1.5000    0.0364 
    #R     50	1		1	    9|10|23|26	     2			6		1		..    .. 	...	  ...     ...    ...    ...   ...
    my $sigma1sq = $sigma1*$sigma1;
    my $sigma2sq = $sigma2*$sigma2;


    my $factorOne = 2*$sigma1sq*$sigma2sq;

    my $a = $p1/$p2 * $sigma2/$sigma1;
    my $b = ($sigma2sq - $sigma1sq)/$factorOne;

    my $factorTwo = ($sigma2sq*$mu1 - $sigma1sq*$mu2);
    my $c = $factorTwo/($sigma2sq - $sigma1sq);

    my $d = $factorTwo*$factorTwo/($factorOne*($sigma2sq - $sigma1sq)) - $factorTwo/$factorOne;

    my $newline = sprintf("R 50 1 1 %3.0f 2 6 1 %5.0f %5.0f %8.4f %8.4f %8.4f %8.4f\n",$rsrType,$atom1,$atom2,$a, $b, $c, $d);			
    
    return $newline;	
}


sub CalcRSRPrev{	
    my $probi = $_[0];
    my $probj = $_[1];
    my $sim = $_[2];
    my $dtemplate = $_[3];
    my $rsrType	 = $_[4];
    my $atom1 = $_[5];
    my $atom2 = $_[6];
    
    #distance Template < cutoff
    #         	avPI   	avSIM  	avDELTA
    # CACAlog	0.78950 0.47877 2.29587
    # NOlog 	0.79391 0.47610 2.04427
    # SCMClog	0.81810 0.50638 1.50130
    # SCSClog	0.87671 0.59417 1.46639  
    my @Avec;	#Parametervalues: a0 a1 a2 a3 	b0 b1 b2 b3 b4 b5 b6 b7 b8 b9	c0 c1 c2 c3 c4 c5	d0 d1 d2 d3 d4 d5
    my $PI;															
    my $Sim;
    my $deltatemplate;
    
    if ($rsrType==9){
	@Avec = (-2.6961,-0.7186,0.1326 ,-0.6153,0.6712 ,0.7634 ,0.1174 ,0.0353 ,0.3422 ,-0.2807,0.2447 ,0.0463 ,-0.0171,0.1796 ,-1.0388,-0.1192,-0.3310,0.0900 ,0.2468 ,-0.2058,0.1296 ,-0.1788,-0.1989,0.0080 ,0.3942 ,0.1421);
	$PI= $probi*$probj - 0.7895;
	$Sim = $sim - 0.4788;
	$deltatemplate = log($dtemplate)- 2.29587;
    }
    elsif ($rsrType==10){
	@Avec = (-2.6514,-0.5396,-0.0368,-0.5819,0.6323 ,0.8551 ,0.0992 ,-0.0011,0.4606 ,0.0202 ,0.2238 ,0.1218 ,-0.0396,0.0448 ,-0.9667,-0.2518,-0.4197,0.0653 ,0.2507 ,-0.2392,0.1500 ,-0.0833,-0.2595,-0.0224,0.3706 ,0.1158);	
	$PI= $probi*$probj - 0.7940;
	$Sim = $sim - 0.4763;
	$deltatemplate = log($dtemplate)- 2.04427;				
    }
    elsif($rsrType==23){
	@Avec = (-3.4835,-0.1144,-0.6024,-0.4093,0.5841 ,0.4945 ,0.0102 ,-0.4040,0.3905 ,0.1852 ,0.0844 ,0.0213 ,0.0668 ,-0.3501,-1.6170,-0.8439,-0.2171,0.0006 ,-0.1970,-0.9967,0.0843 ,-0.2868,-0.0493,-0.0108,0.0096 ,0.3162);
	$PI= $probi*$probj - 0.8893;
	$Sim = $sim - 0.6281;
	$deltatemplate = log($dtemplate)- 1.50130;	
    }
    elsif($rsrType==26){
	@Avec = (-2.3377,-0.3537,1.9132 ,-0.1555,0.5179 ,0.6930 ,0.1961 ,0.9687 ,0.3331 ,-0.0203,0.6672 ,-0.1298,-0.0680,1.0575 ,-0.9797,0.2624 ,-0.5030,0.1753 ,1.0932 ,1.6404 ,0.3804 ,-0.2516,-0.5884,-0.0289,0.5928 ,0.3932);
	$PI= $probi*$probj - 0.9373;
	$Sim = $sim - 0.7220;
	$deltatemplate = log($dtemplate)- 1.46639;		
    }
    my $A = $Avec[0] + $Avec[1] *$Sim + $Avec[2] *$deltatemplate;
    my $B = $Avec[3];
    my $logsigma = $A + $B*$PI;
    
    my $p = $Avec[4]+$Avec[5]*$PI+$Avec[6]*$Sim+$Avec[7]*$deltatemplate+$Avec[8]*$PI*$PI+$Avec[9]*$PI*$Sim+$Avec[10]*$PI*$deltatemplate+$Avec[11]*$Sim*$deltatemplate+$Avec[12]*$Sim*$Sim+$Avec[13]*$deltatemplate*$deltatemplate;
    if ($p>1) {$p=1;}elsif ($p<0) {$p=0;}		
    
    my $logsigmaBG = $Avec[14] + $Avec[15]*$deltatemplate + $Avec[16]*$PI + $Avec[17]*$Sim  + $Avec[18]*$PI*$deltatemplate + $Avec[19]*$deltatemplate*$deltatemplate;
    my $logmuBG = $Avec[20] + $Avec[21]*$deltatemplate + $Avec[22]*$PI + $Avec[23]*$Sim  + $Avec[24]*$PI*$deltatemplate + $Avec[25]*$deltatemplate*$deltatemplate;
    
    #for Modeller Restraints:
    my $p1 = sprintf("%.5f" ,$p);
    my $p2 = sprintf("%.5f" ,(1-$p1));

    my $sigma = exp($logsigma);			
    $sigma = sprintf("%.5f" ,$sigma);		

    my $sigmaBG = exp($logsigmaBG);		
    $sigmaBG = sprintf("%.5f" ,$sigmaBG);				

    my $mu = log($dtemplate);			
    $mu = sprintf("%.5f" ,$mu);

    my $muBG = $logmuBG + $mu;		  	
    $muBG = sprintf("%.5f" ,$muBG);					
    
    #      Form	 Modality	Feature		Group		Numb_atoms	Numb_parameters		Numb_Feat	Atom_indices	Parameter(e.g. mean, standard deviation)	
    #R    	3    	1         	1      	  1      		2   		2   			1     		3     2       	1.5000    0.0364 
    #R		50		1			1	  9|10|23|26	    2			6				1			..    .. 		...		  ...     ...    ...    ...   ...
    my $newline = sprintf("R 50 1 1 %3.0f 2 6 1 %5.0f %5.0f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",$rsrType,$atom1,$atom2,$p1,$mu,$sigma,$p2,$muBG,$sigmaBG);			
    
    return $newline;	
}

sub searchAlternatives{
    #alternatives from /cluster/bioprogs/modeller/modlib/atmcls-melo.lib
    my $Residue = $_[0];
    my $Atom = $_[1];
    #print "\s$Residue:$Atom\s\n";

    my @ATMGRP01 = 	(' ALA:CB ', ' ILE:CG2 ', ' ILE:CD1 ', ' ILE:CD ', ' LEU:CD1 ', ' LEU:CD2 ', ' THR:CG2 ', ' VAL:CG1 ', ' VAL:CG2 ');						
    my @ATMGRP02 = 	(' ILE:CB ', ' LEU:CG ', ' VAL:CB ');						
    my @ATMGRP03 = 	(' ARG:CB ', ' ARG:CG ', ' ASN:CB ', ' ASP:CB ', ' GLU:CB ', ' GLU:CG ', ' ILE:CG1 ', ' GLN:CB ', ' GLN:CG ', ' HIS:CB ', ' HSD:CB ', ' LEU:CB ', ' LYS:CB ', ' LYS:CG ', ' LYS:CD ', ' MET:CB ', ' PHE:CB ', ' PRO:CB ', ' PRO:CG ', ' TRP:CB ', ' TYR:CB ');
    my @ATMGRP04 = 	(' PHE:CG ', ' TRP:CD2 ', ' TYR:CG ');
    my @ATMGRP05 = 	(' PHE:CD1 ', ' PHE:CD2 ', ' PHE:CE1 ', ' PHE:CE2 ', ' PHE:CZ ', ' TRP:CE3 ', ' TRP:CZ2 ', ' TRP:CZ3 ', ' TRP:CH2 ', ' TYR:CD1 ', ' TYR:CD2 ', ' TYR:CE1 ', ' TYR:CE2 ');						
    my @ATMGRP06 = 	(' SER:OG ', ' THR:OG1 ');						
    my @ATMGRP07 = 	(' ASN:ND2 ', ' GLN:NE2 ');						
    my @ATMGRP08 = 	(' CYS:SG ', ' CSS:SG ');						
    my @ATMGRP09 = 	(' ARG:NH1 ', ' ARG:NH2 ');
    my @ATMGRP10 = 	(' HIS:CG ', ' HSD:CG ');						
    my @ATMGRP11 = 	(' HIS:CD2 ', ' HSD:CD2 ', ' TRP:CD1 ');						
    my @ATMGRP12 = 	(' HIS:NE2 ', ' HSD:NE2 ');						
    my @ATMGRP13 = 	(' HIS:CE1 ', ' HSD:CE1 ');						
    my @ATMGRP14 = 	(' ASP:CG ', ' GLU:CD ');						
    my @ATMGRP15 = 	(' ALA:OXT ', ' ARG:OXT ', ' ASN:OXT ', ' ASP:OXT ', ' CYS:OXT ', ' CSS:OXT ', ' GLU:OXT ', ' GLN:OXT ', ' GLY:OXT ', ' HIS:OXT ', ' HSD:OXT ', ' ILE:OXT ', ' LEU:OXT ', ' LYS:OXT ', ' MET:OXT ', ' PHE:OXT ', ' PRO:OXT ', ' SER:OXT ', ' THR:OXT ', ' TRP:OXT ', ' TYR:OXT ', ' VAL:OXT ', ' ASP:OD1 ', ' ASP:OD2 ', ' ASP:OT1 ', ' ASP:OT2 ', ' GLU:OE1 ', ' GLU:OE2 ', ' GLU:OT1 ', ' GLU:OT2 ');
    my @ATMGRP16 = 	(' CYS:CB ', ' CSS:CB ', ' MET:CG ');						
    my @ATMGRP17 = 	(' ASN:CG ', ' GLN:CD ');						
    my @ATMGRP18 = 	(' ASN:OD1 ', ' GLN:OE1 ');						
    my @ATMGRP19 = 	(' HIS:ND1 ', ' HSD:ND1 ');
    
    my @alternatives;
    if (grep(/\s$Residue:$Atom\s/,@ATMGRP01)){@alternatives=@ATMGRP01;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP02)){@alternatives=@ATMGRP02;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP03)){@alternatives=@ATMGRP03;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP04)){@alternatives=@ATMGRP04;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP05)){@alternatives=@ATMGRP05;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP06)){@alternatives=@ATMGRP06;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP07)){@alternatives=@ATMGRP07;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP08)){@alternatives=@ATMGRP08;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP09)){@alternatives=@ATMGRP09;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP10)){@alternatives=@ATMGRP10;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP11)){@alternatives=@ATMGRP11;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP12)){@alternatives=@ATMGRP12;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP13)){@alternatives=@ATMGRP13;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP14)){@alternatives=@ATMGRP14;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP15)){@alternatives=@ATMGRP15;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP16)){@alternatives=@ATMGRP16;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP17)){@alternatives=@ATMGRP17;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP18)){@alternatives=@ATMGRP18;}
    elsif (grep(/\s$Residue:$Atom\s/,@ATMGRP19)){@alternatives=@ATMGRP19;}
    
    return @alternatives;
}


sub CalcWnew {
    my $d = shift;
    my $p1 = shift;
    my $mu1 = shift;
    my $sigma1 = shift;
    my $p2 = shift;
    my $mu2 = shift;
    my $sigma2 = shift;

    my $correct = $p1 / $sigma1 * exp( -0.5 * (($d - $mu1)/$sigma1)*(($d - $mu1)/$sigma1) );
    my $bg = $p2 / $sigma2 * exp(-0.5 * (($d - $mu2)/$sigma2)*(($d - $mu2)/$sigma2) );

    return $correct/($correct + $bg);
}


########################################################################
## input: rsrfile: name of restraints file (full path)
##        models: pdb models generated in round before (full path)
##
## output: generates a new restraints file called:
##         OLDRESTRAINTSFILE.iter
##         this file is then to be used for another round of Modeller
##         model generation
##
## description: assume that a previous Modeller run has generated
## e.g. 4 models. Now, one wants to use the information contained in
## the generated models. 
## For all models, there is a single rsr file. The generated pdb-files
## (models) are "realizations" of these contraints, e.g. take a distance
## restraint for CACA between residues 4 and 39 in a target. 
## Each model has now a concrete value for the distance between 4 and 39
## Interpret now w as:
## w = P(correct_iji'j' | d_ij, d'_i'j', A_ij) = wN(log(d_ij)|..) / (wN(log(d_ij)|..) + (1-w)N(log(d'_ij)|..) 
## and calculate w for each model.
## Take as new w the median of all these w.
########################################################################
sub UpdateRsrFromModels {
    my $rsrfile = shift;
    my @models = @_;


    my @pdbs;

    ## copy pdb files and rename the copied files, because next modeller round will
    ## generate new files named as the old ones
    for (my $i=0; $i<@models; $i++) {
	system("cp $models[$i] $models[$i].prev");
	$models[$i] .= ".prev";
    }

    ## read in models
    ## atomnr -> x,y,z
    for (my $i=0; $i<@models; $i++) {

	my %structureInfo;

	if (not defined( open (MH, "< $models[$i]"))) {
	    print "Cant open $models[$i]. Skipping $models[$i]...\n";
	    next;
	}
	while(<MH>) {
	    if (/ATOM.{2}\s*(\d+)\s.([\S|\s]{3}).(\S{3})\s.\s*(\d+)[\s|\D]\s*(-?\d+.\d+)\s*(-?\d+.\d+)\s*(-?\d+.\d+)/) {
		##          $atomnr   x   y   z
		$structureInfo{$1} = [$5, $6, $7];
	    }
	}
	close(MH);

	## save atomnr->[x,y,z] hash for each model
	$pdbs[$i] = \%structureInfo;
    }

    open (RH, "< $rsrfile") or die "Cant open $rsrfile: $!\n";

    my @newRestraints;

    ## go through restraints and for each (CACA,NO,SCMC,SCSC distance) restraint calculate the acutal distance
    ## in each model and calculate new w for the restraint r. Replace r by new restraint (new w and 1-w)
    while (my $rsr = <RH>) {
	##  50 1 1   9 2 6 1   289   305   0.2070   1.8664   0.1161   0.7930   2.0669   0.4316
	##                                                 atom1   atom2   p1      mu1     sigma1  p2     mu2     sigma2
	if ($rsr =~ /^R\s+50\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
	    my $rsrtype = $1;
	    my $atom1 = $2;
	    my $atom2 = $3;
	    my $p1 = $4;
	    my $mu1 = $5;
	    my $sigma1 = $6;
	    my $p2 = $7;
	    my $mu2 = $8;
	    my $sigma2 = $9;


	    my @allWs;
	    my $finalNewW;

	    for (my $i=0; $i<@models; $i++) {

		## coordinates of atom1 and atom2
		my @coAtom1 = @{$pdbs[$i]{$atom1}};
		my @coAtom2 = @{$pdbs[$i]{$atom2}};

		my $distance = sqrt( ($coAtom1[0]-$coAtom2[0])*($coAtom1[0]-$coAtom2[0]) +($coAtom1[1]-$coAtom2[1])*($coAtom1[1]-$coAtom2[1]) +($coAtom1[2]-$coAtom2[2])*($coAtom1[2]-$coAtom2[2]) );

		my $wUpdated = &CalcWnew(log($distance), $p1, $mu1, $sigma1, $p2, $mu2, $sigma2);

		push(@allWs, $wUpdated);
	    }

	    ## calculate the median of all new w:
	    @allWs = sort {$a <=> $b} @allWs;

	    
	    if (@allWs % 2 == 0 and @allWs >= 2) {
		$finalNewW = ($allWs[@allWs/2-1] + $allWs[@allWs/2]) / 2;
	    }
	    else {
		$finalNewW = $allWs[@allWs/2];
	    }


	    ## update restraint and add it to the new restraints list
	    my $newrsr = sprintf("R 50 1 1 %3.0f 2 6 1 %5.0f %5.0f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", $rsrtype, $atom1, $atom2, $finalNewW, $mu1, $sigma1, (1-$finalNewW), $mu2, $sigma2); 

	    push(@newRestraints, $newrsr);
	}
	else {
	    ## keep other restraints
	    push (@newRestraints, $rsr);
	}
    }

    close(RH);

    ## save new restraints in a new file
    my $newRsrFileName = $rsrfile . ".iter";

    open (IH, "> $newRsrFileName") or die ("Cant write $newRsrFileName: $!\n");
    print (IH @newRestraints);
    close(IH);

}


sub GenerateIterativeModels {

    my $myrsr = shift;
    my $queryName = shift;
    my $tmpdir = shift;
    my $inbase = shift;
    my $outbase = shift;


    &System("cp $inbase.B9999000*.pdb $tmpdir");
    &System("mv $tmpdir/$queryName.B99990003.pdb $tmpdir/$queryName.B99990004.pdb");
    &System("mv $tmpdir/$queryName.B99990002.pdb $tmpdir/$queryName.B99990003.pdb");
    &System("mv $tmpdir/$queryName.B99990001.pdb $tmpdir/$queryName.B99990002.pdb");

    &ModellerNew($myrsr, $queryName, $inbase, $outbase, $tmpdir, 1);

    my @models = <$tmpdir/$queryName.B9999000*.pdb>;


    &UpdateRsrFromModels($myrsr, @models);
    &ModellerNew("$myrsr.iter", $queryName, $inbase, $outbase, $tmpdir, 3);

}


##########
## write out "final" pdb file in CASP compatible format
## 
## outFile: filename of pdb file to write
## srcFile: filename of already existing (Modeller generated) pdb file
##
sub WriteCASPpdbFile {
    my $outFile = shift @_;
    my $srcFile = shift @_;
    my $seqname = shift @_;
    my $server = shift @_;
    my $tList = shift @_;

    my $templates = "";
    my $predTMs = "";
    for (my $i=0; $i<$tList->size(); $i++) {
	$templates .= $tList->get($i)->get_Hit() . " ";
	$predTMs .= sprintf("%.2f", $tList->get($i)->get_predTM()) . " ";
    }

    open(SRC, "< $srcFile") or die ("WriteCASPpdbFile: Cant open $srcFile: $!\n");
    my @mod = <SRC>;
    close(SRC);

    my @coordinates;

    ## skip all lines until first ATOM record is read
    for (my $i=0; $i<@mod; $i++) {
	if ($mod[$i] =~ /^ATOM/) {
	    push (@coordinates, $mod[$i]);
	}
    }

    open(PDB, ">$outFile") or die ("WriteCASPpdbFile: Cant open $outFile: $!\n");
    print(PDB "PFRMAT TS\n");
    print(PDB "TARGET $seqname\n");
    print(PDB "AUTHOR $server\n");
    print(PDB "METHOD $server (http://hhpred.tuebingen.mpg.de)\n");
    print(PDB "MODEL 1\n");
    print(PDB "PARENT $templates\n");
    print(PDB "REMARK PredTM $predTMs\n");
    print(PDB @coordinates);
    print(PDB "TER\n");
    print(PDB "END\n");
    close(PDB);

    FixOccupancy::FixOccupancy($outFile);
}

1;
