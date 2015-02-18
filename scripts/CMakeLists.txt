
INSTALL(PROGRAMS
	Align.pm  
	HHPaths.pm  
	PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
	DESTINATION scripts
)


INSTALL(PROGRAMS 
	addss.pl  
	checkA3M.pl	  
	create_profile_from_hhm.pl	
	hhmakemodel.pl	
	hhsuitedb.pl  
	mergeali.pl
	multithread.pl  
	pdbfilter.pl  
	renumberpdb.pl
	create_profile_from_hmmer.pl	
	pdb2fasta.pl    
	reformat.pl   
	splitfasta.pl
	DESTINATION scripts
)
