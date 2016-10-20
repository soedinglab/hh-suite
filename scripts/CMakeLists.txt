
INSTALL(PROGRAMS
	Align.pm  
	HHPaths.pm  
	PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
	DESTINATION scripts
)


INSTALL(PROGRAMS 
	addss.pl  
	create_profile_from_hhm.pl	
	hhmakemodel.pl	
	hhsuitedb.py  
	mergeali.pl
	multithread.pl  
	pdbfilter.pl  
	renumberpdb.pl
	create_profile_from_hmmer.pl	
	pdb2fasta.pl    
	reformat.pl   
	splitfasta.pl
	check_a3m.py	  
	ffindex.py
	a3m.py
	get_a3m_size.py
	pdbfilter.py
	cif2fasta.py
	hhmakemodel.py
  hh_reader.py
	DESTINATION scripts
)
