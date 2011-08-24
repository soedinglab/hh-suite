// hhfunc.C


// Calculate secondary structure prediction with PSIpred
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile)
{
  // Initialize
  std::string command;
  char line[LINELEN]=""; 
  char filename[NAMELEN];

  strcpy(ss_pred,"-");
  strcpy(ss_conf,"-");
  
  // Run PSIpred
  
  // Check for Psipred ver >= 3.0 (weights.dat4 doesn't exists an more)
  strcpy(filename,par.psipred_data);
  strcat(filename,"/weights.dat4");
  FILE* check_exists = fopen(filename,"r");
  if (check_exists) {  // Psipred version < 3.0
    command = (std::string)par.psipred + "/psipred " + (std::string)tmpfile + ".mtx " + (std::string)par.psipred_data + "/weights.dat " + (std::string)par.psipred_data + "/weights.dat2 " + (std::string)par.psipred_data + "/weights.dat3 " + (std::string)par.psipred_data + "/weights.dat4 > " + (std::string)tmpfile + ".ss";
  } else {
    command = (std::string)par.psipred + "/psipred " + (std::string)tmpfile + ".mtx " + (std::string)par.psipred_data + "/weights.dat " + (std::string)par.psipred_data + "/weights.dat2 " + (std::string)par.psipred_data + "/weights.dat3 > " + (std::string)tmpfile + ".ss";
  }
  runSystem(command,v);
  command = (std::string)par.psipred + "/psipass2 " + (std::string)par.psipred_data + "/weights_p2.dat 1 0.98 1.09 " + (std::string)tmpfile + ".ss2 " + (std::string)tmpfile + ".ss > " + (std::string)tmpfile + ".horiz";
  runSystem(command,v);

  // Read results
  strcpy(filename,tmpfile);
  strcat(filename,".horiz");
  FILE* horizf = fopen(filename,"r");
  if (!horizf) return;

  while (fgets(line,LINELEN,horizf))
    {
      char tmp_seq[NAMELEN]="";
      char* ptr=line;
      if (!strncmp(line,"Conf:",5))
	{
	  ptr+=5;
	  strwrd(tmp_seq,ptr);
	  strcat(ss_conf,tmp_seq);
	}
      if (!strncmp(line,"Pred:",5))
	{
	  ptr+=5;
	  strwrd(tmp_seq,ptr);
	  strcat(ss_pred,tmp_seq);
	}
    }
  fclose(horizf);

  if (v>3)
    {
      printf("SS-pred: %s\n",ss_pred);
      printf("SS-conf: %s\n",ss_conf);
    }
}

// Calculate secondary structure for given HMM and return prediction
void CalculateSS(HMM& q, char *ss_pred, char *ss_conf)
{
  char tmpfile[]="/tmp/hhCalcSSXXXXXX";
  if (mkstemp(tmpfile) == -1) {
    cerr << "ERROR! Could not create tmp-file!\n"; 
    exit(4);
  }
  
  // Write log-odds matrix from q to tmpfile.mtx
  char filename[NAMELEN];
  FILE* mtxf = NULL;

  strcpy(filename,tmpfile);
  strcat(filename,".mtx");
  mtxf = fopen(filename,"w");
  if (!mtxf) OpenFileError(filename);

  fprintf(mtxf,"%i\n",q.L);
  fprintf(mtxf,"%s\n",q.seq[q.nfirst]+1);
  fprintf(mtxf,"2.670000e-03\n4.100000e-02\n-3.194183e+00\n1.400000e-01\n2.670000e-03\n4.420198e-02\n-3.118986e+00\n1.400000e-01\n3.176060e-03\n1.339561e-01\n-2.010243e+00\n4.012145e-01\n");
  
  for (int i = 1; i <= q.L; ++i) 
    {
      fprintf(mtxf,"-32768 ");
      for (int a = 0; a < 20; ++a)
	{
	  int tmp = iround(50*flog2(q.p[i][s2a[a]]/pb[s2a[a]]));
	  fprintf(mtxf,"%5i ",tmp);
	  if (a == 0) {   // insert logodds value for B
	    fprintf(mtxf,"%5i ",-32768);
	  } else if (a == 18) {   // insert logodds value for X
	    fprintf(mtxf,"%5i ",-100);
	  } else if (a == 19) {   // insert logodds value for Z
	    fprintf(mtxf,"%5i ",-32768);
	  }
	}
      fprintf(mtxf,"-32768 -400\n");
    }
  fclose(mtxf);

  // Calculate secondary structure
  CalculateSS(ss_pred, ss_conf, tmpfile);
  
  q.AddSSPrediction(ss_pred, ss_conf);

  // Remove temp-files
  std::string command = "rm " + (std::string)tmpfile + "*";
  runSystem(command,v);
}

// Calculate secondary structure for given HMM
void CalculateSS(HMM& q)
{
  char ss_pred[MAXRES];
  char ss_conf[MAXRES];

  CalculateSS(q, ss_pred, ss_conf);
}


// Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
void ReadInput(char* infile, HMM& q, Alignment* qali=NULL)
{
    char path[NAMELEN];

    // Open query file and determine file type
    char line[LINELEN]=""; // input line
    FILE* inf=NULL;
    if (strcmp(infile,"stdin"))
    {
        inf = fopen(infile, "r");
        if (!inf) OpenFileError(infile);
        Pathname(path,infile);
    }
    else
    {
        inf = stdin;
        if (v>=2) printf("Reading HMM / multiple alignment from standard input ...\n(To get a help list instead, quit and type %s -h.)\n",program_name);
        *path='\0';
    }

    fgetline(line,LINELEN-1,inf);

    // Is infile a HMMER3 file?
    if (!strncmp(line,"HMMER3",6))
    {
        if (v>=2) cout<<"Query file is in HMMER3 format\n";
	cout<<"WARNING! Use of HMMER3 format as input results in dramatically loss of sensitivity!\n";

        // Read 'query HMMER file
        rewind(inf);
        q.ReadHMMer3(inf,path);
    }

    // ... or is infile an old HMMER file?
    else if (!strncmp(line,"HMMER",5))
    {
        if (v>=2) cout<<"Query file is in HMMER format\n";
	cout<<"WARNING! Use of HMMER format as input results in dramatically loss of sensitivity!\n";

        // Read 'query HMMER file
        rewind(inf);
        q.ReadHMMer(inf,path);
    }

    // ... or is it an hhm file?
    else if (!strncmp(line,"NAME",4) || !strncmp(line,"HH",2))
    {
        if (v>=2) cout<<"Query file is in HHM format\n";

        // Rewind to beginning of line and read query hhm file
        rewind(inf);
        q.Read(inf,path);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);

    }
    // ... or is it an alignment file
    else
    {
        Alignment* pali;
        if (qali==NULL) pali=new(Alignment); else pali=qali;
        if (par.calibrate) {
            printf("\nError in %s: only HHM files can be calibrated.\n",program_name);
            printf("Build an HHM file from your alignment with 'hhmake -i %s' and rerun hhsearch with the hhm file\n\n",infile);
            exit(1);
        }

        if (v>=2 && strcmp(infile,"stdin")) cout<<infile<<" is in A2M, A3M or FASTA format\n";

        // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
        pali->Read(inf,infile,line);

        // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
        // and store marked sequences in name[k] and seq[k]
        pali->Compress(infile);

        // Sort out the nseqdis most dissimilar sequences for display in the output alignments
        pali->FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);

        // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
        pali->N_filtered = pali->Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);

	if (par.Neff>=0.999) 
	  pali->FilterNeff();

        // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
        pali->FrequenciesAndTransitions(q);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: alignment %s looks too diverse (Neff=%.1f>11). Better check it with an alignment viewer... \n",q.name,q.Neff_HMM);

        if (qali==NULL) delete(pali);
    }
    fclose(inf);

    return;
}

// Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
void ReadAndPrepare(char* infile, HMM& q, Alignment* qali=NULL)
{
    char path[NAMELEN];

    // Open query file and determine file type
    char line[LINELEN]=""; // input line
    FILE* inf=NULL;
    if (strcmp(infile,"stdin"))
    {
        inf = fopen(infile, "r");
        if (!inf) OpenFileError(infile);
        Pathname(path,infile);
    }
    else
    {
        inf = stdin;
        if (v>=2) printf("Reading HMM / multiple alignment from standard input ...\n(To get a help list instead, quit and type %s -h.)\n",program_name);
        *path='\0';
    }

    fgetline(line,LINELEN-1,inf);

    // Is it an hhm file?
    if (!strncmp(line,"NAME",4) || !strncmp(line,"HH",2))
    {
        if (v>=2) cout<<"Query file is in HHM format\n";

        // Rewind to beginning of line and read query hhm file
        rewind(inf);
        q.Read(inf,path);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);

        // Add transition pseudocounts to query -> q.p[i][a]
        q.AddTransitionPseudocounts();

        if (!*par.clusterfile) { //compute context-specific pseudocounts?
	  // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	  q.PreparePseudocounts();
	  // Add amino acid pseudocounts to query:  q.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	  q.AddAminoAcidPseudocounts(q.has_pseudocounts ? 0:par.pcm, par.pca, par.pcb, par.pcc);;
        } else {
	  // Add context specific pseudocount to query
	  q.AddContextSpecificPseudocounts(q.has_pseudocounts ? 0:par.pcm);
        }
        
        q.CalculateAminoAcidBackground();
    }

    // ... or is it an a2m/a3m alignment file
    else if (line[0]=='#' || line[0]=='>')
    {
        Alignment* pali;
        if (qali==NULL) pali=new(Alignment); else pali=qali;
        if (par.calibrate) {
            printf("\nError in %s: only HHM files can be calibrated.\n",program_name);
            printf("Build an HHM file from your alignment with 'hhmake -i %s' and rerun hhsearch with the hhm file\n\n",infile);
            exit(1);
        }

        if (v>=2 && strcmp(infile,"stdin")) cout<<infile<<" is in A2M, A3M or FASTA format\n";

        // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
        pali->Read(inf,infile,line);

        // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
        // and store marked sequences in name[k] and seq[k]
        pali->Compress(infile);

        // Sort out the nseqdis most dissimilar sequences for display in the output alignments
        pali->FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);

        // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
        pali->N_filtered = pali->Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);

 	if (par.Neff>=0.999) 
	  pali->FilterNeff();

	// Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
        pali->FrequenciesAndTransitions(q);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: alignment %s looks too diverse (Neff=%.1f>11). Better check it with an alignment viewer... \n",q.name,q.Neff_HMM);

        // Add transition pseudocounts to query -> p[i][a]
        q.AddTransitionPseudocounts();

        if (!*par.clusterfile) { //compute context-specific pseudocounts?
	  // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	  q.PreparePseudocounts();
	  // Add amino acid pseudocounts to query:  p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	  q.AddAminoAcidPseudocounts(q.has_pseudocounts ? 0:par.pcm, par.pca, par.pcb, par.pcc);
        } else {
	  // Add context specific pseudocount to query
	  q.AddContextSpecificPseudocounts(q.has_pseudocounts ? 0:par.pcm);
        }

        q.CalculateAminoAcidBackground();

        if (qali==NULL) delete(pali);
    
    } else if (!strncmp(line,"HMMER",5)) {

        ///////////////////////////////////////////////////////////////////////////////////////
        // Don't allow HMMER format as input due to the dramatically loss of sensitivity!!!! (only allowed in HHmake)
        if (strncmp(program_name,"hhmake",6)) {
	  cerr<<endl<<"Error in "<<program_name<<": HMMER format not allowed as input due to the dramatically loss of sensitivity!\n";
	  exit(1);
        }
      
        // Is infile a HMMER3 file?
	if (!strncmp(line,"HMMER3",6))
	  {
	    if (v>=2) cout<<"Query file is in HMMER3 format\n";
	    
	    // Read 'query HMMER file
	    rewind(inf);
	    q.ReadHMMer3(inf,path);
	    
	    // Don't add transition pseudocounts to query!!
	    // DON'T ADD amino acid pseudocounts to query: pcm=0!  q.p[i][a] = f[i][a]
	    q.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
	    q.CalculateAminoAcidBackground();
	  }
	
	// ... or is infile an old HMMER file?
	else if (!strncmp(line,"HMMER",5))
	  {
	    if (v>=2) cout<<"Query file is in HMMER format\n";
	    
	    // Read 'query HMMER file
	    rewind(inf);
	    q.ReadHMMer(inf,path);
	    
	    // Don't add transition pseudocounts to query!!
	    
	    // NEEDED?????
	    
	    // if (!*par.clusterfile) { //compute context-specific pseudocounts?
	    //     // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	    //     q.PreparePseudocounts();
	    // } else {
	    //     // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
	    //     q.PrepareContextSpecificPseudocounts();
	    // }
	    
	    // DON'T ADD amino acid pseudocounts to query: pcm=0!  q.p[i][a] = f[i][a]
	    q.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
	    q.CalculateAminoAcidBackground();
	  }
	
    } else {
      cerr<<endl<<"Error in "<<program_name<<": unrecognized input file format in \'"<<infile<<"\'\n";
      cerr<<"line = "<<line<<"\n";
      exit(1);
    }
    fclose(inf);

    if (par.addss==1)
      CalculateSS(q);

    if (par.forward>=1) q.Log2LinTransitionProbs(1.0);
    return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Do precalculations for q and t to prepare comparison
/////////////////////////////////////////////////////////////////////////////////////
void PrepareTemplate(HMM& q, HMM& t, int format)
{
    if (format==0) // HHM format
    {
        // Add transition pseudocounts to template
        t.AddTransitionPseudocounts();

	// Don't use CS-pseudocounts because of runtime!!!
	// Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	t.PreparePseudocounts();
	
	// Add amino acid pseudocounts to query:  p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	t.AddAminoAcidPseudocounts(t.has_pseudocounts ? 0:par.pcm, par.pca, par.pcb, par.pcc);

        t.CalculateAminoAcidBackground();
    }
    else // HHMER format
    {
        // Don't add transition pseudocounts to template
        // t.AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg, par.gaph, par.gapi, 0.0);

        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        // t.PreparePseudocounts();

        // DON'T ADD amino acid pseudocounts to temlate: pcm=0!  t.p[i][a] = t.f[i][a]
        t.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
        t.CalculateAminoAcidBackground();
    }

    if (par.forward>=1) t.Log2LinTransitionProbs(1.0);

    // Factor Null model into HMM t
    // ATTENTION! t.p[i][a] is divided by pnul[a] (for reasons of efficiency) => do not reuse t.p
    t.IncludeNullModelInHMM(q,t);  // Can go BEFORE the loop if not dependent on template

    return;
}

void WriteToAlifile(FILE* alitabf, Hit* hit)
    {
      if (par.forward==2 || par.realign) 
	{
	  if (hit->nss_dssp >= 0)
	    {
	    // secondary structure determined by dssp 0:-  1:H  2:E  3:C  4:S  5:T  6:G  7:B
	      fprintf(alitabf,"    i     j  score     SS  probab  dssp\n");
	      for (int step=hit->nsteps; step>=1; step--)
		if (hit->states[step]>=MM) 
		  fprintf(alitabf,"%5i %5i %6.2f %6.2f %7.4f %5c\n",hit->i[step],hit->j[step],hit->S[step],hit->S_ss[step],hit->P_posterior[step],hit->seq[hit->nss_dssp][hit->j[step]]);
	    }
	  else 
	    {
	      fprintf(alitabf, "missing dssp\n");
	      fprintf(alitabf,"    i     j  score     SS  probab\n");
	      for (int step=hit->nsteps; step>=1; step--)
	    	if (hit->states[step]>=MM) 
		  fprintf(alitabf,"%5i %5i %6.2f %6.2f %7.4f\n",hit->i[step],hit->j[step],hit->S[step],hit->S_ss[step],hit->P_posterior[step]);
	    }
	} 
      else 
	{
	  fprintf(alitabf,"    i     j  score     SS\n");
	  for (int step=hit->nsteps; step>=1; step--)
	    if (hit->states[step]>=MM) 
	      fprintf(alitabf,"%5i %5i %6.2f %6.2f\n",hit->i[step],hit->j[step],hit->S[step],hit->S_ss[step]);
	}
      return;
   }
