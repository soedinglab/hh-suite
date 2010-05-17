// hhfunc.C


// Calculate secondary structure prediction with PSIpred
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile)
{
  // Initialize
  std::string command;
  std::string tmpinfile = (std::string)tmpfile + ".fas";
  std::string tmppsifile = (std::string)tmpfile + ".psi";
  char line[LINELEN]=""; 

  strcpy(ss_pred," ");
  strcpy(ss_conf," ");
  
  // Create tmp-file
  char rootname[NAMELEN];
  RemovePath(rootname,tmpfile);

  // Create dummy-DB if not exists
  //if (!*dummydb)
  if( access( par.dummydb, R_OK ) == -1 )
    {
      strcpy(par.dummydb,tmpfile);
      strcat(par.dummydb,"_dummy_db");
      command = "cp " + tmpinfile + " " + (std::string)par.dummydb;
      runSystem(command,v);
      command = (std::string)par.blast + "/formatdb -i " + (std::string)par.dummydb + " -l /dev/null > /dev/null";
      runSystem(command,v);
    }

  // Create BLAST checkpoint file
  command = (std::string)par.blast + "/blastpgp -b 1 -j 1 -h 0.001 -d " + (std::string)par.dummydb + " -i " + tmpinfile + " -B " + tmppsifile + " -C " + (std::string)tmpfile + ".chk 1> /dev/null 2> /dev/null";
  runSystem(command,v);
  command = "echo " + (std::string)rootname + ".chk > " + (std::string)tmpfile + ".pn";
  runSystem(command,v);
  command = "echo " + (std::string)rootname + ".fas > " + (std::string)tmpfile + ".sn";
  runSystem(command,v);
  command =  (std::string)par.blast + "/makemat -P " + (std::string)tmpfile;
  runSystem(command,v);

  // Run PSIpred
  command = (std::string)par.psipred + "/psipred " + (std::string)tmpfile + ".mtx " + (std::string)par.psipred_data + "/weights.dat " + (std::string)par.psipred_data + "/weights.dat2 " + (std::string)par.psipred_data + "/weights.dat3 " + (std::string)par.psipred_data + "/weights.dat4 > " + (std::string)tmpfile + ".ss";
  runSystem(command,v);
  command = (std::string)par.psipred + "/psipass2 " + (std::string)par.psipred_data + "/weights_p2.dat 1 0.98 1.09 " + (std::string)tmpfile + ".ss2 " + (std::string)tmpfile + ".ss > " + (std::string)tmpfile + ".horiz";
  runSystem(command,v);

  // Read results
  char filename[NAMELEN];
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

// Calculate secondary structure for given HMM
void CalculateSS(HMM& q)
{
  char ss_pred[MAXRES];
  char ss_conf[MAXRES];

  char tmpfile[]="/tmp/hhCalcSSXXXXXX";
  if (mkstemp(tmpfile) == -1) {
    cerr << "ERROR! Could not create tmp-file!\n"; 
    exit(4);
  }
  
  char queryfile[NAMELEN];
  strcpy(queryfile,tmpfile);
  strcat(queryfile,".fas");
  char alifile[NAMELEN];
  strcpy(alifile,tmpfile);
  strcat(alifile,".psi");

  // Write query-file
  FILE* outf = fopen(queryfile,"w");
  if (!outf) OpenFileError(queryfile);
  fprintf(outf,">%s\n%s\n",q.longname,q.seq[q.nfirst]+1);
  fclose(outf);

  // Write ali-file
  HalfAlignment qa;
  int n = imin(q.n_display,par.nseqdis+(q.nss_dssp>=0)+(q.nss_pred>=0)+(q.nss_conf>=0)+(q.ncons>=0));
  qa.Set(q.name,q.seq,q.sname,n,q.L,q.nss_dssp,q.nss_pred,q.nss_conf,q.nsa_dssp,q.ncons);
  qa.BuildA3M();
  qa.Print(alifile,NULL,"psi");   // print alignment to outfile

  // Calculate secondary structure
  CalculateSS(ss_pred, ss_conf, tmpfile);
  
  q.AddSSPrediction(ss_pred, ss_conf);
  
  // Remove temp-files
  std::string command = "rm " + (std::string)tmpfile + "*";
  runSystem(command,v);
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

    // Is it an hhm file?
    if (!strncmp(line,"NAME",4) || !strncmp(line,"HH",2))
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

    if (par.addss==1)
      CalculateSS(q);

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

    // Is infile a HMMER file?
    if (!strncmp(line,"HMMER",5))
    {
        if (v>=2) cout<<"Query file is in HMMER format\n";

        // Read 'query HMMER file
        q.ReadHMMer(inf,path);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);

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

    // ... or is it an hhm file?
    else if (!strncmp(line,"NAME",4) || !strncmp(line,"HH",2))
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
	  fprintf(alitabf,"    i     j  score     SS  probab\n");
	  for (int step=hit->nsteps; step>=1; step--)
	    if (hit->states[step]>=MM) 
	      fprintf(alitabf,"%5i %5i %6.2f %6.2f %7.4f\n",hit->i[step],hit->j[step],hit->S[step],hit->S_ss[step],hit->P_posterior[step]);
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
