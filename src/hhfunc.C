// hhfunc.C

// Read input file (HMM, HHM, or alignment format)
void ReadQueryFile(char* infile, char& input_format, HMM* q, Alignment* qali =
    NULL);

// Add transition and amino acid pseudocounts to query HMM, calculate aa background etc.
void PrepareQueryHMM(char& input_format, HMM* q);

// Do precalculations for q and t to prepare comparison
void PrepareTemplateHMM(HMM* q, HMM* t, int format);

// Calculate secondary structure prediction with PSIPRED
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile);

// Calculate secondary structure for given HMM and return prediction
void CalculateSS(HMM* q, char *ss_pred, char *ss_conf);

// Calculate secondary structure for given HMM
void CalculateSS(HMM* q);

// Write alignment in tab format (option -atab)
void WriteToAlifile(FILE* alitabf, Hit* hit);

/////////////////////////////////////////////////////////////////////////////////////
// Read input file (HMM, HHM, or alignment format)
/////////////////////////////////////////////////////////////////////////////////////
void ReadQueryFile(char* infile, char& input_format, HMM* q, Alignment* qali) {
  // Open query file and determine file type
  char path[NAMELEN]; // path of input file (is needed to write full path and file name to HMM FILE record)
  char line[LINELEN] = "";   // input line
  FILE* inf = NULL;
  if (strcmp(infile, "stdin")) {
    inf = fopen(infile, "r");
    if (!inf)
      OpenFileError(infile);
    Pathname(path, infile);
  }
  else {
    inf = stdin;
    if (v >= 2)
      printf(
          "Reading HMM / multiple alignment from standard input ...\n(To get a help list instead, quit and type %s -h.)\n",
          program_name);
    *path = '\0';
  }
  
  if (!fgetline(line, LINELEN, inf)) {
    cerr << endl << "Error in " << program_name << ": " << infile
        << " is empty!\n";
    exit(4);
  }
  while (strscn(line) == NULL)
    fgetline(line, LINELEN, inf); // skip lines that contain only white space

  // Is infile a HMMER file?
  if (!strncmp(line, "HMMER", 5)) {
    // Uncomment this line to allow HMMER2/HMMER3 models as queries:
    cerr
        << "Error: Use of HMMER format as input will result in severe loss of sensitivity!\n";
    exit(4);

    // Read query HMMER file
    rewind(inf);
    if (!strncmp(line, "HMMER3", 6)) {
      if (v >= 3)
        cout << "Query file is in HMMER3 format\n";
      q->ReadHMMer3(inf, path);
    }
    else {
      if (v >= 3)
        cout << "Query file is in HMMER format\n";
      q->ReadHMMer(inf, path);
    }
    input_format = 1;
    par.hmmer_used = true;
  }
  
  // ... or is it an hhm file?
  else if (!strncmp(line, "NAME", 4) || !strncmp(line, "HH", 2)) {
    if (v >= 2)
      cout << "Query file is in HHM format\n";

    // Rewind to beginning of line and read query hhm file
    rewind(inf);
    q->Read(inf, path);
    input_format = 0;
  }
  // ... or is it an alignment file
  else if (line[0] == '#' || line[0] == '>')         // read sequence/alignment
      {
    Alignment* pali;
    if (qali == NULL)
      pali = new (Alignment);
    else
      pali = qali;
    if (par.calibrate) {
      printf("\nError in %s: only HHM files can be calibrated.\n",
          program_name);
      printf(
          "Build an HHM file from your alignment with 'hhmake -i %s' and rerun hhsearch with the hhm file\n\n",
          infile);
      exit(1);
    }

    if (v >= 2 && strcmp(infile, "stdin"))
      cout << infile << " is in A2M, A3M or FASTA format\n";

    // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
    pali->Read(inf, infile, line);

    // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
    // and store marked sequences in name[k] and seq[k]
    pali->Compress(infile);

    // Sort out the nseqdis most dissimilar sequences for display in the output alignments
    pali->FilterForDisplay(par.max_seqid, par.coverage, par.qid, par.qsc,
        par.nseqdis);

    // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
    pali->N_filtered = pali->Filter(par.max_seqid, par.coverage, par.qid,
        par.qsc, par.Ndiff);

    if (par.Neff >= 0.999)
      pali->FilterNeff();

    // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
    pali->FrequenciesAndTransitions(q);
    input_format = 0;
    if (qali == NULL)
      delete (pali);
  }
  else {
    cerr << endl << "Error in " << program_name
        << ": unrecognized input file format in \'" << infile << "\'\n";
    cerr << "line = " << line << "\n";
    exit(1);
  }

  if (v >= 2 && input_format == 0 && q->Neff_HMM > 11.0)
    fprintf(stderr,
        "WARNING: MSA %s looks too diverse (Neff=%.1f>11). Better check it with an alignment viewer for non-homologous segments. Also consider building the MSA with hhblits using the - option to limit MSA diversity.\n",
        q->name, q->Neff_HMM);

  fclose(inf);
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Add transition and amino acid pseudocounts to query HMM, calculate aa background etc.
/////////////////////////////////////////////////////////////////////////////////////
void PrepareQueryHMM(char& input_format, HMM* q) {
  // Was query an HHsearch formatted file or MSA (no pseudocounts added yet)?
  if (input_format == 0) {
    // Add transition pseudocounts to query -> q->p[i][a]
    q->AddTransitionPseudocounts();

    // Compute substitutino matrix pseudocounts?
    if (par.nocontxt) {
      // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
      q->PreparePseudocounts();
      // Add amino acid pseudocounts to query:  q->p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
      q->AddAminoAcidPseudocounts(par.aa_pcm, par.aa_pca, par.aa_pcb,
          par.aa_pcc);
    }
    else {
      // Add context specific pseudocount to query
      q->AddContextSpecificPseudocounts(pc, pc_admix);
    }

  }
  
  // or was query a HMMER file? (pseudocounts already added!)
  else if (input_format == 1) {
    // Don't add transition pseudocounts to query!!
    // DON'T ADD amino acid pseudocounts to query: pcm=0!  q->p[i][a] = f[i][a]
    q->AddAminoAcidPseudocounts(0, par.aa_pca, par.aa_pcb, par.aa_pcc);
  }
  
  q->CalculateAminoAcidBackground();
  
  // if (par.addss==1) CalculateSS(q);
  
  if (par.columnscore == 5 && !q->divided_by_local_bg_freqs)
    q->DivideBySqrtOfLocalBackgroundFreqs(
        par.half_window_size_local_aa_bg_freqs);
  
  if (par.forward >= 1)
    q->Log2LinTransitionProbs(1.0);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Do precalculations for q and t to prepare comparison
/////////////////////////////////////////////////////////////////////////////////////
void PrepareTemplateHMM(HMM* q, HMM* t, int format) {
  if (format == 0) // HHM format
      {
    // Add transition pseudocounts to template
    t->AddTransitionPseudocounts();

    // Don't use CS-pseudocounts because of runtime!!!
    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
    t->PreparePseudocounts();

    // Add amino acid pseudocounts to query:  p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    t->AddAminoAcidPseudocounts(par.aa_pcm, par.aa_pca, par.aa_pcb, par.aa_pcc);
  }
  else // HHMER format
  {
    // Don't add transition pseudocounts to template
    // t->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg, par.gaph, par.gapi, 0.0);

    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
    // t->PreparePseudocounts();
    
    // DON'T ADD amino acid pseudocounts to temlate: pcm=0!  t->p[i][a] = t->f[i][a]
    t->AddAminoAcidPseudocounts(0, par.aa_pca, par.aa_pcb, par.aa_pcc);
  }
  t->CalculateAminoAcidBackground();

  if (par.forward >= 1)
    t->Log2LinTransitionProbs(1.0);

  // Factor Null model into HMM t
  // ATTENTION! t->p[i][a] is divided by pnul[a] (for reasons of efficiency) => do not reuse t->p
  t->IncludeNullModelInHMM(q, t); // Can go BEFORE the loop if not dependent on template

  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate secondary structure prediction with PSIPRED
/////////////////////////////////////////////////////////////////////////////////////
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile) {
  // Initialize
  std::string command;
  char line[LINELEN] = "";
  char filename[NAMELEN];
  
  strcpy(ss_pred, "-");
  strcpy(ss_conf, "-");

  // Run PSIPRED
  
  // Check for PSIPRED ver >= 3.0 (weights.dat4 doesn't exists anymore)
  strcpy(filename, par.psipred_data);
  strcat(filename, "/weights.dat4");
  FILE* check_exists = fopen(filename, "r");
  if (check_exists) {  // Psipred version < 3.0
    command = (std::string) par.psipred + "/psipred " + (std::string) tmpfile
        + ".mtx " + (std::string) par.psipred_data + "/weights.dat "
        + (std::string) par.psipred_data + "/weights.dat2 "
        + (std::string) par.psipred_data + "/weights.dat3 "
        + (std::string) par.psipred_data + "/weights.dat4 > "
        + (std::string) tmpfile + ".ss";
  }
  else {
    command = (std::string) par.psipred + "/psipred " + (std::string) tmpfile
        + ".mtx " + (std::string) par.psipred_data + "/weights.dat "
        + (std::string) par.psipred_data + "/weights.dat2 "
        + (std::string) par.psipred_data + "/weights.dat3 > "
        + (std::string) tmpfile + ".ss";
  }
  runSystem(command, v);
  command = (std::string) par.psipred + "/psipass2 "
      + (std::string) par.psipred_data + "/weights_p2.dat 1 0.98 1.09 "
      + (std::string) tmpfile + ".ss2 " + (std::string) tmpfile + ".ss > "
      + (std::string) tmpfile + ".horiz";
  runSystem(command, v);

  // Read results
  strcpy(filename, tmpfile);
  strcat(filename, ".horiz");
  FILE* horizf = fopen(filename, "r");
  if (!horizf)
    return;

  while (fgets(line, LINELEN, horizf)) {
    char tmp_seq[NAMELEN] = "";
    char* ptr = line;
    if (!strncmp(line, "Conf:", 5)) {
      ptr += 5;
      strwrd(tmp_seq, ptr);
      strcat(ss_conf, tmp_seq);
    }
    if (!strncmp(line, "Pred:", 5)) {
      ptr += 5;
      strwrd(tmp_seq, ptr);
      strcat(ss_pred, tmp_seq);
    }
  }
  fclose(horizf);

  if (v > 3) {
    printf("SS-pred: %s\n", ss_pred);
    printf("SS-conf: %s\n", ss_conf);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate secondary structure for given HMM and return prediction
/////////////////////////////////////////////////////////////////////////////////////
void CalculateSS(HMM* q, char *ss_pred, char *ss_conf) {

  if (q->divided_by_local_bg_freqs) {
    cerr
        << "WARNING: Can not add predicted secondary structure when using column score 5!\n";
    return;
  }

  char tmpfile[] = "/tmp/HHsuite_CaluclateSS_XXXXXX";
  if (mkstemp(tmpfile) == -1) {
    cerr << "Error: Could not create tmp file " << tmpfile << "!\n";
    exit(4);
  }
  
  // Write log-odds matrix from q to tmpfile.mtx
  char filename[NAMELEN];
  FILE* mtxf = NULL;
  
  strcpy(filename, tmpfile);
  strcat(filename, ".mtx");
  mtxf = fopen(filename, "w");
  if (!mtxf)
    OpenFileError(filename);

  fprintf(mtxf, "%i\n", q->L);
  fprintf(mtxf, "%s\n", q->seq[q->nfirst] + 1);
  fprintf(mtxf,
      "2.670000e-03\n4.100000e-02\n-3.194183e+00\n1.400000e-01\n2.670000e-03\n4.420198e-02\n-3.118986e+00\n1.400000e-01\n3.176060e-03\n1.339561e-01\n-2.010243e+00\n4.012145e-01\n");

  for (int i = 1; i <= q->L; ++i) {
    fprintf(mtxf, "-32768 ");
    for (int a = 0; a < 20; ++a) {
      int tmp = iround(50 * flog2(q->p[i][s2a[a]] / pb[s2a[a]]));
      fprintf(mtxf, "%5i ", tmp);
      if (a == 0) {   // insert logodds value for B
        fprintf(mtxf, "%5i ", -32768);
      }
      else if (a == 18) {   // insert logodds value for X
        fprintf(mtxf, "%5i ", -100);
      }
      else if (a == 19) {   // insert logodds value for Z
        fprintf(mtxf, "%5i ", -32768);
      }
    }
    fprintf(mtxf, "-32768 -400\n");
  }
  fclose(mtxf);

  // Calculate secondary structure
  CalculateSS(ss_pred, ss_conf, tmpfile);
  
  q->AddSSPrediction(ss_pred, ss_conf);

  // Remove temp-files
  std::string command = "rm " + (std::string) tmpfile + "*";
  runSystem(command, v);
}

// Calculate secondary structure for given HMM
void CalculateSS(HMM* q) {
  char ss_pred[par.maxres];
  char ss_conf[par.maxres];

  CalculateSS(q, ss_pred, ss_conf);
}

/////////////////////////////////////////////////////////////////////////////////////
// Write alignment in tab format (option -atab)
/////////////////////////////////////////////////////////////////////////////////////
void WriteToAlifile(FILE* alitabf, Hit* hit) {
      if (hit->P_posterior != NULL && (par.forward==2 || par.realign)) 
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

/////////////////////////////////////////////////////////////////////////////////////
// Read number of sequences in annotation, after second '|'
/////////////////////////////////////////////////////////////////////////////////////
int SequencesInCluster(char* name) {
  int num = 1;
  char *ptr = strchr(name, '|');
  if (!strncmp(name, "cl|", 3) || !strncmp(name, "UP20|", 5)
      || !strncmp(name, "NR20|", 5))   // kClust formatted database (NR20, ...)
          {
    if (*ptr == '|')
      ptr = strchr(ptr, '|');
    if (*ptr == '|') {
      num = strint(ptr);
      if (num < 0)
        num = 1;
    }
  }
  return num;
}
