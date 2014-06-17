// hhfunc.C

#include "hhfunc.h"

#include <sstream>

/////////////////////////////////////////////////////////////////////////////////////
// Read input file (HMM, HHM, or alignment format)
/////////////////////////////////////////////////////////////////////////////////////
void ReadQueryFile(Parameters& par, FILE* inf, char& input_format,
    char use_global_weights, HMM* q, Alignment& qali, char infile[], float* pb,
    const float S[20][20], const float Sim[20][20]) {
  char line[LINELEN];

  if (!fgetline(line, LINELEN, inf)) {
	HH_LOG(LogLevel::ERROR) << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	HH_LOG(LogLevel::ERROR) << "\t" << infile << " is empty!\n";
    exit(4);
  }
  while (strscn(line) == NULL)
    fgetline(line, LINELEN, inf); // skip lines that contain only white space

  // Is infile a HMMER file?
  if (!strncmp(line, "HMMER", 5)) {
    // Uncomment this line to allow HMMER2/HMMER3 models as queries:
    HH_LOG(LogLevel::ERROR) << "Error: Use of HMMER format as input will result in severe loss of sensitivity!\n";
  }
  // ... or is it an hhm file?
  else if (!strncmp(line, "NAME", 4) || !strncmp(line, "HH", 2)) {
    char path[NAMELEN];
    Pathname(path, infile);

    HH_LOG(LogLevel::INFO) << "Query file is in HHM format\n";

    // Rewind to beginning of line and read query hhm file
    rewind(inf);
    q->Read(inf, par.maxcol, par.nseqdis, pb, path);
    input_format = 0;

    // HHM format
    if (input_format == 0) {
      HH_LOG(LogLevel::INFO) << "Extracting representative sequences from " << infile << " to merge later with matched database sequences\n";
    }

    Alignment ali_tmp;
    ali_tmp.GetSeqsFromHMM(q);
    ali_tmp.Compress(infile, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
    qali = ali_tmp;
  }
  // ... or is it an alignment file
  else if (line[0] == '#' || line[0] == '>') {
    if (par.calibrate) {
      HH_LOG(LogLevel::ERROR) << "Error in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      HH_LOG(LogLevel::ERROR) << "\tonly HHM files can be calibrated.\n";
      HH_LOG(LogLevel::ERROR) << "\tBuild an HHM file from your alignment with hhmake and rerun with the hhm file" << std::endl;
      exit(1);
    }

    if (strcmp(infile, "stdin")) {
    	HH_LOG(LogLevel::INFO) << infile << " is in A2M, A3M or FASTA format\n";
    }

    Alignment ali_tmp;

    // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
    ali_tmp.Read(inf, infile, par.mark, par.maxcol, par.nseqdis, line);

    // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
    // and store marked sequences in name[k] and seq[k]
    ali_tmp.Compress(infile, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);

    // Sort out the nseqdis most dissimilar sequences for display in the output alignments
    ali_tmp.FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid,
        par.qsc, par.nseqdis);

    // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
    ali_tmp.N_filtered = ali_tmp.Filter(par.max_seqid, S, par.coverage, par.qid,
        par.qsc, par.Ndiff);

    if (par.Neff >= 0.999)
    	ali_tmp.FilterNeff(use_global_weights, par.mark, par.cons, par.showcons,
          par.maxres, par.max_seqid, par.coverage, par.Neff, pb, S, Sim);

    // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
    ali_tmp.FrequenciesAndTransitions(q, use_global_weights, par.mark, par.cons,
        par.showcons, par.maxres, pb, Sim);

    qali = ali_tmp;
    input_format = 0;
  }
  else {
	HH_LOG(LogLevel::ERROR) << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	HH_LOG(LogLevel::ERROR) << "\tunrecognized input file format in \'" << infile << "\'\n";
	HH_LOG(LogLevel::ERROR) << "\tline = " << line << "\n";
    exit(1);
  }

  if (input_format == 0 && q->Neff_HMM > 11.0) {
	  HH_LOG(LogLevel::WARNING) << "WARNING: MSA " << q->name << " looks too diverse (Neff=" << q->Neff_HMM << ">11). Better check it with an alignment viewer for non-homologous segments. Also consider building the MSA with hhblits using the - option to limit MSA diversity.\n";
  }
}

void ReadQueryFile(Parameters& par, char* infile, char& input_format,
    char use_global_weights, HMM* q, Alignment& qali, float* pb,
    const float S[20][20], const float Sim[20][20]) {
  // Open query file and determine file type
  char path[NAMELEN]; // path of input file (is needed to write full path and file name to HMM FILE record)
  FILE* inf = NULL;
  if (strcmp(infile, "stdin") == 0) {
    inf = stdin;
    HH_LOG(LogLevel::INFO) << "Reading HMM / multiple alignment from standard input ...\n";
    path[0] = '\0';
  }
  else {
    inf = fopen(infile, "r");
    if (!inf)
      OpenFileError(infile, __FILE__, __LINE__, __func__);
    Pathname(path, infile);
  }
  
  ReadQueryFile(par, inf, input_format, use_global_weights, q, qali, infile, pb,
      S, Sim);

  fclose(inf);
}

/////////////////////////////////////////////////////////////////////////////////////
// Add transition and amino acid pseudocounts to query HMM, calculate aa background etc.
/////////////////////////////////////////////////////////////////////////////////////
void PrepareQueryHMM(Parameters& par, char& input_format, HMM* q,
    cs::Pseudocounts<cs::AA>* pc_hhm_context_engine,
    cs::Admix* pc_hhm_context_mode, const float* pb, const float R[20][20]) {
  // Was query an HHsearch formatted file or MSA (no pseudocounts added yet)?
  if (input_format == 0) {
    // Add transition pseudocounts to query -> q->p[i][a]
    q->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg,
        par.gaph, par.gapi, par.gapb, par.gapb);

    // Compute substitutino matrix pseudocounts?
    if (par.nocontxt) {
      // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
      q->PreparePseudocounts(R);
      // Add amino acid pseudocounts to query:  q->p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
      q->AddAminoAcidPseudocounts(par.pc_hhm_nocontext_mode,
          par.pc_hhm_nocontext_a, par.pc_hhm_nocontext_b,
          par.pc_hhm_nocontext_c);
    }
    else {
      // Add context specific pseudocount to query
      q->AddContextSpecificPseudocounts(pc_hhm_context_engine,
          pc_hhm_context_mode);
    }
  }
  // or was query a HMMER file? (pseudocounts already added!)
  else if (input_format == 1) {
    // Don't add transition pseudocounts to query!!
    // DON'T ADD amino acid pseudocounts to query: pcm=0!  q->p[i][a] = f[i][a]
    q->AddAminoAcidPseudocounts(0, par.pc_hhm_nocontext_a,
        par.pc_hhm_nocontext_b, par.pc_hhm_nocontext_c);
  }
  
  q->CalculateAminoAcidBackground(pb);
  
  // if (par.addss==1) CalculateSS(q);
  
  if (par.columnscore == 5 && !q->divided_by_local_bg_freqs)
    q->DivideBySqrtOfLocalBackgroundFreqs(
        par.half_window_size_local_aa_bg_freqs, pb);
  
  if (par.forward >= 1)
    q->Log2LinTransitionProbs(1.0);
}

/////////////////////////////////////////////////////////////////////////////////////
// Do precalculations for q and t to prepare comparison
/////////////////////////////////////////////////////////////////////////////////////
void PrepareTemplateHMM(Parameters& par, HMM* q, HMM* t, int format,
    const float* pb, const float R[20][20]) {
  // HHM format
  if (format == 0) {
    // Add transition pseudocounts to template
    t->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg,
        par.gaph, par.gapi, par.gapb, par.gapb);

    // Don't use CS-pseudocounts because of runtime!!!
    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
    t->PreparePseudocounts(R);

    // Add amino acid pseudocounts to query:  p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    t->AddAminoAcidPseudocounts(par.pc_hhm_nocontext_mode,
        par.pc_hhm_nocontext_a, par.pc_hhm_nocontext_b, par.pc_hhm_nocontext_c);
  }
  // HHMER format
  else {
    // Don't add transition pseudocounts to template
    // t->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg, par.gaph, par.gapi, 0.0);

    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
    // t->PreparePseudocounts();
    
    // DON'T ADD amino acid pseudocounts to temlate: pcm=0!  t->p[i][a] = t->f[i][a]
    t->AddAminoAcidPseudocounts(0, par.pc_hhm_nocontext_a,
        par.pc_hhm_nocontext_b, par.pc_hhm_nocontext_c);
  }
  t->CalculateAminoAcidBackground(pb);

  if (par.forward >= 1)
    t->Log2LinTransitionProbs(1.0);

  // Factor Null model into HMM t
  // ATTENTION! t->p[i][a] is divided by pnul[a] (for reasons of efficiency) => do not reuse t->p
  t->IncludeNullModelInHMM(q, t, par.columnscore,
      par.half_window_size_local_aa_bg_freqs, pb); // Can go BEFORE the loop if not dependent on template

  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate secondary structure prediction with PSIPRED
/////////////////////////////////////////////////////////////////////////////////////
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile,
    const char* psipred_data, const char* psipred) {
  // Initialize
  std::string command;
  char line[LINELEN] = "";
  char filename[NAMELEN];
  
  strcpy(ss_pred, "-");
  strcpy(ss_conf, "-");

  // Run PSIPRED
  
  // Check for PSIPRED ver >= 3.0 (weights.dat4 doesn't exists anymore)
  strcpy(filename, psipred_data);
  strcat(filename, "/weights.dat4");
  FILE* check_exists = fopen(filename, "r");
  if (check_exists) {  // Psipred version < 3.0
    command = (std::string) psipred + "/psipred " + (std::string) tmpfile
        + ".mtx " + (std::string) psipred_data + "/weights.dat "
        + (std::string) psipred_data + "/weights.dat2 "
        + (std::string) psipred_data + "/weights.dat3 "
        + (std::string) psipred_data + "/weights.dat4 > "
        + (std::string) tmpfile + ".ss";
  }
  else {
    command = (std::string) psipred + "/psipred " + (std::string) tmpfile
        + ".mtx " + (std::string) psipred_data + "/weights.dat "
        + (std::string) psipred_data + "/weights.dat2 "
        + (std::string) psipred_data + "/weights.dat3 > "
        + (std::string) tmpfile + ".ss";
  }
  runSystem(command);
  command = (std::string) psipred + "/psipass2 " + (std::string) psipred_data
      + "/weights_p2.dat 1 0.98 1.09 " + (std::string) tmpfile + ".ss2 "
      + (std::string) tmpfile + ".ss > " + (std::string) tmpfile + ".horiz";
  runSystem(command);

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

  HH_LOG(LogLevel::DEBUG1) << "SS-pred: " << ss_pred << std::endl;
  HH_LOG(LogLevel::DEBUG1) << "SS-conf: " << ss_conf << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate secondary structure for given HMM and return prediction
/////////////////////////////////////////////////////////////////////////////////////
void CalculateSS(HMM* q, char *ss_pred, char *ss_conf, const char* psipred_data,
    const char* psipred, const float* pb) {

  if (q->divided_by_local_bg_freqs) {
    std::cerr
        << "WARNING: Can not add predicted secondary structure when using column score 5!\n";
    return;
  }

  char tmpfile[] = "/tmp/HHsuite_CaluclateSS_XXXXXX";
  if (mkstemp(tmpfile) == -1) {
    std::cerr << "Error: Could not create tmp file " << tmpfile << "!\n";
    exit(4);
  }
  
  // Write log-odds matrix from q to tmpfile.mtx
  char filename[NAMELEN];
  FILE* mtxf = NULL;
  
  strcpy(filename, tmpfile);
  strcat(filename, ".mtx");
  mtxf = fopen(filename, "w");
  if (!mtxf)
    OpenFileError(filename, __FILE__, __LINE__, __func__);

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
  CalculateSS(ss_pred, ss_conf, tmpfile, psipred_data, psipred);
  
  q->AddSSPrediction(ss_pred, ss_conf);

  // Remove temp-files
  std::string command = "rm " + (std::string) tmpfile + "*";
  runSystem(command);
}

// Calculate secondary structure for given HMM
void CalculateSS(HMM* q, const int maxres, const char* psipred_data,
    const char* psipred, const float* pb) {
  char ss_pred[maxres];
  char ss_conf[maxres];

  CalculateSS(q, ss_pred, ss_conf, psipred_data, psipred, pb);
}

/////////////////////////////////////////////////////////////////////////////////////
// Write alignment in tab format (option -atab)
/////////////////////////////////////////////////////////////////////////////////////
void WriteToAlifile(FILE* alitabf, Hit* hit, const char forward,
    const char realign) {
  if (hit->P_posterior != NULL && (forward == 2 || realign)) {
    if (hit->nss_dssp >= 0) {
      // secondary structure determined by dssp 0:-  1:H  2:E  3:C  4:S  5:T  6:G  7:B
      fprintf(alitabf, "    i     j  score     SS  probab  dssp\n");
      for (int step = hit->nsteps; step >= 1; step--)
        if (hit->states[step] >= MM)
          fprintf(alitabf, "%5i %5i %6.2f %6.2f %7.4f %5c\n", hit->i[step],
              hit->j[step], hit->S[step], hit->S_ss[step],
              hit->P_posterior[step], hit->seq[hit->nss_dssp][hit->j[step]]);
    }
    else {
      fprintf(alitabf, "missing dssp\n");
      fprintf(alitabf, "    i     j  score     SS  probab\n");
      for (int step = hit->nsteps; step >= 1; step--)
        if (hit->states[step] >= MM)
          fprintf(alitabf, "%5i %5i %6.2f %6.2f %7.4f\n", hit->i[step],
              hit->j[step], hit->S[step], hit->S_ss[step],
              hit->P_posterior[step]);
    }
  }
  else {
    fprintf(alitabf, "    i     j  score     SS\n");
    for (int step = hit->nsteps; step >= 1; step--)
      if (hit->states[step] >= MM)
        fprintf(alitabf, "%5i %5i %6.2f %6.2f\n", hit->i[step], hit->j[step],
            hit->S[step], hit->S_ss[step]);
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

void InitializePseudocountsEngine(Parameters& par,
    cs::ContextLibrary<cs::AA>* context_lib, cs::Crf<cs::AA>* crf,
    cs::Pseudocounts<cs::AA>* pc_hhm_context_engine,
    cs::Admix* pc_hhm_context_mode,
    cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine,
    cs::Admix* pc_prefilter_context_mode) {
  // Prepare pseudocounts engine
  FILE* fin = fopen(par.clusterfile, "r");
  if (!fin) {
    std::cerr << std::endl << "Error in " << par.argv[0]
        << ": could not open file \'" << par.clusterfile << "\'\n";
    exit(2);
  }
  char ext[100];
  Extension(ext, par.clusterfile);
  if (strcmp(ext, "crf") == 0) {
    crf = new cs::Crf<cs::AA>(fin);
    pc_hhm_context_engine = new cs::CrfPseudocounts<cs::AA>(*crf);
    pc_prefilter_context_engine = new cs::CrfPseudocounts<cs::AA>(*crf);
  }
  else {
    context_lib = new cs::ContextLibrary<cs::AA>(fin);
    cs::TransformToLog(*context_lib);
    pc_hhm_context_engine = new cs::LibraryPseudocounts<cs::AA>(*context_lib,
        par.csw, par.csb);
    pc_prefilter_context_engine = new cs::LibraryPseudocounts<cs::AA>(
        *context_lib, par.csw, par.csb);
  }
  fclose(fin);
  pc_hhm_context_engine->SetTargetNeff(par.pc_hhm_context_engine.target_neff);
  pc_prefilter_context_engine->SetTargetNeff(
      par.pc_prefilter_context_engine.target_neff);

  // Prepare pseudocounts admixture method
  pc_hhm_context_mode = par.pc_hhm_context_engine.CreateAdmix();
  pc_prefilter_context_mode = par.pc_prefilter_context_engine.CreateAdmix();
}

void DeletePseudocountsEngine(cs::ContextLibrary<cs::AA>* context_lib,
    cs::Crf<cs::AA>* crf, cs::Pseudocounts<cs::AA>* pc_hhm_context_engine,
    cs::Admix* pc_hhm_context_mode,
    cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine,
    cs::Admix* pc_prefilter_context_mode) {

  if (context_lib != NULL)
    delete context_lib;
  if (crf != NULL)
    delete crf;
  if (pc_hhm_context_engine != NULL)
    delete pc_hhm_context_engine;
  if (pc_hhm_context_mode != NULL)
    delete pc_hhm_context_mode;
  if (pc_prefilter_context_engine != NULL)
    delete pc_prefilter_context_engine;
  if (pc_prefilter_context_mode != NULL)
    delete pc_prefilter_context_mode;
}




