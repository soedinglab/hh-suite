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
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << std::endl;
    std::cerr << "\t" << infile << " is empty!\n";
    exit(4);
  }
  while (strscn(line) == NULL)
    fgetline(line, LINELEN, inf); // skip lines that contain only white space

  // Is infile a HMMER file?
  if (!strncmp(line, "HMMER", 5)) {
    // Uncomment this line to allow HMMER2/HMMER3 models as queries:
    std::cerr
        << "Error: Use of HMMER format as input will result in severe loss of sensitivity!\n";
    exit(4);
  }
  // ... or is it an hhm file?
  else if (!strncmp(line, "NAME", 4) || !strncmp(line, "HH", 2)) {
    char path[NAMELEN];
    Pathname(path, infile);

    if (v >= 2)
      std::cout << "Query file is in HHM format\n";

    // Rewind to beginning of line and read query hhm file
    rewind(inf);
    q->Read(inf, par.maxcol, par.nseqdis, pb, path);
    input_format = 0;

    if (v >= 1 && input_format == 0)  // HHM format
      printf(
          "Extracting representative sequences from %s to merge later with matched database sequences\n",
          infile);

    Alignment ali_tmp;
    ali_tmp.GetSeqsFromHMM(q);
    ali_tmp.Compress(infile, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
    qali = ali_tmp;
  }
  // ... or is it an alignment file
  else if (line[0] == '#' || line[0] == '>') {
    if (par.calibrate) {
      std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      printf("\tonly HHM files can be calibrated.\n");
      printf(
          "\tBuild an HHM file from your alignment with 'hhmake -i %s' and rerun hhsearch with the hhm file\n\n",
          infile);
      exit(1);
    }

    if (v >= 2 && strcmp(infile, "stdin"))
      std::cout << infile << " is in A2M, A3M or FASTA format\n";

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
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << std::endl;
    std::cerr << "\tunrecognized input file format in \'" << infile << "\'\n";
    std::cerr << "\tline = " << line << "\n";
    exit(1);
  }

  if (v >= 2 && input_format == 0 && q->Neff_HMM > 11.0)
    fprintf(stderr,
        "WARNING: MSA %s looks too diverse (Neff=%.1f>11). Better check it with an alignment viewer for non-homologous segments. Also consider building the MSA with hhblits using the - option to limit MSA diversity.\n",
        q->name, q->Neff_HMM);
}

void ReadQueryFile(Parameters& par, char* infile, char& input_format,
    char use_global_weights, HMM* q, Alignment& qali, float* pb,
    const float S[20][20], const float Sim[20][20]) {
  // Open query file and determine file type
  char path[NAMELEN]; // path of input file (is needed to write full path and file name to HMM FILE record)
  FILE* inf = NULL;
  if (strcmp(infile, "stdin") == 0) {
    inf = stdin;
    if (v >= 2) {
      printf("Reading HMM / multiple alignment from standard input ...\n");
    }
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
  runSystem(command, v);
  command = (std::string) psipred + "/psipass2 " + (std::string) psipred_data
      + "/weights_p2.dat 1 0.98 1.09 " + (std::string) tmpfile + ".ss2 "
      + (std::string) tmpfile + ".ss > " + (std::string) tmpfile + ".horiz";
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
  runSystem(command, v);
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

void AlignByWorker(Parameters& par, Hit* hit, HMM* t, HMM* q, const int format,
    const float* pb, const float R[20][20],
    const float S73[NDSSP][NSSPRED][MAXCF],
    const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF], HitList& hitlist) {
  // Prepare q ant t and compare
  PrepareTemplateHMM(par, q, t, format, pb, R);

  // Do HMM-HMM comparison, store results if score>SMIN, and try next best alignment
  for (hit->irep = 1; hit->irep <= par.altali; hit->irep++) {
    if (par.forward == 0) {
      hit->Viterbi(q, t, par.loc, par.ssm, par.maxres, par.min_overlap,
          par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);
      if (hit->irep > 1 && hit->score <= SMIN)
        break;
      hit->Backtrace(q, t, par.corr, par.ssw, S73, S33);
    }
    else if (par.forward == 2) {
      hit->Forward(q, t, par.ssm, par.min_overlap, par.loc, par.shift, par.ssw,
          par.exclstr, S73, S33);
      hit->Backward(q, t, par.loc, par.shift, par.ssw, S73, S33);
      hit->MACAlignment(q, t, par.loc, par.mact, par.macins);
      hit->BacktraceMAC(q, t, par.corr, par.ssw, S73, S33);
    }
    hit->score_sort = hit->score_aass;
    if (hit->score <= SMIN)
      hit->lastrep = 1;
    else
      hit->lastrep = 0;


    if (par.early_stopping_filter) {
      // Calculate Evalue
      float q_len = log(q->L) / LOG1000;
      float hit_len = log(hit->L) / LOG1000;
      float q_neff = q->Neff_HMM / 10.0;
      float hit_neff = hit->Neff_HMM / 10.0;
      float lamda = lamda_NN(q_len, hit_len, q_neff, hit_neff);
      float mu = mu_NN(q_len, hit_len, q_neff, hit_neff);
      hit->logPval = logPvalue(hit->score, lamda, mu);

      float alpha = 0;
      float log_Pcut = log(par.prefilter_evalue_thresh / par.dbsize);
      float log_dbsize = log(par.dbsize);

      if (par.prefilter)
        alpha = par.alphaa
            + par.alphab * (hit_neff - 1) * (1 - par.alphac * (q_neff - 1));

      hit->Eval = exp(hit->logPval + log_dbsize + (alpha * log_Pcut));
      hit->logEval = hit->logPval + log_dbsize + (alpha * log_Pcut);

      #pragma omp critical
      {
        par.filter_sum -= par.filter_evals[par.filter_counter];
        par.filter_evals[par.filter_counter] = 1.0 / (1.0 + hit->Eval);
        par.filter_sum += par.filter_evals[par.filter_counter];

        par.filter_counter++;
        if (par.filter_counter == par.filter_length) {
          par.filter_counter = 0;
        }
      }
    }

    #pragma omp critical
    {
      hitlist.Push(*(hit)); // insert hit at beginning of list (last repeats first!)
    }

    // find only best alignment for forward algorithm and stochastic sampling
    if (par.forward > 0)
      break;

    // break if score for previous hit is already worse than SMIN
    if (hit->score <= SMIN)
      break;
  }
}

void PerformViterbiByWorker(Parameters& par, Hit* hit, HMM* t, HMM* q,
    const int format, const float* pb, const float R[20][20],
    const float S73[NDSSP][NSSPRED][MAXCF],
    const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF], HitList& hitlist,
    Hash<Hit>* previous_hits) {
  PrepareTemplateHMM(par, q, t, format, pb, R);

  for (hit->irep = 1; hit->irep <= par.altali; hit->irep++) {
    // Break, if no previous_hit with irep is found
    hit->Viterbi(q, t, par.loc, par.ssm, par.maxres, par.min_overlap, par.shift,
        par.egt, par.egq, par.ssw, par.exclstr, S73, S33);
    if (hit->irep > 1 && hit->score <= SMIN)
      break;
    hit->Backtrace(q, t, par.corr, par.ssw, S73, S33);

    hit->score_sort = hit->score_aass;

#pragma omp critical
    {
      std::stringstream ss_tmp;
      ss_tmp << hit->file << "__" << hit->irep;

      if (previous_hits
          && previous_hits->Contains((char*) ss_tmp.str().c_str())) {
        //printf("Previous hits contains %s!\n",(char*)ss_tmp.str().c_str());
        Hit hit_cur = previous_hits->Remove((char*) ss_tmp.str().c_str());
        previous_hits->Add((char*) ss_tmp.str().c_str(), *(hit));

        // Overwrite *hit[bin] with alignment, etc. of hit_cur
        hit_cur.score = hit->score;
        hit_cur.score_aass = hit->score_aass;
        hit_cur.score_ss = hit->score_ss;
        hit_cur.Pval = hit->Pval;
        hit_cur.Pvalt = hit->Pvalt;
        hit_cur.logPval = hit->logPval;
        hit_cur.logPvalt = hit->logPvalt;
        hit_cur.Eval = hit->Eval;
        hit_cur.logEval = hit->logEval;
        hit_cur.Probab = hit->Probab;

        hitlist.Push(hit_cur); // insert hit at beginning of list (last repeats first!)
      }
    }

    // break if score for first hit is already worse than SMIN
    if (hit->score <= SMIN)
      break;
  }
}

///////////////////////////////////////////////////////////////////////////////////////
//// Realign q and with t[bin] in all hits from same tempate using  MAC algorithm
//////////////////////////////////////////////////////////////////////////////////////
void RealignByWorker(Parameters& par, Hit* hit, HMM* q, HMM* t,
    const int format, const float* pb, const float R[20][20],
    const float S73[NDSSP][NSSPRED][MAXCF],
    const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]) {
  // Realign all hits with same template, pointed to by list List<void*>* hit[bin]->plist_phits;
  // This list is set up in HHseach and HHblits at the beginning of perform_realign()

  // Prepare MAC comparison(s)
  PrepareTemplateHMM(par, q, t, format, pb, R);
  t->Log2LinTransitionProbs(1.0);

  hit->irep = 1;
  hit->plist_phits->Reset();
  while (!hit->plist_phits->End()) {
    // Set pointer hit_cur to next hit to be realigned
    Hit* hit_cur = (Hit*) hit->plist_phits->ReadNext();
    // printf("Realigning %s, irep=%i\n",hit_cur->name,hit_cur->irep);  //?????????

    // Realign only around previous Viterbi hit
    // hit[bin] = *hit_cur; is not possible because the pointers to the DP matrices would be overwritten
    hit->i1 = hit_cur->i1;
    hit->i2 = hit_cur->i2;
    hit->j1 = hit_cur->j1;
    hit->j2 = hit_cur->j2;
    hit->nsteps = hit_cur->nsteps;
    hit->i = hit_cur->i;
    hit->j = hit_cur->j;
    hit->realign_around_viterbi = true;

    // Align q to template in *hit[bin]
    hit->Forward(q, t, par.ssm, par.min_overlap, par.loc, par.shift, par.ssw,
        par.exclstr, S73, S33);
    hit->Backward(q, t, par.loc, par.shift, par.ssw, S73, S33);
    hit->MACAlignment(q, t, par.loc, par.mact, par.macins);
    hit->BacktraceMAC(q, t, par.corr, par.ssw, S73, S33);

    // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
    hit->score = hit_cur->score;
    hit->score_ss = hit_cur->score_ss;
    hit->score_aass = hit_cur->score_aass;
    hit->score_sort = hit_cur->score_sort;
    hit->Pval = hit_cur->Pval;
    hit->Pvalt = hit_cur->Pvalt;
    hit->logPval = hit_cur->logPval;
    hit->logPvalt = hit_cur->logPvalt;
    hit->Eval = hit_cur->Eval;
    hit->logEval = hit_cur->logEval;
    hit->Probab = hit_cur->Probab;


    // Replace original hit in hitlist with realigned hit
    hit_cur->Delete(); // delete content of pointers etc. of hit_cur (but not DP matrices)
    *hit_cur = *hit; // copy all variables and pointers from *hit[bin] into hitlist

    hit->irep++;
  }
}

