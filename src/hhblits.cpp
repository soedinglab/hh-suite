/*
 * HHblits.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include "hhblits.h"

//TODO: get rid of exit(1)... throw errors
//TODO: get more flexible database reader
//TODO: add length to HHDatabaseEntry... for Martin

HHblits::HHblits(Parameters& parameters) {
  par = parameters;

  q = NULL;
  q_tmp = NULL;

  v = par.v;
  v1 = 0;

  N_searched = 0;

  //TODO: multiple databases
  for(size_t i = 0; i < par.db_bases.size(); i++) {
	  HHblitsDatabase* db = new HHblitsDatabase(par.db_bases[i].c_str());
	  dbs.push_back(db);
  }

  par.dbsize = 0;
  for (size_t i = 0; i < dbs.size(); i++) {
    par.dbsize += dbs[i]->cs219_database->db_index->n_entries;
  }

  if (par.prefilter) {
    for (size_t i = 0; i < dbs.size(); i++) {
      dbs[i]->initPrefilter(par.cs_library);
    }
  }

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix(par.matrix, pb, P, R, S, Sim);

  // Set secondary structure substitution matrix
  if (par.ssm)
    SetSecStrucSubstitutionMatrix(par.ssa, S73, S33);

  // Prepare pseudocounts
  if (!par.nocontxt && *par.clusterfile) {
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

  // Prepare multi-threading - reserve memory for threads, intialize, etc.
  omp_set_num_threads(par.threads);
  for (int bin = 0; bin < par.threads; bin++) {
    t[bin] = new HMM; // Each bin has a template HMM allocated that was read from the database file
    // Each bin has an object of type Hit allocated ...
    hit[bin] = new Hit;
    // ...with a separate dynamic programming matrix (memory!!)
    hit[bin]->AllocateBacktraceMatrix(par.maxres, par.maxres);
    // Allocate memory for matrix and set to 0
    hit[bin]->AllocateForwardMatrix(par.maxres, par.maxres);
  }
  format = new int[par.threads];
}

HHblits::~HHblits() {
  Reset();

  for (size_t i = 0; i < dbs.size(); i++) {
    delete dbs[i];
  }
  dbs.clear();

  for (int bin = 0; bin < par.threads; bin++) {
    hit[bin]->DeleteBacktraceMatrix(par.maxres);
    hit[bin]->DeleteForwardMatrix(par.maxres);
  }

  for (int bin = 0; bin < par.threads; bin++) {
    delete hit[bin];
    delete t[bin];
  }

  delete[] format;

  DeletePseudocountsEngine(context_lib, crf, pc_hhm_context_engine,
      pc_hhm_context_mode, pc_prefilter_context_engine,
      pc_prefilter_context_mode);
}

void HHblits::ProcessAllArguments(int argc, char** argv, Parameters& par) {
  par.argv = argv;
  par.argc = argc;

  par.premerge = 3;
  par.Ndiff = 1000;
  par.prefilter = true;

  //TODO
  par.early_stopping_filter = true;
  par.filter_thresh = 0.01;

  // Enable changing verbose mode before command line are processed
  for (int i = 1; i < argc; i++) {
    if (argc > 1 && !strcmp(argv[i], "-v0"))
      v = 0;
    else if (argc > 1 && !strcmp(argv[i], "-v1"))
      v = 1;
    else if (argc > 2 && !strcmp(argv[i], "-v"))
      v = atoi(argv[i + 1]);
  }

  par.v = v;

  par.SetDefaultPaths();

  // Process default otpions from .hhdefaults file
  char* argv_conf[MAXOPT];
  int argc_conf = 0;

  ReadDefaultsFile(argc_conf, argv_conf, argv[0]);
  ProcessArguments(argc_conf, argv_conf, par);

  for (int n = 1; n < argc_conf; n++)
    delete[] argv_conf[n];

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc, argv, par);

  // Check needed files
  if (!*par.infile || !strcmp(par.infile, "")) {
    help(par);
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << std::endl;
    std::cerr << "\tinput file missing!" << std::endl;
    exit(4);
  }
  if (par.db_bases.size() == 0) {
    help(par);
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    std::cerr << "\tdatabase missing (see -d)\n";
    exit(4);
  }
  if (par.addss == 1 && (!*par.psipred || !*par.psipred_data)) {
    help(par);
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << std::endl;
    std::cerr
        << "\tmissing PSIPRED directory (see -psipred and -psipred_data).\n"
        << std::endl;
    std::cerr
        << "\tIf you don't need the predicted secondary structure, don't use the -addss option!"
        << std::endl;
    exit(4);
  }
  if (!par.nocontxt) {
    if (!strcmp(par.clusterfile, "")) {
      help(par);
      std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      std::cerr << "\tcontext-specific library missing (see -contxt)\n";
      exit(4);
    }
    if (!strcmp(par.cs_library, "")) {
      help(par);
      std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      std::cerr << "\tcolumn state library (see -cslib)\n";
      exit(4);
    }
  }
  if (par.loc == 0 && par.num_rounds >= 2 && v >= 1) {
    std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": "
        << __func__ << ":" << std::endl;
    std::cerr
        << "\tusing -global alignment for iterative searches is deprecated since non-homologous sequence segments can easily enter the MSA and corrupt it.\n";
  }

  if (par.num_rounds < 1)
    par.num_rounds = 1;
  else if (par.num_rounds > 8) {
    if (v >= 1) {
      std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      std::cerr << "\tNumber of iterations (" << par.num_rounds
          << ") to large => Set to 8 iterations\n";
    }
    par.num_rounds = 8;
  }

  // Premerging can be very time-consuming on large database a3ms, such as from pdb70.
  // Hence it is only done when iteratively searching against uniprot20 or nr20 with their much smaller MSAs:
  if (!(par.num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile
      || *par.alisbasename))
    par.premerge = 0;

  // No outfile given? Name it basename.hhm
  //TODO check if no output at all is specified
//  if (!*par.outfile) {
//    RemoveExtension(par.outfile, par.infile);
//    strcat(par.outfile, ".hhr");
//    if (v >= 2)
//      cout << "Search results will be written to " << par.outfile << "\n";
//  }

  // Check option compatibilities
  if (par.nseqdis > MAXSEQDIS - 3 - par.showcons)
    par.nseqdis = MAXSEQDIS - 3 - par.showcons; //3 reserved for secondary structure
  if (par.aliwidth < 20)
    par.aliwidth = 20;
  if (par.pc_hhm_context_engine.pca < 0.001)
    par.pc_hhm_context_engine.pca = 0.001; // to avoid log(0)
  if (par.pc_prefilter_context_engine.pca < 0.001)
    par.pc_prefilter_context_engine.pca = 0.001; // to avoid log(0)
  if (par.b > par.B)
    par.B = par.b;
  if (par.z > par.Z)
    par.Z = par.z;
  if (par.maxmem < 1.0) {
    cerr << "WARNING: setting -maxmem to its minimum allowed value of 1.0\n";
    par.maxmem = 1.0;
  }
  if (par.mact >= 1.0)
    par.mact = 0.999;
  else if (par.mact < 0)
    par.mact = 0.0;
  if (par.macins >= 1.0)
    par.macins = 0.999;
  else if (par.macins < 0)
    par.macins = 0.0;
}

void HHblits::Reset() {
  if (q) {
    delete q;
    q = NULL;
  }

  if (q_tmp) {
    delete q_tmp;
    q_tmp = NULL;
  }

  hitlist.Reset();
  while (!hitlist.End())
    hitlist.Delete().Delete();

  reducedHitlist.Reset();
  while (!reducedHitlist.End())
    reducedHitlist.Delete().Delete();

  alis.clear();
}

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void HHblits::help(Parameters& par, char all) {
  char program_name[NAMELEN];
  RemovePathAndExtension(program_name, par.argv[0]);

  printf("\n");
  printf(
      "HHblits %s:\nHMM-HMM-based lightning-fast iterative sequence search\n",
      VERSION_AND_DATE);
  printf(
      "HHblits is a sensitive, general-purpose, iterative sequence search tool that represents\n");
  printf(
      "both query and database sequences by HMMs. You can search HHblits databases starting\n");
  printf(
      "with a single query sequence, a multiple sequence alignment (MSA), or an HMM. HHblits\n");
  printf(
      "prints out a ranked list of database HMMs/MSAs and can also generate an MSA by merging\n");
  printf("the significant database HMMs/MSAs onto the query MSA.\n");
  printf("\n");
  printf("%s", HHBLITS_REFERENCE);
  printf("%s", COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query [options] \n", program_name);
  printf(
      " -i <file>      input/query: single sequence or multiple sequence alignment (MSA)\n");
  printf(
      "                in a3m, a2m, or FASTA format, or HMM in hhm format\n");
  if (all) {
    printf("\n");
    printf("<file> may be 'stdin' or 'stdout' throughout.\n");
  }
  printf("\n");
  printf(
      "Options:                                                                       \n");
  printf(
      " -d <name>      database name (e.g. uniprot20_29Feb2012)                       \n");
  printf(
      " -n     [1,8]   number of iterations (default=%i)                              \n",
      par.num_rounds);
  printf(
      " -e     [0,1]   E-value cutoff for inclusion in result alignment (def=%G)      \n",
      par.e);
  printf("\n");
  printf(
      "Input alignment format:                                                       \n");
  printf(
      " -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf(
      "               ' -' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(
      " -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(
      " -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");
  printf("Output options: \n");
  printf(
      " -o <file>      write results in standard format to file (default=<infile.hhr>)\n");
  printf(
      " -ored <file>   write filtered and wiggled alignments in standard format to file\n");
  printf(
      " -oa3m <file>   write result MSA with significant matches in a3m format\n");
  if (!all) {
    printf("                Analogous for -opsi and -ohhm\n");
  }
  if (all) {
    printf(
        " -opsi <file>   write result MSA of significant matches in PSI-BLAST format\n");
    printf(
        " -ohhm <file>   write HHM file for result MSA of significant matches\n");
  }
  printf(" -oalis <name>  write MSAs in A3M format after each iteration\n");
  if (all) {
    printf(
        " -Ofas <file>   write pairwise alignments of significant matches in FASTA format\n");
    printf(
        "                Analogous for output in a3m and a2m format (e.g. -Oa3m)\n");
    printf(
        " -qhhm <file>   write query input HHM file of last iteration (default=off)      \n");
    printf(
        " -seq <int>     max. number of query/template sequences displayed (default=%i)  \n",
        par.nseqdis);
    printf(
        " -aliw <int>    number of columns per line in alignment list (default=%i)       \n",
        par.aliwidth);
    printf(
        " -p [0,100]     minimum probability in summary and alignment list (default=%G)  \n",
        par.p);
    printf(
        " -E [0,inf[     maximum E-value in summary and alignment list (default=%G)      \n",
        par.E);
    printf(
        " -Z <int>       maximum number of lines in summary hit list (default=%zu)        \n",
        par.Z);
    printf(
        " -z <int>       minimum number of lines in summary hit list (default=%zu)        \n",
        par.z);
    printf(
        " -B <int>       maximum number of alignments in alignment list (default=%zu)     \n",
        par.B);
    printf(
        " -b <int>       minimum number of alignments in alignment list (default=%zu)     \n",
        par.b);
    printf("\n");
    printf(
        "Prefilter options                                                               \n");
    printf(
        " -noprefilt                disable all filter steps                                        \n");
    printf(
        " -noaddfilter              disable all filter steps (except for fast prefiltering)         \n");
    printf(
        " -nodbfilter               disable additional filtering of prefiltered HMMs                \n");
    printf(
        " -noblockfilter            search complete matrix in Viterbi                               \n");
    printf(
        " -maxfilt                  max number of hits allowed to pass 2nd prefilter (default=%i)   \n",
        par.maxnumdb);
    printf(
        " -min_prefilter_hits       min number of hits to pass prefilter (default=%i)               \n",
        par.min_prefilter_hits);
    printf("\n");
  }
  printf(
      "Filter options applied to query MSA, database MSAs, and result MSA              \n");
  printf(
      " -all           show all sequences in result MSA; do not filter result MSA      \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (def=%i)\n",
      par.max_seqid);
  printf(
      " -diff [0,inf[  filter MSAs by selecting most diverse set of sequences, keeping \n");
  printf(
      "                at least this many seqs in each MSA block of length 50 (def=%i) \n",
      par.Ndiff);
  printf(
      " -cov  [0,100]  minimum coverage with master sequence (%%) (def=%i)             \n",
      par.coverage);
  printf(
      " -qid  [0,100]  minimum sequence identity with master sequence (%%) (def=%i)    \n",
      par.qid);
  printf(
      " -qsc  [0,100]  minimum score per column with master sequence (default=%.1f)    \n",
      par.qsc);
  printf(
      " -neff [1,inf]  target diversity of multiple sequence alignment (default=off)   \n");
  printf("\n");
  printf(
      "HMM-HMM alignment options:                                                       \n");
  printf(
      " -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)   \n");
  printf(
      " -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf(
      "                ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n",
      par.mact);
  printf(
      " -macins [0,1[  controls the cost of internal gap positions in the MAC algorithm.\n");
  printf(
      "                0:dense alignments  1:gappy alignments (default=%.2f)\n",
      par.macins);
  printf(
      " -glob/-loc     use global/local alignment mode for searching/ranking (def=local)\n");
  if (all) {
    printf(
        " -realign_max <int>  realign max. <int> hits (default=%zu)                        \n",
        par.realign_max);
    printf(
        " -alt <int>     show up to this many significant alternative alignments(def=%i)  \n",
        par.altali);
    printf(
        " -premerge <int> merge <int> hits to query MSA before aligning remaining hits (def=%i)\n",
        par.premerge);
    printf(
        " -shift [-1,1]  profile-profile score offset (def=%-.2f)                         \n",
        par.shift);
    printf(
        " -corr [0,1]    weight of term for pair correlations (def=%.2f)                \n",
        par.corr);
    printf(
        " -ssm {0,..,4}  0:   no ss scoring                                             \n");
    printf(
        "                1,2: ss scoring after or during alignment  [default=%1i]         \n",
        par.ssm);
    printf(
        "                3,4: ss scoring after or during alignment, predicted vs. predicted\n");
    printf(
        " -ssw [0,1]     weight of ss score  (def=%-.2f)                                  \n",
        par.ssw);
    printf(
        " -wg            use global sequence weighting for realignment!                   \n");
    printf("\n");
    printf(
        "Gap cost options:                                                                \n");
    printf(
        " -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                     \n",
        par.gapb);
    printf(
        " -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)    \n",
        par.gapd);
    printf(
        " -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)      \n",
        par.gape);
    printf(
        " -gapf ]0,inf]  factor to increase/reduce gap open penalty for deletes (def=%-.2f) \n",
        par.gapf);
    printf(
        " -gapg ]0,inf]  factor to increase/reduce gap open penalty for inserts (def=%-.2f) \n",
        par.gapg);
    printf(
        " -gaph ]0,inf]  factor to increase/reduce gap extend penalty for deletes(def=%-.2f)\n",
        par.gaph);
    printf(
        " -gapi ]0,inf]  factor to increase/reduce gap extend penalty for inserts(def=%-.2f)\n",
        par.gapi);
    printf(
        " -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f) \n",
        par.egq);
    printf(
        " -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)\n",
        par.egt);
    printf("\n");
    printf(
        "Pseudocount (pc) options:                                                        \n");
    printf(" Context specific hhm pseudocounts:\n");
    printf(
        "  -pc_hhm_contxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",
        par.pc_hhm_context_engine.admix);
    printf(
        "               0: no pseudo counts:    tau = 0                                  \n");
    printf(
        "               1: constant             tau = a                                  \n");
    printf(
        "               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf(
        "               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf(
        "               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf(
        "  -pc_hhm_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",
        par.pc_hhm_context_engine.pca);
    printf(
        "  -pc_hhm_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",
        par.pc_hhm_context_engine.pcb);
    printf(
        "  -pc_hhm_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",
        par.pc_hhm_context_engine.pcc);

    printf(
        " Context independent hhm pseudocounts (used for templates; used for query if contxt file is not available):\n");
    printf(
        "  -pc_hhm_nocontxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",
        par.pc_hhm_nocontext_mode);
    printf(
        "               0: no pseudo counts:    tau = 0                                  \n");
    printf(
        "               1: constant             tau = a                                  \n");
    printf(
        "               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    //  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf(
        "               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf(
        "  -pc_hhm_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",
        par.pc_hhm_nocontext_a);
    printf(
        "  -pc_hhm_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",
        par.pc_hhm_nocontext_b);
    printf(
        "  -pc_hhm_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",
        par.pc_hhm_nocontext_c);

    printf(" Context specific prefilter pseudocounts:\n");
    printf(
        "  -pc_prefilter_contxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",
        par.pc_prefilter_context_engine.admix);
    printf(
        "               0: no pseudo counts:    tau = 0                                  \n");
    printf(
        "               1: constant             tau = a                                  \n");
    printf(
        "               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf(
        "               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf(
        "               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf(
        "  -pc_prefilter_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",
        par.pc_prefilter_context_engine.pca);
    printf(
        "  -pc_prefilter_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",
        par.pc_prefilter_context_engine.pcb);
    printf(
        "  -pc_prefilter_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",
        par.pc_prefilter_context_engine.pcc);

    printf(
        " Context independent prefilter pseudocounts (used if context file is not available):\n");
    printf(
        "  -pc_prefilter_nocontxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",
        par.pc_prefilter_nocontext_mode);
    printf(
        "               0: no pseudo counts:    tau = 0                                  \n");
    printf(
        "               1: constant             tau = a                                  \n");
    printf(
        "               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    //  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf(
        "               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf(
        "  -pc_prefilter_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",
        par.pc_prefilter_nocontext_a);
    printf(
        "  -pc_prefilter_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",
        par.pc_prefilter_nocontext_b);
    printf(
        "  -pc_prefilter_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",
        par.pc_prefilter_nocontext_c);

    printf("\n");
    printf(
        " Context-specific pseudo-counts:                                                  \n");
    printf(
        "  -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
    printf(
        "  -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",
        par.clusterfile);
    //should not be in the section of pseudocounts ... associated to prefiltering
    printf("\n");
    printf("Predict secondary structure\n");
    printf(
        " -addss         add 2ndary structure predicted with PSIPRED to result MSA \n");
    printf(
        " -psipred <dir> directory with PSIPRED executables (default=%s)  \n",
        par.psipred);
    printf(" -psipred_data <dir>  directory with PSIPRED data (default=%s) \n",
        par.psipred_data);
    printf("\n");
  }
  printf(
      "Other options:                                                                   \n");
  printf(
      " -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose (def=%i)\n",
      v);
  printf(
      " -neffmax ]1,20] skip further search iterations when diversity Neff of query MSA \n");
  printf("                becomes larger than neffmax (default=%.1f)\n",
      par.neffmax);
  printf(
      " -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=%i)      \n",
      par.threads);
  if (all) {
    printf(
        " -scores <file> write scores for all pairwise comparisions to file               \n");
    printf(
        " -atab   <file> write all alignments in tabular layout to file                   \n");
    printf(" -maxres <int>  max number of HMM columns (def=%5i)             \n",
        par.maxres);
    printf(
        " -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n",
        par.maxmem);
  }
  printf("\n");
  if (!all) {
    printf(
        "An extended list of options can be obtained by calling 'hhblits -help'\n");
  }
  printf("\n");
  printf("Example: %s -i query.fas -oa3m query.a3m -n 1  \n", program_name);
  cout << endl;
}

void HHblits::ProcessArguments(int argc, char** argv, Parameters& par) {
  char program_name[NAMELEN];
  RemovePathAndExtension(program_name, par.argv[0]);

  //Processing command line input
  for (int i = 1; i < argc; i++) {
    if (v >= 4)
      cout << i << "  " << argv[i] << endl; //PRINT
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no query file following -i\n";
        exit(4);
      }
      else
        strcpy(par.infile, argv[i]);
    }
    else if (!strcmp(argv[i], "-d")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no database basename following -d\n";
        exit(4);
      }
      else {
    	  std::string db(argv[i]);
    	  par.db_bases.push_back(db);
      }
    }
    else if (!strcmp(argv[i], "-contxt") || !strcmp(argv[i], "-context_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no lib following -contxt\n";
        exit(4);
      }
      else
        strcpy(par.clusterfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-cslib") || !strcmp(argv[i], "-cs_lib")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no lib following -cslib\n";
        exit(4);
      }
      else
        strcpy(par.cs_library, argv[i]);
    }
    else if (!strcmp(argv[i], "-psipred")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no directory following -psipred\n";
        exit(4);
      }
      else
        strcpy(par.psipred, argv[i]);
    }
    else if (!strcmp(argv[i], "-psipred_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no database directory following -psipred_data\n";
        exit(4);
      }
      else
        strcpy(par.psipred_data, argv[i]);
    }
    else if (!strcmp(argv[i], "-o")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-ored")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.reduced_outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa3m")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -oa3m\n";
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-ohhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -ohhm\n";
        exit(4);
      }
      else
        strcpy(par.hhmfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-opsi")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -opsi\n";
        exit(4);
      }
      else
        strcpy(par.psifile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oalis")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no file basename following -oalis\n";
        exit(4);
      }
      else
        strcpy(par.alisbasename, argv[i]);
    }
    else if (!strcmp(argv[i], "-Ofas")) {
      par.append = 0;
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Oa2m")) {
      par.append = 0;
      par.outformat = 2;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Oa3m")) {
      par.append = 0;
      par.outformat = 3;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-qhhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no filename following -qhhm\n";
        exit(4);
      }
      else
        strcpy(par.query_hhmfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-scores")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no file following -scores\n";
        exit(4);
      }
      else {
        strcpy(par.scorefile, argv[i]);
      }
    }
    else if (!strcmp(argv[i], "-atab")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
            << ": no file following -atab\n";
        exit(4);
      }
      else {
        strcpy(par.alitabfile, argv[i]);
      }
    }
    else if (!strcmp(argv[i], "-atab_scop"))
      par.alitab_scop = true;
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
      help(par, 1);
      exit(0);
    }
    else if (!strcmp(argv[i], "-v") && (i < argc - 1) && argv[i + 1][0] != '-')
      par.v = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-v"))
      par.v = 2;
    else if (!strcmp(argv[i], "-v0"))
      par.v = 0;
    else if (!strcmp(argv[i], "-v1"))
      par.v = 1;
    else if (!strcmp(argv[i], "-n") && (i < argc - 1))
      par.num_rounds = atoi(argv[++i]);
    else if (!strncmp(argv[i], "-BLOSUM", 7)
        || !strncmp(argv[i], "-Blosum", 7)) {
      if (!strcmp(argv[i] + 7, "30"))
        par.matrix = 30;
      else if (!strcmp(argv[i] + 7, "40"))
        par.matrix = 40;
      else if (!strcmp(argv[i] + 7, "50"))
        par.matrix = 50;
      else if (!strcmp(argv[i] + 7, "62"))
        par.matrix = 62;
      else if (!strcmp(argv[i] + 7, "65"))
        par.matrix = 65;
      else if (!strcmp(argv[i] + 7, "80"))
        par.matrix = 80;
      else
        cerr << endl << "WARNING: Ignoring unknown option " << argv[i]
            << " ...\n";
    }
    else if (!strcmp(argv[i], "-M") && (i < argc - 1))
      if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m"))
        par.M = 1;
      else if (!strcmp(argv[i], "first"))
        par.M = 3;
      else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
        par.Mgaps = atoi(argv[i]);
        par.M = 2;
      }
      else
        cerr << endl << "WARNING: Ignoring unknown argument: -M " << argv[i]
            << "\n";
    else if (!strcmp(argv[i], "-p") && (i < argc - 1))
      par.p = atof(argv[++i]);
    else if (!strcmp(argv[i], "-P") && (i < argc - 1))
      par.p = atof(argv[++i]);
    else if (!strcmp(argv[i], "-E") && (i < argc - 1))
      par.E = atof(argv[++i]);
    else if (!strcmp(argv[i], "-b") && (i < argc - 1))
      par.b = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-B") && (i < argc - 1))
      par.B = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-z") && (i < argc - 1))
      par.z = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-Z") && (i < argc - 1))
      par.Z = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-realign_max") && (i < argc - 1))
      par.realign_max = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-e") && (i < argc - 1))
      par.e = atof(argv[++i]);
    else if (!strncmp(argv[i], "-nopred", 7) || !strncmp(argv[i], "-noss", 5))
      par.showpred = 0;
    else if (!strncmp(argv[i], "-addss", 6))
      par.addss = 1;
    else if (!strcmp(argv[i], "-seq") && (i < argc - 1))
      par.nseqdis = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-aliw") && (i < argc - 1))
      par.aliwidth = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-id") && (i < argc - 1))
      par.max_seqid = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-qid") && (i < argc - 1))
      par.qid = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-qsc") && (i < argc - 1))
      par.qsc = atof(argv[++i]);
    else if (!strcmp(argv[i], "-cov") && (i < argc - 1))
      par.coverage = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-diff") && (i < argc - 1))
      par.Ndiff = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-all") || !strcmp(argv[i], "-nodiff")) {
      par.allseqs = true;
    }
    else if (!strcmp(argv[i], "-neffmax") && (i < argc - 1))
      par.neffmax = atof(argv[++i]);
    else if ((!strcmp(argv[i], "-neff") || !strcmp(argv[i], "-Neff"))
        && (i < argc - 1))
      par.Neff = atof(argv[++i]);
    //pc hhm context variables
    else if (!strcmp(argv[i], "-pc_hhm_contxt_mode") && (i < argc - 1))
      par.pc_hhm_context_engine.admix = (Pseudocounts::Admix) atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_contxt_a") && (i < argc - 1))
      par.pc_hhm_context_engine.pca = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_contxt_b") && (i < argc - 1))
      par.pc_hhm_context_engine.pcb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_contxt_c") && (i < argc - 1))
      par.pc_hhm_context_engine.pcc = atof(argv[++i]);
    //pc prefilter context variables
    else if (!strcmp(argv[i], "-pc_prefilter_contxt_mode") && (i < argc - 1))
      par.pc_prefilter_context_engine.admix = (Pseudocounts::Admix) atoi(
          argv[++i]);
    else if (!strcmp(argv[i], "-pc_prefilter_contxt_a") && (i < argc - 1))
      par.pc_prefilter_context_engine.pca = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_prefilter_contxt_b") && (i < argc - 1))
      par.pc_prefilter_context_engine.pcb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_prefilter_contxt_c") && (i < argc - 1))
      par.pc_prefilter_context_engine.pcc = atof(argv[++i]);
    //pc hhm nocontext variables
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_mode") && (i < argc - 1))
      par.pc_hhm_nocontext_mode = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_a") && (i < argc - 1))
      par.pc_hhm_nocontext_a = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_b") && (i < argc - 1))
      par.pc_hhm_nocontext_b = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_c") && (i < argc - 1))
      par.pc_hhm_nocontext_c = atof(argv[++i]);
    //pc prefilter nocontext variables
    else if (!strcmp(argv[i], "-pc_prefilter_nocontxt_mode") && (i < argc - 1))
      par.pc_prefilter_nocontext_mode = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pre_pca") && (i < argc - 1))
      par.pc_hhm_nocontext_a = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pre_pcb") && (i < argc - 1))
      par.pc_hhm_nocontext_b = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pre_pcc") && (i < argc - 1))
      par.pc_hhm_nocontext_c = atof(argv[++i]);

    else if (!strcmp(argv[i], "-gapb") && (i < argc - 1)) {
      par.gapb = atof(argv[++i]);
      if (par.gapb <= 0.01)
        par.gapb = 0.01;
    }
    else if (!strcmp(argv[i], "-gapd") && (i < argc - 1))
      par.gapd = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gape") && (i < argc - 1))
      par.gape = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapf") && (i < argc - 1))
      par.gapf = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapg") && (i < argc - 1))
      par.gapg = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gaph") && (i < argc - 1))
      par.gaph = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapi") && (i < argc - 1))
      par.gapi = atof(argv[++i]);
    else if (!strcmp(argv[i], "-egq") && (i < argc - 1))
      par.egq = atof(argv[++i]);
    else if (!strcmp(argv[i], "-egt") && (i < argc - 1))
      par.egt = atof(argv[++i]);
    else if (!strcmp(argv[i], "-alphaa") && (i < argc - 1))
      par.alphaa = atof(argv[++i]);
    else if (!strcmp(argv[i], "-alphab") && (i < argc - 1))
      par.alphab = atof(argv[++i]);
    else if (!strcmp(argv[i], "-alphac") && (i < argc - 1))
      par.alphac = atof(argv[++i]);
    else if (!strcmp(argv[i], "-noprefilt")) {
      par.prefilter = false;
      par.already_seen_filter = false;
    }
    else if (!strcmp(argv[i], "-noaddfilter")) {
      par.already_seen_filter = false;
    }
    else if (!strcmp(argv[i], "-maxfilt") && (i < argc - 1))
      par.maxnumdb = par.maxnumdb_no_prefilter = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-min_prefilter_hits") && (i < argc - 1))
      par.min_prefilter_hits = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-prepre_smax_thresh") && (i < argc - 1))
      par.preprefilter_smax_thresh = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pre_evalue_thresh") && (i < argc - 1))
      par.prefilter_evalue_thresh = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pre_bitfactor") && (i < argc - 1))
      par.prefilter_bit_factor = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pre_gap_open") && (i < argc - 1))
      par.prefilter_gap_open = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pre_gap_extend") && (i < argc - 1))
      par.prefilter_gap_extend = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pre_score_offset") && (i < argc - 1))
      par.prefilter_score_offset = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-realignoldhits"))
      par.realign_old_hits = true;
    else if (!strcmp(argv[i], "-realign"))
      par.realign = 1;
    else if (!strcmp(argv[i], "-norealign"))
      par.realign = 0;
    else if (!strcmp(argv[i], "-ssm") && (i < argc - 1))
      par.ssm = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ssw") && (i < argc - 1))
      par.ssw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-wg")) {
      par.wg = 1;
    }
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = 2 * par.maxres;
    }
    else if (!strncmp(argv[i], "-glo", 3)) {
      par.loc = 0;
      if (par.mact > 0.35 && par.mact < 0.3502) {
        par.mact = 0;
      }
    }
    else if (!strncmp(argv[i], "-loc", 4))
      par.loc = 1;
    else if (!strncmp(argv[i], "-alt", 4) && (i < argc - 1))
      par.altali = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-shift") && (i < argc - 1))
      par.shift = atof(argv[++i]);
    else if ((!strcmp(argv[i], "-mact") || !strcmp(argv[i], "-mapt"))
        && (i < argc - 1))
      par.mact = atof(argv[++i]);
    else if (!strcmp(argv[i], "-macins") && (i < argc - 1))
      par.macins = atof(argv[++i]);
    else if (!strcmp(argv[i], "-scwin") && (i < argc - 1)) {
      par.columnscore = 5;
      par.half_window_size_local_aa_bg_freqs = std::max(1, atoi(argv[++i]));
    }
    else if (!strncmp(argv[i], "-cpu", 4) && (i < argc - 1)) {
      par.threads = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i], "-maxmem") && (i < argc - 1)) {
      par.maxmem = atof(argv[++i]);
    }
    else if (!strncmp(argv[i], "-premerge", 9) && (i < argc - 1))
      par.premerge = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-nocontxt"))
      par.nocontxt = 1;
    else if (!strcmp(argv[i], "-csb") && (i < argc - 1))
      par.csb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-csw") && (i < argc - 1))
      par.csw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-corr") && (i < argc - 1))
      par.corr = atof(argv[++i]);
    else
      cerr << endl << "WARNING: Ignoring unknown option " << argv[i]
          << " ...\n";
    if (v >= 4)
      cout << i << "  " << argv[i] << endl; //PRINT
  } // end of for-loop for command line input
}

void HHblits::getTemplateA3M(HHblitsDatabase* db, char* entry_name,
    long& ftellpos, Alignment& tali) {
  if (db->use_compressed) {
    ffindex_entry_t* entry = ffindex_get_entry_by_name(
        db->ca3m_database->db_index, entry_name);

    char* data = ffindex_get_data_by_entry(db->ca3m_database->db_data, entry);

    if (data == NULL) {
      std::cerr << "Could not fetch data for a3m " << entry_name << "!"
          << std::endl;
      exit(4);
    }

    ftellpos = entry->offset;
    tali.ReadCompressed(entry, data, db->sequence_database->db_index,
        db->sequence_database->db_data, db->header_database->db_index,
        db->header_database->db_data, par.mark, par.maxcol);
  }
  else {
    FILE* dbf = ffindex_fopen_by_name(db->a3m_database->db_data,
        db->a3m_database->db_index, entry_name);

    if (dbf == NULL) {
      cerr << endl << "Error: opening A3M " << entry_name << std::endl;
      exit(4);
    }

    ftellpos = ftell(dbf);

    char line[LINELEN];
    if (!fgetline(line, LINELEN, dbf)) {
      std::cerr << "this should not happen!" << std::endl;
      //TODO: throw error
    }

    while (strscn(line) == NULL)
      fgetline(line, LINELEN, dbf); // skip lines that contain only white space

    tali.Read(dbf, entry_name, par.mark, par.maxcol, par.nseqdis, line);
    fclose(dbf);
  }

  tali.Compress(entry_name, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
}

HHblitsDatabase* HHblits::getHHblitsDatabase(HHDatabaseEntry& entry,
    std::vector<HHblitsDatabase*>& dbs) {
  for (std::vector<HHblitsDatabase*>::size_type i = 0; i < dbs.size(); i++) {
    if (dbs[i]->id == entry.ffdatabase->superId) {
      return dbs[i];
    }
  }

  return NULL;
}

void HHblits::getTemplateHMM(HHDatabaseEntry& entry, char use_global_weights,
    long& ftellpos, int& format, HMM* t) {
  if (entry.ffdatabase->isCompressed) {
    Alignment tali;

    char* data = ffindex_get_data_by_entry(entry.ffdatabase->db_data,
        entry.entry);

    if (data == NULL) {
      std::cerr << "Could not fetch data for a3m " << entry.entry->name << "!"
          << std::endl;
      exit(4);
    }

    //TODO: get rid of ftellpos
    ftellpos = entry.entry->offset;

    HHblitsDatabase* db = getHHblitsDatabase(entry, dbs);
    if (db == NULL) {
      //TODO throw error
      std::cerr << "this should not happen!!!!" << std::endl;
      exit(0);
    }

    tali.ReadCompressed(entry.entry, data, db->sequence_database->db_index,
        db->sequence_database->db_data, db->header_database->db_index,
        db->header_database->db_data, par.mark, par.maxcol);

    tali.Compress(entry.entry->name, par.cons, par.maxres, par.maxcol, par.M,
        par.Mgaps);

    tali.N_filtered = tali.Filter(par.max_seqid_db, S, par.coverage_db,
        par.qid_db, par.qsc_db, par.Ndiff_db);
    t->name[0] = t->longname[0] = t->fam[0] = '\0';
    tali.FrequenciesAndTransitions(t, use_global_weights, par.mark, par.cons,
        par.showcons, par.maxres, pb, Sim);

    format = 0;
  }
  else {
    FILE* dbf = ffindex_fopen_by_entry(entry.ffdatabase->db_data, entry.entry);

    if (dbf != NULL) {
      ftellpos = ftell(dbf); // record position in dbfile of next HMM to be read

      char line[LINELEN];
      if (!fgetline(line, LINELEN, dbf)) {
        std::cerr << "this should not happen!" << std::endl;
        //TODO: throw error
      }

      while (strscn(line) == NULL)
        fgetline(line, LINELEN, dbf); // skip lines that contain only white space

      if (!strncmp(line, "HMMER3", 6))      // read HMMER3 format
          {
        format = 1;
        t->ReadHMMer3(dbf, par.showcons, pb, entry.entry->name);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HMMER", 5))      // read HMMER format
          {
        format = 1;
        t->ReadHMMer(dbf, par.showcons, pb, entry.entry->name);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HH", 2))    // read HHM format
          {
        char path[NAMELEN];
        Pathname(path, entry.entry->name);

        format = 0;
        t->Read(dbf, par.maxcol, par.nseqdis, pb, path);

      }
      //TODO: old hhm format discarded
      // read a3m alignment
      else if (line[0] == '#' || line[0] == '>') {
        Alignment tali;
        tali.Read(dbf, entry.entry->name, par.mark, par.maxcol, par.nseqdis,
            line);
        tali.Compress(entry.entry->name, par.cons, par.maxres, par.maxcol,
            par.M, par.Mgaps);
        //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
        tali.N_filtered = tali.Filter(par.max_seqid_db, S, par.coverage_db,
            par.qid_db, par.qsc_db, par.Ndiff_db);
        t->name[0] = t->longname[0] = t->fam[0] = '\0';
        tali.FrequenciesAndTransitions(t, use_global_weights, par.mark,
            par.cons, par.showcons, par.maxres, pb, Sim);
        format = 0;
      }
      else {
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tunrecognized HMM file format in \'" << entry.entry->name
            << "\'. \n";
        cerr << "Context:\n'" << line << "\n";
        fgetline(line, LINELEN, dbf);
        cerr << line << "\n";
        fgetline(line, LINELEN, dbf);
        cerr << line << "'\n";
        exit(1);
      }

      fclose(dbf);
      return;
    }
  }
}

void HHblits::DoViterbiSearch(std::vector<HHDatabaseEntry*>& prefiltered_hits,
    Hash<Hit>* previous_hits, bool alignByWorker) {

  // Search databases
  for (int bin = 0; bin < par.threads; bin++) {
    hit[bin]->realign_around_viterbi = false;
  }

  double filter_cutoff = par.filter_length * par.filter_thresh;
  par.filter_sum = par.filter_length;
  par.filter_counter = 0;
  par.filter_evals = new double[par.filter_length];

  for (int a = 0; a < par.filter_length; a++) {
    par.filter_evals[a] = 1;
  }

  // For all the databases comming through prefilter
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t idb = 0; idb < prefiltered_hits.size(); idb++) {
    // Allocate free bin (no need to lock, since slave processes cannot change FREE to other status)
    int bin = omp_get_thread_num();

    if (!(par.early_stopping_filter && par.filter_sum < filter_cutoff)) {
      #pragma omp critical
      {
        hit[bin]->index = N_searched; // give hit a unique index for HMM
        ++N_searched;
      }

      getTemplateHMM(*prefiltered_hits[idb], 1, hit[bin]->ftellpos, format[bin],
          t[bin]);

      if (v >= 4)
        printf("Aligning with %s\n", t[bin]->name);

      hit[bin]->dbfile =
          new char[strlen(prefiltered_hits[idb]->entry->name) + 1];
      strcpy(hit[bin]->dbfile, prefiltered_hits[idb]->entry->name);
      hit[bin]->entry = prefiltered_hits[idb];

      if (alignByWorker)
        AlignByWorker(par, hit[bin], t[bin], q, format[bin], pb, R, S73, S33, hitlist);
      else
        PerformViterbiByWorker(par, hit[bin], t[bin], q, format[bin], pb, R,
            S73, S33, hitlist, previous_hits);
    }
  }

  delete[] par.filter_evals;
}

void HHblits::ViterbiSearch(std::vector<HHDatabaseEntry*>& prefiltered_hits,
    Hash<Hit>* previous_hits, int db_size) {
  // Initialize
  v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2;

  //hand over number of HMMs scanned to hitlist (for E-value calculation)
  hitlist.N_searched = db_size;

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  DoViterbiSearch(prefiltered_hits, previous_hits);

  if (v1 >= 2)
    cout << "\n";
  v = v1;

  // Sort list according to sortscore
  if (v >= 3)
    printf("Sorting hit list ...\n");
  hitlist.SortList();

  // Use NN prediction of lamda and mu
  hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);

  // Calculate E-values as combination of P-value for Viterbi HMM-HMM comparison and prefilter E-value: E = Ndb P (Epre/Ndb)^alpha
  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa, par.alphab,
        par.alphac, par.prefilter_evalue_thresh);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant of ViterbiSearch() function for rescoring previously found HMMs with Viterbi algorithm.
// Perform Viterbi search on each hit object in global hash previous_hits, but keep old alignment
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HHblits::RescoreWithViterbiKeepAlignment(int db_size,
    Hash<Hit>* previous_hits) {
  // Initialize
  v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2;

  std::vector<HHDatabaseEntry*> hits_to_rescore;

  // Get dbfiles of previous hits
  previous_hits->Reset();
  while (!previous_hits->End()) {
    Hit hit_cur = previous_hits->ReadNext();
    if (hit_cur.irep == 1) {
      hits_to_rescore.push_back(hit_cur.entry);
    }
  }

  hitlist.N_searched = db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  DoViterbiSearch(hits_to_rescore, previous_hits, false);

  if (v1 >= 2)
    cout << "\n";
  v = v1;

  // Sort list according to sortscore
  if (v >= 3)
    printf("Sorting hit list ...\n");
  hitlist.SortList();

  hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw); // Use NN prediction of lamda and mu

  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa, par.alphab,
        par.alphac, par.prefilter_evalue_thresh);
}

void HHblits::perform_realign(std::vector<HHDatabaseEntry*>& hits_to_realign, const size_t premerge,
    Hash<char>* premerged_hits) {
  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
  size_t nhits = 0;
  size_t N_aligned = 0;

  // Longest allowable length of database HMM (backtrace: 5 chars, fwd, bwd: 1 double
  long int Lmaxmem = (par.maxmem * 1024 * 1024 * 1024) / sizeof(double) / q->L
      / par.threads;
  long int Lmax = 0;      // length of longest HMM to be realigned

  // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
  // This list can be sorted by ftellpos to access one template after the other efficiently during realignment
  Hash<List<Realign_hitpos>*>* phash_plist_realignhitpos;
  phash_plist_realignhitpos = new Hash<List<Realign_hitpos>*>(100031, NULL);

  // Some templates have several (suboptimal) alignments in hitlist. For realignment, we need to efficiently
  // access all hit objects in hitlist belonging to one template (because we don't want to read templates twice)
  // We therefore need for each template (identified by its index between 0 and N_searched-1) a list of elements
  // in hitlist that store the alignments with the template of that index.
  // This list is pointed to by array_plist_phits[index].
  List<void*>** array_plist_phits;
  array_plist_phits = new List<void*>*[N_searched];
  for (int index = 0; index < N_searched; index++)
    array_plist_phits[index] = NULL; // initialize

  // Store all dbfiles and ftell positions of templates to be displayed and realigned
  hitlist.Reset();
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();
    if (nhits >= par.realign_max && nhits >= std::max(par.B, par.Z))
      break;
    if (hit_cur.Eval > par.e) {
      if (nhits >= std::max(par.B, par.Z))
        continue;
      if (nhits >= std::max(par.b, par.z) && hit_cur.Probab < par.p)
        continue;
      if (nhits >= std::max(par.b, par.z) && hit_cur.Eval > par.E)
        continue;
    }

    if (hit_cur.L > Lmax)
      Lmax = hit_cur.L;
    if (hit_cur.L > Lmaxmem) {
      nhits++;
      continue;
    }

    //fprintf(stderr,"hit.name=%-15.15s  hit.index=%-5i hit.ftellpos=%-8i  hit.dbfile=%s\n",hit_cur.name,hit_cur.index,(unsigned int)hit_cur.ftellpos,hit_cur.dbfile);

    // realign the first premerge hits consecutively to query profile
    if (nhits >= premerge) {
      if (hit_cur.irep == 1) {
        // For each template (therefore irep==1), store template index and position on disk in a list
        Realign_hitpos realign_hitpos;
        realign_hitpos.ftellpos = hit_cur.ftellpos; // stores position on disk of template for current hit
        realign_hitpos.index = hit_cur.index; // stores index of template of current hit
        realign_hitpos.entry = hit_cur.entry;
        if (!phash_plist_realignhitpos->Contains(hit_cur.dbfile)) {
          List<Realign_hitpos>* newlist = new List<Realign_hitpos>;
          phash_plist_realignhitpos->Add(hit_cur.dbfile, newlist);
        }
        // Add template index and ftellpos to list which belongs to key dbfile in hash
        phash_plist_realignhitpos->Show(hit_cur.dbfile)->Push(realign_hitpos);
      }
      // pointer at index is still NULL
      if (!array_plist_phits[hit_cur.index]) {
        List<void*>* newlist = new List<void*>; // create new list of pointers to all aligments of a template
        array_plist_phits[hit_cur.index] = newlist; // set array[index] to newlist
      }
      // Push(hitlist.ReadCurrentAddress()) :  Add address of current hit in hitlist to list...
      // array_plist_phits[hit_cur.index]-> :  pointed to by hit_cur.index'th element of array_plist_phits
      array_plist_phits[hit_cur.index]->Push(hitlist.ReadCurrentAddress());
    }

    nhits++;
  }
  if (v >= 2)
    printf(
        "Realigning %zu HMM-HMM alignments using Maximum Accuracy algorithm\n",
        nhits);

  if (Lmax > Lmaxmem) {
    Lmax = Lmaxmem;
    if (v >= 1) {
      cerr << "WARNING: Realigning sequences only up to length " << Lmaxmem
          << "." << endl;
      cerr
          << "This is genarally unproboblematic but may lead to slightly sub-optimal alignments for longer sequences."
          << endl;
      cerr
          << "You can increase available memory using the -maxmem <GB> option (currently "
          << par.maxmem << " GB)." << endl; // still to be implemented
      cerr
          << "The maximum length realignable is approximately maxmem/query_length/(cpus+1)/8B."
          << endl;
    }
  }

  v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2; // Supress verbose output during iterative realignment and realignment

  //////////////////////////////////////////////////////////////////////////////////
  // start premerge:
  // Align the first premerge templates
  if (premerge > 0) {
    if (v >= 2)
      printf("Merging %zu best hits to query alignment ...\n", premerge);

    int bin = 0;
    nhits = 0;
    hitlist.Reset();

    while (!hitlist.End() && nhits < premerge) {
      Hit hit_cur = hitlist.ReadNext();
      // JS: removed bug on 13 Feb 13 due to which premerged hits with E-value > par.e were not realigned
      if (hit_cur.Eval > par.e) {
        if (nhits >= std::max(par.B, par.Z))
          break;
        if (nhits >= std::max(par.b, par.z) && hit_cur.Probab < par.p)
          break;
        if (nhits >= std::max(par.b, par.z) && hit_cur.Eval > par.E)
          continue;
      }
      nhits++;

      if (hit_cur.L > Lmaxmem)
        continue;  // Don't align too long sequences due to memory limit

      // Forward stream position to start of next database HMM to be realigned
      hit[bin]->index = hit_cur.index; // give hit a unique index for HMM
      hit[bin]->dbfile = new char[strlen(hit_cur.dbfile) + 1];
      strcpy(hit[bin]->dbfile, hit_cur.dbfile); // record db file name from which next HMM is read
      hit[bin]->irep = 1; // Needed for min_overlap calculation in InitializeForAlignment in hhhit.C

      getTemplateHMM(*hit_cur.entry, par.wg, hit[bin]->ftellpos, format[bin],
          t[bin]);

      if (v >= 2)
        fprintf(stderr, "Realigning with %s ***** \n", t[bin]->name);
      ///////////////////////////////////////////////////

      N_aligned++;
      if (v1 >= 2 && !(N_aligned % 10)) {
        cout << ".";
        if (!(N_aligned % 500))
          printf(" %-4zu HMMs aligned\n", N_aligned);
        cout.flush();
      }

      // Prepare MAC comparison(s)
      PrepareTemplateHMM(par, q, t[bin], format[bin], pb, R);
      t[bin]->Log2LinTransitionProbs(1.0);

      // Realign only around previous Viterbi hit
      hit[bin]->i1 = hit_cur.i1;
      hit[bin]->i2 = hit_cur.i2;
      hit[bin]->j1 = hit_cur.j1;
      hit[bin]->j2 = hit_cur.j2;
      hit[bin]->nsteps = hit_cur.nsteps;
      hit[bin]->i = hit_cur.i;
      hit[bin]->j = hit_cur.j;
      hit[bin]->realign_around_viterbi = true;

      // Align q to template in *hit[bin]
      hit[bin]->Forward(q, t[bin], par.ssm, par.min_overlap, par.loc, par.shift,
          par.ssw, par.exclstr, S73, S33);
      hit[bin]->Backward(q, t[bin], par.loc, par.shift, par.ssw, S73, S33);
      hit[bin]->MACAlignment(q, t[bin], par.loc, par.mact, par.macins);
      hit[bin]->BacktraceMAC(q, t[bin], par.corr, par.ssw, S73, S33);

      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
      hit[bin]->score = hit_cur.score;
      hit[bin]->score_ss = hit_cur.score_ss;
      hit[bin]->score_aass = hit_cur.score_aass;
      hit[bin]->score_sort = hit_cur.score_sort;
      hit[bin]->Pval = hit_cur.Pval;
      hit[bin]->Pvalt = hit_cur.Pvalt;
      hit[bin]->logPval = hit_cur.logPval;
      hit[bin]->logPvalt = hit_cur.logPvalt;
      hit[bin]->Eval = hit_cur.Eval;
      hit[bin]->logEval = hit_cur.logEval;
      hit[bin]->Probab = hit_cur.Probab;
      hit[bin]->irep = hit_cur.irep;

      hit[bin]->entry = hit_cur.entry;

      // Replace original hit in hitlist with realigned hit
      //hitlist.ReadCurrent().Delete();
      hitlist.Delete().Delete();  // delete the list record and hit object
      hitlist.Insert(*hit[bin]);

      // merge only when hit length > MINCOLS_REALIGN (don't merge 1 column matches)
      if (hit[bin]->matched_cols < MINCOLS_REALIGN)
        continue;

      // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
      // Reading in next db MSA and merging it onto Qali

      HHblitsDatabase* db = getHHblitsDatabase(*hit_cur.entry, dbs);
      if (db == NULL) {
        std::cerr << "Could not find database for premerge!" << std::endl;
        continue;
      }

      Alignment Tali;
      long ftellpos;
      getTemplateA3M(db, hit[bin]->entry->entry->name, ftellpos, Tali);

      if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
        Qali_allseqs.MergeMasterSlave(*hit[bin], Tali, hit[bin]->dbfile,
            par.maxcol);

      Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
          par.qid_db, par.qsc_db, par.Ndiff_db);

      Qali.MergeMasterSlave(*hit[bin], Tali, hit[bin]->dbfile, par.maxcol);

      // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
      Qali.Compress("merged A3M file", par.cons, par.maxres, par.maxcol, par.M,
          par.Mgaps);

      // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
      Qali.N_filtered = Qali.Filter(par.max_seqid, S, par.coverage, par.qid,
          par.qsc, par.Ndiff);

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons,
          par.showcons, par.maxres, pb, Sim);

      stringstream ss_tmp;
      ss_tmp << hit[bin]->file << "__" << hit[bin]->irep;
      premerged_hits->Add((char*) ss_tmp.str().c_str());

      if (par.notags)
        q->NeutralizeTags(pb);

      // Compute substitution matrix pseudocounts?
      if (par.nocontxt) {
        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        q->PreparePseudocounts(R);
        // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        q->AddAminoAcidPseudocounts(par.pc_hhm_nocontext_mode,
            par.pc_hhm_nocontext_a, par.pc_hhm_nocontext_b,
            par.pc_hhm_nocontext_c);
      }
      else {
        // Add full context specific pseudocounts to query
        q->AddContextSpecificPseudocounts(pc_hhm_context_engine,
            pc_hhm_context_mode);
      }

      q->CalculateAminoAcidBackground(pb);
      if (par.columnscore == 5 && !q->divided_by_local_bg_freqs)
        q->DivideBySqrtOfLocalBackgroundFreqs(
            par.half_window_size_local_aa_bg_freqs, pb);

      // Transform transition freqs to lin space if not already done
      q->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg,
          par.gaph, par.gapi, par.gapb, par.gapb);
      q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
    }
  }
  // end premerge
  //////////////////////////////////////////////////////////////////////////////////

  // Read all HMMs whose position is given in phash_plist_realignhitpos
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t idb = 0; idb < hits_to_realign.size(); idb++) {
    // Can we skip dbfiles[idb] because it contains no template to be realigned?
    if (!phash_plist_realignhitpos->Contains(hits_to_realign[idb]->entry->name))
      continue;

    // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
    // This list is now sorted by ftellpos in ascending order to access one template after the other efficiently
    phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->SortList();

    ///////////////////////////////////////////////////////////////////////////////////////
    // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
    phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->Reset();
    while (!phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->End()) {
      // Submit jobs until no bin is free anymore
      while (!phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->End()) {
        // Allocate free bin
        int bin = omp_get_thread_num();

        // Forward stream position to start of next database HMM to be realigned
        Realign_hitpos hitpos_curr = phash_plist_realignhitpos->Show(
            hits_to_realign[idb]->entry->name)->ReadNext();
        hit[bin]->index = hitpos_curr.index; // give hit[bin] a unique index for HMM
        hit[bin]->entry = hitpos_curr.entry;

        // Give hit[bin] the pointer to the list of pointers to hitlist elements of same template (for realignment)
        hit[bin]->plist_phits = array_plist_phits[hitpos_curr.index];

        // record db file name from which next HMM is read
        hit[bin]->dbfile = new char[strlen(hits_to_realign[idb]->entry->name)
            + 1];
        strcpy(hit[bin]->dbfile, hits_to_realign[idb]->entry->name);


        getTemplateHMM(*hitpos_curr.entry, par.wg, hit[bin]->ftellpos,
            format[bin], t[bin]);

        if (v >= 2)
          fprintf(stderr, "Realigning with %s\n", t[bin]->name);

#pragma omp critical
        {
          N_aligned++;
        }

        RealignByWorker(par, hit[bin], q, t[bin], format[bin], pb, R, S73, S33);

        //TODO???
        break;
      }
    }
  }

  if (v1 >= 2)
    cout << "\n";
  v = v1;

  // Delete all hitlist entries with too short alignments
  nhits = 0;
  hitlist.Reset();
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();

    if (nhits > par.realign_max && nhits >= std::max(par.B, par.Z))
      break;
    if (hit_cur.Eval > par.e) {
      if (nhits >= std::max(par.B, par.Z))
        continue;
      if (nhits >= std::max(par.b, par.z) && hit_cur.Probab < par.p)
        continue;
      if (nhits >= std::max(par.b, par.z) && hit_cur.Eval > par.E)
        continue;
    }

    if (hit_cur.matched_cols < MINCOLS_REALIGN) {
      if (v >= 3)
        printf("Deleting alignment of %s with length %i\n", hit_cur.name,
            hit_cur.matched_cols);
      // delete the list record and hit object
      hitlist.Delete().Delete();
    }
    nhits++;
  }

  // Delete hash phash_plist_realignhitpos with lists
  phash_plist_realignhitpos->Reset();
  while (!phash_plist_realignhitpos->End())
    delete (phash_plist_realignhitpos->ReadNext()); // delete list to which phash_plist_realignhitpos->ReadNext() points
  delete (phash_plist_realignhitpos);

  // Delete array_plist_phits with lists
  for (int index = 0; index < N_searched; index++)
    if (array_plist_phits[index])
      delete (array_plist_phits[index]); // delete list to which array[index] points
  delete[] (array_plist_phits);
}

void HHblits::reduceRedundancyOfHitList(int n_redundancy, int query_length,
    HitList& hitlist, HitList& reducedHitList) {
  int* coverage = new int[query_length + 1];
  for (int i = 0; i <= query_length; i++) {
    coverage[i] = 0;
  }

  int total = 0;

  hitlist.Reset();
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();

    if (!hit_cur.P_posterior) {
      continue;
    }
    total++;

    //calculate actual coverage
    int length_actual_coverage = 0;
    for (int alignment_index = 1; alignment_index <= hit_cur.nsteps;
        alignment_index++) {
      length_actual_coverage +=
          (hit_cur.P_posterior[alignment_index] > 0.5) ? 1 : 0;
    }

    int length_actual_contribution = 0;
    for (int alignment_index = 1; alignment_index <= hit_cur.nsteps;
        alignment_index++) {
      int i = hit_cur.i[alignment_index];
      if (hit_cur.P_posterior[alignment_index] > 0.5
          && coverage[i] < n_redundancy) {
        length_actual_contribution++;
      }
    }

    if (length_actual_contribution > 0.5 * length_actual_coverage) {
      reducedHitList.Insert(hit_cur);

      //update coverage
      for (int alignment_index = 1; alignment_index <= hit_cur.nsteps;
          alignment_index++) {
        int i = hit_cur.i[alignment_index];

        coverage[i] += (hit_cur.P_posterior[alignment_index] > 0.5) ? 1 : 0;
      }
    }
  }

  delete[] coverage;
}

void HHblits::recalculateAlignmentsForDifferentQSC(HitList& hitlist,
    Alignment& Qali, char inputformat, float* qsc, size_t nqsc,
    HitList& recalculated_hitlist) {
  char v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2; // Supress verbose output during iterative realignment and realignment

  const int COV_ABS = 25;
  int cov_tot = std::max(std::min((int) (COV_ABS / Qali.L * 100 + 0.5), 70),
      par.coverage);

  Alignment qali;
  qali = Qali;
  HMM* q = new HMM();

  HitList realigned_viterbi_hitlist;

  for (size_t qsc_index = 0; qsc_index < nqsc; qsc_index++) {
    float actual_qsc = qsc[qsc_index];

    qali.Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M,
        par.Mgaps);
    qali.N_filtered = qali.Filter(par.max_seqid, S, cov_tot, par.qid,
        actual_qsc, par.Ndiff);
    qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons,
        par.maxres, pb, Sim, NULL, false);
    PrepareQueryHMM(par, inputformat, q, pc_hhm_context_engine,
        pc_hhm_context_mode, pb, R);

    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_ref = hitlist.ReadNext();

      HMM* t = new HMM();

      int format;
      long ftellpos;
      getTemplateHMM(*hit_ref.entry, 1, ftellpos, format, t);

      PrepareTemplateHMM(par, q, t, format, pb, R);

      Hit hit;
      hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
      hit.self = 0;
      hit.realign_around_viterbi = false;

      hit.dbfile = new char[strlen(hit_ref.dbfile) + 1];
      strcpy(hit.dbfile, hit_ref.dbfile);

      for (int irep = 1; irep <= par.altali; irep++) {
        hit.irep = irep;
        hit.Viterbi(q, t, par.loc, par.ssm, par.maxres, par.min_overlap,
            par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);

        if (hit.irep > 1 && hit.score <= SMIN)
          break;

        hit.Backtrace(q, t, par.corr, par.ssw, S73, S33);
        realigned_viterbi_hitlist.Push(hit);
      }

      hit.DeleteBacktraceMatrix(q->L + 2);

      delete t;
    }

    realigned_viterbi_hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);
    realigned_viterbi_hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa,
        par.alphab, par.alphac, par.prefilter_evalue_thresh);

    q->Log2LinTransitionProbs(1.0);

    realigned_viterbi_hitlist.Reset();
    while (!realigned_viterbi_hitlist.End()) {
      Hit hit_ref = realigned_viterbi_hitlist.ReadNext();

      HMM* t = new HMM();

      int format;
      long ftellpos;
      getTemplateHMM(*hit_ref.entry, par.wg, ftellpos, format, t);

      PrepareTemplateHMM(par, q, t, format, pb, R);
      t->Log2LinTransitionProbs(1.0);

      Hit hit;
      hit.AllocateForwardMatrix(q->L + 2, par.maxres + 1);
      hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
      hit.irep = 1;
      hit.self = 0;

      hit.dbfile = new char[strlen(hit_ref.dbfile) + 1];
      strcpy(hit.dbfile, hit_ref.dbfile); // record db file name from which next HMM is read

      hit.i1 = hit_ref.i1;
      hit.i2 = hit_ref.i2;
      hit.j1 = hit_ref.j1;
      hit.j2 = hit_ref.j2;
      hit.nsteps = hit_ref.nsteps;
      hit.i = hit_ref.i;
      hit.j = hit_ref.j;
      hit.realign_around_viterbi = false;

      // Align q to template in *hit[bin]
      hit.Forward(q, t, par.ssm, par.min_overlap, par.loc, par.shift, par.ssw,
          par.exclstr, S73, S33);
      hit.Backward(q, t, par.loc, par.shift, par.ssw, S73, S33);
      hit.MACAlignment(q, t, par.loc, par.mact, par.macins);
      hit.BacktraceMAC(q, t, par.corr, par.ssw, S73, S33);

      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
      hit.score = hit_ref.score;
      hit.score_ss = hit_ref.score_ss;
      hit.score_aass = hit_ref.score_aass;
      hit.score_sort = hit_ref.score_sort;
      hit.Pval = hit_ref.Pval;
      hit.Pvalt = hit_ref.Pvalt;
      hit.logPval = hit_ref.logPval;
      hit.logPvalt = hit_ref.logPvalt;
      hit.Eval = hit_ref.Eval;
      hit.logEval = hit_ref.logEval;
      hit.Probab = hit_ref.Probab;

      hit.DeleteForwardMatrix(q->L + 2);
      hit.DeleteBacktraceMatrix(q->L + 2);

      if (hit.matched_cols >= MINCOLS_REALIGN) {
        recalculated_hitlist.Insert(hit);
      }

      delete t;
    }

    realigned_viterbi_hitlist.Reset();
    while (!realigned_viterbi_hitlist.End()) {
      realigned_viterbi_hitlist.Delete().Delete();
    }
  }

  delete q;

  v = v1;
}

void HHblits::uniqueHitlist(HitList& hitlist) {
  std::set<std::string> ids;

  hitlist.Reset();
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();

    std::string id = std::string(hit_cur.name);

    if (ids.find(id) != ids.end()) {
      hitlist.Delete();
    }
    else {
      ids.insert(id);
    }
  }
}

void HHblits::wiggleQSC(HitList& hitlist, int n_redundancy, Alignment& Qali,
    char inputformat, float* qsc, size_t nqsc, HitList& reducedFinalHitList) {
  int query_length = Qali.L;
  HitList* reducedHitList = new HitList();
  //  HitList* wiggledHitList = new HitList();

  //filter by 2*n_redundancy
  reduceRedundancyOfHitList(2 * n_redundancy, query_length, hitlist,
      *reducedHitList);

  //recalculate alignments for different qsc and choose best alignment for each template
  recalculateAlignmentsForDifferentQSC(*reducedHitList, Qali, inputformat, qsc,
      nqsc, reducedFinalHitList);  //*wiggledHitList);
  //  uniqueHitlist(*wiggledHitList);

  //filter by n_redundancy
  //  reduceRedundancyOfHitList(n_redundancy, query_length, *wiggledHitList, reducedFinalHitList);
}

void HHblits::wiggleQSC(int n_redundancy, float* qsc, size_t nqsc,
    HitList& reducedFinalHitList) {
  wiggleQSC(hitlist, n_redundancy, Qali, input_format, qsc, nqsc,
      reducedFinalHitList);
  reducedFinalHitList.N_searched = hitlist.N_searched;
}

void HHblits::run(FILE* query_fh, char* query_path) {
  int cluster_found = 0;
  int seqs_found = 0;
  size_t premerge = par.premerge;

  v1 = v;
  if (v > 0 && v <= 2)
    v = 1;
  else
    v--;

  Hit hit_cur;
  Hash<Hit>* previous_hits = new Hash<Hit>(1631, hit_cur);
  Hash<char>* premerged_hits = new Hash<char>(1631);

  q = new HMM;
  q_tmp = new HMM;

  // Read query input file (HHM or alignment format) without adding pseudocounts
  Qali.N_in = 0;
  ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb, S, Sim);

  if (Qali.N_in - Qali.N_ss > 1)
    premerge = 0;

  if (par.allseqs) {
    Qali_allseqs = Qali; // make a *deep* copy of Qali!
    for (int k = 0; k < Qali_allseqs.N_in; ++k)
      Qali_allseqs.keep[k] = 1; // keep *all* sequences (reset filtering in Qali)
  }

  v = v1;

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags)
    q->NeutralizeTags(pb);

  //save all entries pointer in this vector to delete, when it's safe
  std::vector<HHDatabaseEntry*> all_entries;
  std::vector<HHDatabaseEntry*> new_entries;
  std::vector<HHDatabaseEntry*> old_entries;

  if (!par.prefilter) {
    for (size_t i = 0; i < dbs.size(); i++) {
      dbs[0]->initNoPrefilter(new_entries);
    }
    all_entries.insert(all_entries.end(), new_entries.begin(),
        new_entries.end());
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  for (int round = 1; round <= par.num_rounds; round++) {
    if (v >= 2)
      printf("\nIteration %i\n", round);

    // Settings for different rounds
    if (premerge > 0 && round > 1
        && previous_hits->Size() >= premerge) {
      if (v > 3)
        printf(
            "Set premerge to 0! (premerge: %zu   iteration: %i   hits.Size: %zu)\n",
            premerge, round, previous_hits->Size());
      premerge = 0;
    }
    else {
      premerge -= previous_hits->Size();
    }

    // Save HMM without pseudocounts for prefilter query-profile
    *q_tmp = *q;

    //TODO: Write query HHM file? (not the final HMM, which will be written to par.hhmfile)
//    if (*query_hhmfile) {
//      v1 = v;
//      if (v > 0 && v <= 3)
//        v = 1;
//      else
//        v -= 2;
//
//      // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
//      q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
//      q->CalculateAminoAcidBackground();
//
//      q->WriteToFile(query_hhmfile);
//
//      v = v1;
//    }

    PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine,
        pc_hhm_context_mode, pb, R);

    ////////////////////////////////////////////
    // Prefiltering
    ////////////////////////////////////////////

    if (par.prefilter) {
      if (v >= 2)
        printf("Prefiltering database\n");

      new_entries.clear();
      old_entries.clear();

      // Add Pseudocounts to q_tmp
      if (par.nocontxt) {
        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        q_tmp->PreparePseudocounts(R);
        // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        q_tmp->AddAminoAcidPseudocounts(par.pc_prefilter_nocontext_mode,
            par.pc_prefilter_nocontext_a, par.pc_prefilter_nocontext_b,
            par.pc_prefilter_nocontext_c);
      }
      else {
        // Add context specific pseudocounts (now always used, because clusterfile is necessary)
        q_tmp->AddContextSpecificPseudocounts(pc_prefilter_context_engine,
            pc_prefilter_context_mode);
      }

      q_tmp->CalculateAminoAcidBackground(pb);

      for (size_t i = 0; i < dbs.size(); i++) {
        dbs[i]->prefilter_db(q_tmp, previous_hits, par.threads,
            par.prefilter_gap_open, par.prefilter_gap_extend,
            par.prefilter_score_offset, par.prefilter_bit_factor,
            par.prefilter_evalue_thresh, par.prefilter_evalue_coarse_thresh,
            par.preprefilter_smax_thresh, par.min_prefilter_hits, R,
            new_entries, old_entries);
      }
      all_entries.insert(all_entries.end(), new_entries.begin(),
          new_entries.end());
      all_entries.insert(all_entries.end(), old_entries.begin(),
          old_entries.end());
    }

    if (v >= 2 && new_entries.size() == 0) {
      printf("No HMMs pass prefilter => Stop searching!\n");
      break;
    }

    // Search datbases
    if (v >= 2) {
      printf(
          "HMMs passed 2nd prefilter (gapped profile-profile alignment)   : %6zu\n",
          (new_entries.size() + old_entries.size()));
      printf(
          "HMMs passed 2nd prefilter and not found in previous iterations : %6zu\n",
          new_entries.size());
      printf("Scoring %zu HMMs using HMM-HMM Viterbi alignment\n",
          new_entries.size());
    }

    // Main Viterbi HMM-HMM search
    // Starts with empty hitlist (hits of previous iterations were deleted) and creates a hitlist with the hits of this iteration
    ViterbiSearch(new_entries, previous_hits,
        (new_entries.size() + old_entries.size()));

    // check for new hits or end with iteration
    int new_hits = 0;
    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();
      if (hit_cur.Eval > 100.0 * par.e)
        break; // E-value much too large
      if (hit_cur.Eval > par.e)
        continue; // E-value too large
      new_hits++;
    }

    if (new_hits == 0 || round == par.num_rounds) {
      if (round < par.num_rounds && v >= 2)
        printf("No new hits found in iteration %i => Stop searching\n", round);

      if (old_entries.size() > 0 && par.realign_old_hits) {
        if (v > 0) {
          printf("Rescoring previously found HMMs with Viterbi algorithm\n");
        }
        ViterbiSearch(old_entries, previous_hits,
            (new_entries.size() + old_entries.size()));

        // Add dbfiles_old to dbfiles_new for realign
        for (size_t a = 0; a < old_entries.size(); a++) {
          new_entries.push_back(old_entries[a]);
        }
      }
      else if (!par.realign_old_hits && previous_hits->Size() > 0) {
        if (v > 0) {
          printf("Rescoring previously found HMMs with Viterbi algorithm\n");
        }

        RescoreWithViterbiKeepAlignment(
            new_entries.size() + previous_hits->Size(), previous_hits);
      }
    }

    // Realign hits with MAC algorithm
    if (par.realign)
      perform_realign(new_entries, premerge, premerged_hits);

    // Generate alignment for next iteration
    if (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile || *par.alisbasename) {
      v1 = v;
      if (v > 0 && v <= 3)
        v = 1;
      else
        v -= 2;

      // If new hits found, merge hits to query alignment
      if (new_hits > 0) {
        if (v >= 1)
          printf("Merging hits to query profile\n");

        // For each template below threshold
        hitlist.Reset();
        while (!hitlist.End()) {
          Hit hit_cur = hitlist.ReadNext();

          if (hit_cur.Eval > 100.0 * par.e)
            break; // E-value much too large
          if (hit_cur.Eval > par.e)
            continue; // E-value too large
          if (hit_cur.matched_cols < MINCOLS_REALIGN)
            continue; // leave out too short alignments

          // Already in alignment
          stringstream ss_tmp;
          ss_tmp << hit_cur.file << "__" << hit_cur.irep;
          if (previous_hits->Contains((char*) ss_tmp.str().c_str()))
            continue;

          // Add number of sequences in this cluster to total found
          seqs_found += SequencesInCluster(hit_cur.name); // read number after second '|'
          cluster_found++;

          // Skip merging this hit if hit alignment was already merged during premerging
          if (premerged_hits->Contains((char*) ss_tmp.str().c_str()))
            continue;

          // Read a3m alignment of hit from <file>.a3m file
          // Reading in next db MSA and merging it onto Qali
          HHblitsDatabase* db = getHHblitsDatabase(*hit_cur.entry, dbs);
          if (db == NULL) {
            std::cerr << "Could not find database for merging!" << std::endl;
            continue;
          }

          Alignment Tali;
          long ftellpos;
          getTemplateA3M(db, hit_cur.entry->entry->name, ftellpos, Tali);

          if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
            Qali_allseqs.MergeMasterSlave(hit_cur, Tali, hit_cur.dbfile,
                par.maxcol);
          Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
              par.qid_db, par.qsc_db, par.Ndiff_db);
          Qali.MergeMasterSlave(hit_cur, Tali, hit_cur.dbfile, par.maxcol);

          if (Qali.N_in >= MAXSEQ)
            break; // Maximum number of sequences reached
        }

        // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
        Qali.Compress("merged A3M file", par.cons, par.maxres, par.maxcol,
            par.M, par.Mgaps);

        // Sort out the nseqdis most dissimilacd r sequences for display in the result alignments
        Qali.FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid,
            par.qsc, par.nseqdis);

        // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
        const float COV_ABS = 25;     // min. number of aligned residues
        int cov_tot = std::max(std::min((int) (COV_ABS / Qali.L * 100 + 0.5), 70),
            par.coverage);
        if (v > 2)
          printf("Filter new alignment with cov %3i%%\n", cov_tot);
        Qali.N_filtered = Qali.Filter(par.max_seqid, S, cov_tot, par.qid,
            par.qsc, par.Ndiff);
      }

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons,
          par.showcons, par.maxres, pb, Sim, NULL, true);

      if (par.notags)
        q->NeutralizeTags(pb);

      // Calculate SSpred if we need to print out alis after each iteration or if last iteration
      // TODO: we should get rid of this... since it calls psipred on the command line and is untested
      if (par.addss && (*par.alisbasename || round == par.num_rounds || new_hits == 0)) {
        char ss_pred[par.maxres];
        char ss_conf[par.maxres];

        CalculateSS(q, ss_pred, ss_conf, par.psipred_data, par.psipred, pb);

        Qali.AddSSPrediction(ss_pred, ss_conf);
      }

      if (*par.alisbasename) {
        if (par.allseqs) {
          alis[round] = Qali_allseqs;
        }
        else {
          alis[round] = Qali;
        }
      }

      v = v1;
    }
    // Update counts for log
    else if (round == par.num_rounds) {
      hitlist.Reset();
      while (!hitlist.End()) {
        Hit hit_cur = hitlist.ReadNext();

        // E-value much too large
        if (hit_cur.Eval > 100.0 * par.e)
          break;

        // E-value too large
        if (hit_cur.Eval > par.e)
          continue;

        stringstream ss_tmp;
        ss_tmp << hit_cur.file << "__" << hit_cur.irep;
        // Already in alignment?
        if (previous_hits->Contains((char*) ss_tmp.str().c_str()))
          continue;

        // Add number of sequences in this cluster to total found
        // read number after second '|'
        seqs_found += SequencesInCluster(hit_cur.name);
        cluster_found++;
      }
    }

    if (v >= 2)
      printf(
          "%i sequences belonging to %i database HMMs found with an E-value < %-12.4g\n",
          seqs_found, cluster_found, par.e);

    if (v >= 2
        && (round < par.num_rounds || *par.alnfile || *par.psifile
            || *par.hhmfile || *par.alisbasename))
      printf(
          "Number of effective sequences of resulting query HMM: Neff = %4.2f\n",
          q->Neff_HMM);

    if (v >= 2 && q->Neff_HMM > par.neffmax && round < par.num_rounds) {
      printf(
          "Diversity is above threshold (%4.2f). Stop searching! (Change threshold using -neffmax <float>.)\n",
          par.neffmax);
    }

    if (v >= 2 && Qali.N_in >= MAXSEQ)
      printf(
          "Maximun number of sequences in query alignment reached (%i). Stop searching!\n",
          MAXSEQ);

    if (new_hits == 0 || round == par.num_rounds || q->Neff_HMM > par.neffmax
        || Qali.N_in >= MAXSEQ)
      break;

    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();
      char strtmp[NAMELEN + 6];
      sprintf(strtmp, "%s__%i%c", hit_cur.file, hit_cur.irep, '\0');
      if (!par.already_seen_filter || hit_cur.Eval > par.e
          || previous_hits->Contains(strtmp))
        hit_cur.Delete(); // Delete hit object (deep delete with Hit::Delete())
      else
        previous_hits->Add(strtmp, hit_cur);

      hitlist.Delete(); // Delete list record (flat delete)
    }
  }

  // Warn, if HMMER files were used
  if (par.hmmer_used && v >= 1)
    fprintf(stderr,
        "WARNING: Using HMMER files results in a drastically reduced sensitivity (>10%%).\nWe recommend to use HHMs build by hhmake.\n");

  if (*par.reduced_outfile) {
    float qscs[] = { -20, 0, 0.1, 0.2 };
    wiggleQSC(par.n_redundancy, qscs, 4, reducedHitlist);
  }

  for (size_t i = 0; i < all_entries.size(); i++) {
    delete all_entries[i];
  }
  all_entries.clear();

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete(); // Delete hit object
  delete previous_hits;

  delete premerged_hits;
}

void HHblits::writeHHRFile(char* hhrFile) {
  if (*hhrFile) {
    hitlist.PrintHHR(q_tmp, hhrFile, par.maxdbstrlen, par.showconf,
        par.showcons, par.showdssp, par.showpred, par.b, par.B, par.z, par.Z,
        par.aliwidth, par.nseqdis, par.p, par.E, par.argc, par.argv, S);
  }
}

void HHblits::writeAlisFile(char* basename) {
  if (*basename) {
    std::map<int, Alignment>::iterator it;
    for (it = alis.begin(); it != alis.end(); it++) {
      stringstream ss_tmp;
      ss_tmp << basename << "_" << (*it).first << ".a3m";
      std::string id = ss_tmp.str();

      (*it).second.WriteToFile(id.c_str(), par.append, "a3m");
    }
  }
}

void HHblits::writeScoresFile(char* scoresFile) {
  if (*scoresFile) {
    hitlist.PrintScoreFile(q, scoresFile);
  }
}

void HHblits::writePairwiseAlisFile(char* pairwiseAlisFile, char outformat) {
  if (*pairwiseAlisFile) {
    hitlist.PrintAlignments(q, pairwiseAlisFile, par.showconf, par.showcons,
        par.showdssp, par.showpred, par.p, par.aliwidth, par.nseqdis, par.b,
        par.B, par.E, S, outformat);
  }
}

void HHblits::writeAlitabFile(char* alitabFile) {
  if (*alitabFile) {
    hitlist.WriteToAlifile(q, alitabFile, par.alitab_scop);
  }
}

void HHblits::writeReducedHHRFile(char* reducedHHRFile) {
  if (*reducedHHRFile) {
    reducedHitlist.PrintHHR(q_tmp, reducedHHRFile, par.maxdbstrlen,
        par.showconf, par.showcons, par.showdssp, par.showpred, par.b, par.B,
        par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E, par.argc,
        par.argv, S);
  }
}

void HHblits::writePsiFile(char* psiFile) {
  // Write output PSI-BLAST-formatted alignment?
  if (*psiFile) {
    if (par.allseqs)
      Qali_allseqs.WriteToFile(psiFile, par.append, "psi");
    else
      Qali.WriteToFile(psiFile, par.append, "psi");
  }
}

void HHblits::writeHMMFile(char* HMMFile) {
  // Write output HHM file?
  if (*HMMFile) {
    // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
    q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
    q->CalculateAminoAcidBackground(pb);

    q->WriteToFile(HMMFile, par.append, par.max_seqid, par.coverage, par.qid,
        par.Ndiff, par.qsc, par.argc, par.argv, pb);
  }
}

void HHblits::writeA3MFile(char* A3MFile) {
  // Write output A3M alignment?
  if (*A3MFile) {
    if (par.allseqs)
      Qali_allseqs.WriteToFile(A3MFile, par.append, "a3m");
    else
      Qali.WriteToFile(A3MFile, par.append, "a3m");
  }
}

std::map<int, Alignment>& HHblits::getAlis() {
  return alis;
}

void HHblits::writeHHRFile(std::stringstream& out) {
  hitlist.PrintHHR(q_tmp, out, par.maxdbstrlen, par.showconf, par.showcons,
      par.showdssp, par.showpred, par.b, par.B, par.z, par.Z, par.aliwidth,
      par.nseqdis, par.p, par.E, par.argc, par.argv, S);
}

void HHblits::writeScoresFile(std::stringstream& out) {
  hitlist.PrintScoreFile(q, out);
}

void HHblits::writePairwiseAlisFile(char outformat, std::stringstream& out) {
  hitlist.PrintAlignments(q, out, par.showconf, par.showcons, par.showdssp,
      par.showpred, par.p, par.aliwidth, par.nseqdis, par.b, par.B, par.E, S,
      outformat);
}

void HHblits::writeAlitabFile(std::stringstream& out) {
  hitlist.WriteToAlifile(q, out, par.alitab_scop);
}

void HHblits::writeReducedHHRFile(std::stringstream& out) {
  reducedHitlist.PrintHHR(q_tmp, out, par.maxdbstrlen, par.showconf,
      par.showcons, par.showdssp, par.showpred, par.b, par.B, par.z, par.Z,
      par.aliwidth, par.nseqdis, par.p, par.E, par.argc, par.argv, S);
}

void HHblits::writePsiFile(std::stringstream& out) {
  if (par.allseqs)
    Qali_allseqs.WriteToFile(out, "psi");
  else
    Qali.WriteToFile(out, "psi");
}

void HHblits::writeHMMFile(std::stringstream& out) {
  // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
  q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
  q->CalculateAminoAcidBackground(pb);

  q->WriteToFile(out, par.max_seqid, par.coverage, par.qid, par.Ndiff, par.qsc,
      par.argc, par.argv, pb);
}

void HHblits::writeA3MFile(std::stringstream& out) {
  if (par.allseqs)
    Qali_allseqs.WriteToFile(out, "a3m");
  else
    Qali.WriteToFile(out, "a3m");
}
