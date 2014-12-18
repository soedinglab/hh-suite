/*
 * HHblits.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include "hhblits.h"

//TODO: get rid of exit(1)... throw errors
//TODO: get more flexible database reader

HHblits::HHblits(Parameters& parameters,
                 std::vector<HHblitsDatabase*>& databases) {
  par = parameters;
  dbs = databases;

  Qali = NULL;
  Qali_allseqs = NULL;

  q = NULL;
  q_tmp = NULL;

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
    } else {
      context_lib = new cs::ContextLibrary<cs::AA>(fin);
      cs::TransformToLog(*context_lib);
      pc_hhm_context_engine = new cs::LibraryPseudocounts<cs::AA>(*context_lib,
                                                                  par.csw,
                                                                  par.csb);
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
  for (int bin = 0; bin < par.threads; bin++) {
    viterbiMatrices[bin] = new ViterbiMatrix();
    posteriorMatrices[bin] = new PosteriorMatrix();
  }
}

HHblits::~HHblits() {
  Reset();

  for (int bin = 0; bin < par.threads; bin++) {
    delete viterbiMatrices[bin];
    delete posteriorMatrices[bin];
  }

  DeletePseudocountsEngine(context_lib, crf, pc_hhm_context_engine,
                           pc_hhm_context_mode, pc_prefilter_context_engine,
                           pc_prefilter_context_mode);
}

void HHblits::prepareDatabases(Parameters& par,
                               std::vector<HHblitsDatabase*>& databases) {
  for (size_t i = 0; i < par.db_bases.size(); i++) {
    HHblitsDatabase* db = new HHblitsDatabase(par.db_bases[i].c_str());
    databases.push_back(db);
  }

  par.dbsize = 0;
  for (size_t i = 0; i < databases.size(); i++) {
    par.dbsize += databases[i]->cs219_database->db_index->n_entries;
  }

  if (par.prefilter) {
    for (size_t i = 0; i < databases.size(); i++) {
      databases[i]->initPrefilter(par.cs_library);
    }
  }
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
  int v = 2;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      v = atoi(argv[i + 1]);
      break;
    }
  }
  par.v = Log::from_int(v);
  Log::reporting_level() = par.v;

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
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
              << ":" << std::endl;
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
  if (par.loc == 0 && par.num_rounds >= 2) {
    HH_LOG(WARNING)
        << "Warning in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << "\n"
        << "\tusing -global alignment for iterative searches is deprecated "
        "since non-homologous sequence segments can easily enter the "
        "MSA and corrupt it."
        << std::endl;
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
    par.nseqdis = MAXSEQDIS - 3 - par.showcons;  //3 reserved for secondary structure
  if (par.aliwidth < 20)
    par.aliwidth = 20;
  if (par.pc_hhm_context_engine.pca < 0.001)
    par.pc_hhm_context_engine.pca = 0.001;  // to avoid log(0)
  if (par.pc_prefilter_context_engine.pca < 0.001)
    par.pc_prefilter_context_engine.pca = 0.001;  // to avoid log(0)
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

  if (Qali) {
    delete Qali;
    Qali = NULL;
  }

  if (Qali_allseqs) {
    delete Qali_allseqs;
    Qali_allseqs = NULL;
  }

  hitlist.Reset();
  while (!hitlist.End())
    hitlist.Delete().Delete();

  optimized_hitlist.Reset();
  while (!optimized_hitlist.End())
    optimized_hitlist.Delete().Delete();

  std::map<int, Alignment*>::iterator it;
  for (it = alis.begin(); it != alis.end(); it++) {
    delete (*it).second;
  }
  alis.clear();

  for (int bin = 0; bin < par.threads; bin++) {
    viterbiMatrices[bin]->DeleteBacktraceMatrix();
    posteriorMatrices[bin]->DeleteProbabilityMatrix();
  }
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
        " -Z <int>       maximum number of lines in summary hit list (default=%i)        \n",
        par.Z);
    printf(
        " -z <int>       minimum number of lines in summary hit list (default=%i)        \n",
        par.z);
    printf(
        " -B <int>       maximum number of alignments in alignment list (default=%i)     \n",
        par.B);
    printf(
        " -b <int>       minimum number of alignments in alignment list (default=%i)     \n",
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
        " -realign_max <int>  realign max. <int> hits (default=%i)                        \n",
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
        " -sc   <int>    amino acid score         (tja: template HMM at column j) (def=%i)\n",
        par.columnscore);
    printf(
        "        0       = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
    printf(
        "        1       = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
    printf(
        "        2       = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
    printf(
        "        3       = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
    printf(
        "        5       local amino acid composition correction                     \n");
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
      par.v);
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
    HH_LOG(DEBUG1) << i << "  " << argv[i] << endl;
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no query file following -i\n";
        exit(4);
      } else
        strcpy(par.infile, argv[i]);
    } else if (!strcmp(argv[i], "-d")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no database basename following -d\n";
        exit(4);
      } else {
        std::string db(argv[i]);
        par.db_bases.push_back(db);
      }
    } else if (!strcmp(argv[i], "-contxt")
        || !strcmp(argv[i], "-context_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no lib following -contxt\n";
        exit(4);
      } else
        strcpy(par.clusterfile, argv[i]);
    } else if (!strcmp(argv[i], "-cslib") || !strcmp(argv[i], "-cs_lib")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no lib following -cslib\n";
        exit(4);
      } else
        strcpy(par.cs_library, argv[i]);
    } else if (!strcmp(argv[i], "-psipred")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no directory following -psipred\n";
        exit(4);
      } else
        strcpy(par.psipred, argv[i]);
    } else if (!strcmp(argv[i], "-psipred_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no database directory following -psipred_data\n";
        exit(4);
      } else
        strcpy(par.psipred_data, argv[i]);
    } else if (!strcmp(argv[i], "-o")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -o\n";
        exit(4);
      } else
        strcpy(par.outfile, argv[i]);
    } else if (!strcmp(argv[i], "-omat")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -omat\n";
        exit(4);
      } else
        strcpy(par.matrices_output_file, argv[i]);
    } else if (!strcmp(argv[i], "-oopt")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -o\n";
        exit(4);
      } else
        strcpy(par.opt_outfile, argv[i]);
    } else if (!strcmp(argv[i], "-oa3m")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -oa3m\n";
        exit(4);
      } else
        strcpy(par.alnfile, argv[i]);
    } else if (!strcmp(argv[i], "-opt")) {
      par.optimize_qsc = true;
    } else if (!strcmp(argv[i], "-ohhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -ohhm\n";
        exit(4);
      } else
        strcpy(par.hhmfile, argv[i]);
    } else if (!strcmp(argv[i], "-opsi")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -opsi\n";
        exit(4);
      } else
        strcpy(par.psifile, argv[i]);
    } else if (!strcmp(argv[i], "-oalis")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no file basename following -oalis\n";
        exit(4);
      } else
        strcpy(par.alisbasename, argv[i]);
    } else if (!strcmp(argv[i], "-Ofas")) {
      par.append = 0;
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -o\n";
        exit(4);
      } else
        strcpy(par.pairwisealisfile, argv[i]);
    } else if (!strcmp(argv[i], "-Oa2m")) {
      par.append = 0;
      par.outformat = 2;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -o\n";
        exit(4);
      } else
        strcpy(par.pairwisealisfile, argv[i]);
    } else if (!strcmp(argv[i], "-Oa3m")) {
      par.append = 0;
      par.outformat = 3;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no output file following -o\n";
        exit(4);
      } else
        strcpy(par.pairwisealisfile, argv[i]);
    } else if (!strcmp(argv[i], "-qhhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no filename following -qhhm\n";
        exit(4);
      } else
        strcpy(par.query_hhmfile, argv[i]);
    } else if (!strcmp(argv[i], "-scores")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no file following -scores\n";
        exit(4);
      } else {
        strcpy(par.scorefile, argv[i]);
      }
    } else if (!strcmp(argv[i], "-atab")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        cerr << endl << "Error in " << program_name
             << ": no file following -atab\n";
        exit(4);
      } else {
        strcpy(par.alitabfile, argv[i]);
      }
    } else if (!strcmp(argv[i], "-atab_scop"))
      par.alitab_scop = true;
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
      help(par, 1);
      exit(0);
    } else if (!strcmp(argv[i], "-v") && (i < argc - 1)
        && argv[i + 1][0] != '-') {
      int v = atoi(argv[++i]);
      par.v = Log::from_int(v);
      Log::reporting_level() = par.v;
    } else if (!strcmp(argv[i], "-n") && (i < argc - 1))
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
    } else if (!strcmp(argv[i], "-M") && (i < argc - 1))
      if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m"))
        par.M = 1;
      else if (!strcmp(argv[i], "first"))
        par.M = 3;
      else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
        par.Mgaps = atoi(argv[i]);
        par.M = 2;
      } else
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
    } else if (!strcmp(argv[i], "-neffmax") && (i < argc - 1))
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
    } else if (!strcmp(argv[i], "-gapd") && (i < argc - 1))
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
    } else if (!strcmp(argv[i], "-noaddfilter")) {
      par.already_seen_filter = false;
    } else if (!strcmp(argv[i], "-maxfilt") && (i < argc - 1))
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
    } else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = 2 * par.maxres;
    } else if (!strncmp(argv[i], "-glo", 3)) {
      par.loc = 0;
      if (par.mact > 0.35 && par.mact < 0.3502) {
        par.mact = 0;
      }
    } else if (!strncmp(argv[i], "-loc", 4))
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
    else if (!strcmp(argv[i], "-sc") && (i < argc - 1))
      par.columnscore = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-scwin") && (i < argc - 1)) {
      par.columnscore = 5;
      par.half_window_size_local_aa_bg_freqs = std::max(1, atoi(argv[++i]));
    } else if (!strncmp(argv[i], "-cpu", 4) && (i < argc - 1)) {
      par.threads = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-maxmem") && (i < argc - 1)) {
      par.maxmem = atof(argv[++i]);
    } else if (!strncmp(argv[i], "-premerge", 9) && (i < argc - 1))
      par.premerge = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-nocontxt"))
      par.nocontxt = 1;
    else if (!strcmp(argv[i], "-csb") && (i < argc - 1))
      par.csb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-csw") && (i < argc - 1))
      par.csw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-corr") && (i < argc - 1))
      par.corr = atof(argv[++i]);
    else {
      HH_LOG(WARNING) << endl << "WARNING: Ignoring unknown option "
                                << argv[i] << " ...\n";
    }

    HH_LOG(DEBUG1) << i << "  " << argv[i] << endl;
  }  // end of for-loop for command line input
}

void HHblits::mergeHitsToQuery(Hash<Hit>* previous_hits,
                               Hash<char>* premerged_hits, int& seqs_found,
                               int& cluster_found) {
  // For each template below threshold
  hitlist.Reset();
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();

    if (hit_cur.Eval > 100.0 * par.e)
      break;  // E-value much too large
    if (hit_cur.Eval > par.e)
      continue;  // E-value too large
    if (hit_cur.matched_cols < MINCOLS_REALIGN)
      continue;  // leave out too short alignments

    // Already in alignment
    stringstream ss_tmp;
    ss_tmp << hit_cur.file << "__" << hit_cur.irep;
    if (previous_hits->Contains((char*) ss_tmp.str().c_str()))
      continue;

    // Add number of sequences in this cluster to total found
    seqs_found += SequencesInCluster(hit_cur.name);  // read number after second '|'
    cluster_found++;

    // Skip merging this hit if hit alignment was already merged during premerging
    if (premerged_hits->Contains((char*) ss_tmp.str().c_str()))
      continue;

    // Read a3m alignment of hit from <file>.a3m file
    // Reading in next db MSA and merging it onto Qali

    Alignment Tali;
    hit_cur.entry->getTemplateA3M(par, pb, S, Sim, Tali);

    if (par.allseqs)  // need to keep *all* sequences in Qali_allseqs? => merge before filtering
      Qali_allseqs->MergeMasterSlave(hit_cur, Tali, hit_cur.name, par.maxcol);
    Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
                                  par.qid_db, par.qsc_db, par.Ndiff_db);
    Qali->MergeMasterSlave(hit_cur, Tali, hit_cur.name, par.maxcol);

    if (Qali->N_in >= MAXSEQ)
      break;  // Maximum number of sequences reached
  }

  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
  Qali->Compress("merged A3M file", par.cons, par.maxres, par.maxcol, par.M,
                 par.Mgaps);

  // Sort out the nseqdis most dissimilacd r sequences for display in the result alignments
  Qali->FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid,
                         par.qsc, par.nseqdis);

  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
  const float COV_ABS = 25;     // min. number of aligned residues
  int cov_tot = std::max(std::min((int) (COV_ABS / Qali->L * 100 + 0.5), 70),
                         par.coverage);

  HH_LOG(DEBUG) << "Filter new alignment with cov " << cov_tot
                          << std::endl;

  Qali->N_filtered = Qali->Filter(par.max_seqid, S, cov_tot, par.qid, par.qsc,
                                  par.Ndiff);
}

void HHblits::add_hits_to_hitlist(std::vector<Hit>& hits, HitList& hitlist) {
  for (std::vector<Hit>::size_type i = 0; i != hits.size(); i++) {
    hitlist.Push(hits[i]);
  }

  // Sort list according to sortscore
  HH_LOG(DEBUG) << "Sorting hit list ...\n";
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
void HHblits::RescoreWithViterbiKeepAlignment(HMMSimd& q_vec,
                                              Hash<Hit>* previous_hits) {
  // Initialize
  std::vector<HHEntry*> hits_to_rescore;

  // Get dbfiles of previous hits
  previous_hits->Reset();
  while (!previous_hits->End()) {
    Hit hit_cur = previous_hits->ReadNext();
    if (hit_cur.irep == 1) {
      hits_to_rescore.push_back(hit_cur.entry);
    }
  }

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
//	DoViterbiSearch(hits_to_rescore, previous_hits, false);

  ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
  std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec,
                                                         hits_to_rescore,
                                                         par.qsc_db, pb, S, Sim,
                                                         R);

  for (std::vector<Hit>::size_type i = 0; i != hits_to_add.size(); i++) {
    stringstream ss_tmp;
    ss_tmp << hits_to_add[i].name << "__" << hits_to_add[i].irep;
    if (previous_hits->Contains((char*) ss_tmp.str().c_str())) {
      Hit hit_cur = previous_hits->Remove((char*) ss_tmp.str().c_str());
      previous_hits->Add((char*) ss_tmp.str().c_str(), hits_to_add[i]);
      // Overwrite *hit[bin] with alignment, etc. of hit_cur
      hit_cur.score = hits_to_add[i].score;
      hit_cur.score_aass = hits_to_add[i].score_aass;
      hit_cur.score_ss = hits_to_add[i].score_ss;
      hit_cur.Pval = hits_to_add[i].Pval;
      hit_cur.Pvalt = hits_to_add[i].Pvalt;
      hit_cur.logPval = hits_to_add[i].logPval;
      hit_cur.logPvalt = hits_to_add[i].logPvalt;
      hit_cur.Eval = hits_to_add[i].Eval;
      hit_cur.logEval = hits_to_add[i].logEval;
      hit_cur.Probab = hits_to_add[i].Probab;
      hitlist.Push(hit_cur);  // insert hit at beginning of list (last repeats first!)
    } else {
      hits_to_add[i].Delete();
    }
  }

  // Sort list according to sortscore
  HH_LOG(DEBUG) << "Sorting hit list ..." << std::endl;

  hitlist.SortList();

  hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);  // Use NN prediction of lamda and mu

  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa, par.alphab,
                                    par.alphac, par.prefilter_evalue_thresh);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Realign hits with MAC algorithm
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HHblits::perform_realign(HMMSimd& q_vec, const char input_format,
                              std::vector<HHEntry*>& hits_to_realign,
                              const int premerge, Hash<char>* premerged_hits) {

  if (par.optimize_qsc) {
    realignWithOptimalQSC(hitlist, q_vec, input_format);
  } else {
    // 19/02/2014: F/B-algos are calculated in log-space
    //  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
    int nhits = 0;

    // Longest allowable length of database HMM (backtrace: 5 chars, fwd: 1 double, bwd: 1 double
    long int Lmaxmem = ((par.maxmem - 0.5) * 1024 * 1024 * 1024)
        / (2 * sizeof(double) + 8) / q->L / par.threads;
    int Lmax = 0;      // length of longest HMM to be realigned

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Categorize the hits: first irep, then dbfile
    int n_realignments = 0;
    std::map<short int, std::vector<Hit *> > alignments;
    hitlist.Reset();

    std::vector<Hit *> hit_vector;

    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();
      if (n_realignments >= par.realign_max
          && n_realignments >= imax(par.B, par.Z))
        break;
      if (hit_cur.Eval > par.e) {
        if (n_realignments >= imax(par.B, par.Z))
          continue;
        if (n_realignments >= imax(par.b, par.z) && hit_cur.Probab < par.p)
          continue;
        if (n_realignments >= imax(par.b, par.z) && hit_cur.Eval > par.E)
          continue;
      }

      if (hit_cur.L > Lmax) {
        Lmax = hit_cur.L;
      }
      if (hit_cur.L > Lmaxmem) {
        nhits++;
        continue;
      }
      //TODO:
      if (true) {  //nhits >= par.premerge) {
        hit_vector.push_back(hitlist.ReadCurrentAddress());
        n_realignments++;
      }
      nhits++;
    }

    int t_maxres = Lmax + 2;
    for (int i = 0; i < par.threads; i++) {
      posteriorMatrices[i]->allocateMatrix(q->L, t_maxres);
    }

    // Sort hits in descending order
    std::qsort(&hit_vector[0], hit_vector.size(), sizeof(Hit*),
               compareHitLengths);

    for (int elem = 0; elem < (int) hit_vector.size(); elem++) {
      alignments[hit_vector.at(elem)->irep].push_back(hit_vector.at(elem));
    }

    PosteriorDecoderRunnerInputData input_data(dbs, hits_to_realign, alignments,
                                               n_realignments, t_maxres);

    // Initialize a Null-value as a return value if not items are available anymore
    PosteriorDecoderRunner runner(input_data, q_vec, posteriorMatrices,
                                  viterbiMatrices, par.threads);

    HH_LOG(INFO)
        << "Realigning " << nhits
        << " HMM-HMM alignments using Maximum Accuracy algorithm" << std::endl;

    runner.executeComputation(par, par.qsc_db, pb, S, Sim, R);

  }

//  // Delete all hitlist entries with too short alignments
//  nhits = 0;
//  hitlist.Reset();
//  while (!hitlist.End()) {
//    Hit hit_cur = hitlist.ReadNext();
//    //printf("Deleting alignment of %s with length %i? irep=%i nhits=%-2i  par.B=%-3i  par.Z=%-3i par.e=%.2g par.b=%-3i  par.z=%-3i par.p=%.2g\n",hit_cur.name,hit_cur.matched_cols,hit_cur.irep,nhits,par.B,par.Z,par.e,par.b,par.z,par.p);
//
//    if (nhits > par.realign_max && nhits >= imax(par.B, par.Z))
//      break;
//    if (hit_cur.Eval > par.e) {
//      if (nhits >= imax(par.B, par.Z))
//        continue;
//      if (nhits >= imax(par.b, par.z) && hit_cur.Probab < par.p)
//        continue;
//      if (nhits >= imax(par.b, par.z) && hit_cur.Eval > par.E)
//        continue;
//    }
//
//    if (hit_cur.matched_cols < MINCOLS_REALIGN) {
//      HH_LOG(LogLevel::DEBUG) << "Deleting alignment of " << hit_cur.name
//          << " with length " << hit_cur.matched_cols << std::endl;
//      hitlist.Delete().Delete();        // delete the list record and hit object
//      // // Make sure only realigned alignments get displayed! JS: Why? better unrealigned than none.
//      // if (last_round)
//      // if (par.B>par.Z) par.B--; else if (par.B==par.Z) {par.B--; par.Z--;} else par.Z--;
//    }
//    nhits++;
//  }
}

void HHblits::realignWithOptimalQSC(HitList& hitlist, HMMSimd& q_vec,
                                    char query_input_format) {

  const int COV_ABS = 25;
  const int cov_tot = std::max(std::min((int) (COV_ABS / q->L * 100 + 0.5), 70),
                               par.coverage);

  // Longest allowable length of database HMM (backtrace: 5 chars, fwd: 1 double, bwd: 1 double
  long int Lmaxmem = ((par.maxmem - 0.5) * 1024 * 1024 * 1024)
      / (2 * sizeof(double) + 8) / q->L / par.threads;

  const int nqsc = 6;
  float qscs[nqsc] = { -20, 0, 0.1, 0.2, 0.3, 0.4 };

  int n_realignments = 0;
  int Lmax = 0;
  int nhits = 0;

  //find interesting hits
  HitList interesting_list;
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();

    if (n_realignments >= par.realign_max
        && n_realignments >= imax(par.B, par.Z))
      break;

    if (hit_cur.Eval > par.e) {
      if (n_realignments >= imax(par.B, par.Z))
        continue;
      if (n_realignments >= imax(par.b, par.z) && hit_cur.Probab < par.p)
        continue;
      if (n_realignments >= imax(par.b, par.z) && hit_cur.Eval > par.E)
        continue;
    }

    if (hit_cur.L > Lmax) {
      Lmax = hit_cur.L;
    }

    if (hit_cur.L > Lmaxmem) {
      nhits++;
      continue;
    }

    //move relevant hits
    Hit hit_new;
    hit_new.initHitFromHit(hit_cur, q->L);
    interesting_list.Insert(hit_new);
    hitlist.Delete().Delete();

    n_realignments++;
    nhits++;
  }

  //init alignment matrices
  int t_maxres = Lmax + 2;
  for (int i = 0; i < par.threads; i++) {
    posteriorMatrices[i]->allocateMatrix(q->L, t_maxres);
  }

  HitList all_list;

  HitList curr_list;
  for (int i = 0; i < nqsc; i++) {
    float actual_qsc = qscs[i];

    Qali->Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M,
                   par.Mgaps);
    Qali->N_filtered = Qali->Filter(par.max_seqid, S, cov_tot, par.qid,
                                    actual_qsc, par.Ndiff);
    Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons,
                                    par.maxres, pb, Sim, NULL, false);
    PrepareQueryHMM(par, query_input_format, q, pc_hhm_context_engine,
                    pc_hhm_context_mode, pb, R);

    q_vec.MapOneHMM(q);

    std::vector<Hit *> hit_vector;
    std::vector<HHEntry*> hits_to_realign;

    //copy interesting hits to curr_list
    while (!interesting_list.End()) {
      Hit hit_cur = interesting_list.ReadNext();

      //copy relevant hits
      Hit hit_new;
      hit_new.initHitFromHit(hit_cur, q->L);

      curr_list.Insert(hit_new);
    }

    //prepare data-structures for realignment
    curr_list.Reset();
    while (!curr_list.End()) {
      Hit hit_cur = curr_list.ReadNext();
      hit_vector.push_back(curr_list.ReadCurrentAddress());
      hits_to_realign.push_back(hit_cur.entry);
    }

    std::qsort(&hit_vector[0], hit_vector.size(), sizeof(Hit*),
               compareHitLengths);

    std::map<short int, std::vector<Hit *> > alignments;
    for (int elem = 0; elem < (int) hit_vector.size(); elem++) {
      alignments[hit_vector.at(elem)->irep].push_back(hit_vector.at(elem));
    }

    PosteriorDecoderRunnerInputData input_data(dbs, hits_to_realign, alignments,
                                               n_realignments, t_maxres);

    // Initialize a Null-value as a return value if not items are available anymore
    PosteriorDecoderRunner runner(input_data, q_vec, posteriorMatrices,
                                  viterbiMatrices, par.threads);

    runner.executeComputation(par, actual_qsc, pb, S, Sim, R);

    //copy hits from curr_list to all_list
    curr_list.Reset();
    while (!curr_list.End()) {
      Hit hit_cur = curr_list.ReadNext();
      hit_cur.predicted_alignment_quality = hit_cur.estimateAlignmentQuality(q);
      hit_cur.qsc = actual_qsc;
      all_list.Push(hit_cur);

      curr_list.Delete();
    }
  }

  interesting_list.Reset();
  while (!interesting_list.End()) {
    interesting_list.Delete().Delete();
  }

  //select of each hit_irep the alignment with the best predicted alignment quality
  std::set<std::string> output_set;
  //sort all_list by predicted alignment quality
  all_list.SortList(&Hit::compare_predicted_alignment_quality);
  all_list.Reset();

  while (!all_list.End()) {
    Hit hit_cur = all_list.ReadNext();

    stringstream ss_tmp;
    ss_tmp << hit_cur.name << "__" << hit_cur.irep;

    if (output_set.find(ss_tmp.str()) != output_set.end()) {
      all_list.Delete().Delete();
    } else {
      output_set.insert(ss_tmp.str());

      hitlist.Push(hit_cur);
      all_list.Delete();
    }
  }

  hitlist.SortList();
}

void HHblits::optimizeQSC(std::vector<HHEntry*>& selected_entries,
                          const int N_searched, HMMSimd& q_vec,
                          char query_input_format, HitList& output_list) {
  const int COV_ABS = 25;
  const int cov_tot = std::max(std::min((int) (COV_ABS / q->L * 100 + 0.5), 70),
                               par.coverage);

  const int nqsc = 6;
  float qscs[nqsc] = { -20, 0, 0.1, 0.2, 0.3, 0.4 };

  HitList all_list;

  for (int i = 0; i < nqsc; i++) {
    float actual_qsc = qscs[i];

    HitList tmp_list;
    tmp_list.N_searched = N_searched;

    Qali->Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M,
                   par.Mgaps);
    Qali->N_filtered = Qali->Filter(par.max_seqid, S, cov_tot, par.qid,
                                    actual_qsc, par.Ndiff);
    Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons,
                                    par.maxres, pb, Sim, NULL, false);
    PrepareQueryHMM(par, query_input_format, q, pc_hhm_context_engine,
                    pc_hhm_context_mode, pb, R);

    q_vec.MapOneHMM(q);

    ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
    std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec,
                                                           selected_entries,
                                                           actual_qsc, pb, S,
                                                           Sim, R);

    add_hits_to_hitlist(hits_to_add, tmp_list);

    std::vector<Hit *> hit_vector;
    std::vector<HHEntry*> hits_to_realign;

    int t_maxres = 0;
    int n_realignments = 0;
    tmp_list.Reset();
    while (!tmp_list.End()) {
      Hit hit_cur = tmp_list.ReadNext();
      t_maxres = std::max(t_maxres, hit_cur.L + 2);
      hits_to_realign.push_back(hit_cur.entry);
      hit_vector.push_back(tmp_list.ReadCurrentAddress());
      n_realignments++;
    }

    for (int i = 0; i < par.threads; i++) {
      posteriorMatrices[i]->allocateMatrix(q->L, t_maxres);
    }

    std::qsort(&hit_vector[0], hit_vector.size(), sizeof(Hit*),
               compareHitLengths);

    std::map<short int, std::vector<Hit *> > alignments;
    for (int elem = 0; elem < (int) hit_vector.size(); elem++) {
      alignments[hit_vector.at(elem)->irep].push_back(hit_vector.at(elem));
    }

    PosteriorDecoderRunnerInputData input_data(dbs, hits_to_realign, alignments,
                                               n_realignments, t_maxres);

    // Initialize a Null-value as a return value if not items are available anymore
    PosteriorDecoderRunner runner(input_data, q_vec, posteriorMatrices,
                                  viterbiMatrices, par.threads);

    runner.executeComputation(par, actual_qsc, pb, S, Sim, R);

    tmp_list.Reset();
    while (!tmp_list.End()) {
      Hit hit_cur = tmp_list.ReadNext();
      hit_cur.predicted_alignment_quality = hit_cur.estimateAlignmentQuality(q);
      hit_cur.qsc = actual_qsc;
      all_list.Push(hit_cur);

      tmp_list.Delete();
    }
  }

//TODO: more sophisticated stuff
//	std::map<std::string, ViterbiScores> best_evalue_scores;
//
//	//select for each hit_irep the best evalue
//	all_list.SortList(&Hit::compare_evalue);
//  all_list.Reset();
//  while (!all_list.End()) {
//    Hit hit_cur = all_list.ReadNext();
//
//    stringstream ss_tmp;
//    ss_tmp << hit_cur.name << "__" << hit_cur.irep;
//
//    if(best_evalue_scores.find(ss_tmp.str()) == best_evalue_scores.end()) {
//      ViterbiScores scores(hit_cur);
//      best_evalue_scores.insert(std::pair<std::string, ViterbiScores>(ss_tmp.str(), scores));
//    }
//  }

  //select of each hit_irep the alignment with the best predicted alignment quality
  std::set<std::string> output_set;
  //sort all_list by predicted alignment quality
  all_list.SortList(&Hit::compare_predicted_alignment_quality);
  all_list.Reset();
  while (!all_list.End()) {
    Hit hit_cur = all_list.ReadNext();

    stringstream ss_tmp;
    ss_tmp << hit_cur.name << "__" << hit_cur.irep;
//		std::cout << ss_tmp.str() << "\t" << hit_cur.predicted_alignment_quality << "\t" << hit_cur.score_sort << std::endl;

    if (output_set.find(ss_tmp.str()) != output_set.end()) {
      all_list.Delete().Delete();
    } else {
      output_set.insert(ss_tmp.str());

//      ViterbiScores best_scores = best_evalue_scores[ss_tmp.str()];
//      hit_cur.score = best_scores.score;
//      hit_cur.score_aass = best_scores.score_aass;
//      hit_cur.score_ss = best_scores.score_ss;
//      hit_cur.Pval = best_scores.Pval;
//      hit_cur.Pvalt = best_scores.Pvalt;
//      hit_cur.logPval = best_scores.logPval;
//      hit_cur.logPvalt = best_scores.logPvalt;
//      hit_cur.Eval = best_scores.Eval;
//      hit_cur.logEval = best_scores.logEval;
//      hit_cur.Probab = best_scores.Probab;

      output_list.Push(hit_cur);
      all_list.Delete();
    }
  }

  output_list.N_searched = N_searched;
  output_list.SortList();
}

void HHblits::get_entries_of_selected_hits(
    HitList& input, std::vector<HHEntry*>& selected_entries) {
  std::set<std::string> output_set;

  int n_realignments = 0;

  input.Reset();
  while (!input.End()) {
    Hit hit_cur = input.ReadNext();
    if (n_realignments >= par.realign_max
        && n_realignments >= imax(par.B, par.Z))
      break;

    if (hit_cur.Eval > par.e) {
      if (n_realignments >= imax(par.B, par.Z))
        continue;
      if (n_realignments >= imax(par.b, par.z) && hit_cur.Probab < par.p)
        continue;
      if (n_realignments >= imax(par.b, par.z) && hit_cur.Eval > par.E)
        continue;
    }

    std::string ss_tmp = std::string(hit_cur.entry->getName());

    if (output_set.find(ss_tmp) == output_set.end()) {
      output_set.insert(ss_tmp);
      selected_entries.push_back(hit_cur.entry);
      n_realignments++;
    }
  }
}

void HHblits::get_entries_of_all_hits(HitList& input,
                                      std::vector<HHEntry*>& selected_entries) {
  std::set<std::string> output_set;

  input.Reset();
  while (!input.End()) {
    Hit hit_cur = input.ReadNext();

    std::string ss_tmp = std::string(hit_cur.entry->getName());

    if (output_set.find(ss_tmp) == output_set.end()) {
      output_set.insert(ss_tmp);
      selected_entries.push_back(hit_cur.entry);
    }
  }
}

//void HHblits::perform_realign(std::vector<HHDatabaseEntry*>& hits_to_realign,
//    const int premerge, Hash<char>* premerged_hits) {
//  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
//  int nhits = 0;
//  size_t N_aligned = 0;
//
//  // Longest allowable length of database HMM (backtrace: 5 chars, fwd, bwd: 1 double
//  long int Lmaxmem = (par.maxmem * 1024 * 1024 * 1024) / sizeof(double) / q->L
//      / par.threads;
//  long int Lmax = 0;      // length of longest HMM to be realigned
//
//  // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
//  // This list can be sorted by ftellpos to access one template after the other efficiently during realignment
//  Hash<List<Realign_hitpos>*>* phash_plist_realignhitpos;
//  phash_plist_realignhitpos = new Hash<List<Realign_hitpos>*>(100031, NULL);
//
//  // Some templates have several (suboptimal) alignments in hitlist. For realignment, we need to efficiently
//  // access all hit objects in hitlist belonging to one template (because we don't want to read templates twice)
//  // We therefore need for each template (identified by its index between 0 and N_searched-1) a list of elements
//  // in hitlist that store the alignments with the template of that index.
//  // This list is pointed to by array_plist_phits[index].
//  List<void*>** array_plist_phits;
//  array_plist_phits = new List<void*>*[N_searched];
//  for (int index = 0; index < N_searched; index++)
//    array_plist_phits[index] = NULL; // initialize
//
//  // Store all dbfiles and ftell positions of templates to be displayed and realigned
//  hitlist.Reset();
//  while (!hitlist.End()) {
//    Hit hit_cur = hitlist.ReadNext();
//    if (nhits >= par.realign_max && nhits >= std::max(par.B, par.Z))
//      break;
//    if (hit_cur.Eval > par.e) {
//      if (nhits >= std::max(par.B, par.Z))
//        continue;
//      if (nhits >= std::max(par.b, par.z) && hit_cur.Probab < par.p)
//        continue;
//      if (nhits >= std::max(par.b, par.z) && hit_cur.Eval > par.E)
//        continue;
//    }
//
//    if (hit_cur.L > Lmax)
//      Lmax = hit_cur.L;
//    if (hit_cur.L > Lmaxmem) {
//      nhits++;
//      continue;
//    }
//
//    //fprintf(stderr,"hit.name=%-15.15s  hit.index=%-5i hit.ftellpos=%-8i  hit.dbfile=%s\n",hit_cur.name,hit_cur.index,(unsigned int)hit_cur.ftellpos,hit_cur.dbfile);
//
//    // realign the first premerge hits consecutively to query profile
//    if (nhits >= premerge) {
//      if (hit_cur.irep == 1) {
//        // For each template (therefore irep==1), store template index and position on disk in a list
//        Realign_hitpos realign_hitpos;
//        realign_hitpos.ftellpos = hit_cur.ftellpos; // stores position on disk of template for current hit
//        realign_hitpos.index = hit_cur.index; // stores index of template of current hit
//        realign_hitpos.entry = hit_cur.entry;
//        if (!phash_plist_realignhitpos->Contains(hit_cur.entry->entry->name)) {
//          List<Realign_hitpos>* newlist = new List<Realign_hitpos>;
//          phash_plist_realignhitpos->Add(hit_cur.entry->entry->name, newlist);
//        }
//        // Add template index and ftellpos to list which belongs to key dbfile in hash
//        phash_plist_realignhitpos->Show(hit_cur.entry->entry->name)->Push(
//            realign_hitpos);
//      }
//      // pointer at index is still NULL
//      if (!array_plist_phits[hit_cur.index]) {
//        List<void*>* newlist = new List<void*>; // create new list of pointers to all aligments of a template
//        array_plist_phits[hit_cur.index] = newlist; // set array[index] to newlist
//      }
//      // Push(hitlist.ReadCurrentAddress()) :  Add address of current hit in hitlist to list...
//      // array_plist_phits[hit_cur.index]-> :  pointed to by hit_cur.index'th element of array_plist_phits
//      array_plist_phits[hit_cur.index]->Push(hitlist.ReadCurrentAddress());
//    }
//
//    nhits++;
//  }
//
//  HH_LOG(LogLevel::INFO) << "Realigning " << nhits
//      << " HMM-HMM alignments using Maximum Accuracy algorithm" << std::endl;
//
//  if (Lmax > Lmaxmem) {
//    Lmax = Lmaxmem;
//
//    HH_LOG(LogLevel::WARNING)
//        << "WARNING: Realigning sequences only up to length " << Lmaxmem << "."
//        << endl;
//    HH_LOG(LogLevel::WARNING)
//        << "This is genarally unproboblematic but may lead to slightly sub-optimal alignments for longer sequences."
//        << endl;
//    HH_LOG(LogLevel::WARNING)
//        << "You can increase available memory using the -maxmem <GB> option (currently "
//        << par.maxmem << " GB)." << endl; // TODO: still to be implemented
//    HH_LOG(LogLevel::WARNING)
//        << "The maximum length realignable is approximately maxmem/query_length/(cpus+1)/8B."
//        << endl;
//  }
//
//  //////////////////////////////////////////////////////////////////////////////////
//  // start premerge:
//  // Align the first premerge templates
//  if (premerge > 0) {
//    HH_LOG(LogLevel::INFO) << "Merging " << premerge
//        << " best hits to query alignment ..." << std::endl;
//
//    int bin = 0;
//    nhits = 0;
//    hitlist.Reset();
//
//    while (!hitlist.End() && nhits < premerge) {
//      Hit hit_cur = hitlist.ReadNext();
//      // JS: removed bug on 13 Feb 13 due to which premerged hits with E-value > par.e were not realigned
//      if (hit_cur.Eval > par.e) {
//        if (nhits >= std::max(par.B, par.Z))
//          break;
//        if (nhits >= std::max(par.b, par.z) && hit_cur.Probab < par.p)
//          break;
//        if (nhits >= std::max(par.b, par.z) && hit_cur.Eval > par.E)
//          continue;
//      }
//      nhits++;
//
//      if (hit_cur.L > Lmaxmem)
//        continue;  // Don't align too long sequences due to memory limit
//
//      // Forward stream position to start of next database HMM to be realigned
//      hit[bin]->index = hit_cur.index; // give hit a unique index for HMM
//      hit[bin]->irep = 1; // Needed for min_overlap calculation in InitializeForAlignment in hhhit.C
//
//      getTemplateHMM(par, *hit_cur.entry, dbs, par.wg, format[bin], pb, S, Sim,
//          t[bin]);
//
//      HH_LOG(LogLevel::DEBUG1) << "Realigning with " << t[bin]->name
//          << std::endl;
//
//      N_aligned++;
//
//      // Prepare MAC comparison(s)
//      PrepareTemplateHMM(par, q, t[bin], format[bin], pb, R);
//      t[bin]->Log2LinTransitionProbs(1.0);
//
//      // Realign only around previous Viterbi hit
//      hit[bin]->i1 = hit_cur.i1;
//      hit[bin]->i2 = hit_cur.i2;
//      hit[bin]->j1 = hit_cur.j1;
//      hit[bin]->j2 = hit_cur.j2;
//      hit[bin]->nsteps = hit_cur.nsteps;
//      hit[bin]->i = hit_cur.i;
//      hit[bin]->j = hit_cur.j;
//      hit[bin]->realign_around_viterbi = true;
//
//      // Align q to template in *hit[bin]
//      hit[bin]->Forward(q, t[bin], par.ssm, par.min_overlap, par.loc, par.shift,
//          par.ssw, par.exclstr, S73, S33);
//      hit[bin]->Backward(q, t[bin], par.loc, par.shift, par.ssw, S73, S33);
//      hit[bin]->MACAlignment(q, t[bin], par.loc, par.mact, par.macins);
//      hit[bin]->BacktraceMAC(q, t[bin], par.corr, par.ssw, S73, S33);
//
//      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
//      hit[bin]->score = hit_cur.score;
//      hit[bin]->score_ss = hit_cur.score_ss;
//      hit[bin]->score_aass = hit_cur.score_aass;
//      hit[bin]->score_sort = hit_cur.score_sort;
//      hit[bin]->Pval = hit_cur.Pval;
//      hit[bin]->Pvalt = hit_cur.Pvalt;
//      hit[bin]->logPval = hit_cur.logPval;
//      hit[bin]->logPvalt = hit_cur.logPvalt;
//      hit[bin]->Eval = hit_cur.Eval;
//      hit[bin]->logEval = hit_cur.logEval;
//      hit[bin]->Probab = hit_cur.Probab;
//      hit[bin]->irep = hit_cur.irep;
//
//      hit[bin]->entry = hit_cur.entry;
//
//      // Replace original hit in hitlist with realigned hit
//      //hitlist.ReadCurrent().Delete();
//      hitlist.Delete().Delete();  // delete the list record and hit object
//      hitlist.Insert(*hit[bin]);
//
//      // merge only when hit length > MINCOLS_REALIGN (don't merge 1 column matches)
//      if (hit[bin]->matched_cols < MINCOLS_REALIGN)
//        continue;
//
//      // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
//      // Reading in next db MSA and merging it onto Qali
//
//      HHblitsDatabase* db = getHHblitsDatabase(*hit_cur.entry, dbs);
//      if (db == NULL) {
//        std::cerr << "Could not find database for premerge!" << std::endl;
//        continue;
//      }
//
//      Alignment Tali;
//      long ftellpos;
//      getTemplateA3M(db, hit[bin]->entry->entry->name, ftellpos, Tali);
//
//      if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
//        Qali_allseqs.MergeMasterSlave(*hit[bin], Tali, hit[bin]->name,
//            par.maxcol);
//
//      Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
//          par.qid_db, par.qsc_db, par.Ndiff_db);
//
//      Qali.MergeMasterSlave(*hit[bin], Tali, hit[bin]->name, par.maxcol);
//
//      // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
//      Qali.Compress("merged A3M file", par.cons, par.maxres, par.maxcol, par.M,
//          par.Mgaps);
//
//      // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
//      Qali.N_filtered = Qali.Filter(par.max_seqid, S, par.coverage, par.qid,
//          par.qsc, par.Ndiff);
//
//      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
//      Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons,
//          par.showcons, par.maxres, pb, Sim);
//
//      stringstream ss_tmp;
//      ss_tmp << hit[bin]->file << "__" << hit[bin]->irep;
//      premerged_hits->Add((char*) ss_tmp.str().c_str());
//
//      if (par.notags)
//        q->NeutralizeTags(pb);
//
//      // Compute substitution matrix pseudocounts?
//      if (par.nocontxt) {
//        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
//        q->PreparePseudocounts(R);
//        // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
//        q->AddAminoAcidPseudocounts(par.pc_hhm_nocontext_mode,
//            par.pc_hhm_nocontext_a, par.pc_hhm_nocontext_b,
//            par.pc_hhm_nocontext_c);
//      }
//      else {
//        // Add full context specific pseudocounts to query
//        q->AddContextSpecificPseudocounts(pc_hhm_context_engine,
//            pc_hhm_context_mode);
//      }
//
//      q->CalculateAminoAcidBackground(pb);
//      if (par.columnscore == 5 && !q->divided_by_local_bg_freqs)
//        q->DivideBySqrtOfLocalBackgroundFreqs(
//            par.half_window_size_local_aa_bg_freqs, pb);
//
//      // Transform transition freqs to lin space if not already done
//      q->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg,
//          par.gaph, par.gapi, par.gapb, par.gapb);
//      q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
//    }
//  }
//  // end premerge
//  //////////////////////////////////////////////////////////////////////////////////
//
//  // Read all HMMs whose position is given in phash_plist_realignhitpos
//#pragma omp parallel for schedule(dynamic, 1)
//  for (size_t idb = 0; idb < hits_to_realign.size(); idb++) {
//    // Can we skip dbfiles[idb] because it contains no template to be realigned?
//    if (!phash_plist_realignhitpos->Contains(hits_to_realign[idb]->entry->name))
//      continue;
//
//    // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
//    // This list is now sorted by ftellpos in ascending order to access one template after the other efficiently
//    phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->SortList();
//
//    ///////////////////////////////////////////////////////////////////////////////////////
//    // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
//    phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->Reset();
//    while (!phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->End()) {
//      // Submit jobs until no bin is free anymore
//      while (!phash_plist_realignhitpos->Show(hits_to_realign[idb]->entry->name)->End()) {
//        // Allocate free bin
//        int bin = omp_get_thread_num();
//
//        // Forward stream position to start of next database HMM to be realigned
//        Realign_hitpos hitpos_curr = phash_plist_realignhitpos->Show(
//            hits_to_realign[idb]->entry->name)->ReadNext();
//        hit[bin]->index = hitpos_curr.index; // give hit[bin] a unique index for HMM
//        hit[bin]->entry = hitpos_curr.entry;
//
//        // Give hit[bin] the pointer to the list of pointers to hitlist elements of same template (for realignment)
//        hit[bin]->plist_phits = array_plist_phits[hitpos_curr.index];
//
//        getTemplateHMM(par, *hitpos_curr.entry, dbs, par.wg, format[bin], pb, S,
//            Sim, t[bin]);
//
//        HH_LOG(LogLevel::DEBUG) << "Realigning with " << t[bin]->name
//            << std::endl;
//
//#pragma omp critical
//        {
//          N_aligned++;
//        }
//
//        RealignByWorker(par, hit[bin], q, t[bin], format[bin], pb, R, S73, S33);
//
//        //TODO???
//        break;
//      }
//    }
//  }
//
//  // Delete all hitlist entries with too short alignments
//  nhits = 0;
//  hitlist.Reset();
//  while (!hitlist.End()) {
//    Hit hit_cur = hitlist.ReadNext();
//
//    if (nhits > par.realign_max && nhits >= std::max(par.B, par.Z))
//      break;
//    if (hit_cur.Eval > par.e) {
//      if (nhits >= std::max(par.B, par.Z))
//        continue;
//      if (nhits >= std::max(par.b, par.z) && hit_cur.Probab < par.p)
//        continue;
//      if (nhits >= std::max(par.b, par.z) && hit_cur.Eval > par.E)
//        continue;
//    }
//
//    if (hit_cur.matched_cols < MINCOLS_REALIGN) {
//      HH_LOG(LogLevel::DEBUG) << "Deleting alignment of " << hit_cur.name
//          << " with length " << hit_cur.matched_cols << std::endl;
//      // delete the list record and hit object
//      hitlist.Delete().Delete();
//    }
//    nhits++;
//  }
//
//  // Delete hash phash_plist_realignhitpos with lists
//  phash_plist_realignhitpos->Reset();
//  while (!phash_plist_realignhitpos->End())
//    delete (phash_plist_realignhitpos->ReadNext()); // delete list to which phash_plist_realignhitpos->ReadNext() points
//  delete (phash_plist_realignhitpos);
//
//  // Delete array_plist_phits with lists
//  for (int index = 0; index < N_searched; index++)
//    if (array_plist_phits[index])
//      delete (array_plist_phits[index]); // delete list to which array[index] points
//  delete[] (array_plist_phits);
//}

//void HHblits::reduceRedundancyOfHitList(int n_redundancy, int query_length,
//    HitList& hitlist, HitList& reducedHitList) {
//  int* coverage = new int[query_length + 1];
//  for (int i = 0; i <= query_length; i++) {
//    coverage[i] = 0;
//  }
//
//  int total = 0;
//
//  hitlist.Reset();
//  while (!hitlist.End()) {
//    Hit hit_cur = hitlist.ReadNext();
//
//    if (!hit_cur.P_posterior) {
//      continue;
//    }
//    total++;
//
//    //calculate actual coverage
//    int length_actual_coverage = 0;
//    for (int alignment_index = 1; alignment_index <= hit_cur.nsteps;
//        alignment_index++) {
//      length_actual_coverage +=
//          (hit_cur.P_posterior[alignment_index] > 0.5) ? 1 : 0;
//    }
//
//    int length_actual_contribution = 0;
//    for (int alignment_index = 1; alignment_index <= hit_cur.nsteps;
//        alignment_index++) {
//      int i = hit_cur.i[alignment_index];
//      if (hit_cur.P_posterior[alignment_index] > 0.5
//          && coverage[i] < n_redundancy) {
//        length_actual_contribution++;
//      }
//    }
//
//    if (length_actual_contribution > 0.5 * length_actual_coverage) {
//      reducedHitList.Insert(hit_cur);
//
//      //update coverage
//      for (int alignment_index = 1; alignment_index <= hit_cur.nsteps;
//          alignment_index++) {
//        int i = hit_cur.i[alignment_index];
//
//        coverage[i] += (hit_cur.P_posterior[alignment_index] > 0.5) ? 1 : 0;
//      }
//    }
//  }
//
//  delete[] coverage;
//}

//void HHblits::recalculateAlignmentsForDifferentQSC(HitList& hitlist,
//    Alignment& Qali, char inputformat, float* qsc, size_t nqsc,
//    HitList& recalculated_hitlist) {
//  const int COV_ABS = 25;
//  int cov_tot = std::max(std::min((int) (COV_ABS / Qali.L * 100 + 0.5), 70),
//      par.coverage);
//
//  Alignment qali;
//  qali = Qali;
//  HMM* q = new HMM();
//
//  HitList realigned_viterbi_hitlist;
//
//  for (size_t qsc_index = 0; qsc_index < nqsc; qsc_index++) {
//    float actual_qsc = qsc[qsc_index];
//
//    qali.Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M,
//        par.Mgaps);
//    qali.N_filtered = qali.Filter(par.max_seqid, S, cov_tot, par.qid,
//        actual_qsc, par.Ndiff);
//    qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons,
//        par.maxres, pb, Sim, NULL, false);
//    PrepareQueryHMM(par, inputformat, q, pc_hhm_context_engine,
//        pc_hhm_context_mode, pb, R);
//
//    hitlist.Reset();
//    while (!hitlist.End()) {
//      Hit hit_ref = hitlist.ReadNext();
//
//      HMM* t = new HMM();
//
//      int format;
//      long ftellpos;
//      getTemplateHMM(par, *hit_ref.entry, dbs, 1, format, pb, S, Sim, t);
//
//      PrepareTemplateHMM(par, q, t, format, pb, R);
//
//      Hit hit;
//      hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
//      hit.self = 0;
//      hit.realign_around_viterbi = false;
//
//      for (int irep = 1; irep <= par.altali; irep++) {
//        hit.irep = irep;
//        hit.Viterbi(q, t, par.loc, par.ssm, par.maxres, par.min_overlap,
//            par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);
//
//        if (hit.irep > 1 && hit.score <= SMIN)
//          break;
//
//        hit.Backtrace(q, t, par.corr, par.ssw, S73, S33);
//        realigned_viterbi_hitlist.Push(hit);
//      }
//
//      hit.DeleteBacktraceMatrix(q->L + 2);
//
//      delete t;
//    }
//
//    realigned_viterbi_hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);
//    realigned_viterbi_hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa,
//        par.alphab, par.alphac, par.prefilter_evalue_thresh);
//
//    q->Log2LinTransitionProbs(1.0);
//
//    realigned_viterbi_hitlist.Reset();
//    while (!realigned_viterbi_hitlist.End()) {
//      Hit hit_ref = realigned_viterbi_hitlist.ReadNext();
//
//      HMM* t = new HMM();
//
//      int format;
//      long ftellpos;
//      getTemplateHMM(par, *hit_ref.entry, dbs, par.wg, format, pb, S, Sim, t);
//
//      PrepareTemplateHMM(par, q, t, format, pb, R);
//      t->Log2LinTransitionProbs(1.0);
//
//      Hit hit;
//      hit.AllocateForwardMatrix(q->L + 2, par.maxres + 1);
//      hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
//      hit.irep = 1;
//      hit.self = 0;
//
//      hit.i1 = hit_ref.i1;
//      hit.i2 = hit_ref.i2;
//      hit.j1 = hit_ref.j1;
//      hit.j2 = hit_ref.j2;
//      hit.nsteps = hit_ref.nsteps;
//      hit.i = hit_ref.i;
//      hit.j = hit_ref.j;
//      hit.realign_around_viterbi = false;
//
//      // Align q to template in *hit[bin]
//      hit.Forward(q, t, par.ssm, par.min_overlap, par.loc, par.shift, par.ssw,
//          par.exclstr, S73, S33);
//      hit.Backward(q, t, par.loc, par.shift, par.ssw, S73, S33);
//      hit.MACAlignment(q, t, par.loc, par.mact, par.macins);
//      hit.BacktraceMAC(q, t, par.corr, par.ssw, S73, S33);
//
//      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
//      hit.score = hit_ref.score;
//      hit.score_ss = hit_ref.score_ss;
//      hit.score_aass = hit_ref.score_aass;
//      hit.score_sort = hit_ref.score_sort;
//      hit.Pval = hit_ref.Pval;
//      hit.Pvalt = hit_ref.Pvalt;
//      hit.logPval = hit_ref.logPval;
//      hit.logPvalt = hit_ref.logPvalt;
//      hit.Eval = hit_ref.Eval;
//      hit.logEval = hit_ref.logEval;
//      hit.Probab = hit_ref.Probab;
//
//      hit.DeleteForwardMatrix(q->L + 2);
//      hit.DeleteBacktraceMatrix(q->L + 2);
//
//      if (hit.matched_cols >= MINCOLS_REALIGN) {
//        recalculated_hitlist.Insert(hit);
//      }
//
//      delete t;
//    }
//
//    realigned_viterbi_hitlist.Reset();
//    while (!realigned_viterbi_hitlist.End()) {
//      realigned_viterbi_hitlist.Delete().Delete();
//    }
//  }
//
//  delete q;
//}

//void HHblits::uniqueHitlist(HitList& hitlist) {
//  std::set<std::string> ids;
//
//  hitlist.Reset();
//  while (!hitlist.End()) {
//    Hit hit_cur = hitlist.ReadNext();
//
//    std::string id = std::string(hit_cur.name);
//
//    if (ids.find(id) != ids.end()) {
//      hitlist.Delete();
//    }
//    else {
//      ids.insert(id);
//    }
//  }
//}

//void HHblits::wiggleQSC(HitList& hitlist, int n_redundancy, Alignment& Qali,
//    char inputformat, float* qsc, size_t nqsc, HitList& reducedFinalHitList) {
//  int query_length = Qali.L;
//  HitList* reducedHitList = new HitList();
//  //  HitList* wiggledHitList = new HitList();
//
//  //filter by 2*n_redundancy
//  reduceRedundancyOfHitList(2 * n_redundancy, query_length, hitlist,
//      *reducedHitList);
//
//  //recalculate alignments for different qsc and choose best alignment for each template
//  recalculateAlignmentsForDifferentQSC(*reducedHitList, Qali, inputformat, qsc,
//      nqsc, reducedFinalHitList);  //*wiggledHitList);
//  //  uniqueHitlist(*wiggledHitList);
//
//  //filter by n_redundancy
//  //  reduceRedundancyOfHitList(n_redundancy, query_length, *wiggledHitList, reducedFinalHitList);
//}
//
//void HHblits::wiggleQSC(int n_redundancy, float* qsc, size_t nqsc,
//    HitList& reducedFinalHitList) {
//  //TODO: parameter
//  char input_format = 1;
//  wiggleQSC(hitlist, n_redundancy, Qali, input_format, qsc, nqsc,
//      reducedFinalHitList);
//  reducedFinalHitList.N_searched = hitlist.N_searched;
//}

void HHblits::run(FILE* query_fh, char* query_path) {
  int cluster_found = 0;
  int seqs_found = 0;
  int premerge = par.premerge;

  SearchCounter search_counter;

  Hit hit_cur;
  Hash<Hit>* previous_hits = new Hash<Hit>(1631, hit_cur);
  Hash<char>* premerged_hits = new Hash<char>(1631);

  Qali = new Alignment();
  Qali_allseqs = new Alignment();

  q = new HMM(MAXSEQDIS, par.maxres);
  HMMSimd q_vec(par.maxres);
  q_tmp = new HMM(MAXSEQDIS, par.maxres);

  // Read query input file (HHM or alignment format) without adding pseudocounts
  Qali->N_in = 0;
  char input_format;
  ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb, S,
                Sim);

  if (Qali->N_in - Qali->N_ss > 1)
    premerge = 0;

  if (par.allseqs) {
    *Qali_allseqs = *Qali;  // make a *deep* copy of Qali!
    for (int k = 0; k < Qali_allseqs->N_in; ++k)
      Qali_allseqs->keep[k] = 1;  // keep *all* sequences (reset filtering in Qali)
  }

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags)
    q->NeutralizeTags(pb);

  //save all entries pointer in this vector to delete, when it's safe
  std::vector<HHEntry*> all_entries;
  std::vector<HHEntry*> new_entries;
  std::vector<HHEntry*> old_entries;

  if (!par.prefilter) {
    for (size_t i = 0; i < dbs.size(); i++) {
      dbs[i]->initNoPrefilter(new_entries);
    }
    all_entries.insert(all_entries.end(), new_entries.begin(),
                       new_entries.end());

    for (size_t i = 0; i < all_entries.size(); i++) {
      search_counter.append(std::string(all_entries[i]->getName()));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  for (int round = 1; round <= par.num_rounds; round++) {
    HH_LOG(INFO) << "Iteration " << round << std::endl;

    // Settings for different rounds
    if (premerge > 0 && round > 1 && previous_hits->Size() >= premerge) {
      HH_LOG(DEBUG1) << "Set premerge to 0! (premerge: " << premerge
                               << "   iteration: " << round << "   hits.Size: "
                               << previous_hits->Size() << ")" << std::endl;
      premerge = 0;
    } else {
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
    q_vec.MapOneHMM(q);

    ////////////////////////////////////////////
    // Prefiltering
    ////////////////////////////////////////////

    if (par.prefilter) {
      HH_LOG(INFO) << "Prefiltering database" << std::endl;

      new_entries.clear();
      old_entries.clear();

      // Add Pseudocounts to q_tmp
      if (par.nocontxt) {
        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        q_tmp->PreparePseudocounts(R);
        // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        q_tmp->AddAminoAcidPseudocounts(par.pc_prefilter_nocontext_mode,
                                        par.pc_prefilter_nocontext_a,
                                        par.pc_prefilter_nocontext_b,
                                        par.pc_prefilter_nocontext_c);
      } else {
        // Add context specific pseudocounts (now always used, because clusterfile is necessary)
        q_tmp->AddContextSpecificPseudocounts(pc_prefilter_context_engine,
                                              pc_prefilter_context_mode);
      }

      q_tmp->CalculateAminoAcidBackground(pb);

      for (size_t i = 0; i < dbs.size(); i++) {
        dbs[i]->prefilter_db(q_tmp, previous_hits, par.threads,
                             par.prefilter_gap_open, par.prefilter_gap_extend,
                             par.prefilter_score_offset,
                             par.prefilter_bit_factor,
                             par.prefilter_evalue_thresh,
                             par.prefilter_evalue_coarse_thresh,
                             par.preprefilter_smax_thresh,
                             par.min_prefilter_hits, par.maxnumdb, R,
                             new_entries, old_entries);
      }

      for (size_t i = 0; i < new_entries.size(); i++) {
        search_counter.append(std::string(new_entries[i]->getName()));
      }

      all_entries.insert(all_entries.end(), new_entries.begin(),
                         new_entries.end());
      all_entries.insert(all_entries.end(), old_entries.begin(),
                         old_entries.end());
    }

    int max_template_length = getMaxTemplateLength(new_entries);
    if (max_template_length > par.maxres) {
      HH_LOG(WARNING)
          << "database contains sequnces that exceeds maximum allowed size (maxres = "
          << par.maxres << "). Maxres can be increased with parameter -maxres."
          << std::endl;
    }
    max_template_length = std::min(max_template_length, par.maxres);
    for (int i = 0; i < par.threads; i++) {
      viterbiMatrices[i]->AllocateBacktraceMatrix(q->L, max_template_length);
    }

    hitlist.N_searched = search_counter.getCounter();

    if (new_entries.size() == 0) {
      HH_LOG(INFO) << "No HMMs pass prefilter => Stop searching!"
                             << std::endl;
      break;
    }

    // Search datbases
    HH_LOG(INFO)
        << "HMMs passed 2nd prefilter (gapped profile-profile alignment)   : "
        << new_entries.size() + old_entries.size() << std::endl;
    HH_LOG(INFO)
        << "HMMs passed 2nd prefilter and not found in previous iterations : "
        << new_entries.size() << std::endl;
    HH_LOG(INFO) << "Scoring " << new_entries.size()
                           << " HMMs using HMM-HMM Viterbi alignment"
                           << std::endl;

    // Main Viterbi HMM-HMM search
    ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
    std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec,
                                                           new_entries,
                                                           par.qsc_db, pb, S,
                                                           Sim, R);

    add_hits_to_hitlist(hits_to_add, hitlist);

    // check for new hits or end with iteration
    int new_hits = 0;
    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();
      if (hit_cur.Eval > 100.0 * par.e)
        break;  // E-value much too large
      if (hit_cur.Eval > par.e)
        continue;  // E-value too large
      new_hits++;
    }

    if (new_hits == 0 || round == par.num_rounds) {
      if (round < par.num_rounds) {
        HH_LOG(INFO) << "No new hits found in iteration " << round
                               << " => Stop searching" << std::endl;
      }

      if (old_entries.size() > 0 && par.realign_old_hits) {
        HH_LOG(INFO)
            << "Rescoring previously found HMMs with Viterbi algorithm"
            << std::endl;

        ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
        std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec,
                                                               old_entries,
                                                               par.qsc_db, pb,
                                                               S, Sim, R);

        add_hits_to_hitlist(hits_to_add, hitlist);

//				// Add dbfiles_old to dbfiles_new for realign
        new_entries.insert(new_entries.end(), old_entries.begin(),
                           old_entries.end());
      } else if (!par.realign_old_hits && previous_hits->Size() > 0) {
        HH_LOG(INFO)
            << "Rescoring previously found HMMs with Viterbi algorithm"
            << std::endl;
        RescoreWithViterbiKeepAlignment(q_vec, previous_hits);
      }
    }

    // Realign hits with MAC algorithm
    if (par.realign)
      perform_realign(q_vec, input_format, new_entries, premerge,
                      premerged_hits);

    // Generate alignment for next iteration
    if (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile
        || *par.alisbasename) {
      if (new_hits > 0) {
        mergeHitsToQuery(previous_hits, premerged_hits, seqs_found,
                         cluster_found);
      }

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons,
                                      par.showcons, par.maxres, pb, Sim, NULL,
                                      true);

      if (par.notags)
        q->NeutralizeTags(pb);

      // Calculate SSpred if we need to print out alis after each iteration or if last iteration
      // TODO: we should get rid of this... since it calls psipred on the command line and is untested
      if (par.addss
          && (*par.alisbasename || round == par.num_rounds || new_hits == 0)) {
        char ss_pred[par.maxres];
        char ss_conf[par.maxres];

        CalculateSS(q, ss_pred, ss_conf, par.psipred_data, par.psipred, pb);

        Qali->AddSSPrediction(ss_pred, ss_conf);
      }

      if (*par.alisbasename) {
        Alignment* tmp = new Alignment();
        if (par.allseqs) {
          (*tmp) = (*Qali_allseqs);
        } else {
          (*tmp) = (*Qali);
        }

        alis[round] = tmp;
      }
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

    HH_LOG(INFO) << seqs_found << " sequences belonging to "
                           << cluster_found
                           << " database HMMs found with an E-value < " << par.e
                           << std::endl;

    if (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile
        || *par.alisbasename) {
      HH_LOG(INFO)
          << "Number of effective sequences of resulting query HMM: Neff = "
          << q->Neff_HMM << std::endl;
    }

    if (q->Neff_HMM > par.neffmax && round < par.num_rounds) {
      HH_LOG(INFO)
          << "Diversity is above threshold (" << par.neffmax
          << "). Stop searching! (Change threshold using -neffmax <float>.)"
          << std::endl;
    }

    if (Qali->N_in >= MAXSEQ) {
      HH_LOG(INFO)
          << "Maximun number of sequences in query alignment reached ("
          << MAXSEQ << "). Stop searching!" << std::endl;
    }

    if (new_hits == 0 || round == par.num_rounds || q->Neff_HMM > par.neffmax
        || Qali->N_in >= MAXSEQ)
      break;

    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();
      char strtmp[NAMELEN + 6];
      sprintf(strtmp, "%s__%i%c", hit_cur.file, hit_cur.irep, '\0');
      if (!par.already_seen_filter || hit_cur.Eval > par.e
          || previous_hits->Contains(strtmp))
        hit_cur.Delete();  // Delete hit object (deep delete with Hit::Delete())
      else {
        previous_hits->Add(strtmp, hit_cur);
      }

      hitlist.Delete();  // Delete list record (flat delete)
    }
  }

  if (*par.opt_outfile) {
    std::vector<HHEntry*> selected_entries;
    get_entries_of_selected_hits(hitlist, selected_entries);
    optimizeQSC(selected_entries, hitlist.N_searched, q_vec, input_format,
                optimized_hitlist);
  }

  // Warn, if HMMER files were used
  if (par.hmmer_used) {
    HH_LOG(WARNING)
        << "WARNING: Using HMMER files results in a drastically reduced sensitivity (>10%%).\n"
        << " We recommend to use HHMs build by hhmake." << std::endl;
  }

//  if (*par.reduced_outfile) {
//    float qscs[] = { -20, 0, 0.1, 0.2 };
//    wiggleQSC(par.n_redundancy, qscs, 4, reducedHitlist);
//  }

  for (size_t i = 0; i < all_entries.size(); i++) {
    delete all_entries[i];
  }
  all_entries.clear();

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete();  // Delete hit object
  delete previous_hits;

  delete premerged_hits;
}

void HHblits::writeHHRFile(char* hhrFile) {
  if (*hhrFile) {
    hitlist.PrintHHR(q_tmp, hhrFile, par.maxdbstrlen, par.showconf,
                     par.showcons, par.showdssp, par.showpred, par.b, par.B,
                     par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E,
                     par.argc, par.argv, S);
  }
}

void HHblits::writeAlisFile(char* basename) {
  if (*basename) {
    std::map<int, Alignment*>::iterator it;
    for (it = alis.begin(); it != alis.end(); it++) {
      stringstream ss_tmp;
      ss_tmp << basename << "_" << (*it).first << ".a3m";
      std::string id = ss_tmp.str();

      (*it).second->WriteToFile(id.c_str(), par.append, "a3m");
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
                            par.showdssp, par.showpred, par.p, par.aliwidth,
                            par.nseqdis, par.b, par.B, par.E, S, outformat);
  }
}

void HHblits::writeAlitabFile(char* alitabFile) {
  if (*alitabFile) {
    hitlist.WriteToAlifile(q, alitabFile, par.b, par.B, par.z, par.Z, par.p,
                           par.E);
  }
}

void HHblits::writeOptimizedHHRFile(char* reducedHHRFile) {
  if (*reducedHHRFile) {
    optimized_hitlist.PrintHHR(q_tmp, reducedHHRFile, par.maxdbstrlen,
                               par.showconf, par.showcons, par.showdssp,
                               par.showpred, par.b, par.B, par.z, par.Z,
                               par.aliwidth, par.nseqdis, par.p, par.E,
                               par.argc, par.argv, S);
  }
}

void HHblits::writePsiFile(char* psiFile) {
  // Write output PSI-BLAST-formatted alignment?
  if (*psiFile) {
    if (par.allseqs)
      Qali_allseqs->WriteToFile(psiFile, par.append, "psi");
    else
      Qali->WriteToFile(psiFile, par.append, "psi");
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
      Qali_allseqs->WriteToFile(A3MFile, par.append, "a3m");
    else
      Qali->WriteToFile(A3MFile, par.append, "a3m");
  }
}

std::map<int, Alignment*>& HHblits::getAlis() {
  return alis;
}

void HHblits::writeHHRFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintHHR(hhblits.q_tmp, out, hhblits.par.maxdbstrlen,
                           hhblits.par.showconf, hhblits.par.showcons,
                           hhblits.par.showdssp, hhblits.par.showpred,
                           hhblits.par.b, hhblits.par.B, hhblits.par.z,
                           hhblits.par.Z, hhblits.par.aliwidth,
                           hhblits.par.nseqdis, hhblits.par.p, hhblits.par.E,
                           hhblits.par.argc, hhblits.par.argv, hhblits.S);
}

void HHblits::writeScoresFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintScoreFile(hhblits.q, out);
}

void HHblits::writePairwiseAlisFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintAlignments(hhblits.q, out, hhblits.par.showconf,
                                  hhblits.par.showcons, hhblits.par.showdssp,
                                  hhblits.par.showpred, hhblits.par.p,
                                  hhblits.par.aliwidth, hhblits.par.nseqdis,
                                  hhblits.par.b, hhblits.par.B, hhblits.par.E,
                                  hhblits.S, hhblits.par.outformat);
}

void HHblits::writeAlitabFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.WriteToAlifile(hhblits.q, out, hhblits.par.b, hhblits.par.B,
                                 hhblits.par.z, hhblits.par.Z, hhblits.par.p,
                                 hhblits.par.E);
}

void HHblits::writeOptimizedHHRFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.optimized_hitlist.PrintHHR(hhblits.q_tmp, out,
                                     hhblits.par.maxdbstrlen,
                                     hhblits.par.showconf, hhblits.par.showcons,
                                     hhblits.par.showdssp, hhblits.par.showpred,
                                     hhblits.par.b, hhblits.par.B,
                                     hhblits.par.z, hhblits.par.Z,
                                     hhblits.par.aliwidth, hhblits.par.nseqdis,
                                     hhblits.par.p, hhblits.par.E,
                                     hhblits.par.argc, hhblits.par.argv,
                                     hhblits.S);
}

void HHblits::writePsiFile(HHblits& hhblits, std::stringstream& out) {
  if (hhblits.par.allseqs)
    hhblits.Qali_allseqs->WriteToFile(out, "psi");
  else
    hhblits.Qali->WriteToFile(out, "psi");
}

void HHblits::writeHMMFile(HHblits& hhblits, std::stringstream& out) {
  // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
  hhblits.q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
  hhblits.q->CalculateAminoAcidBackground(hhblits.pb);

  hhblits.q->WriteToFile(out, hhblits.par.max_seqid, hhblits.par.coverage,
                         hhblits.par.qid, hhblits.par.Ndiff, hhblits.par.qsc,
                         hhblits.par.argc, hhblits.par.argv, hhblits.pb);
}

void HHblits::writeA3MFile(HHblits& hhblits, std::stringstream& out) {
  if (hhblits.par.allseqs)
    hhblits.Qali_allseqs->WriteToFile(out, "a3m");
  else
    hhblits.Qali->WriteToFile(out, "a3m");
}

void HHblits::writeMatricesFile(char* matricesOutputFileName) {
  if (*matricesOutputFileName) {
    hitlist.PrintMatrices(q, matricesOutputFileName, par.max_number_matrices,
                          S);
  }
}

void HHblits::writeMatricesFile(HHblits& hhblits, stringstream& out) {
  hhblits.hitlist.PrintMatrices(hhblits.q, out, hhblits.par.max_number_matrices,
                                hhblits.S);
}

