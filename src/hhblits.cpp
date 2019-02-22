/*
 * HHblits.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include "hhblits.h"
#include "hhsuite_config.h"

HHblits::HHblits(Parameters& par, std::vector<HHblitsDatabase*>& databases) : par(par) {
  dbs = databases;

  context_lib = NULL;
  crf = NULL;
  pc_hhm_context_engine = NULL;
  pc_hhm_context_mode = NULL;
  pc_prefilter_context_engine = NULL;
  pc_prefilter_context_mode = NULL;

  Qali = NULL;
  Qali_allseqs = NULL;

  q = NULL;
  q_tmp = NULL;

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix(par.matrix, pb, P, R, S, Sim);

  // Set secondary structure substitution matrix
  if (par.ssm)
      SetSecStrucSubstitutionMatrix(par.ssa, S73, S37, S33);

  // Prepare pseudocounts
  if (!par.nocontxt) {
    InitializePseudocountsEngine(par, context_lib, crf, pc_hhm_context_engine, pc_hhm_context_mode, pc_prefilter_context_engine, pc_prefilter_context_mode);
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

  DeletePseudocountsEngine(context_lib, crf, pc_hhm_context_engine, pc_hhm_context_mode, pc_prefilter_context_engine, pc_prefilter_context_mode);
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

void HHblits::ProcessAllArguments(Parameters& par) {
  const int argc = par.argc;
  const char** argv = par.argv;

  par.Ndiff = 1000;
  par.prefilter = true;

  par.early_stopping_filter = true;
  par.filter_thresh = 0.01;

  ProcessArguments(par);

  // Check needed files
  if (!*par.infile || !strcmp(par.infile, "")) {
    help(par);
    HH_LOG(ERROR) << "Input file is missing! (see -i)" << std::endl;
    exit(4);
  }
  if (par.db_bases.size() == 0) {
    help(par);
    HH_LOG(ERROR) << "Database is missing! (see -d)" << std::endl;
    exit(4);
  }

  if (par.loc == 0 && par.num_rounds >= 2) {
    HH_LOG(WARNING)
        << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << "\n"
        << "\tusing -global alignment for iterative searches is deprecated "
        "since non-homologous sequence segments can easily enter the "
        "MSA and corrupt it."
        << std::endl;
  }

  if (par.num_rounds < 1)
    par.num_rounds = 1;
  else if (par.num_rounds > 8) {
    HH_LOG(WARNING) << "Number of iterations (" << par.num_rounds << ") too large => Set to 8 iterations\n";
    par.num_rounds = 8;
  }

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
    HH_LOG(WARNING) << "Setting -maxmem to its minimum allowed value of 1.0" << std::endl;
    par.maxmem = 1.0;
  }
  if (par.mact >= 1.0)
    par.mact = 0.999;
  else if (par.mact < 0)
    par.mact = 0.0;
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
  hitlist.Reset();

  std::map<int, Alignment*>::iterator it;
  for (it = alis.begin(); it != alis.end(); it++) {
    delete (*it).second;
  }
  alis.clear();

  // Prepare multi-threading - reserve memory for threads, intialize, etc.
  for (int bin = 0; bin < par.threads; bin++) {
    viterbiMatrices[bin]->DeleteBacktraceMatrix();
    posteriorMatrices[bin]->DeleteProbabilityMatrix();
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void HHblits::help(Parameters& par, char all) {
  printf("HHblits %i.%i.%i:\nHMM-HMM-based lightning-fast iterative sequence search\n", HHSUITE_VERSION_MAJOR, HHSUITE_VERSION_MINOR, HHSUITE_VERSION_PATCH);
  printf("HHblits is a sensitive, general-purpose, iterative sequence search tool that represents\n");
  printf("both query and database sequences by HMMs. You can search HHblits databases starting\n");
  printf("with a single query sequence, a multiple sequence alignment (MSA), or an HMM. HHblits\n");
  printf("prints out a ranked list of database HMMs/MSAs and can also generate an MSA by merging\n");
  printf("the significant database HMMs/MSAs onto the query MSA.\n");
  printf("\n");
  printf("%s", HHBLITS_REFERENCE);
  printf("%s", COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query [options] \n", par.program_name);
  printf(" -i <file>      input/query: single sequence or multiple sequence alignment (MSA)\n");
  printf("                in a3m, a2m, or FASTA format, or HMM in hhm format\n");

  if (all) {
    printf("\n");
    printf("<file> may be 'stdin' or 'stdout' throughout.\n");
  }
  printf("\n");

  printf("Options:                                                                        \n");
  printf(" -d <name>      database name (e.g. uniprot20_29Feb2012)                        \n");
  printf("                Multiple databases may be specified with '-d <db1> -d <db2> ...'\n");
  printf(" -n     [1,8]   number of iterations (default=%i)                               \n", par.num_rounds);
  printf(" -e     [0,1]   E-value cutoff for inclusion in result alignment (def=%G)       \n", par.e);
  printf("\n");

  printf("Input alignment format:                                                       \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               ' -' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf(" -tags/-notags  do NOT / do neutralize His-, C-myc-, FLAG-tags, and trypsin \n");
  printf("                recognition sequence to background distribution (def=-notags)  \n");
  printf("\n");

  printf("Output options: \n");
  printf(" -o <file>      write results in standard format to file (default=<infile.hhr>)\n");
  printf(" -oa3m <file>   write result MSA with significant matches in a3m format\n");
  if (!all) {
    printf("                Analogous for -opsi and -ohhm\n");
  }
  if (all) {
    printf(" -opsi <file>   write result MSA of significant matches in PSI-BLAST format\n");
    printf(" -ohhm <file>   write HHM file for result MSA of significant matches\n");
  }
  printf(" -oalis <name>  write MSAs in A3M format after each iteration\n");
  printf(" -add_cons      generate consensus sequence as master sequence of query MSA (default=don't)\n");
  printf(" -hide_cons     don't show consensus sequence in alignments (default=show)     \n");
  printf(" -hide_pred     don't show predicted 2ndary structure in alignments (default=show)\n");
  printf(" -hide_dssp     don't show DSSP 2ndary structure in alignments (default=show)  \n");
  printf(" -show_ssconf   show confidences for predicted 2ndary structure in alignments\n");
  if (all) {
    printf(" -Ofas <file>   write pairwise alignments in FASTA xor A2M (-Oa2m) xor A3M (-Oa3m) format   \n");
    printf(" -seq <int>     max. number of query/template sequences displayed (default=%i)  \n", par.nseqdis);
    printf(" -aliw <int>    number of columns per line in alignment list (default=%i)       \n", par.aliwidth);
    printf(" -p [0,100]     minimum probability in summary and alignment list (default=%G)  \n", par.p);
    printf(" -E [0,inf[     maximum E-value in summary and alignment list (default=%G)      \n", par.E);
    printf(" -Z <int>       maximum number of lines in summary hit list (default=%i)        \n", par.Z);
    printf(" -z <int>       minimum number of lines in summary hit list (default=%i)        \n", par.z);
    printf(" -B <int>       maximum number of alignments in alignment list (default=%i)     \n", par.B);
    printf(" -b <int>       minimum number of alignments in alignment list (default=%i)     \n", par.b);
    printf("\n");
  }


  if(all) {
    printf("Prefilter options                                                               \n");
    printf(" -noprefilt                disable all filter steps                                        \n");
    printf(" -noaddfilter              disable all filter steps (except for fast prefiltering)         \n");
    printf(" -maxfilt                  max number of hits allowed to pass 2nd prefilter (default=%i)   \n", par.maxnumdb);
    printf(" -min_prefilter_hits       min number of hits to pass prefilter (default=%i)               \n", par.min_prefilter_hits);
    printf(" -prepre_smax_thresh       min score threshold of ungapped prefilter (default=%i)               \n", par.preprefilter_smax_thresh);
    printf(" -pre_evalue_thresh        max E-value threshold of Smith-Waterman prefilter score (default=%.1f)\n", par.prefilter_evalue_thresh);
    printf(" -pre_bitfactor            prefilter scores are in units of 1 bit / pre_bitfactor (default=%i)\n", par.prefilter_bit_factor);
    printf(" -pre_gap_open             gap open penalty in prefilter Smith-Waterman alignment (default=%i)\n", par.prefilter_gap_open);
    printf(" -pre_gap_extend           gap extend penalty in prefilter Smith-Waterman alignment (default=%i)\n", par.prefilter_gap_extend);
    printf(" -pre_score_offset         offset on sequence profile scores in prefilter S-W alignment (default=%i)\n", par.prefilter_score_offset);

    printf("\n");
  }

  printf("Filter options applied to query MSA, database MSAs, and result MSA              \n");
  printf(" -all           show all sequences in result MSA; do not filter result MSA      \n");
  printf(" -interim_filter NONE|FULL  \n");
  printf("                filter sequences of query MSA during merging to avoid early stop (default: FULL)\n");
  printf("                  NONE: disables the intermediate filter \n");
  printf("                  FULL: if an early stop occurs compare filter seqs in an all vs. all comparison\n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (def=%i)\n", par.max_seqid);
  printf(" -diff [0,inf[  filter MSAs by selecting most diverse set of sequences, keeping \n");
  printf("                at least this many seqs in each MSA block of length 50 \n");
  printf("                Zero and non-numerical values turn off the filtering. (def=%i) \n", par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with master sequence (%%) (def=%i)             \n", par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with master sequence (%%) (def=%i)    \n", par.qid);
  printf(" -qsc  [0,100]  minimum score per column with master sequence (default=%.1f)    \n", par.qsc);
  printf(" -neff [1,inf]  target diversity of multiple sequence alignment (default=off)   \n");
  printf(" -mark          do not filter out sequences marked by \">@\"in their name line  \n");
  printf("\n");

  printf("HMM-HMM alignment options:                                                       \n");
  printf(" -norealign           do NOT realign displayed hits with MAC algorithm (def=realign)   \n");
  printf(" -realign_old_hits    realign hits from previous iterations                          \n");
  printf(" -mact [0,1[          posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                      ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n", par.mact);
  printf(" -glob/-loc           use global/local alignment mode for searching/ranking (def=local)\n");
  if (all) {
    printf(" -realign             realign displayed hits with max. accuracy (MAC) algorithm \n");
    printf(" -realign_max <int>   realign max. <int> hits (default=%i)                        \n", par.realign_max);
    printf(" -ovlp <int>          banded alignment: forbid <ovlp> largest diagonals |i-j| of DP matrix (def=%i)\n", par.min_overlap);
    printf(" -alt <int>           show up to this many alternative alignments with raw score > smin(def=%i)  \n", par.altali);
    printf(" -smin <float>        minimum raw score for alternative alignments (def=%.1f)  \n", par.smin);
    printf(" -shift [-1,1]        profile-profile score offset (def=%-.2f)                         \n", par.shift);
    printf(" -corr [0,1]          weight of term for pair correlations (def=%.2f)                \n", par.corr);
    printf(" -sc   <int>          amino acid score         (tja: template HMM at column j) (def=%i)\n", par.columnscore);
    printf("              0       = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
    printf("              1       = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
    printf("              2       = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
    printf("              3       = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
    printf("              5       local amino acid composition correction                     \n");
    printf(" -ssm {0,..,4}        0:   no ss scoring                                             \n");
    printf("                      1,2: ss scoring after or during alignment  [default=%1i]         \n", par.ssm);
    printf("                      3,4: ss scoring after or during alignment, predicted vs. predicted\n");
    printf(" -ssw [0,1]           weight of ss score  (def=%-.2f)                                  \n", par.ssw);
    printf(" -ssa [0,1]           ss confusion matrix = (1-ssa)*I + ssa*psipred-confusion-matrix [def=%-.2f)\n", par.ssa);
    printf(" -wg                  use global sequence weighting for realignment!                   \n");
    printf("\n");

    printf("Gap cost options:                                                                \n");
    printf(" -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                     \n", par.gapb);
    printf(" -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)    \n", par.gapd);
    printf(" -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)      \n", par.gape);
    printf(" -gapf ]0,inf]  factor to increase/reduce gap open penalty for deletes (def=%-.2f) \n", par.gapf);
    printf(" -gapg ]0,inf]  factor to increase/reduce gap open penalty for inserts (def=%-.2f) \n", par.gapg);
    printf(" -gaph ]0,inf]  factor to increase/reduce gap extend penalty for deletes(def=%-.2f)\n", par.gaph);
    printf(" -gapi ]0,inf]  factor to increase/reduce gap extend penalty for inserts(def=%-.2f)\n", par.gapi);
    printf(" -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f) \n",  par.egq);
    printf(" -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)\n", par.egt);
    printf("\n");

    printf("Pseudocount (pc) options:                                                        \n");
    printf(" Context specific hhm pseudocounts:\n");
    printf("  -pc_hhm_contxt_mode {0,..,3}   position dependence of pc admixture 'tau' (pc mode, default=%-i) \n", par.pc_hhm_context_engine.admix);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_hhm_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n", par.pc_hhm_context_engine.pca);
    printf("  -pc_hhm_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n", par.pc_hhm_context_engine.pcb);
    printf("  -pc_hhm_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n", par.pc_hhm_context_engine.pcc);
    printf("\n");

    printf(" Context independent hhm pseudocounts (used for templates; used for query if contxt file is not available):\n");
    printf("  -pc_hhm_nocontxt_mode {0,..,3}   position dependence of pc admixture 'tau' (pc mode, default=%-i) \n", par.pc_hhm_nocontext_mode);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_hhm_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n", par.pc_hhm_nocontext_a);
    printf("  -pc_hhm_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n", par.pc_hhm_nocontext_b);
    printf("  -pc_hhm_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n", par.pc_hhm_nocontext_c);
    printf("\n");

    printf(" Context specific prefilter pseudocounts:\n");
    printf("  -pc_prefilter_contxt_mode {0,..,3}   position dependence of pc admixture 'tau' (pc mode, default=%-i) \n", par.pc_prefilter_context_engine.admix);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_prefilter_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n", par.pc_prefilter_context_engine.pca);
    printf("  -pc_prefilter_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n", par.pc_prefilter_context_engine.pcb);
    printf("  -pc_prefilter_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n", par.pc_prefilter_context_engine.pcc);
    printf("\n");

    printf(" Context independent prefilter pseudocounts (used if context file is not available):\n");
    printf("  -pc_prefilter_nocontxt_mode {0,..,3}   position dependence of pc admixture 'tau' (pc mode, default=%-i) \n", par.pc_prefilter_nocontext_mode);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_prefilter_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n", par.pc_prefilter_nocontext_a);
    printf("  -pc_prefilter_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n", par.pc_prefilter_nocontext_b);
    printf("  -pc_prefilter_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n", par.pc_prefilter_nocontext_c);
    printf("\n");

    printf(" Context-specific pseudo-counts:                                                  \n");
    printf("  -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
    printf("  -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n", par.clusterfile.c_str());
    printf("  -csw  [0,inf]  weight of central position in cs pseudocount mode (def=%.1f)\n", par.csw);
    printf("  -csb  [0,1]    weight decay parameter for positions in cs pc mode (def=%.1f)\n", par.csb);
    printf("\n");
  }

  printf("Other options:                                                                   \n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose (def=%i)\n", par.v);
  printf(" -neffmax ]1,20] skip further search iterations when diversity Neff of query MSA \n");
  printf("                becomes larger than neffmax (default=%.1f)\n", par.neffmax);
  printf(" -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=%i)      \n", par.threads);
  if (all) {
	printf(" -scores <file> write scores for all pairwise comparisons to file               \n");
	  printf(" -filter_matrices filter matrices for similarity to output at most 100 matrices\n");
    printf(" -atab   <file> write all alignments in tabular layout to file                   \n");
    printf(" -maxseq <int>  max number of input rows (def=%5i)\n", par.maxseq);
    printf(" -maxres <int>  max number of HMM columns (def=%5i)\n", par.maxres);
    printf(" -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n", par.maxmem);
  }
  printf("\n");

  if (!all) {
    printf("An extended list of options can be obtained by calling 'hhblits -h all'\n");
    printf("\n");
  }
  printf("Examples:\n");
  printf("hhblits -i query.fas -o query.hhr -d ./uniclust30\n");
  printf("\n");
  printf("hhblits -i query.fas -o query.hhr -oa3m query.a3m -n 1 -d ./uniclust30\n");
  printf("\n");
  printf("Download databases from <http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/>.\n");
}


void HHblits::ProcessArguments(Parameters& par) {
  const int argc = par.argc;
  const char** argv = par.argv;

  //Processing command line input
  for (int i = 1; i < argc; i++) {
    HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No query file following -i" << std::endl;
        exit(4);
      } else
        strcpy(par.infile, argv[i]);
    } else if (!strcmp(argv[i], "-d")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No database basename following -d" << std::endl;
        exit(4);
      } else {
        std::string db(argv[i]);
        if (HHDatabase::checkDatabaseConflicts(argv[i])) {
          HH_LOG(ERROR) << "Ambiguous database basename. Choose either a A3M or CA3M database." << std::endl;
          exit(4);
        }
        par.db_bases.push_back(db);
      }
    } else if (!strcmp(argv[i], "-contxt")
        || !strcmp(argv[i], "-context_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No lib following -contxt" << std::endl;
        exit(4);
      } else
        par.clusterfile = argv[i];
    } else if (!strcmp(argv[i], "-cslib") || !strcmp(argv[i], "-cs_lib")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No lib following -cslib" << std::endl;
        exit(4);
      } else
        par.cs_library = argv[i];
    } else if (!strcmp(argv[i], "-o")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      } else
        strcpy(par.outfile, argv[i]);
    }
    //no help required
    else if (!strcmp(argv[i], "-omat")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -omat" << std::endl;
        exit(4);
      } else
        strcpy(par.matrices_output_file, argv[i]);
    } else if (!strcmp(argv[i], "-oa3m")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -oa3m" << std::endl;
        exit(4);
      } else
        strcpy(par.alnfile, argv[i]);
    } else if (!strcmp(argv[i], "-filter_matrices")) {
        par.filter_matrices = true;
    } else if (!strcmp(argv[i], "-ohhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -ohhm" << std::endl;
        exit(4);
      } else
        strcpy(par.hhmfile, argv[i]);
    } else if (!strcmp(argv[i], "-opsi")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -opsi" << std::endl;
        exit(4);
      } else
        strcpy(par.psifile, argv[i]);
    } else if (!strcmp(argv[i], "-oalis")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No file basename following -oalis" << std::endl;
        exit(4);
      } else
        strcpy(par.alisbasename, argv[i]);
    }
    else if (!strcmp(argv[i], "-Ofas")) {
      par.append = 0;
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      } else
        strcpy(par.pairwisealisfile, argv[i]);
    } else if (!strcmp(argv[i], "-Oa2m")) {
      par.append = 0;
      par.outformat = 2;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      } else
        strcpy(par.pairwisealisfile, argv[i]);
    } else if (!strcmp(argv[i], "-Oa3m")) {
      par.append = 0;
      par.outformat = 3;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      } else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-scores")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No file following -scores" << std::endl;
        exit(4);
      } else {
        strcpy(par.scorefile, argv[i]);
      }
    }
    else if (!strcmp(argv[i], "-blasttab")) {
        if (++i >= argc || argv[i][0] == '-') {
            help(par);
            HH_LOG(ERROR) << "No file following -blasttab" << std::endl;
            exit(4);
        } else {
            strcpy(par.m8file, argv[i]);
        }
    }
    else if (!strcmp(argv[i], "-atab")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No file following -atab" << std::endl;
        exit(4);
      } else {
        strcpy(par.alitabfile, argv[i]);
      }
    }
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
    //no help required
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
        HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << std::endl;
    } else if (!strcmp(argv[i], "-M") && (i < argc - 1)) {
      if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m"))
        par.M = 1;
      else if (!strcmp(argv[i], "first"))
        par.M = 3;
      else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
        par.Mgaps = atoi(argv[i]);
        par.M = 2;
      } else
        HH_LOG(WARNING) << "Ignoring unknown argument: -M " << argv[i] << std::endl;
    }
    else if (!strcmp(argv[i], "-p") && (i < argc - 1))
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
    else if (!strncmp(argv[i], "-hide_cons", 7))
      par.showcons = 0;
    else if (!strncmp(argv[i], "-hide_pred", 7))
      par.showpred = 0;
    else if (!strncmp(argv[i], "-hide_dssp", 7))
      par.showdssp = 0;
    else if (!strncmp(argv[i], "-show_ssconf", 7))
      par.showconf = 1;
    else if (!strncmp(argv[i], "-mark", 7))
      par.mark = 1;
    else if (!strncmp(argv[i], "-add_cons", 5))
      par.cons = 1;
    else if (!strcmp(argv[i], "-realign_max") && (i < argc - 1))
      par.realign_max = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-e") && (i < argc - 1))
      par.e = atof(argv[++i]);
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
      //Will be zero if "Inf" or other non-numeric values are passed to atoi
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
    else if (!strcmp(argv[i], "-pc_prefilter_nocontxt_a") && (i < argc - 1))
      par.pc_hhm_nocontext_a = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_prefilter_nocontxt_b") && (i < argc - 1))
      par.pc_hhm_nocontext_b = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_prefilter_nocontxt_c") && (i < argc - 1))
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
    //no help required
    else if (!strcmp(argv[i], "-alphaa") && (i < argc - 1))
      par.alphaa = atof(argv[++i]);
    //no help required
    else if (!strcmp(argv[i], "-alphab") && (i < argc - 1))
      par.alphab = atof(argv[++i]);
    //no help required
    else if (!strcmp(argv[i], "-alphac") && (i < argc - 1))
      par.alphac = atof(argv[++i]);
    else if (!strcmp(argv[i], "-noprefilt")) {
      par.prefilter = false;
      par.already_seen_filter = false;
    } else if (!strcmp(argv[i], "-noaddfilter")) {
      par.already_seen_filter = false;
    } else if (!strcmp(argv[i], "-min_prefilter_hits") && (i < argc - 1))
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
    else if (!strcmp(argv[i], "-realign_old_hits"))
      par.realign_old_hits = true;
    else if (!strcmp(argv[i], "-realign"))
      par.realign = 1;
    else if (!strcmp(argv[i], "-norealign"))
      par.realign = 0;
    else if (!strcmp(argv[i], "-ssm") && (i < argc - 1))
      par.ssm = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ssw") && (i < argc - 1))
      par.ssw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ssa") && (i < argc - 1))
      par.ssa = atof(argv[++i]);
    else if (!strcmp(argv[i], "-wg"))
      par.wg = 1;
    else if (!strcmp(argv[i], "-maxseq") && (i < argc - 1))
        par.maxseq = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = 2 * par.maxres;
    } else if (!strncmp(argv[i], "-glob", 5)) {
      par.loc = 0;
      if (par.mact > 0.35 && par.mact < 0.3502) {
        par.mact = 0;
      }
    } else if (!strncmp(argv[i], "-loc", 4))
      par.loc = 1;
    else if (!strncmp(argv[i], "-alt", 4) && (i < argc - 1))
      par.altali = atoi(argv[++i]);
    else if (!strncmp(argv[i], "-smin", 4) && (i < argc - 1))
      par.smin = atof(argv[++i]);
    else if (!strcmp(argv[i], "-shift") && (i < argc - 1))
      par.shift = atof(argv[++i]);
    else if ((!strcmp(argv[i], "-mact") || !strcmp(argv[i], "-mapt"))
        && (i < argc - 1))
      par.mact = atof(argv[++i]);
    else if (!strcmp(argv[i], "-sc") && (i < argc - 1))
      par.columnscore = atoi(argv[++i]);
    //no help required
    else if (!strcmp(argv[i], "-scwin") && (i < argc - 1)) {
      par.columnscore = 5;
      par.half_window_size_local_aa_bg_freqs = std::max(1, atoi(argv[++i]));
    } else if (!strncmp(argv[i], "-cpu", 4) && (i < argc - 1)) {
      par.threads = atoi(argv[++i]);

      if(par.threads > MAXBINS) {
		  HH_LOG(WARNING) << "Reduced threads to max. allowed threads " << MAXBINS << "!" << std::endl;
		  par.threads = MAXBINS;
      }
    } else if (!strcmp(argv[i], "-maxmem") && (i < argc - 1)) {
      par.maxmem = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-nocontxt"))
      par.nocontxt = 1;
    else if (!strcmp(argv[i], "-csb") && (i < argc - 1))
      par.csb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-csw") && (i < argc - 1))
      par.csw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-corr") && (i < argc - 1))
      par.corr = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ovlp") && (i < argc - 1))
      par.min_overlap = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-tags"))
      par.notags = 0;
    else if (!strcmp(argv[i], "-notags"))
      par.notags = 1;
    else if (!strcmp(argv[i], "-maxfilt") && (i < argc - 1))
      par.maxnumdb = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-interim_filter")){
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No state following -interim_filter" << std::endl;
        exit(4);
      } else {
        if(!strcmp(argv[i], "NONE")) {
          par.interim_filter = INTERIM_FILTER_NONE;
        } else if(!strcmp(argv[i], "FULL")) {
          par.interim_filter = INTERIM_FILTER_FULL;
        } else {
          help(par);
          HH_LOG(ERROR) << "No state out of NONE|FULL following -interim_filter" << std::endl;
          exit(4);
        }
      }
    }
    else {
      HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << std::endl;
    }

    HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
  }  // end of for-loop for command line input
}

void HHblits::mergeHitsToQuery(Hash<Hit>* previous_hits,
                               int& seqs_found, int& cluster_found) {

  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
  const float COV_ABS = 25;     // min. number of aligned residues
  int cov_tot = std::max(std::min((int) (COV_ABS / Qali->L * 100 + 0.5), 70), par.coverage);

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
    std::stringstream ss_tmp;
    ss_tmp << hit_cur.file << "__" << hit_cur.irep;
    if (previous_hits->Contains((char*) ss_tmp.str().c_str()))
      continue;

    // Add number of sequences in this cluster to total found
    seqs_found += SequencesInCluster(hit_cur.name);  // read number after second '|'
    cluster_found++;

    // Read a3m alignment of hit from <file>.a3m file
    // Reading in next db MSA and merging it onto Qali
    Alignment Tali(par.maxseq, par.maxres);
    hit_cur.entry->getTemplateA3M(par, pb, S, Sim, Tali);

    if (par.allseqs)  // need to keep *all* sequences in Qali_allseqs? => merge before filtering
      Qali_allseqs->MergeMasterSlave(hit_cur, Tali, hit_cur.name, par.maxcol);

    Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
                                  par.qid_db, par.qsc_db, par.Ndiff_db);

    if(par.interim_filter == INTERIM_FILTER_FULL && Tali.N_filtered + Qali->N_in >= par.maxseq) {
  	  Qali->N_filtered = Qali->Filter(par.max_seqid, S, cov_tot, par.qid, par.qsc, par.Ndiff);
      Qali->Shrink();
    }

    Qali->MergeMasterSlave(hit_cur, Tali, hit_cur.name, par.maxcol);

    if (Qali->N_in >= par.maxseq)
      break;  // Maximum number of sequences reached
  }

  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
  Qali->Compress("merged A3M file", par.cons, par.maxcol, 1, par.Mgaps);

  // Sort out the nseqdis most dissimilacd r sequences for display in the result alignments
  Qali->FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid, par.qsc, par.nseqdis);


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
                                                         R, par.ssm, S73, S33, S37);

  for (std::vector<Hit>::size_type i = 0; i != hits_to_add.size(); i++) {
    std::stringstream ss_tmp;
    ss_tmp << hits_to_add[i].file << "__" << hits_to_add[i].irep;
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
                              std::vector<HHEntry*>& hits_to_realign) {

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

    hit_vector.push_back(hitlist.ReadCurrentAddress());
    n_realignments++;
    nhits++;
  }

  int t_maxres = Lmax + 2;
  for (int i = 0; i < par.threads; i++) {
    posteriorMatrices[i]->allocateMatrix(q->L, t_maxres);
  }

  // Initialize a Null-value as a return value if not items are available anymore
  PosteriorDecoderRunner runner(posteriorMatrices, viterbiMatrices, par.threads, par.ssw, S73, S33, S37);

  HH_LOG(INFO)
      << "Realigning " << nhits
      << " HMM-HMM alignments using Maximum Accuracy algorithm" << std::endl;

    runner.executeComputation(*q, hit_vector, par, par.qsc_db, pb, S, Sim, R);


  // Delete all hitlist entries with too short alignments
  nhits = 0;
  hitlist.Reset();
  while (!hitlist.End()) {
    Hit hit_cur = hitlist.ReadNext();
    //printf("Deleting alignment of %s with length %i? irep=%i nhits=%-2i  par.B=%-3i  par.Z=%-3i par.e=%.2g par.b=%-3i  par.z=%-3i par.p=%.2g\n",hit_cur.name,hit_cur.matched_cols,hit_cur.irep,nhits,par.B,par.Z,par.e,par.b,par.z,par.p);

    if (nhits > par.realign_max && nhits >= imax(par.B, par.Z))
      break;
    if (hit_cur.Eval > par.e) {
      if (nhits >= imax(par.B, par.Z))
        continue;
      if (nhits >= imax(par.b, par.z) && hit_cur.Probab < par.p)
        continue;
      if (nhits >= imax(par.b, par.z) && hit_cur.Eval > par.E)
        continue;
    }

    if (hit_cur.matched_cols < MINCOLS_REALIGN) {
      HH_LOG(DEBUG) << "Deleting alignment of " << hit_cur.name
          << " with length " << hit_cur.matched_cols << std::endl;
      hitlist.Delete().Delete();        // delete the list record and hit object
      // // Make sure only realigned alignments get displayed! JS: Why? better unrealigned than none.
      // if (last_round)
      // if (par.B>par.Z) par.B--; else if (par.B==par.Z) {par.B--; par.Z--;} else par.Z--;
    }
    nhits++;
  }
}

void HHblits::run(FILE* query_fh, char* query_path) {
  int cluster_found = 0;
  int seqs_found = 0;

  std::set<std::string> search_counter;

  Hit hit_cur;
  Hash<Hit>* previous_hits = new Hash<Hit>(1631, hit_cur);

  Qali = new Alignment(par.maxseq, par.maxres);
  Qali_allseqs = new Alignment(par.maxseq, par.maxres);

  q = new HMM(MAXSEQDIS, par.maxres);
  HMMSimd q_vec(par.maxres);
  q_tmp = new HMM(MAXSEQDIS, par.maxres);

  // Read query input file (HHM or alignment format) without adding pseudocounts
  Qali->N_in = 0;
  char input_format;
  ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb, S,
                Sim);

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
      search_counter.insert(all_entries[i]->getName());
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  for (int round = 1; round <= par.num_rounds; round++) {
    HH_LOG(INFO) << "Iteration " << round << std::endl;

    // Save HMM without pseudocounts for prefilter query-profile
    *q_tmp = *q;
    HMM* q_rescore = new HMM(MAXSEQDIS, par.maxres);


    PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine,
                    pc_hhm_context_mode, pb, R);
    //deep copy for rescoring of q
    *q_rescore = *q;
    q_vec.MapOneHMM(q_rescore);

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
        search_counter.insert(new_entries[i]->getName());
      }

      all_entries.insert(all_entries.end(), new_entries.begin(),
                         new_entries.end());
      all_entries.insert(all_entries.end(), old_entries.begin(),
                         old_entries.end());
    }

    int max_template_length = getMaxTemplateLength(new_entries);
    if (max_template_length > par.maxres) {
      HH_LOG(WARNING)
          << "database contains sequences that exceeds maximum allowed size (maxres = "
          << par.maxres << "). Maxres can be increased with parameter -maxres."
          << std::endl;
    }
    max_template_length = std::min(max_template_length, par.maxres);
    for (int i = 0; i < par.threads; i++) {
      viterbiMatrices[i]->AllocateBacktraceMatrix(q->L, max_template_length);
    }

    hitlist.N_searched = search_counter.size();

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
                                                           Sim, R, par.ssm, S73, S33, S37);

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
                                                                  S, Sim, R, par.ssm, S73, S33, S37);
        add_hits_to_hitlist(hits_to_add, hitlist);
        // Add dbfiles_old to dbfiles_new for realign
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
      perform_realign(q_vec, input_format, new_entries);

    // Generate alignment for next iteration
    if (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile || *par.alisbasename) {
      if (new_hits > 0) {
        mergeHitsToQuery(previous_hits, seqs_found, cluster_found);
      }

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, pb, Sim, NULL, true);

      if (par.notags)
        q->NeutralizeTags(pb);

      if (*par.alisbasename) {
        Alignment* tmp = new Alignment(par.maxseq, par.maxres);
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

        std::stringstream ss_tmp;
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

    if (Qali->N_in >= par.maxseq) {
      HH_LOG(INFO)
          << "Maximun number of sequences in query alignment reached (" << par.maxseq << "). Stop searching!" << std::endl;
    }

    if (new_hits == 0 || round == par.num_rounds || q->Neff_HMM > par.neffmax || Qali->N_in >= par.maxseq) {
//      if (new_hits == 0 && round < par.num_rounds) {
//        HH_LOG(INFO) << "No new hits found in iteration " << round
//                               << " => Stop searching" << std::endl;
//      }
//
//      if (old_entries.size() > 0 && par.realign_old_hits) {
//        HH_LOG(INFO)
//            << "Recalculating previously found HMMs with Viterbi algorithm"
//            << std::endl;
//
//        ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
//        std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec,
//                                                               old_entries,
//                                                               par.qsc_db, pb,
//                                                               S, Sim, R, par.ssm, S73, S33, S37);
//
//        add_hits_to_hitlist(hits_to_add, hitlist);
//
//        if (par.realign)
//          perform_realign(q_vec, input_format, old_entries);
//      } else if (!par.realign_old_hits && previous_hits->Size() > 0) {
//        HH_LOG(INFO)
//            << "Rescoring previously found HMMs with Viterbi algorithm"
//            << std::endl;
//        RescoreWithViterbiKeepAlignment(q_vec, previous_hits);
//      }

      delete q_rescore;
      break;
    }


    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();

      std::stringstream ss_tmp;
      ss_tmp << hit_cur.file << "__" << hit_cur.irep;

      if (!par.already_seen_filter || hit_cur.Eval > par.e
          || previous_hits->Contains((char*) ss_tmp.str().c_str()))
        hit_cur.Delete();  // Delete hit object (deep delete with Hit::Delete())
      else {
        previous_hits->Add((char*) ss_tmp.str().c_str(), hit_cur);
      }

      hitlist.Delete();  // Delete list record (flat delete)
    }

    delete q_rescore;
  }

  // Warn, if HMMER files were used
  if (par.hmmer_used) {
    HH_LOG(WARNING)
        << "Using HMMER files results in a drastically reduced sensitivity (>10%%).\n"
        << " We recommend to use HHMs build by hhmake." << std::endl;
  }

  for (size_t i = 0; i < all_entries.size(); i++) {
    delete all_entries[i];
  }
  all_entries.clear();

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete();  // Delete hit object
  delete previous_hits;
}


void HHblits::run(ffindex_entry_t* entry, char* data,
      ffindex_index_t* sequence_index, char* seq,
      ffindex_index_t* header_index, char* header) {
    
    
    q=new HMM(MAXSEQDIS, par.maxres);
    
    Qali = new Alignment(par.maxseq, par.maxres);
    Qali_allseqs = new Alignment(par.maxseq, par.maxres);



    // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
    Qali->ReadCompressed(entry,data,
      sequence_index, seq,
      header_index, header,
      par.mark, par.maxcol);


    // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
    // and store marked sequences in name[k] and seq[k]
    Qali->Compress(par.infile, par.cons, par.maxcol, par.M, par.Mgaps);

    Qali->Shrink();

    // Sort out the nseqdis most dissimilar sequences for display in the output alignments
    Qali->FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid, par.qsc, par.nseqdis);


    // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
    Qali->N_filtered = Qali->Filter(par.max_seqid, S, par.coverage, par.qid, par.qsc, par.Ndiff);

    if (par.Neff >= 0.999)
    	Qali->FilterNeff(par.wg, par.mark, par.cons, par.showcons, par.max_seqid, par.coverage, par.Neff, pb, S, Sim);

    // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
    Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons,par.showcons, pb, Sim);



  int cluster_found = 0;
  int seqs_found = 0;

  std::set<std::string> search_counter;

  Hit hit_cur;
  Hash<Hit>* previous_hits = new Hash<Hit>(1631, hit_cur);



  HMMSimd q_vec(par.maxres);
  q_tmp = new HMM(MAXSEQDIS, par.maxres);

  // Read query input file (HHM or alignment format) without adding pseudocounts
  Qali->N_in = 0;
  char input_format;
  
  /*ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb, S,
                Sim);
                */
                

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
      search_counter.insert(all_entries[i]->getName());
    }
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  for (int round = 1; round <= par.num_rounds; round++) {
    HH_LOG(INFO) << "Iteration " << round << std::endl;

    // Save HMM without pseudocounts for prefilter query-profile
    *q_tmp = *q;
    HMM* q_rescore = new HMM(MAXSEQDIS, par.maxres);


    PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine,
                    pc_hhm_context_mode, pb, R);
    //deep copy for rescoring of q
    *q_rescore = *q;
    q_vec.MapOneHMM(q_rescore);

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
        search_counter.insert(new_entries[i]->getName());
      }

      all_entries.insert(all_entries.end(), new_entries.begin(),
                         new_entries.end());
      all_entries.insert(all_entries.end(), old_entries.begin(),
                         old_entries.end());
    }

    int max_template_length = getMaxTemplateLength(new_entries);
    if (max_template_length > par.maxres) {
      HH_LOG(WARNING)
          << "database contains sequences that exceeds maximum allowed size (maxres = "
          << par.maxres << "). Maxres can be increased with parameter -maxres."
          << std::endl;
    }
    max_template_length = std::min(max_template_length, par.maxres);
    for (int i = 0; i < par.threads; i++) {
      viterbiMatrices[i]->AllocateBacktraceMatrix(q->L, max_template_length);
    }

    hitlist.N_searched = search_counter.size();

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
                                                           Sim, R, par.ssm, S73, S33, S37);

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
                                                                  S, Sim, R, par.ssm, S73, S33, S37);
        add_hits_to_hitlist(hits_to_add, hitlist);
        // Add dbfiles_old to dbfiles_new for realign
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
      perform_realign(q_vec, input_format, new_entries);

    // Generate alignment for next iteration
    if (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile || *par.alisbasename) {
      if (new_hits > 0) {
        mergeHitsToQuery(previous_hits, seqs_found, cluster_found);
      }

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons,par.showcons, pb, Sim, NULL,true);

      if (par.notags)
        q->NeutralizeTags(pb);

      if (*par.alisbasename) {
        Alignment* tmp = new Alignment(par.maxseq, par.maxres);
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

        std::stringstream ss_tmp;
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

    if (Qali->N_in >= par.maxseq) {
      HH_LOG(INFO)
          << "Maximun number of sequences in query alignment reached ("
          << par.maxseq << "). Stop searching!" << std::endl;
    }

    if (new_hits == 0 || round == par.num_rounds || q->Neff_HMM > par.neffmax || Qali->N_in >= par.maxseq) {
//      if (new_hits == 0 && round < par.num_rounds) {
//        HH_LOG(INFO) << "No new hits found in iteration " << round
//                               << " => Stop searching" << std::endl;
//      }
//
//      if (old_entries.size() > 0 && par.realign_old_hits) {
//        HH_LOG(INFO)
//            << "Recalculating previously found HMMs with Viterbi algorithm"
//            << std::endl;
//
//        ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
//        std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec,
//                                                               old_entries,
//                                                               par.qsc_db, pb,
//                                                               S, Sim, R, par.ssm, S73, S33, S37);
//
//        add_hits_to_hitlist(hits_to_add, hitlist);
//
//        if (par.realign)
//          perform_realign(q_vec, input_format, old_entries);
//      } else if (!par.realign_old_hits && previous_hits->Size() > 0) {
//        HH_LOG(INFO)
//            << "Rescoring previously found HMMs with Viterbi algorithm"
//            << std::endl;
//        RescoreWithViterbiKeepAlignment(q_vec, previous_hits);
//      }

      delete q_rescore;
      break;
    }


    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_cur = hitlist.ReadNext();

      std::stringstream ss_tmp;
      ss_tmp << hit_cur.file << "__" << hit_cur.irep;

      if (!par.already_seen_filter || hit_cur.Eval > par.e
          || previous_hits->Contains((char*) ss_tmp.str().c_str()))
        hit_cur.Delete();  // Delete hit object (deep delete with Hit::Delete())
      else {
        previous_hits->Add((char*) ss_tmp.str().c_str(), hit_cur);
      }

      hitlist.Delete();  // Delete list record (flat delete)
    }

    delete q_rescore;
  }

  // Warn, if HMMER files were used
  if (par.hmmer_used) {
    HH_LOG(WARNING)
        << "Using HMMER files results in a drastically reduced sensitivity (>10%%).\n"
        << " We recommend to use HHMs build by hhmake." << std::endl;
  }

  for (size_t i = 0; i < all_entries.size(); i++) {
    delete all_entries[i];
  }
  all_entries.clear();

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete();  // Delete hit object
  delete previous_hits;
}



void HHblits::writeHHRFile(char* hhrFile) {
  if (*hhrFile) {
    hitlist.PrintHHR(q_tmp, hhrFile, par.maxdbstrlen, par.showconf,
                     par.showcons, par.showdssp, par.showpred, par.b, par.B,
                     par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E,
                     par.argc, par.argv, S, par.maxseq);
  }
}

void HHblits::writeAlisFile(char* basename) {
  if (*basename) {
    std::map<int, Alignment*>::iterator it;
    for (it = alis.begin(); it != alis.end(); it++) {
      std::stringstream ss_tmp;
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

void HHblits::writeM8(char* m8File) {
    if (*m8File) {
        hitlist.PrintM8File(q, m8File);
    }
}

void HHblits::writePairwiseAlisFile(char* pairwiseAlisFile, char outformat) {
  if (*pairwiseAlisFile) {
    hitlist.PrintAlignments(q, pairwiseAlisFile, par.showconf, par.showcons,
                            par.showdssp, par.showpred, par.p, par.aliwidth,
                            par.nseqdis, par.b, par.B, par.E, S, par.maxseq, outformat);
  }
}

void HHblits::writeAlitabFile(char* alitabFile) {
  if (*alitabFile) {
    hitlist.WriteToAlifile(q, alitabFile, par.b, par.B, par.z, par.Z, par.p,
                           par.E);
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

void HHblits::writeHHRFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintHHR(hhblits.q_tmp, out, hhblits.par.maxdbstrlen,
                           hhblits.par.showconf, hhblits.par.showcons,
                           hhblits.par.showdssp, hhblits.par.showpred,
                           hhblits.par.b, hhblits.par.B, hhblits.par.z,
                           hhblits.par.Z, hhblits.par.aliwidth,
                           hhblits.par.nseqdis, hhblits.par.p, hhblits.par.E,
                           hhblits.par.argc, hhblits.par.argv, hhblits.S,
                           hhblits.par.maxseq);
}

void HHblits::printHitList() {
  char output[] = {"stdout"};
  hitlist.PrintHitList(q, output, par.maxdbstrlen, par.z, par.Z, par.p, par.E, par.argc, par.argv);
}

void HHblits::printHHRFile() {
  char output[] = {"stdout"};
  hitlist.PrintHHR(q_tmp, output, par.maxdbstrlen,
                           par.showconf, par.showcons,
                           par.showdssp, par.showpred,
                           par.b, par.B, par.z,
                           par.Z, par.aliwidth,
                           par.nseqdis, par.p, par.E,
                           par.argc, par.argv, S, par.maxseq);
}

void HHblits::writeScoresFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintScoreFile(hhblits.q, out);
}

void HHblits::writeM8(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintM8File(hhblits.q, out);
}

void HHblits::writePairwiseAlisFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintAlignments(hhblits.q, out, hhblits.par.showconf,
                                  hhblits.par.showcons, hhblits.par.showdssp,
                                  hhblits.par.showpred, hhblits.par.p,
                                  hhblits.par.aliwidth, hhblits.par.nseqdis,
                                  hhblits.par.b, hhblits.par.B, hhblits.par.E,
                                  hhblits.S, hhblits.par.maxseq, hhblits.par.outformat);
}

void HHblits::writeAlitabFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.WriteToAlifile(hhblits.q, out, hhblits.par.b, hhblits.par.B,
                                 hhblits.par.z, hhblits.par.Z, hhblits.par.p,
                                 hhblits.par.E);
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
    hitlist.PrintMatrices(q, matricesOutputFileName, par.filter_matrices,
                          par.max_number_matrices, S);
  }
}

void HHblits::writeMatricesFile(HHblits& hhblits, std::stringstream& out) {
  hhblits.hitlist.PrintMatrices(hhblits.q, out, hhblits.par.filter_matrices,
                                hhblits.par.max_number_matrices,
                                hhblits.S);
}

