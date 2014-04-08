/*
 * HHblits.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include "HHblits.h"

#define SWAP(tmp, arg1, arg2) tmp = arg1; arg1 = arg2; arg2 = tmp;

HHblits::HHblits(Parameters& parameters) {
  Init(parameters);
}

HHblits::~HHblits() {
  Reset();

  for (int bin = 0; bin < par.threads; bin++) {
    hit[bin]->DeleteBacktraceMatrix(par.maxres);
    hit[bin]->DeleteForwardMatrix(par.maxres);
  }

  for (int bin = 0; bin < par.threads; bin++) {
    delete hit[bin];
    delete t[bin];
  }

  delete[] format;

  if (use_compressed_a3m) {
    fclose(dbca3m_data_file);
    fclose(dbca3m_index_file);
    fclose(dbuniprot_header_data_file);
    fclose(dbuniprot_header_index_file);
    fclose(dbuniprot_sequence_data_file);
    fclose(dbuniprot_sequence_index_file);
    free(dbca3m_index);
    free(dbuniprot_sequence_index);
    free(dbuniprot_header_index);
  }
  else {
    fclose(dbhhm_data_file);
    fclose(dbhhm_index_file);

    if (dba3m_index_file != NULL) {
      fclose(dba3m_data_file);
      fclose(dba3m_index_file);
    }

    free(dbhhm_index);
    free(dba3m_index);
  }

  delete cs_lib;
  DeletePseudocountsEngine(context_lib, crf, pc_hhm_context_engine, pc_hhm_context_mode, pc_prefilter_context_engine, pc_prefilter_context_mode);

  if (par.prefilter) {
    free(length);
    free(first);

    for (size_t n = 0; n < num_dbs; n++)
      delete[] dbnames[n];
    free(dbnames);

    fclose(db_data_file);
  }

  delete[] dbfiles_new;
  delete[] dbfiles_old;
}

void HHblits::Init(Parameters& parameters) {
  par = parameters;

  v = par.v;

  N_searched = 0;
  dbfiles_new = new char*[par.maxnumdb_no_prefilter + 1];
  dbfiles_old = new char*[par.maxnumdb + 1];

  SetDatabase(par.db_base);

  if (par.prefilter) {
    init_prefilter();
  }
  else {
    init_no_prefiltering();
  }

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix(par.matrix, pb, P, R, S, Sim);

  // Set secondary structure substitution matrix
  if (par.ssm)
    SetSecStrucSubstitutionMatrix(par.ssa, S73, S33);

  // Prepare pseudocounts
  if (!par.nocontxt && *par.clusterfile) {
    char ext[100];
    Extension(ext, par.clusterfile);
    InitializePseudocountsEngine(par, context_lib, crf, pc_hhm_context_engine, pc_hhm_context_mode, pc_prefilter_context_engine, pc_prefilter_context_mode);
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

  // Prepare column state lib (context size =1 )
  FILE* fin = fopen(par.cs_library, "r");
  if (!fin)
    OpenFileError(par.cs_library, __FILE__, __LINE__, __func__);
  cs_lib = new cs::ContextLibrary<cs::AA>(fin);
  fclose(fin);
  cs::TransformToLin(*cs_lib);
}

void HHblits::ProcessAllArguments(int argc, char** argv, Parameters& par) {
	  par.argv = argv;
	  par.argc = argc;

	  par.premerge = 3;
	  par.Ndiff = 1000;
	  par.prefilter = true;

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
	    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	    std::cerr << "\tinput file missing!" << std::endl;
	    exit(4);
	  }
	  if (!*par.db_base) {
	    help(par);
	    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	    std::cerr << "\tdatabase missing (see -d)\n";
	    exit(4);
	  }
	  if (par.addss == 1 && (!*par.psipred || !*par.psipred_data)) {
	    help(par);
	    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	    std::cerr << "\tmissing PSIPRED directory (see -psipred and -psipred_data).\n" << std::endl;
	    std::cerr << "\tIf you don't need the predicted secondary structure, don't use the -addss option!" << std::endl;
	    exit(4);
	  }
	  if (!par.nocontxt) {
	    if (!strcmp(par.clusterfile, "")) {
	      help(par);
		  std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	      std::cerr << "\tcontext-specific library missing (see -contxt)\n";
	      exit(4);
	    }
	    if (!strcmp(par.cs_library, "")) {
	      help(par);
          std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
          std::cerr << "\tcolumn state library (see -cslib)\n";
	      exit(4);
	    }
	  }
	  if (par.loc == 0 && par.num_rounds >= 2 && v >= 1) {
		std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	    std::cerr << "\tusing -global alignment for iterative searches is deprecated since non-homologous sequence segments can easily enter the MSA and corrupt it.\n";
	  }

	  if (par.num_rounds < 1)
	    par.num_rounds = 1;
	  else if (par.num_rounds > 8) {
	    if (v >= 1) {
	      std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	      std::cerr << "\tNumber of iterations (" << par.num_rounds << ") to large => Set to 8 iterations\n";
	    }
	    par.num_rounds = 8;
	  }

	  // Premerging can be very time-consuming on large database a3ms, such as from pdb70.
	  // Hence it is only done when iteratively searching against uniprot20 or nr20 with their much smaller MSAs:
	  if (!(par.num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile || *par.alisbasename))
	    par.premerge = 0;

	  // No outfile given? Name it basename.hhm
	  if (!*par.outfile) {     // outfile not given? Name it basename.hhm
	    RemoveExtension(par.outfile, par.infile);
	    strcat(par.outfile, ".hhr");
	    if (v >= 2)
	      cout << "Search results will be written to " << par.outfile << "\n";
	  }

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

void HHblits::SetDatabase(char* db_base) {
  // Set databases
  strcpy(dbcs_base, db_base);
  strcat(dbcs_base, "_cs219");
  strcpy(dbcs_index_filename, dbcs_base);
  strcat(dbcs_index_filename, ".ffindex");
  strcpy(dbcs_data_filename, dbcs_base);
  strcat(dbcs_data_filename, ".ffdata");

  strcpy(dbhhm_base, db_base);
  strcat(dbhhm_base, "_hhm");
  strcpy(dbhhm_index_filename, dbhhm_base);
  strcat(dbhhm_index_filename, ".ffindex");
  strcpy(dbhhm_data_filename, dbhhm_base);
  strcat(dbhhm_data_filename, ".ffdata");

  strcpy(dba3m_base, db_base);
  strcat(dba3m_base, "_a3m");
  strcpy(dba3m_index_filename, dba3m_base);
  strcat(dba3m_index_filename, ".ffindex");
  strcpy(dba3m_data_filename, dba3m_base);
  strcat(dba3m_data_filename, ".ffdata");

  strcpy(dbca3m_base, db_base);
  strcat(dbca3m_base, "_ca3m");
  strcpy(dbca3m_index_filename, dbca3m_base);
  strcat(dbca3m_index_filename, ".ffindex");
  strcpy(dbca3m_data_filename, dbca3m_base);
  strcat(dbca3m_data_filename, ".ffdata");

  strcpy(dbuniprot_header_index_filename, db_base);
  strcat(dbuniprot_header_index_filename, "_header.ffindex");
  strcpy(dbuniprot_header_data_filename, db_base);
  strcat(dbuniprot_header_data_filename, "_header.ffdata");
  strcpy(dbuniprot_sequence_index_filename, db_base);
  strcat(dbuniprot_sequence_index_filename, "_sequence.ffindex");
  strcpy(dbuniprot_sequence_data_filename, db_base);
  strcat(dbuniprot_sequence_data_filename, "_sequence.ffdata");

  if (file_exists(dbca3m_index_filename) && file_exists(dbca3m_data_filename)
      && file_exists(dbuniprot_header_index_filename)
      && file_exists(dbuniprot_header_data_filename)
      && file_exists(dbuniprot_sequence_index_filename)
      && file_exists(dbuniprot_sequence_data_filename)) {
    use_compressed_a3m = true;
  }
  else if (!(file_exists(dba3m_data_filename)
      && file_exists(dba3m_index_filename))) {
    if (par.num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile
        || *par.alisbasename) {

      std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
      std::cerr << "\tCould not open A3M database " << dba3m_data_filename << " (needed to construct result MSA)" << endl;
      exit(4);
    }
    dba3m_data_filename[0] = 0;
  }

  if (use_compressed_a3m) {
    //ca3m database
    dbca3m_data_file = fopen(dbca3m_data_filename, "r");
    dbca3m_index_file = fopen(dbca3m_index_filename, "r");

    if (dbca3m_data_file == NULL) {
      OpenFileError(dbca3m_data_filename, __FILE__, __LINE__, __func__);
    }

    if (dbca3m_index_file == NULL) {
      OpenFileError(dbca3m_index_filename, __FILE__, __LINE__, __func__);
    }

    size_t ca3m_data_size = CountLinesInFile(dbca3m_index_filename);

    dbca3m_index = ffindex_index_parse(dbca3m_index_file, ca3m_data_size);
    if (dbca3m_index == NULL) {
      cerr << "Error in " << par.argv[0] << ": could not read index file"
          << dbca3m_index_filename << ". Is the file empty or corrupted?\n";
      exit(1);
    }

    dbca3m_data = ffindex_mmap_data(dbca3m_data_file, &ca3m_data_offset);

    //uniprot sequences database
    dbuniprot_sequence_data_file = fopen(dbuniprot_sequence_data_filename, "r");
    dbuniprot_sequence_index_file = fopen(dbuniprot_sequence_index_filename,
        "r");

    if (dbuniprot_sequence_data_file == NULL) {
      OpenFileError(dbuniprot_sequence_data_filename, __FILE__, __LINE__, __func__);
    }

    if (dbuniprot_sequence_index_file == NULL) {
      OpenFileError(dbuniprot_sequence_index_filename, __FILE__, __LINE__, __func__);
    }

    size_t uniprot_sequence_data_size = CountLinesInFile(
        dbuniprot_sequence_index_filename);

    dbuniprot_sequence_index = ffindex_index_parse(
        dbuniprot_sequence_index_file, uniprot_sequence_data_size);
    if (dbuniprot_sequence_index == NULL) {
      cerr << "Error in " << par.argv[0] << ": could not read index file"
          << dbuniprot_sequence_index_filename
          << ". Is the file empty or corrupted?\n";
      exit(1);
    }

    dbuniprot_sequence_data = ffindex_mmap_data(dbuniprot_sequence_data_file,
        &uniprot_sequence_data_offset);

    //uniprot sequences database
    dbuniprot_header_data_file = fopen(dbuniprot_header_data_filename, "r");
    dbuniprot_header_index_file = fopen(dbuniprot_header_index_filename, "r");

    if (dbuniprot_header_data_file == NULL) {
      OpenFileError(dbuniprot_header_data_filename, __FILE__, __LINE__, __func__);
    }

    if (dbuniprot_header_index_file == NULL) {
      OpenFileError(dbuniprot_header_index_filename, __FILE__, __LINE__, __func__);
    }

    size_t uniprot_header_data_size = CountLinesInFile(
        dbuniprot_header_index_filename);

    dbuniprot_header_index = ffindex_index_parse(dbuniprot_header_index_file,
        uniprot_header_data_size);
    if (dbuniprot_header_index == NULL) {
      cerr << "Error in " << par.argv[0] << ": could not read index file"
          << dbuniprot_sequence_index_filename
          << ". Is the file empty or corrupted?\n";
      exit(1);
    }

    dbuniprot_header_data = ffindex_mmap_data(dbuniprot_header_data_file,
        &uniprot_header_data_offset);
  }
  else {
    // Prepare index-based databases
    dbhhm_data_file = fopen(dbhhm_data_filename, "r");
    if (!dbhhm_data_file)
      OpenFileError(dbhhm_data_filename, __FILE__, __LINE__, __func__);

    dbhhm_index_file = fopen(dbhhm_index_filename, "r");
    if (!dbhhm_index_file)
      OpenFileError(dbhhm_index_filename, __FILE__, __LINE__, __func__);

    int filesize;
    filesize = CountLinesInFile(dbhhm_index_filename);

    dbhhm_index = ffindex_index_parse(dbhhm_index_file, filesize);
    if (dbhhm_index == NULL) {
      cerr << "Error in " << par.argv[0] << ": could not read index file"
          << dbhhm_index_filename << ". Is the file empty or corrupted?\n";
      exit(1);
    }
    dbhhm_data = ffindex_mmap_data(dbhhm_data_file, &data_size);

    if (!*dba3m_data_filename) {
      dba3m_data_file = dba3m_index_file = NULL;
      dba3m_index = NULL;
    }
    else {
      dba3m_data_file = fopen(dba3m_data_filename, "r");
      if (!dba3m_data_file)
        OpenFileError(dba3m_data_filename, __FILE__, __LINE__, __func__);

      filesize = CountLinesInFile(dba3m_index_filename);

      dba3m_index_file = fopen(dba3m_index_filename, "r");
      if (!dba3m_index_file)
        OpenFileError(dba3m_index_filename, __FILE__, __LINE__, __func__);

      dba3m_index = ffindex_index_parse(dba3m_index_file, filesize);
      if (dba3m_index == NULL) {
        cerr << "Error in " << par.argv[0] << ": could not read index file"
            << dba3m_index_filename << ". Is the file empty or corrupted?\n";
        exit(1);
      }
      dba3m_data = ffindex_mmap_data(dba3m_data_file, &data_size);
    }
  }
}

void HHblits::Reset() {
  if(q) {
    delete q;
    q=NULL;
  }

  if(q_tmp) {
    delete q_tmp;
    q_tmp = NULL;
  }

  for (int idb = 0; idb < ndb_new; idb++)
    delete[] dbfiles_new[idb];

  for (int idb = 0; idb < ndb_old; idb++)
    delete[] dbfiles_old[idb];

  ndb_old = 0;
  ndb_new = 0;

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
      " -d <name>      database name (e.g. uniprot20_29Feb2012) (default=%s)          \n",
      par.db_base);
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
      else
        strcpy(par.db_base, argv[i]);
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
      par.half_window_size_local_aa_bg_freqs = imax(1, atoi(argv[++i]));
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


void HHblits::ReadQueryA3MFile() {
  // Open query a3m MSA
}

void HHblits::getTemplateA3M(char* entry_name, long& ftellpos,
    Alignment& tali) {
  if (use_compressed_a3m) {
    ffindex_entry_t* entry = ffindex_get_entry_by_name(dbca3m_index,
        entry_name);

    if (entry == NULL) {
      std::cerr << "Could not fetch entry for a3m " << entry_name << "!"
          << std::endl;
      exit(4);
    }

    char* data = ffindex_get_data_by_entry(dbca3m_data, entry);

    if (data == NULL) {
      std::cerr << "Could not fetch data for a3m " << entry_name << "!"
          << std::endl;
      exit(4);
    }

    ftellpos = entry->offset;
    tali.ReadCompressed(entry, data, dbuniprot_sequence_index,
        dbuniprot_sequence_data, dbuniprot_header_index, dbuniprot_header_data, par.mark, par.maxcol);
  }
  else {
    FILE* dbf = ffindex_fopen_by_name(dba3m_data, dba3m_index, entry_name);

    if (dbf == NULL) {
      cerr << endl << "Error: opening A3M " << entry_name << std::endl;
      if (dba3m_index_file == NULL) {
        cerr << endl << "Error: A3M database missing" << std::endl;
      }
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

void HHblits::getTemplateHMMFromA3M(char* entry_name, char use_global_weights,
    long& ftellpos, int& format, HMM* t) {
  Alignment tali;
  getTemplateA3M(entry_name, ftellpos, tali);
  tali.N_filtered = tali.Filter(par.max_seqid_db, S, par.coverage_db, par.qid_db,
      par.qsc_db, par.Ndiff_db);
  t->name[0] = t->longname[0] = t->fam[0] = '\0';
  tali.FrequenciesAndTransitions(t, use_global_weights, par.mark, par.cons, par.showcons, par.maxres, pb, Sim);

  format = 0;
}

void HHblits::getTemplateHMM(char* entry_name, char use_global_weights,
    long& ftellpos, int& format, HMM* t) {
  if (!use_compressed_a3m) {
    FILE* dbf = ffindex_fopen_by_name(dbhhm_data, dbhhm_index, entry_name);

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
        t->ReadHMMer3(dbf, par.showcons, pb, entry_name);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HMMER", 5))      // read HMMER format
          {
        format = 1;
        t->ReadHMMer(dbf, par.showcons, pb, entry_name);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HH", 2))    // read HHM format
      {
        char path[NAMELEN];
        Pathname(path, entry_name);

        format = 0;
        t->Read(dbf, par.maxcol, par.nseqdis, pb, path);
      }
      //    else if (!strncmp(line, "NAME", 4)) // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
      //    {
      //      char path[NAMELEN];
      //      Pathname(path, dbfiles[idb]);
      //
      //      fseek(dbf, hit[bin]->ftellpos, SEEK_SET); // rewind to beginning of line
      //      format = 0;
      //      t->Read(dbf, path);
      //    }
      else {
    	std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    	std::cerr << "\tunrecognized HMM file format in \'" << entry_name << "\'. \n";
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

  getTemplateHMMFromA3M(entry_name, use_global_weights, ftellpos, format, t);
}

void HHblits::DoViterbiSearch(char *dbfiles[], int ndb, Hash<Hit>* previous_hits, bool alignByWorker) {
  // Search databases
  for (int bin = 0; bin < par.threads; bin++) {
    hit[bin]->realign_around_viterbi = false;
  }

  // For all the databases comming through prefilter
  #pragma omp parallel for schedule(dynamic, 5)
  for (int idb = 0; idb < ndb; idb++) {
    // Allocate free bin (no need to lock, since slave processes cannot change FREE to other status)
    int bin = omp_get_thread_num();

    #pragma omp critical
    {
      hit[bin]->index = N_searched;     // give hit a unique index for HMM
      ++N_searched;
    }

    getTemplateHMM(dbfiles[idb], 1, hit[bin]->ftellpos, format[bin], t[bin]);

    if (v >= 4)
      printf("Aligning with %s\n", t[bin]->name);

    hit[bin]->dbfile = new char[strlen(dbfiles[idb]) + 1];
    strcpy(hit[bin]->dbfile, dbfiles[idb]); // record db file name from which next HMM is read

    if (alignByWorker)
      AlignByWorker(par, hit[bin], t[bin], q, format[bin], pb, R, S73, S33, hitlist);
    else
      PerformViterbiByWorker(par, hit[bin], t[bin], q, format[bin], pb, R, S73, S33, hitlist, previous_hits);
  }
}

void HHblits::ViterbiSearch(char *dbfiles[], int ndb, Hash<Hit>* previous_hits, int db_size) {
  // Initialize
  v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2;

  hitlist.N_searched = db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  DoViterbiSearch(dbfiles, ndb, previous_hits);

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
    hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa, par.alphab, par.alphac, par.prefilter_evalue_thresh);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant of ViterbiSearch() function for rescoring previously found HMMs with Viterbi algorithm.
// Perform Viterbi search on each hit object in global hash previous_hits, but keep old alignment
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HHblits::RescoreWithViterbiKeepAlignment(int db_size, Hash<Hit>* previous_hits) {
  // Initialize
  v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2;

  char *dbfiles[db_size + 1];
  int ndb = 0;

  // Get dbfiles of previous hits
  previous_hits->Reset();
  while (!previous_hits->End()) {
    Hit hit_cur = previous_hits->ReadNext();
    if (hit_cur.irep == 1) {
      dbfiles[ndb] = new char[strlen(hit_cur.dbfile) + 1];
      strcpy(dbfiles[ndb], hit_cur.dbfile);
      ++ndb;
    }
  }

  hitlist.N_searched = db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  DoViterbiSearch(dbfiles, ndb, previous_hits, false);

  if (v1 >= 2)
    cout << "\n";
  v = v1;

  for (int n = 0; n < ndb; ++n)
    delete[] (dbfiles[n]);

  // Sort list according to sortscore
  if (v >= 3)
    printf("Sorting hit list ...\n");
  hitlist.SortList();

  hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);  // Use NN prediction of lamda and mu

  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa, par.alphab, par.alphac, par.prefilter_evalue_thresh);

}

void HHblits::perform_realign(char *dbfiles[], int ndb, Hash<char>* premerged_hits) {
  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
  int nhits = 0;
  int N_aligned = 0;

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
    if (nhits >= par.realign_max && nhits >= imax(par.B, par.Z))
      break;
    if (hit_cur.Eval > par.e) {
      if (nhits >= imax(par.B, par.Z))
        continue;
      if (nhits >= imax(par.b, par.z) && hit_cur.Probab < par.p)
        continue;
      if (nhits >= imax(par.b, par.z) && hit_cur.Eval > par.E)
        continue;
    }

    if (hit_cur.L > Lmax)
      Lmax = hit_cur.L;
    if (hit_cur.L > Lmaxmem) {
      nhits++;
      continue;
    }

    //fprintf(stderr,"hit.name=%-15.15s  hit.index=%-5i hit.ftellpos=%-8i  hit.dbfile=%s\n",hit_cur.name,hit_cur.index,(unsigned int)hit_cur.ftellpos,hit_cur.dbfile);

    if (nhits >= par.premerge) // realign the first premerge hits consecutively to query profile
        {
      if (hit_cur.irep == 1) {
        // For each template (therefore irep==1), store template index and position on disk in a list
        Realign_hitpos realign_hitpos;
        realign_hitpos.ftellpos = hit_cur.ftellpos; // stores position on disk of template for current hit
        realign_hitpos.index = hit_cur.index; // stores index of template of current hit
        realign_hitpos.entry = NULL;  //TODO
        if (!phash_plist_realignhitpos->Contains(hit_cur.dbfile)) {
          List<Realign_hitpos>* newlist = new List<Realign_hitpos>;
          phash_plist_realignhitpos->Add(hit_cur.dbfile, newlist);
        }
        // Add template index and ftellpos to list which belongs to key dbfile in hash
        phash_plist_realignhitpos->Show(hit_cur.dbfile)->Push(realign_hitpos);
      }
      if (!array_plist_phits[hit_cur.index]) // pointer at index is still NULL
      {
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
        "Realigning %i HMM-HMM alignments using Maximum Accuracy algorithm\n",
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
  // Align the first par.premerge templates
  if (par.premerge > 0) {
    if (v >= 2)
      printf("Merging %i best hits to query alignment ...\n", par.premerge);

    int bin = 0;
    nhits = 0;
    hitlist.Reset();

    while (!hitlist.End() && nhits < par.premerge) {
      Hit hit_cur = hitlist.ReadNext();
      if (hit_cur.Eval > par.e) // JS: removed bug on 13 Feb 13 due to which premerged hits with E-value > par.e were not realigned
          {
        if (nhits >= imax(par.B, par.Z))
          break;
        if (nhits >= imax(par.b, par.z) && hit_cur.Probab < par.p)
          break;
        if (nhits >= imax(par.b, par.z) && hit_cur.Eval > par.E)
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

      getTemplateHMM(hit_cur.dbfile, par.wg, hit[bin]->ftellpos, format[bin],
          t[bin]);

      if (v >= 2)
        fprintf(stderr, "Realigning with %s ***** \n", t[bin]->name);
      ///////////////////////////////////////////////////

      N_aligned++;
      if (v1 >= 2 && !(N_aligned % 10)) {
        cout << ".";
        if (!(N_aligned % 500))
          printf(" %-4i HMMs aligned\n", N_aligned);
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
      hit[bin]->Forward(q, t[bin], par.ssm, par.min_overlap, par.loc, par.shift, par.ssw, par.exclstr, S73, S33);
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

      // Replace original hit in hitlist with realigned hit
      //hitlist.ReadCurrent().Delete();
      hitlist.Delete().Delete();  // delete the list record and hit object
      hitlist.Insert(*hit[bin]);

      // merge only when hit length > MINCOLS_REALIGN (don't merge 1 column matches)
      if (hit[bin]->matched_cols < MINCOLS_REALIGN)
        continue;

      // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
      // Reading in next db MSA and merging it onto Qali
      Alignment Tali;
      long ftellpos;
      getTemplateA3M(hit[bin]->dbfile, ftellpos, Tali);

      if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
        Qali_allseqs.MergeMasterSlave(*hit[bin], Tali, hit[bin]->dbfile, par.maxcol);

      Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
          par.qid_db, par.qsc_db, par.Ndiff_db);

      Qali.MergeMasterSlave(*hit[bin], Tali, hit[bin]->dbfile, par.maxcol);

      // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
      Qali.Compress("merged A3M file", par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);

      // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
      Qali.N_filtered = Qali.Filter(par.max_seqid, S, par.coverage, par.qid,
          par.qsc, par.Ndiff);

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, par.maxres, pb, Sim);

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
  #pragma omp parallel for schedule(dynamic, 5)
  for (int idb = 0; idb < ndb; idb++) {
    // Can we skip dbfiles[idb] because it contains no template to be realigned?
    if (!phash_plist_realignhitpos->Contains(dbfiles[idb]))
      continue;

    // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
    // This list is now sorted by ftellpos in ascending order to access one template after the other efficiently
    phash_plist_realignhitpos->Show(dbfiles[idb])->SortList();

    ///////////////////////////////////////////////////////////////////////////////////////
    // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
    phash_plist_realignhitpos->Show(dbfiles[idb])->Reset();
    while (!phash_plist_realignhitpos->Show(dbfiles[idb])->End()) {
      // Submit jobs until no bin is free anymore
      while (!phash_plist_realignhitpos->Show(dbfiles[idb])->End()) {
        // Allocate free bin
        int bin = omp_get_thread_num();

        // Forward stream position to start of next database HMM to be realigned
        Realign_hitpos hitpos_curr = phash_plist_realignhitpos->Show(
            dbfiles[idb])->ReadNext();
        hit[bin]->index = hitpos_curr.index; // give hit[bin] a unique index for HMM

        // Give hit[bin] the pointer to the list of pointers to hitlist elements of same template (for realignment)
        hit[bin]->plist_phits = array_plist_phits[hitpos_curr.index];

        hit[bin]->dbfile = new char[strlen(dbfiles[idb]) + 1];
        strcpy(hit[bin]->dbfile, dbfiles[idb]); // record db file name from which next HMM is read
        getTemplateHMM(dbfiles[idb], par.wg, hit[bin]->ftellpos, format[bin],
            t[bin]);

        if (v >= 2)
          fprintf(stderr, "Realigning with %s\n", t[bin]->name);

        #pragma omp critical
        {
          N_aligned++;
        }

        RealignByWorker(par, hit[bin], q, t[bin], format[bin], pb, R, S73, S33);
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
  int cov_tot = imax(imin((int) (COV_ABS / Qali.L * 100 + 0.5), 70),
      par.coverage);

  Alignment qali;
  qali = Qali;
  HMM* q = new HMM();

  HitList realigned_viterbi_hitlist;

  for (size_t qsc_index = 0; qsc_index < nqsc; qsc_index++) {
    float actual_qsc = qsc[qsc_index];

    qali.Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
    qali.N_filtered = qali.Filter(par.max_seqid, S, cov_tot, par.qid, actual_qsc, par.Ndiff);
    qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, par.maxres, pb, Sim, NULL, false);
    PrepareQueryHMM(par, inputformat, q, pc_hhm_context_engine, pc_hhm_context_mode, pb, R);

    hitlist.Reset();
    while (!hitlist.End()) {
      Hit hit_ref = hitlist.ReadNext();

      HMM* t = new HMM();

      int format;
      long ftellpos;
      getTemplateHMM(hit_ref.dbfile, 1, ftellpos, format, t);

      PrepareTemplateHMM(par, q, t, format, pb, R);

      Hit hit;
      hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
      hit.self = 0;
      hit.realign_around_viterbi = false;

      hit.dbfile = new char[strlen(hit_ref.dbfile) + 1];
      strcpy(hit.dbfile, hit_ref.dbfile);

      for (int irep = 1; irep <= par.altali; irep++) {
        hit.irep = irep;
        hit.Viterbi(q, t, par.loc, par.ssm, par.maxres, par.min_overlap, par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);

        if (hit.irep > 1 && hit.score <= SMIN)
          break;

        hit.Backtrace(q, t, par.corr, par.ssw, S73, S33);
        realigned_viterbi_hitlist.Push(hit);
      }

      hit.DeleteBacktraceMatrix(q->L + 2);

      delete t;
    }

    realigned_viterbi_hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);
    realigned_viterbi_hitlist.CalculateHHblitsEvalues(q, par.dbsize, par.alphaa, par.alphab, par.alphac, par.prefilter_evalue_thresh);

    q->Log2LinTransitionProbs(1.0);

    realigned_viterbi_hitlist.Reset();
    while (!realigned_viterbi_hitlist.End()) {
      Hit hit_ref = realigned_viterbi_hitlist.ReadNext();

      HMM* t = new HMM();

      int format;
      long ftellpos;
      getTemplateHMM(hit_ref.dbfile, par.wg, ftellpos, format, t);

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
      hit.Forward(q, t, par.ssm, par.min_overlap, par.loc, par.shift, par.ssw, par.exclstr, S73, S33);
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
    char inputformat, float* qsc, size_t nqsc,
    HitList& reducedFinalHitList) {
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

int HHblits::swStripedByte(unsigned char *querySeq, int queryLength,
    unsigned char *dbSeq, int dbLength, unsigned short gapOpen,
    unsigned short gapExtend, __m128i *pvHLoad, __m128i *pvHStore, __m128i *pvE,
    unsigned short bias) {
  int i, j;
  int score;

  int cmp;
  int iter = (queryLength + 15) / 16;

  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore;
  __m128i vBias;
  __m128i vGapOpen;
  __m128i vGapExtend;

  __m128i vTemp;
  __m128i vZero;

  __m128i *pvScore;

  __m128i *pvQueryProf = (__m128i *) querySeq;

  /* Load the bias to all elements of a constant */
  vBias = _mm_set1_epi8(bias);

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi8(gapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi8(gapExtend);

  vMaxScore = _mm_setzero_si128();
  vZero = _mm_setzero_si128();

  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i) {
    _mm_store_si128(pvE + i, vMaxScore);
    _mm_store_si128(pvHStore + i, vMaxScore);
  }

  for (i = 0; i < dbLength; ++i) {
    /* fetch first data asap. */
    pvScore = pvQueryProf + dbSeq[i] * iter;

    /* zero out F. */
    vF = _mm_setzero_si128();

    /* load the next h value */
    vH = _mm_load_si128(pvHStore + iter - 1);
    vH = _mm_slli_si128(vH, 1);

    pv = pvHLoad;
    pvHLoad = pvHStore;
    pvHStore = pv;

    for (j = 0; j < iter; ++j) {
      /* load values of vF and vH from previous row (one unit up) */
      vE = _mm_load_si128(pvE + j);

      /* add score to vH */
      vH = _mm_adds_epu8(vH, *(pvScore++));
      vH = _mm_subs_epu8(vH, vBias);

      /* Update highest score encountered this far */
      vMaxScore = _mm_max_epu8(vMaxScore, vH);

      /* get max from vH, vE and vF */
      vH = _mm_max_epu8(vH, vE);
      vH = _mm_max_epu8(vH, vF);

      /* save vH values */
      _mm_store_si128(pvHStore + j, vH);

      /* update vE value */
      vH = _mm_subs_epu8(vH, vGapOpen);
      vE = _mm_subs_epu8(vE, vGapExtend);
      vE = _mm_max_epu8(vE, vH);

      /* update vF value */
      vF = _mm_subs_epu8(vF, vGapExtend);
      vF = _mm_max_epu8(vF, vH);

      /* save vE values */
      _mm_store_si128(pvE + j, vE);

      /* load the next h value */
      vH = _mm_load_si128(pvHLoad + j);
    }

    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm_load_si128(pvHStore);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm_slli_si128(vF, 1);
    vTemp = _mm_subs_epu8(vH, vGapOpen);
    vTemp = _mm_subs_epu8(vF, vTemp);
    vTemp = _mm_cmpeq_epi8(vTemp, vZero);
    cmp = _mm_movemask_epi8(vTemp);

    while (cmp != 0xffff) {
      vE = _mm_load_si128(pvE + j);

      vH = _mm_max_epu8(vH, vF);

      /* save vH values */
      _mm_store_si128(pvHStore + j, vH);

      /*  update vE incase the new vH value would change it */
      vH = _mm_subs_epu8(vH, vGapOpen);
      vE = _mm_max_epu8(vE, vH);
      _mm_store_si128(pvE + j, vE);

      /* update vF value */
      vF = _mm_subs_epu8(vF, vGapExtend);

      ++j;
      if (j >= iter) {
        j = 0;
        vF = _mm_slli_si128(vF, 1);
      }

      vH = _mm_load_si128(pvHStore + j);

      vTemp = _mm_subs_epu8(vH, vGapOpen);
      vTemp = _mm_subs_epu8(vF, vTemp);
      vTemp = _mm_cmpeq_epi8(vTemp, vZero);
      cmp = _mm_movemask_epi8(vTemp);
    }
  }

  /* find largest score in the vMaxScore vector */
  vTemp = _mm_srli_si128(vMaxScore, 8);
  vMaxScore = _mm_max_epu8(vMaxScore, vTemp);
  vTemp = _mm_srli_si128(vMaxScore, 4);
  vMaxScore = _mm_max_epu8(vMaxScore, vTemp);
  vTemp = _mm_srli_si128(vMaxScore, 2);
  vMaxScore = _mm_max_epu8(vMaxScore, vTemp);
  vTemp = _mm_srli_si128(vMaxScore, 1);
  vMaxScore = _mm_max_epu8(vMaxScore, vTemp);

  /* store in temporary variable */
  score = _mm_extract_epi16(vMaxScore, 0);
  score = score & 0x00ff;

  /* return largest score */
  return score;
}

// d = i-j+LT-1 is index of diagonal
int HHblits::ungapped_sse_score(const unsigned char* query_profile,
    const int query_length, const unsigned char* db_sequence,
    const int dbseq_length, const unsigned char score_offset, __m128i* workspace)
{
      int i; // position in query bands (0,..,W-1)
    int j;// position in db sequence (0,..,dbseq_length-1)
    int W = (query_length + 15) / 16;// width of bands in query and score matrix = hochgerundetes LQ/16

    __m128i *p;
    __m128i S;// 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15)
    __m128i Smax = _mm_setzero_si128();
    __m128i Soffset;// all scores in query profile are shifted up by Soffset to obtain pos values
    __m128i *s_prev, *s_curr;// pointers to Score(i-1,j-1) and Score(i,j), resp.
    __m128i *qji;// query profile score in row j (for residue x_j)
    __m128i *s_prev_it, *s_curr_it;
    __m128i *query_profile_it = (__m128i *) query_profile;
    __m128i Zero = _mm_setzero_si128();

    // Load the score offset to all 16 unsigned byte elements of Soffset
    Soffset = _mm_set1_epi8(score_offset);

    // Initialize  workspace to zero
    for (i=0, p=workspace; i < 2*W; ++i)
    _mm_store_si128(p++, Zero);

    s_curr = workspace;
    s_prev = workspace + W;

    for (j=0; j<dbseq_length; ++j)// loop over db sequence positions
    {
      // Get address of query scores for row j
      qji = query_profile_it + db_sequence[j]*W;

      // Load the next S value
      S = _mm_load_si128(s_curr + W - 1);
      S = _mm_slli_si128(S, 1);

      // Swap s_prev and s_curr, smax_prev and smax_curr
      SWAP(p,s_prev,s_curr);

      s_curr_it = s_curr;
      s_prev_it = s_prev;

      for (i=0; i<W; ++i)// loop over query band positions
      {
        // Saturated addition and subtraction to score S(i,j)
        S = _mm_adds_epu8(S, *(qji++));// S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset)
        S = _mm_subs_epu8(S, Soffset);// S(i,j) = max(0, S(i,j) - Soffset)
        _mm_store_si128(s_curr_it++, S);// store S to s_curr[i]
        Smax = _mm_max_epu8(Smax, S);// Smax(i,j) = max(Smax(i,j), S(i,j))

        // Load the next S and Smax values
        S = _mm_load_si128(s_prev_it++);
      }
    }

    /* find largest score in the Smax vector */
    S = _mm_srli_si128 (Smax, 8);
    Smax = _mm_max_epu8 (Smax, S);
    S = _mm_srli_si128 (Smax, 4);
    Smax = _mm_max_epu8 (Smax, S);
    S = _mm_srli_si128 (Smax, 2);
    Smax = _mm_max_epu8 (Smax, S);
    S = _mm_srli_si128 (Smax, 1);
    Smax = _mm_max_epu8 (Smax, S);

    /* store in temporary variable */
    int score = _mm_extract_epi16 (Smax, 0);
    score = score & 0x00ff;

    /* return largest score */
    return score;
  }

///////////////////////////////////////////////////////////////////////////////////////////////////
// Pull out all names from prefilter db file and copy into dbfiles_new for full HMM-HMM comparison
///////////////////////////////////////////////////////////////////////////////////////////////////
void HHblits::init_no_prefiltering() {
  FILE* db_index_file = fopen(dbcs_index_filename, "r");
  if (db_index_file == NULL)
    OpenFileError(dbcs_index_filename, __FILE__, __LINE__, __func__);
  ffindex_index_t* db_index = ffindex_index_parse(db_index_file, 0);
  fclose(db_index_file);

  num_dbs = db_index->n_entries;
  par.dbsize = db_index->n_entries;

  if (par.dbsize > par.maxnumdb_no_prefilter) {
	std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
	std::cerr << "\tWithout prefiltering, the max. number of database HHMs is "
        << par.maxnumdb_no_prefilter << " (actual: " << par.dbsize << ")\n";
	std::cerr << "\tYou can increase the allowed maximum using the -maxfilt <max> option.\n";
    exit(4);
  }

  for (size_t n = 0; n < num_dbs; n++) {
    ffindex_entry_t* entry = ffindex_get_entry_by_index(db_index, n);
    dbfiles_new[n] = new char[strlen(entry->name) + 1];
    strcpy(dbfiles_new[n], entry->name);
  }
  ndb_new = num_dbs;

  if (v >= 2)
    cout << "Searching " << ndb_new << " database HHMs without prefiltering"
        << endl;
}

void HHblits::checkCSFormat(size_t nr_checks) {
  for (size_t n = 0; n < std::min(nr_checks, num_dbs); n++) {
    if (first[n][0] == '>') {
      nr_checks--;
    }
  }

  if (nr_checks == 0) {
	std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    std::cerr << "\tYour cs database is in an old format!" << std::endl;
    std::cerr << "\tThis format is no longer supportet!" << std::endl;
    std::cerr << "\tCorrespond to the user manual!" << std::endl;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////
// Reading in column state sequences for prefiltering
//////////////////////////////////////////////////////////////
void HHblits::init_prefilter() {
  // Map data file into memory
  db_data_file = fopen(dbcs_data_filename, "rb");
  if (db_data_file == NULL)
    OpenFileError(dbcs_data_filename, __FILE__, __LINE__, __func__);
  size_t db_data_size;
  db_data = (unsigned char*) ffindex_mmap_data(db_data_file, &db_data_size);

  // Read index
  FILE* db_index_file = fopen(dbcs_index_filename, "r");
  if (db_index_file == NULL)
    OpenFileError(dbcs_index_filename, __FILE__, __LINE__, __func__);
  ffindex_index_t* db_index = ffindex_index_parse(db_index_file, 0);
  fclose(db_index_file);

  // Set up variables for prefiltering
  num_dbs = db_index->n_entries;
  par.dbsize = db_index->n_entries;
  first = (unsigned char**) memalign(16, num_dbs * sizeof(unsigned char*));
  length = (int*) memalign(16, num_dbs * sizeof(int));
  dbnames = (char**) memalign(16, num_dbs * sizeof(char*));
  for (size_t n = 0; n < num_dbs; n++) {
    ffindex_entry_t* entry = ffindex_get_entry_by_index(db_index, n);
    first[n] = (unsigned char*) ffindex_get_data_by_entry((char*) db_data,
        entry);
    length[n] = entry->length - 1;
    dbnames[n] = new char[strlen(entry->name) + 1];
    strcpy(dbnames[n], entry->name);
  }

  free(db_index);

  //check if cs219 format is new binary format
  checkCSFormat(5);

  if (v >= 2) {
    printf("Searching %zu column state sequences.\n", num_dbs);
  }
}

////////////////////////////////////////////////////////////////////////
// Prepare query profile for prefitering
////////////////////////////////////////////////////////////////////////
void HHblits::stripe_query_profile() {
  int LQ = q_tmp->L;
  float** query_profile = NULL;
  int a, h, i, j, k;

  // Add Pseudocounts
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

  // Build query profile with 219 column states
  query_profile = new float*[LQ + 1];
  for (i = 0; i < LQ + 1; ++i)
    query_profile[i] = (float*) memalign(16, NUMCOLSTATES * sizeof(float),
        "the query profile during prefiltering");

  const cs::ContextLibrary<cs::AA>& lib = *cs_lib;

  // log S(i,k) = log ( SUM_a p(i,a) * p(k,a) / f(a) )   k: column state, i: pos in ali, a: amino acid
  for (i = 0; i < LQ; ++i)
    for (k = 0; k < NUMCOLSTATES; ++k) {
      float sum = 0;
      for (a = 0; a < 20; ++a)
        sum += (q_tmp->p[i][a] * lib[k].probs[0][a]) / q_tmp->pav[a];
      query_profile[i + 1][k] = sum;
    }

  /////////////////////////////////////////
  // Stripe query profile with chars
  qc = (unsigned char*) memalign(16,
      (NUMCOLSTATES + 1) * (LQ + 15) * sizeof(unsigned char),
      "the striped query profile during prefiltering"); // query profile (states + 1 because of ANY char)
  W = (LQ + 15) / 16;   // band width = hochgerundetes LQ/16

  for (a = 0; a < NUMCOLSTATES; ++a) {
    h = a * W * 16;
    for (i = 0; i < W; ++i) {
      j = i;
      for (k = 0; k < 16; ++k) {
        if (j >= LQ)
          qc[h] = (unsigned char) par.prefilter_score_offset;
        else {
          float dummy = flog2(query_profile[j + 1][a])
              * par.prefilter_bit_factor + par.prefilter_score_offset + 0.5;
          // if (dummy>255.0) qc[h] = 255;
          // else if (dummy<0) qc[h] = 0;
          // else qc[h] = (unsigned char) dummy;  // 1/3 bits & make scores >=0 everywhere
          qc[h] = (unsigned char) fmax(0.0, fmin(255.0, dummy));
        }
        ++h;
        j += W;
      }
    }
  }

  // Add extra ANY-state (220'th state)
  h = NUMCOLSTATES * W * 16;
  for (i = 0; i < W; ++i) {
    j = i;
    for (k = 0; k < 16; ++k) {
      if (j >= LQ)
        qc[h] = (unsigned char) par.prefilter_score_offset;
      else
        qc[h] = (unsigned char) (par.prefilter_score_offset - 1);
      h++;
      j += W;
    }
  }

  //////////////////////////////////////////////+
  // Stripe query profile with shorts
  qw = (unsigned short*) memalign(16,
      (NUMCOLSTATES + 1) * (LQ + 7) * sizeof(unsigned short),
      "the striped 2B query profile during prefiltering"); // query profile (states + 1 because of ANY char)
  Ww = (LQ + 7) / 8;

  /////////////////////////////////////////
  // Stripe query profile
  for (a = 0; a < NUMCOLSTATES; ++a) {
    h = a * Ww * 8;
    for (i = 0; i < Ww; ++i) {
      j = i;
      for (k = 0; k < 8; ++k) {
        if (j >= LQ)
          qw[h] = 0;
        else {
          float dummy = flog2(query_profile[j + 1][a])
              * par.prefilter_bit_factor;
          qw[h] = (unsigned short) dummy; // 1/3 bits & make scores >=0 everywhere
        }
        ++h;
        j += Ww;
      }
    }
  }

  // Add extra ANY-state
  h = NUMCOLSTATES * Ww * 8;
  for (i = 0; i < Ww; ++i) {
    j = i;
    for (k = 0; k < 8; ++k) {
      if (j >= LQ)
        qw[h] = 0;
      else
        qw[h] = (unsigned short) -1;
      h++;
      j += W;
    }
  }

  free(qw);

  for (i = 0; i < LQ + 1; ++i)
    free(query_profile[i]);
  delete[] query_profile;
}
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Main prefilter function
////////////////////////////////////////////////////////////////////////
void HHblits::prefilter_db(Hash<Hit>* previous_hits) {
  Hash<char>* doubled = new Hash<char>;
  doubled->New(16381, 0);
  for (int idb = 0; idb < ndb_new; idb++)
    delete[] (dbfiles_new[idb]);
  for (int idb = 0; idb < ndb_old; idb++)
    delete[] (dbfiles_old[idb]);
  ndb_new = ndb_old = 0;

  // Prefilter with SW evalue preprefilter backtrace

  stripe_query_profile();

//  int* prefiltered_hits = new int[par.dbsize+1];
  int* backtrace_hits = new int[par.maxnumdb + 1];

  __m128i ** workspace = new __m128i *[par.threads];

  int score;
  double evalue;
  vector<pair<double, int> > first_prefilter;

  vector<pair<double, int> > hits;

  int thread_id = 0;
  int count_dbs = 0;
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;
  int LQ = q_tmp->L;
  const float log_qlen = flog2(LQ);
  const double factor = (double) par.dbsize * LQ;

  for (int i = 0; i < par.threads; i++)
    workspace[i] = (__m128i *) memalign(16, 3 * (LQ + 15) * sizeof(char),
        "the dynamic programming workspace during prefiltering");

#pragma omp parallel for schedule(static) private(score, thread_id)
  // Loop over all database sequences
  for (size_t n = 0; n < num_dbs; n++) {
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif

    // Perform search step
    score = ungapped_sse_score(qc, LQ, first[n], length[n],
        par.prefilter_score_offset, workspace[thread_id]);

    score = score
        - (int) (par.prefilter_bit_factor * (log_qlen + flog2(length[n])));

#pragma omp critical
    first_prefilter.push_back(pair<double, int>(score, n));

    if (v >= 2 && !(n % 100000)) {
      cout << ".";
      cout.flush();
    }
  }

  //filter after calculation of ungapped sse score to include at least par.min_prefilter_hits
  vector<pair<double, int> >::iterator it;

  sort(first_prefilter.begin(), first_prefilter.end());
  std::reverse(first_prefilter.begin(), first_prefilter.end());

  vector<pair<double, int> >::iterator first_prefilter_begin_erase =
      first_prefilter.end();
  vector<pair<double, int> >::iterator first_prefilter_end_erase =
      first_prefilter.end();
  count_dbs = 0;
  for (it = first_prefilter.begin(); it < first_prefilter.end(); it++) {
    if (count_dbs >= par.min_prefilter_hits
        && (*it).first < par.preprefilter_smax_thresh) {
      first_prefilter_begin_erase = it;
      break;
    }
    else {
      count_dbs++;
    }
  }

  first_prefilter.erase(first_prefilter_begin_erase, first_prefilter_end_erase);

  if (v >= 2) {
    printf(
        "\nHMMs passed 1st prefilter (gapless profile-profile alignment)  : %6i\n",
        count_dbs);
  }

#pragma omp parallel for schedule(static) private(evalue, score, thread_id)
  // Loop over all database sequences
//  for (int n = 0; n < count_dbs; n++) {
  for (it = first_prefilter.begin(); it < first_prefilter.end(); it++) {
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#endif

    int n = (*it).second;

    // Perform search step
    score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend,
        workspace[thread_id], workspace[thread_id] + W,
        workspace[thread_id] + 2 * W, par.prefilter_score_offset);

    evalue = factor * length[n] * fpow2(-score / par.prefilter_bit_factor);

    if (evalue < par.prefilter_evalue_coarse_thresh) {
#pragma omp critical
      hits.push_back(pair<double, int>(evalue, n));
    }
  }

  //filter after calculation of evalues to include at least par.min_prefilter_hits
  sort(hits.begin(), hits.end());

  vector<pair<double, int> >::iterator second_prefilter_begin_erase =
      hits.end();
  vector<pair<double, int> >::iterator second_prefilter_end_erase = hits.end();
  count_dbs = 0;
  for (it = hits.begin(); it < hits.end(); it++) {
    if (count_dbs >= par.min_prefilter_hits
        && (*it).first > par.prefilter_evalue_thresh) {
      second_prefilter_begin_erase = it;
      break;
    }
    else {
      count_dbs++;
    }
  }

  hits.erase(second_prefilter_begin_erase, second_prefilter_end_erase);

  count_dbs = 0;

  for (it = hits.begin(); it < hits.end(); it++) {
    backtrace_hits[count_dbs++] = (*it).second;

    // Add hit to dbfiles
    char name[NAMELEN];
    strcpy(name, dbnames[(*it).second]);

    char db_name[NAMELEN];
    strcpy(db_name, name);

    if (!doubled->Contains(db_name)) {
      doubled->Add(db_name);
      // check, if DB was searched in previous rounds
      strcat(name, "__1");  // irep=1

      if (previous_hits->Contains(name)) {
        dbfiles_old[ndb_old] = new char[strlen(db_name) + 1];
        strcpy(dbfiles_old[ndb_old], db_name);
        ndb_old++;
      }
      else {
        dbfiles_new[ndb_new] = new char[strlen(db_name) + 1];
        strcpy(dbfiles_new[ndb_new], db_name);
        ndb_new++;
      }
    }

    if (count_dbs >= par.maxnumdb) {
      if (v >= 2) {
        fprintf(stderr,
            "WARNING: Number of hits passing 2nd prefilter reduced from %6i to allowed maximum of %i!\n",
            (int) hits.size(), par.maxnumdb);
        fprintf(stderr,
            "You can increase the allowed maximum using the -maxfilt <max> option.\n\n");
      }
      break;
    }
  }

  // Free memory
  free(qc);
  for (int i = 0; i < par.threads; i++)
    free(workspace[i]);
  delete[] workspace;
  delete[] backtrace_hits;
  if (doubled)
    delete doubled;
}


//void HHblits::AlignByWorker(int bin) {
//  // Prepare q ant t and compare
//  PrepareTemplateHMM(par, q, t[bin], format[bin], pb, R);
//
//  // Do HMM-HMM comparison, store results if score>SMIN, and try next best alignment
//  for (hit[bin]->irep = 1; hit[bin]->irep <= par.altali; hit[bin]->irep++) {
//    if (par.forward == 0) {
//      hit[bin]->Viterbi(q, t[bin], par.loc, par.ssm, par.maxres, par.min_overlap, par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);
//      if (hit[bin]->irep > 1 && hit[bin]->score <= SMIN)
//        break;
//      hit[bin]->Backtrace(q, t[bin], par.corr, par.ssw, S73, S33);
//    }
//    else if (par.forward == 2) {
//      hit[bin]->Forward(q, t[bin], par.ssm, par.min_overlap, par.loc, par.shift, par.ssw, par.exclstr, S73, S33);
//      hit[bin]->Backward(q, t[bin], par.loc, par.shift, par.ssw, S73, S33);
//      hit[bin]->MACAlignment(q, t[bin], par.loc, par.mact, par.macins);
//      hit[bin]->BacktraceMAC(q, t[bin], par.corr, par.ssw, S73, S33);
//    }
//    hit[bin]->score_sort = hit[bin]->score_aass;
//    if (hit[bin]->score <= SMIN)
//      hit[bin]->lastrep = 1;
//    else
//      hit[bin]->lastrep = 0;
//
//    #pragma omp critical
//    {
//      hitlist.Push(*(hit[bin])); // insert hit at beginning of list (last repeats first!)
//    }
//
//    // find only best alignment for forward algorithm and stochastic sampling
//    if (par.forward > 0)
//      break;
//
//    // break if score for previous hit is already worse than SMIN
//    if (hit[bin]->score <= SMIN)
//      break;
//  }
//}
//
//
//void HHblits::PerformViterbiByWorker(int bin, Hash<Hit>* previous_hits) {
//  PrepareTemplateHMM(par, q, t[bin], format[bin], pb, R);
//
//  for (hit[bin]->irep = 1; hit[bin]->irep <= par.altali; hit[bin]->irep++) {
//    // Break, if no previous_hit with irep is found
//    hit[bin]->Viterbi(q, t[bin], par.loc, par.ssm, par.maxres, par.min_overlap, par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);
//    if (hit[bin]->irep > 1 && hit[bin]->score <= SMIN)
//      break;
//    hit[bin]->Backtrace(q, t[bin], par.corr, par.ssw, S73, S33);
//
//    hit[bin]->score_sort = hit[bin]->score_aass;
//
//    #pragma omp critical
//    {
//      stringstream ss_tmp;
//      ss_tmp << hit[bin]->file << "__" << hit[bin]->irep;
//
//      if (previous_hits->Contains((char*) ss_tmp.str().c_str())) {
//        //printf("Previous hits contains %s!\n",(char*)ss_tmp.str().c_str());
//        Hit hit_cur = previous_hits->Remove((char*) ss_tmp.str().c_str());
//        previous_hits->Add((char*) ss_tmp.str().c_str(), *(hit[bin]));
//
//        // Overwrite *hit[bin] with alignment, etc. of hit_cur
//        hit_cur.score = hit[bin]->score;
//        hit_cur.score_aass = hit[bin]->score_aass;
//        hit_cur.score_ss = hit[bin]->score_ss;
//        hit_cur.Pval = hit[bin]->Pval;
//        hit_cur.Pvalt = hit[bin]->Pvalt;
//        hit_cur.logPval = hit[bin]->logPval;
//        hit_cur.logPvalt = hit[bin]->logPvalt;
//        hit_cur.Eval = hit[bin]->Eval;
//        hit_cur.logEval = hit[bin]->logEval;
//        hit_cur.Probab = hit[bin]->Probab;
//
//        hitlist.Push(hit_cur); // insert hit at beginning of list (last repeats first!)
//      }
//    }
//
//    // break if score for first hit is already worse than SMIN
//    if (hit[bin]->score <= SMIN)
//      break;
//  }
//}
//
//
/////////////////////////////////////////////////////////////////////////////////////////
////// Realign q and with t[bin] in all hits from same tempate using  MAC algorithm
////////////////////////////////////////////////////////////////////////////////////////
//void HHblits::RealignByWorker(int bin) {
//  // Realign all hits with same template, pointed to by list List<void*>* hit[bin]->plist_phits;
//  // This list is set up in HHseach and HHblits at the beginning of perform_realign()
//  Hit* hit_cur;
//
//  // Prepare MAC comparison(s)
//  PrepareTemplateHMM(par, q, t[bin], format[bin], pb, R);
//  t[bin]->Log2LinTransitionProbs(1.0);
//
//  hit[bin]->irep = 1;
//  hit[bin]->plist_phits->Reset();
//  while (!hit[bin]->plist_phits->End()) {
//    // Set pointer hit_cur to next hit to be realigned
//    hit_cur = (Hit*) hit[bin]->plist_phits->ReadNext();
//    // printf("Realigning %s, irep=%i\n",hit_cur->name,hit_cur->irep);  //?????????
//
//    // Realign only around previous Viterbi hit
//    // hit[bin] = *hit_cur; is not possible because the pointers to the DP matrices would be overwritten
//    hit[bin]->i1 = hit_cur->i1;
//    hit[bin]->i2 = hit_cur->i2;
//    hit[bin]->j1 = hit_cur->j1;
//    hit[bin]->j2 = hit_cur->j2;
//    hit[bin]->nsteps = hit_cur->nsteps;
//    hit[bin]->i = hit_cur->i;
//    hit[bin]->j = hit_cur->j;
//    hit[bin]->realign_around_viterbi = true;
//
//    // Align q to template in *hit[bin]
//    hit[bin]->Forward(q, t[bin], par.ssm, par.min_overlap, par.loc, par.shift, par.ssw, par.exclstr, S73, S33);
//    hit[bin]->Backward(q, t[bin], par.loc, par.shift, par.ssw, S73, S33);
//    hit[bin]->MACAlignment(q, t[bin], par.loc, par.mact, par.macins);
//    hit[bin]->BacktraceMAC(q, t[bin], par.corr, par.ssw, S73, S33);
//
//    // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
//    hit[bin]->score = hit_cur->score;
//    hit[bin]->score_ss = hit_cur->score_ss;
//    hit[bin]->score_aass = hit_cur->score_aass;
//    hit[bin]->score_sort = hit_cur->score_sort;
//    hit[bin]->Pval = hit_cur->Pval;
//    hit[bin]->Pvalt = hit_cur->Pvalt;
//    hit[bin]->logPval = hit_cur->logPval;
//    hit[bin]->logPvalt = hit_cur->logPvalt;
//    hit[bin]->Eval = hit_cur->Eval;
//    hit[bin]->logEval = hit_cur->logEval;
//    hit[bin]->Probab = hit_cur->Probab;
//
//    // Replace original hit in hitlist with realigned hit
//    hit_cur->Delete(); // delete content of pointers etc. of hit_cur (but not DP matrices)
//    *hit_cur = *hit[bin]; // copy all variables and pointers from *hit[bin] into hitlist
//
//    hit[bin]->irep++;
//  }
//
//  return;
//}

void HHblits::wiggleQSC(int n_redundancy, float* qsc, size_t nqsc, HitList& reducedFinalHitList) {
	wiggleQSC(hitlist, n_redundancy, Qali, input_format, qsc, nqsc, reducedFinalHitList);
	reducedFinalHitList.N_searched = hitlist.N_searched;
}

void HHblits::run(FILE* query_fh, char* query_path) {
  int cluster_found = 0;
  int seqs_found = 0;

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

  // Read query input file (HHM, HMMER, or alignment format) without adding pseudocounts
  Qali.N_in = 0;
  ReadQueryFile(par, query_fh, input_format, par.wg, q, &Qali, query_path, pb, S, Sim);

  if (Qali.N_in - Qali.N_ss > 1)
    par.premerge = 0;

  if (par.allseqs) {
    Qali_allseqs = Qali; // make a *deep* copy of Qali!
    for (int k = 0; k < Qali_allseqs.N_in; ++k)
      Qali_allseqs.keep[k] = 1; // keep *all* sequences (reset filtering in Qali)
  }

  v = v1;

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags)
    q->NeutralizeTags(pb);

  // Input parameters
  if (v >= 3) {
    cout << "Input file       :   " << par.infile << "\n";
    cout << "Output file      :   " << par.outfile << "\n";
    cout << "Prefilter DB     :   " << dbcs_data_filename << " "
        << dbcs_index_filename << "\n";
    cout << "HHM DB           :   " << dbhhm_data_filename << " "
        << dbhhm_index_filename << "\n";
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  for (int round = 1; round <= par.num_rounds; round++) {
    if (v >= 2)
      printf("\nIteration %i\n", round);

    // Settings for different rounds
    if (par.premerge > 0 && round > 1
        && previous_hits->Size() >= par.premerge) {
      if (v > 3)
        printf(
            "Set premerge to 0! (premerge: %i   iteration: %i   hits.Size: %i)\n",
            par.premerge, round, previous_hits->Size());
      par.premerge = 0;
    }
    else {
      par.premerge -= previous_hits->Size();
    }

    // Save HMM without pseudocounts for prefilter query-profile
    *q_tmp = *q;

    // Write query HHM file? (not the final HMM, which will be written to par.hhmfile)
    //TODO
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

    PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine, pc_hhm_context_mode, pb, R);

    ////////////////////////////////////////////
    // Prefiltering
    ////////////////////////////////////////////

    if (par.prefilter) {
      if (v >= 2)
        printf("Prefiltering database\n");
      prefilter_db(previous_hits);  // in hhprefilter.C
    }

    if (v >= 2 && ndb_new == 0) {
      printf("No HMMs pass prefilter => Stop searching!\n");
      break;
    }

    // Search datbases
    if (v >= 2) {
      printf(
          "HMMs passed 2nd prefilter (gapped profile-profile alignment)   : %6i\n",
          (ndb_new + ndb_old));
      printf(
          "HMMs passed 2nd prefilter and not found in previous iterations : %6i\n",
          ndb_new);
      printf("Scoring %i HMMs using HMM-HMM Viterbi alignment\n", ndb_new);
    }

    // Main Viterbi HMM-HMM search
    // Starts with empty hitlist (hits of previous iterations were deleted) and creates a hitlist with the hits of this iteration
    ViterbiSearch(dbfiles_new, ndb_new, previous_hits, (ndb_new + ndb_old));

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
      last_round = true;
      if (round < par.num_rounds && v >= 2)
        printf("No new hits found in iteration %i => Stop searching\n", round);

      if (ndb_old > 0 && par.realign_old_hits) {
        if (v > 0) {
          printf("Rescoring previously found HMMs with Viterbi algorithm\n");
        }
        ViterbiSearch(dbfiles_old, ndb_old, previous_hits, (ndb_new + ndb_old));
        // Add dbfiles_old to dbfiles_new for realign
        for (int a = 0; a < ndb_old; a++) {
          dbfiles_new[ndb_new] = new char[strlen(dbfiles_old[a]) + 1];
          strcpy(dbfiles_new[ndb_new], dbfiles_old[a]);
          ndb_new++;
        }
      }
      else if (!par.realign_old_hits && previous_hits->Size() > 0) {
        if (v > 0) {
          printf("Rescoring previously found HMMs with Viterbi algorithm\n");
        }
        RescoreWithViterbiKeepAlignment(ndb_new + previous_hits->Size(), previous_hits);
      }
    }

    // Realign hits with MAC algorithm
    if (par.realign)
      perform_realign(dbfiles_new, ndb_new, premerged_hits);

    // Generate alignment for next iteration
    if (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile
        || *par.alisbasename) {
      v1 = v;
      if (v > 0 && v <= 3)
        v = 1;
      else
        v -= 2;

      // If new hits found, merge hits to query alignment
      if (new_hits != 0) {
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
          Alignment Tali;
          long ftellpos;
          getTemplateA3M(hit_cur.dbfile, ftellpos, Tali);

          if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
            Qali_allseqs.MergeMasterSlave(hit_cur, Tali, hit_cur.dbfile, par.maxcol);
          Tali.N_filtered = Tali.Filter(par.max_seqid_db, S, par.coverage_db,
              par.qid_db, par.qsc_db, par.Ndiff_db);
          Qali.MergeMasterSlave(hit_cur, Tali, hit_cur.dbfile, par.maxcol);

          if (Qali.N_in >= MAXSEQ)
            break; // Maximum number of sequences reached
        }

        // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
        Qali.Compress("merged A3M file", par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);

        // Sort out the nseqdis most dissimilacd r sequences for display in the result alignments
        Qali.FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid, par.qsc,
            par.nseqdis);

        // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
        const float COV_ABS = 25;     // min. number of aligned residues
        int cov_tot = imax(imin((int) (COV_ABS / Qali.L * 100 + 0.5), 70), par.coverage);
        if (v > 2)
          printf("Filter new alignment with cov %3i%%\n", cov_tot);
        Qali.N_filtered = Qali.Filter(par.max_seqid, S, cov_tot, par.qid, par.qsc, par.Ndiff);
      }

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, par.maxres, pb, Sim, NULL, true);

      if (par.notags)
        q->NeutralizeTags(pb);

      // Calculate SSpred if we need to print out alis after each iteration or if last iteration
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
        && (round < par.num_rounds || *par.alnfile || *par.psifile || *par.hhmfile
            || *par.alisbasename))
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
    float qscs[] = {-20, 0, 0.1, 0.2};
	wiggleQSC(par.n_redundancy, qscs, 4, reducedHitlist);
  }

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete(); // Delete hit object
  delete previous_hits;

  delete premerged_hits;
}


void HHblits::writeHHRFile(char* hhrFile) {
	if(*hhrFile) {
		hitlist.PrintHHR(q_tmp, hhrFile, par.maxdbstrlen, par.showconf, par.showcons, par.showdssp, par.showpred, par.b, par.B, par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E, par.argc, par.argv, S);
	}
}


void HHblits::writeAlisFile(char* basename) {
	if(*basename) {
	  std::map<int, Alignment>::iterator it;
	  for(it = alis.begin(); it != alis.end(); it++) {
		stringstream ss_tmp;
		ss_tmp << basename << "_" << (*it).first << ".a3m";
		std::string id = ss_tmp.str();

        (*it).second.WriteToFile(id.c_str(), par.append, "a3m");
	  }
	}
}


void HHblits::writeScoresFile(char* scoresFile) {
	if(*scoresFile) {
		hitlist.PrintScoreFile(q, scoresFile);
	}
}


void HHblits::writePairwiseAlisFile(char* pairwiseAlisFile, char outformat) {
	if (*pairwiseAlisFile) {
		hitlist.PrintAlignments(q, pairwiseAlisFile, par.showconf, par.showcons, par.showdssp, par.showpred, par.p, par.aliwidth, par.nseqdis, par.b, par.B, par.E, S, outformat);
	}
}


void HHblits::writeAlitabFile(char* alitabFile) {
	if (*alitabFile) {
		hitlist.WriteToAlifile(q, alitabFile, par.alitab_scop);
	}
}


void HHblits::writeReducedHHRFile(char* reducedHHRFile) {
	if(*reducedHHRFile) {
		reducedHitlist.PrintHHR(q_tmp, reducedHHRFile, par.maxdbstrlen, par.showconf, par.showcons, par.showdssp, par.showpred, par.b, par.B, par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E, par.argc, par.argv, S);
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

	    q->WriteToFile(HMMFile, par.append, par.max_seqid, par.coverage, par.qid, par.Ndiff, par.qsc, par.argc, par.argv, pb);
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


std::stringstream* HHblits::writeHHRFile() {
	std::stringstream* out = new std::stringstream();
	hitlist.PrintHHR(q_tmp, *out, par.maxdbstrlen, par.showconf, par.showcons, par.showdssp, par.showpred, par.b, par.B, par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E, par.argc, par.argv, S);
	return out;
}


std::stringstream* HHblits::writeScoresFile() {
	std::stringstream* out = new std::stringstream();
	hitlist.PrintScoreFile(q, *out);
	return out;
}


std::stringstream* HHblits::writePairwiseAlisFile(char outformat) {
	std::stringstream* out = new std::stringstream();
	hitlist.PrintAlignments(q, *out, par.showconf, par.showcons, par.showdssp, par.showpred, par.p, par.aliwidth, par.nseqdis, par.b, par.B, par.E, S, outformat);
	return out;
}


std::stringstream* HHblits::writeAlitabFile() {
	std::stringstream* out = new std::stringstream();
	hitlist.WriteToAlifile(q, *out, par.alitab_scop);
	return out;
}


std::stringstream* HHblits::writeReducedHHRFile() {
	std::stringstream* out = new std::stringstream();
	reducedHitlist.PrintHHR(q_tmp, *out, par.maxdbstrlen, par.showconf, par.showcons, par.showdssp, par.showpred, par.b, par.B, par.z, par.Z, par.aliwidth, par.nseqdis, par.p, par.E, par.argc, par.argv, S);
	return out;
}


std::stringstream* HHblits::writePsiFile() {
	std::stringstream* out = new std::stringstream();
	if (par.allseqs)
      Qali_allseqs.WriteToFile(*out, "psi");
    else
      Qali.WriteToFile(*out, "psi");
	return out;
}


std::stringstream* HHblits::writeHMMFile() {
    // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
    q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
    q->CalculateAminoAcidBackground(pb);

	std::stringstream* out = new std::stringstream();
	q->WriteToFile(*out, par.max_seqid, par.coverage, par.qid, par.Ndiff, par.qsc, par.argc, par.argv, pb);
	return out;
}


std::stringstream* HHblits::writeA3MFile() {
	std::stringstream* out = new std::stringstream();
	if (par.allseqs)
      Qali_allseqs.WriteToFile(*out, "a3m");
    else
      Qali.WriteToFile(*out, "a3m");
	return out;
}
