// hhconsensus.cpp: read A3M/FASTA file and calculate consensus sequence

#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>   // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror()
#include <cassert>
#include <stdexcept>

#include "hhsuite_config.h"
#include "cs.h"          // context-specific pseudocounts
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"

#include "list.h"        // list data structure
#include "hash.h"        // hash data structure
#include "util.h"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhdecl.h"      // Constants, global variables, struct Parameters
#include "hhutil.h"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.h"  // BLOSUM50, GONNET, HSDM
#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfunc.h"      // some functions common to hh programs

// Help functions
void help(Parameters& par) {
  printf("\n");
  printf("HHconsensus %i.%i.%i\n", HHSUITE_VERSION_MAJOR, HHSUITE_VERSION_MINOR, HHSUITE_VERSION_PATCH);
  printf("Calculate the consensus sequence for an A3M/FASTA input file.   \n");
  printf("%s", COPYRIGHT);
  printf("%s", REFERENCE);
  printf("\n");
  printf("Usage: %s -i <file> [options]                           \n", par.program_name);
  printf(
      " -i <file>     query alignment (A2M, A3M, or FASTA), or query HMM          \n");
  printf("\n");
  printf(
      "Output options:                                                            \n");
  printf(
      " -s <file>     append consensus sequence in FASTA (default=<infile.seq>)   \n");
  printf(
      " -o <file>     write alignment with consensus sequence in A3M              \n");
  printf(
      " -oa3m <file>  same                                                        \n");
  printf(
      " -oa2m <file>  write alignment with consensus sequence in A2M              \n");
  printf(
      " -ofas <file>  write alignment with consensus sequence in FASTA            \n");
  printf(
      " -v <int>      verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf("\n");
  printf(
      "Filter input alignment (options can be combined):                         \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",
      par.max_seqid);
  printf(
      " -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf(
      "               many sequences in each block of >50 columns (def=%i)\n",
      par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",
      par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",
      par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",
      par.qsc);
  printf("\n");
  printf(
      "Input alignment format:                                                    \n");
  printf(
      " -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf(
      "               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(
      " -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(
      " -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");
  printf(
      "Other options:                                                               \n");
  printf(
      " -maxres <int> max number of HMM columns (def=%5i)                           \n",
            par.maxres);
  printf("\n");
  printf("Example: %s -i stdin -s stdout\n", par.program_name);
  printf("\n");
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(Parameters& par) {
  const int argc = par.argc;
  const char** argv = par.argv;

  // Read command line options
  for (int i = 1; i <= argc - 1; i++) {
    HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No input file following -i" << std::endl;
        exit(4);
      }
      else
        strcpy(par.infile, argv[i]);
    }
    else if (!strcmp(argv[i], "-s")) {
      if (++i >= argc) {
        help(par);
        HH_LOG(ERROR) << "No output file following -s" << std::endl;
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-o")) {
      par.outformat = 3;
      if (++i >= argc) {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-ofas")) {
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa2m")) {
      par.outformat = 2;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa3m")) {
      par.outformat = 3;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      help(par);
      exit(0);
    }
    else if (!strcmp(argv[i], "-v") && (i < argc - 1) && argv[i + 1][0] != '-') {
		int v = atoi(argv[++i]);
		par.v = Log::from_int(v);
		Log::reporting_level() = par.v;
    }
    else if (!strcmp(argv[i], "-seq") && (i < argc - 1))
      par.nseqdis = atoi(argv[++i]);
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
    else if (!strcmp(argv[i], "-M") && (i < argc - 1)) {
      if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m"))
        par.M = 1;
      else if (!strcmp(argv[i], "first"))
        par.M = 3;
      else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
        par.Mgaps = atoi(argv[i]);
        par.M = 2;
      }
      else
        HH_LOG(WARNING) << "Ignoring unknown argument: -M " << argv[i] << std::endl;
    }
    else if (!strcmp(argv[i], "-Gonnet"))
      par.matrix = 0;
    else if (!strncmp(argv[i], "-BLOSUM", 7)
        || !strncmp(argv[i], "-Blosum", 7)) {
      if (!strcmp(argv[i] + 7, "30"))
        par.matrix = 30;
      else if (!strcmp(argv[i] + 7, "40"))
        par.matrix = 40;
      else if (!strcmp(argv[i] + 7, "50"))
        par.matrix = 50;
      else if (!strcmp(argv[i] + 7, "65"))
        par.matrix = 65;
      else if (!strcmp(argv[i] + 7, "80"))
        par.matrix = 80;
      else
        HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << std::endl;
    }
    else if (!strcmp(argv[i], "-pcm") && (i < argc - 1))
      par.pc_hhm_context_engine.admix = (Pseudocounts::Admix) atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pca") && (i < argc - 1))
      par.pc_hhm_context_engine.pca = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pcb") && (i < argc - 1))
      par.pc_hhm_context_engine.pcb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pcc") && (i < argc - 1))
      par.pc_hhm_context_engine.pcc = atof(argv[++i]);
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
    else if (!strcmp(argv[i], "-maxseq") && (i < argc - 1))
      par.maxseq = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = par.maxres * 2;
    } else if (!strcmp(argv[i], "-nocontxt"))
      par.nocontxt = 1;
    else if (!strcmp(argv[i], "-csb") && (i < argc - 1))
      par.csb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-csw") && (i < argc - 1))
      par.csw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-cs")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No query file following -cs" << std::endl;
        exit(4);
      }
      else
        par.clusterfile = argv[i];
    }
    else if (!strcmp(argv[i], "-name")) {
        // skip this, its handled somewhere else
        ;
    }
    else {
		HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << " ...\n";
    }

    HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
  } // end of for-loop for command line input
}

/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char **argv) {
  Parameters par(argc, argv);

  strcpy(par.infile, "");
  strcpy(par.outfile, "");
  strcpy(par.alnfile, "");

  //Default parameter settings
  par.nseqdis = par.maxseq - 1;        // maximum number of sequences to be written
  par.showcons = 0;
  par.cons = 1;
  par.Ndiff = 0;
  par.max_seqid = 100;
  par.coverage = 0;
  par.pc_hhm_context_engine.pca = 0.0;  // no amino acid pseudocounts
  par.pc_hhm_nocontext_a = 0.0;  // no amino acid pseudocounts
  par.gapb = 0.0; // no transition pseudocounts

  ProcessArguments(par);

  Alignment* qali = new Alignment(par.maxseq, par.maxres);
  HMM* q = new HMM(MAXSEQDIS, par.maxres);        //Create a HMM with maximum of par.maxres match states

  // q is only available after maxres is known, so we had to move this here
  for (int i = 1; i <= argc - 1; i++) {
      if (!strcmp(argv[i], "-name") && (i < argc - 1)) {
          strmcpy(q->name, argv[++i], NAMELEN - 1); //copy longname to name...
          strmcpy(q->longname, argv[i], DESCLEN - 1);   //copy full name to longname
      }
  }

  // Check command line input and default values
  if (!*par.infile) {
    help(par);
    HH_LOG(ERROR) << "Input file is missing!" << std::endl;
    exit(4);
  }

  // Get basename
  RemoveExtension(q->file, par.infile); //Get basename of infile (w/o extension):

  // Outfile not given? Name it basename.hhm
  if (!*par.outfile && !*par.alnfile) {
    RemoveExtension(par.outfile, par.infile);
    strcat(par.outfile, ".seq");
  }

  cs::ContextLibrary<cs::AA>* context_lib = NULL;
  cs::Crf<cs::AA>* crf = NULL;
  cs::Pseudocounts<cs::AA>* pc_hhm_context_engine = NULL;
  cs::Admix* pc_hhm_context_mode = NULL;
  cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine = NULL;
  cs::Admix* pc_prefilter_context_mode = NULL;

  // Prepare CS pseudocounts lib
  if (!par.nocontxt) {
    InitializePseudocountsEngine(par, context_lib, crf, pc_hhm_context_engine,
        pc_hhm_context_mode, pc_prefilter_context_engine,
        pc_prefilter_context_mode);
  }

  // substitution matrix flavours
  float __attribute__((aligned(16))) P[20][20];
  float __attribute__((aligned(16))) R[20][20];
  float __attribute__((aligned(16))) Sim[20][20];
  float __attribute__((aligned(16))) S[20][20];
  float __attribute__((aligned(16))) pb[21];

  // Set substitution matrix; adjust to query aa distribution if par.pcm==3
  SetSubstitutionMatrix(par.matrix, pb, P, R, S, Sim);

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  char input_format = 0;
  ReadQueryFile(par, par.infile, input_format, par.wg, q, qali, pb, S, Sim);

  // Same code as in PrepareQueryHMM(par.infile,input_format,q,qali), except that we add SS prediction
  // Add Pseudocounts, if no HMMER input
  if (input_format == 0) {
    // Transform transition freqs to lin space if not already done
    q->AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg,
        par.gaph, par.gapi, par.gapb, par.gapb);

    // Comput substitution matrix pseudocounts
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
  }
  else {
    q->AddAminoAcidPseudocounts(0, par.pc_hhm_nocontext_a,
        par.pc_hhm_nocontext_b, par.pc_hhm_nocontext_c);
  }

  q->CalculateAminoAcidBackground(pb);

  if (par.columnscore == 5 && !q->divided_by_local_bg_freqs)
    q->DivideBySqrtOfLocalBackgroundFreqs(
        par.half_window_size_local_aa_bg_freqs, pb);

  // Write consensus sequence to sequence file
  // Consensus sequence is calculated in hhalignment.C, Alignment::FrequenciesAndTransitions()
  if (*par.outfile) {
    FILE* outf = NULL;
    if (strcmp(par.outfile, "stdout")) {
      outf = fopen(par.outfile, "a");
      if (!outf)
        OpenFileError(par.outfile, __FILE__, __LINE__, __func__);
    }
    else
      outf = stdout;
    // OLD
    //// ">name_consensus" -> ">name consensus"
    //strsubst(q->sname[q->nfirst],"_consensus"," consensus");
    //fprintf(outf,">%s\n%s\n",q->sname[q->nfirst],q->seq[q->nfirst]+1);
    // NEW (long header needed for NR30cons database)
    fprintf(outf, ">%s\n%s\n", q->longname, q->seq[q->nfirst] + 1);
    fclose(outf);
  }

  // Print A3M/A2M/FASTA output alignment
  if (*par.alnfile) {
    HalfAlignment qa(MAXSEQDIS);
    int n = imin(q->n_display,
        par.nseqdis + (q->nss_dssp >= 0) + (q->nss_pred >= 0)
            + (q->nss_conf >= 0) + (q->ncons >= 0));
    qa.Set(q->name, q->seq, q->sname, n, q->L, q->nss_dssp, q->nss_pred,
        q->nss_conf, q->nsa_dssp, q->ncons);

    if (par.outformat == 1)
      qa.BuildFASTA();
    else if (par.outformat == 2)
      qa.BuildA2M();
    else if (par.outformat == 3)
      qa.BuildA3M();
    if (qali->readCommentLine)
      qa.Print(par.alnfile, par.append, qali->longname); // print alignment to outfile
    else
      qa.Print(par.alnfile, par.append);   // print alignment to outfile
  }

  delete qali;
  delete q;

  DeletePseudocountsEngine(context_lib, crf, pc_hhm_context_engine,
      pc_hhm_context_mode, pc_prefilter_context_engine,
      pc_prefilter_context_mode);
}
