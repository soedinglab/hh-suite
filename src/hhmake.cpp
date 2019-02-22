// hhmake.C: build profile HMM from input alignment for HMM-HMM comparison

//     (C) Johannes Soeding and Michael Remmert 2012

//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

//     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

//     Reference:
//     Remmert M., Biegert A., Hauser A., and Soding J.
//     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
//     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#include <fstream>    // ofstream, ifstream
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <errno.h>    // perror(), strerror(errno)
#include <ctype.h>    // islower, isdigit etc
#include <cassert>

#include "hhsuite_config.h"
#include "cs.h"          // context-specific pseudocounts
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"

#include "list.h"        // list data structure
#include "hash.h"        // hash data structure
#include "util.h"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhdecl.h"      // Constants, global variables, struct Parameters
#include "hhutil.h"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID
#include "hhmatrices.h"  // BLOSUM50, GONNET, HSDM
#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhfunc.h"      // some functions common to hh programs

// Help functions
void help(Parameters& par, char all = 0) {
  printf("HHmake %i.%i.%i\n", HHSUITE_VERSION_MAJOR, HHSUITE_VERSION_MINOR, HHSUITE_VERSION_PATCH);
  printf(
      "Build an HMM from an input alignment in A2M, A3M, or FASTA format,   \n");
  printf(
      "or convert between HMMER format (.hmm) and HHsearch format (.hhm).   \n");
  printf("%s", REFERENCE);
  printf("%s", COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i file [options]                                       \n", par.program_name);
  printf(" -i <file>     query alignment (A2M, A3M, or FASTA), or query HMM         \n");
  if (all) {
    printf("\n");
    printf("<file> may be 'stdin' or 'stdout' throughout.\n");
  }
  printf("\n");
  printf(
      "Output options:                                                           \n");
  printf(
      " -o <file>     HMM file to be written to  (default=<infile.hhm>)          \n");
  printf(
      " -a <file>     HMM file to be appended to                                 \n");
  printf(
      " -v <int>      verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(
      " -seq <int>    max. number of query/template sequences displayed (def=%i)  \n",
      par.nseqdis);
  printf(
      "               Beware of overflows! All these sequences are stored in memory.\n");
  printf(
      " -add_cons         make consensus sequence master sequence of query MSA \n");
  printf(
      " -name <name>  use this name for HMM (default: use name of first sequence)   \n");
  printf("\n");
  printf(
      "Filter query multiple sequence alignment                                     \n");
  printf(
      " -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",
      par.max_seqid);
  printf(
      " -diff [0,inf[  filter MSA by selecting most diverse set of sequences, keeping \n");
  printf(
      "                at least this many seqs in each MSA block of length 50 (def=%i) \n",
      par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",
      par.coverage);
  printf(
      " -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",
      par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",
      par.qsc);
  printf(" -neff [1,inf]  target diversity of alignment (default=off)\n");
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
  if (all) {
    printf("Pseudocount (pc) options:                                                        \n");
    printf(" Context specific hhm pseudocounts:\n");
    printf("  -pc_hhm_contxt_mode {0,..,3}   position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc_hhm_context_engine.admix);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_hhm_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc_hhm_context_engine.pca);
    printf("  -pc_hhm_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",par.pc_hhm_context_engine.pcb);
    printf("  -pc_hhm_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",par.pc_hhm_context_engine.pcc);

    printf(" Context independent hhm pseudocounts (used for templates; used for query if contxt file is not available):\n");
    printf("  -pc_hhm_nocontxt_mode {0,..,3}   position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc_hhm_nocontext_mode);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
  //  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_hhm_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc_hhm_nocontext_a);
    printf("  -pc_hhm_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",par.pc_hhm_nocontext_b);
    printf("  -pc_hhm_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",par.pc_hhm_nocontext_c);
    printf(
        " Context-specific pseudo-counts:                                                  \n");
    printf(
        "  -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
    printf(
        "  -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n\n", par.clusterfile.c_str());
    printf("Other options:                                                                   \n");
    printf(" -maxres <int>  max number of HMM columns (def=%5i)             \n", par.maxres);
    printf("\n");
  }

  printf("Example: %s -i test.a3m \n", par.program_name);
  printf("\n");
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(Parameters& par, std::string& name) {
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
    else if (!strcmp(argv[i], "-o")) {
      par.append = 0;
      if (++i >= argc) {
        help(par);
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-a")) {
      par.append = 1;
      if (++i >= argc) {
        help(par);
        HH_LOG(ERROR) << "No output file following -a" << std::endl;
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      help(par, 1);
      exit(0);
    }
    else if (!strcmp(argv[i], "-v") && (i < argc - 1) && argv[i + 1][0] != '-') {
      int v = atoi(argv[++i]);
      par.v = Log::from_int(v);
      Log::reporting_level() = par.v;
    }
    else if (!strcmp(argv[i], "-seq") && (i < argc - 1))
      par.nseqdis = atoi(argv[++i]);
    else if (!strncmp(argv[i], "-add_cons", 5))
      par.cons = 1;
    else if (!strncmp(argv[i], "-mark", 5))
      par.mark = 1;
    else if (!strcmp(argv[i], "-name") && (i < argc - 1)) {
      name = argv[++i];
    }
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
    else if (!strcmp(argv[i], "-neff") && (i < argc - 1))
      par.Neff = atof(argv[++i]);
    else if (!strcmp(argv[i], "-Neff") && (i < argc - 1))
      par.Neff = atof(argv[++i]);
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
    else if (!strcmp(argv[i], "-wg"))
      par.wg = 1;
    else if (!strcmp(argv[i], "-maxseq") && (i < argc - 1))
      par.maxseq = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = 2 * par.maxres;
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
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1))
      par.maxres = par.maxcol = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-nocontxt")) {
      par.nocontxt = 1;
    }
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
	else {
		HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << std::endl;
	}

	HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;

  } // end of for-loop for command line input
}


int main(int argc, const char **argv) {
  Parameters par(argc, argv);

  strcpy(par.infile, "");
  strcpy(par.outfile, "");
  strcpy(par.alnfile, "");

  //Default parameter settings
  par.showcons = 1;              // write consensus sequence into hhm file
  par.append = 0;                // overwrite output file
  par.nseqdis = 10; // maximum number of query or template sequences to be recoreded in HMM and diplayed in output alignments
  par.mark = 0; // 1: only marked sequences (or first) get displayed; 0: most divergent ones get displayed
  par.max_seqid = 90;         // default for maximum sequence identity threshold
  par.qid = 0;                // default for maximum sequence identity threshold
  par.qsc = -20.0f;           // default for minimum score per column with query
  par.coverage = 0;              // default for minimum coverage threshold
  par.Ndiff = 100;         // pick Ndiff most different sequences from alignment
  par.M = 1;                     // match state assignment is by A2M/A3M
  par.Mgaps = 50; // above this percentage of gaps, columns are assigned to insert states
  par.matrix = 0;    // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50 3: BLOSUM62
  par.pc_hhm_context_engine.admix = (Pseudocounts::Admix) 0; // no amino acid and transition pseudocounts added
  par.gapb = 0.0; // default values for transition pseudocounts; 0.0: add no transition pseudocounts!
  par.wg = 0;               // 0: use local sequence weights   1: use local ones

  std::string name;
  ProcessArguments(par, name);

  // Check command line input and default values
  if (!*par.infile) {
    help(par);
    HH_LOG(ERROR) << "Input file is missing" << std::endl;
    exit(4);
  }
  if (par.nseqdis > MAXSEQDIS - 3) {
    // 3 reserve for secondary structure
    par.nseqdis = MAXSEQDIS - 3;
  }

  HMM* q = new HMM(par.nseqdis, par.maxres);
  strmcpy(q->name, name.c_str(), NAMELEN - 1);
  strmcpy(q->longname, name.c_str(), DESCLEN - 1);


  // Get basename
  RemoveExtension(q->file, par.infile); //Get basename of infile (w/o extension):

  // Outfile not given? Name it basename.hhm
  if (!*par.outfile) {
    RemoveExtension(par.outfile, par.infile);
    strcat(par.outfile, ".hhm");
  }

  cs::ContextLibrary<cs::AA>* context_lib = NULL;
  cs::Crf<cs::AA>* crf = NULL;
  cs::Pseudocounts<cs::AA>* pc_hhm_context_engine = NULL;
  cs::Admix* pc_hhm_context_mode = NULL;
  cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine = NULL;
  cs::Admix* pc_prefilter_context_mode = NULL;

  // Prepare CS pseudocounts lib
  if (!par.nocontxt) {
    InitializePseudocountsEngine(par, context_lib, crf, pc_hhm_context_engine, pc_hhm_context_mode, pc_prefilter_context_engine, pc_prefilter_context_mode);
  }

  // substitution matrix flavours
  float __attribute__((aligned(16))) P[20][20];
  float __attribute__((aligned(16))) R[20][20];
  float __attribute__((aligned(16))) Sim[20][20];
  float __attribute__((aligned(16))) S[20][20];
  float __attribute__((aligned(16))) pb[21];
  float __attribute__((aligned(16))) qav[21];

// secondary structure matrices
  float S73[NDSSP][NSSPRED][MAXCF];
  float S33[NSSPRED][MAXCF][NSSPRED][MAXCF];

  // Set substitution matrix; adjust to query aa distribution if par.pcm==3
  SetSubstitutionMatrix(par.matrix, pb, P, R, S, Sim);

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  char input_format = 0;
  Alignment* Qali = new Alignment(par.maxseq, par.maxres);
  ReadQueryFile(par, par.infile, input_format, par.wg, q, Qali, pb, S, Sim);
  PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine, pc_hhm_context_mode, pb, R);

  // Write HMM to output file in HHsearch format
  q->WriteToFile(par.outfile, par.append, par.max_seqid, par.coverage, par.qid, par.Ndiff, par.qsc, par.argc, par.argv, pb);


  delete q;
  delete Qali;
  DeletePseudocountsEngine(context_lib, crf, pc_hhm_context_engine, pc_hhm_context_mode, pc_prefilter_context_engine, pc_prefilter_context_mode);
}
