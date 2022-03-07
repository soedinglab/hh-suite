// hhfilterC: filter alignment in a2m format with maximum sequence identity of match states and minimum coverage
//
//     (C) Johannes Soeding 2012

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

#include "hhsuite_config.h"

#include "util.h"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhdecl.h"      // Constants, global variables, struct Parameters
#include "hhutil.h"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID
#include "hhmatrices.h"  // BLOSUM50, GONNET, HSDM
#include "hhalignment.h" // class Alignment
#include "hhfunc.h"      // some functions common to hh programs

void help(Parameters& par) {
  printf("HHfilter %i.%i.%i\n", HHSUITE_VERSION_MAJOR, HHSUITE_VERSION_MINOR, HHSUITE_VERSION_PATCH);
  printf("Filter an alignment by maximum pairwise sequence identity, minimum coverage,\n");
  printf("minimum sequence identity, or score per column to the first (seed) sequence.n");
  printf("%s", COPYRIGHT);
  printf("%s", REFERENCE);
  printf("\n");
  printf("Usage: hhfilter -i infile -o outfile [options]\n");
  printf(" -i <file>      read input file in A3M/A2M or FASTA format                 \n");
  printf(" -o <file>      write to output file in A3M format                         \n");
  printf(" -a <file>      append to output file in A3M format                        \n");
  printf("\n");
  printf("Options:                                                                  \n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n", par.max_seqid);
  printf(" -diff [0,inf]  filter MSA by selecting most diverse set of sequences, keeping \n");
  printf("                at least this many seqs in each MSA block of length 50 (def=%i) \n", par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n", par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n", par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n", par.qsc);
  printf(" -neff [1,inf]  target diversity of alignment (default=off)\n");
  printf("\n");
  printf("Input alignment format:                                                    \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("                                                                          \n");
  printf("Other options:\n");
  printf(" -maxseq <int>  max number of input rows (def=%5i)\n", par.maxseq);
  printf(" -maxres <int>  max number of HMM columns (def=%5i)\n", par.maxres);
  printf("Example: hhfilter -id 50 -i d1mvfd_.a2m -o d1mvfd_.fil.a2m          \n\n");
}

//// Processing input options from command line
void ProcessArguments(Parameters& par) {
  const int argc = par.argc;
  const char** argv = par.argv;

  // Read command line options
  for (int i = 1; i <= argc - 1; i++) {
	HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
	if (!strcmp(argv[i], "-i")) {
      if (++i > argc - 1 || argv[i][0] == '-') {
        HH_LOG(ERROR) << "No input file following -f" << std::endl;
        exit(4);
      }
      else {
        strcpy(par.infile, argv[i]);
      }
    } else if (!strcmp(argv[i], "-o")) {
      par.append = 0;
      if (++i > argc - 1) {
        HH_LOG(ERROR) << "No output file following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    } else if (!strcmp(argv[i], "-a")) {
      par.append = 1;
      if (++i > argc - 1) {
        HH_LOG(ERROR) << "No output file following -a" << std::endl;
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    } else if (!strcmp(argv[i], "-v") && (i + 1 < argc) && argv[i + 1][0] != '-') {
		int v = atoi(argv[++i]);
		par.v = Log::from_int(v);
		Log::reporting_level() = par.v;
    }
    else if (!strcmp(argv[i], "-maxseq") && (i < argc - 1))
      par.maxseq = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = par.maxres * 2;
    } else if (!strcmp(argv[i], "-id") && (i < argc - 1))
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
    } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      help(par);
      exit(0);
    } else {
    	HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << std::endl;
    }
  }
}

int main(int argc, const char **argv) {
  Parameters par(argc, argv);

  strcpy(par.infile, "");
  strcpy(par.outfile, "");

  // maximum number of sequences to be written
  par.nseqdis = par.maxseq - 1;
  // no filtering for maximum diversity
  par.Ndiff = 0;

  ProcessArguments(par);

  // Check command line input and default values
  if (!*par.infile) {
    help(par);
    HH_LOG(ERROR) << "Input file is missing!" << std::endl;
    exit(4);
  }
  if (!*par.outfile) {
    help(par);
    HH_LOG(ERROR) << "Output file is missing!" << std::endl;
    exit(4);
  }

  HH_LOG(INFO) << "Input file = " << par.infile << "\n";
  HH_LOG(INFO) << "Output file = " << par.outfile << "\n";

  // Reads in an alignment from par.infile into matrix X[k][l] as ASCII
  FILE* inf = NULL;
  if (strcmp(par.infile, "stdin")) {
    inf = fopen(par.infile, "r");
    if (!inf) {
      OpenFileError(par.infile, __FILE__, __LINE__, __func__);
    }
  }
  else {
    inf = stdin;
  }

  Alignment qali(par.maxseq, par.maxres);
  qali.Read(inf, par.infile, par.mark, par.maxcol, par.nseqdis);
  fclose(inf);

  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
  // and store marked sequences in name[k] and seq[k]
  qali.Compress(par.infile, par.cons, par.maxcol, par.M, par.Mgaps);

  // substitution matrix flavours
  float __attribute__((aligned(16))) P[20][20];
  float __attribute__((aligned(16))) R[20][20];
  float __attribute__((aligned(16))) Sim[20][20];
  float __attribute__((aligned(16))) S[20][20];
  float __attribute__((aligned(16))) pb[21];
  SetSubstitutionMatrix(par.matrix, pb, P, R, S, Sim);

  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
  qali.N_filtered = qali.Filter(par.max_seqid, S, par.coverage, par.qid, par.qsc,par.Ndiff);

  // Atune alignment diversity q.Neff with qsc to value Neff_goal
  if (par.Neff >= 1.0) {
    qali.FilterNeff(par.wg, par.mark, par.cons, par.showcons, par.max_seqid, par.coverage, par.Neff, pb, S, Sim);
  }

  // Write filtered alignment WITH insert states (lower case) to alignment file
  qali.WriteToFile(par.outfile, par.append);
}
