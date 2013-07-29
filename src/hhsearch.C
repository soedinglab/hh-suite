// hhsearch.C:
// Search for a multiple alignment (transformed into HMM) in a profile HMM database

// Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: command line error  6: internal logic error  7: internal numeric error


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

//     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

//     Reference: 
//     Remmert M., Biegert A., Hauser A., and Soding J.
//     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
//     Nat. Methods 9:173-175 (2011); epub Dec 25, doi: 10.1038/NMETH.1818

//
// TODOs
//
// * Allow using the base name (without extension) in hhsearch to make it compatible with the ussage of hhblits
//



////#define WINDOWS
#define PTHREAD
#define MAIN
#include <iostream>   // cin, cout, cerr
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
#include <map>
#include <unistd.h>   // access()

#ifdef PTHREAD
#include <pthread.h>  // POSIX pthread functions and data structures
#endif

#include <sys/time.h>
//#include <new>
//#include "efence.h"
//#include "efence.c"

#ifdef HH_SSE41
#include <tmmintrin.h>   // SSSE3
#include <smmintrin.h>   // SSE4.1
#define HH_SSE3
#endif

#ifdef HH_SSE3
#include <pmmintrin.h>   // SSE3
#define HH_SSE2
#endif

#ifdef HH_SSE2
#ifndef __SUNPRO_C
#include <emmintrin.h>   // SSE2
#else
#include <sunmedia_intrin.h>
#endif
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

#include "cs.h"          // context-specific pseudocounts
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure

#include "hhdecl.C"      // Constants, global variables, struct Parameters

std::map<std::string, unsigned char*> columnStateSequences;
ColumnStateScoring* columnStateScoring;

std::vector<Alignment_Matrices*> matrices;
std::map<std::string, float> ali_probabilities;

#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM

#include "hhhmm.h"       // class Hit
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfullalignment.h" // class FullAlignment
#include "hhhitlist.h"   // class Hit

#include "hhhmm.C"       // class HMM
#include "hhalignment.C" // class Alignment
#include "hhhit.C"       // class Hit
#include "hhhalfalignment.C" // class HalfAlignment
#include "hhfullalignment.C" // class FullAlignment
#include "hhhitlist.C"   // class HitList
#include "hhfunc.C"      // some functions common to hh programs


/////////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////////

const char HHSEARCH_REFERENCE[]="Soding, J. Protein homology detection by HMM-HMM comparison. Bioinformatics 21:951-960 (2005).\n";
const int MAXTHREADS=256; // maximum number of threads (i.e. CPUs) for parallel computation
const int MAXBINS=384;    // maximum number of bins (positions in thread queue)
enum bin_states {FREE=0, SUBMITTED=1, RUNNING=2};
int v1;                   // verbose mode
int threads=0;            // number of threads (apart from the main thread which reads from the databases file) 0:no multithreading
int bins;                 // number of bins; jobs gets allocated to a FREE bin were they are waiting for exection by a thread
char bin_status[MAXBINS]; // The status for each bin is FREE, SUBMITTED, or RUNNING
int jobs_running;         // number of active jobs, i.e. number of bins set to RUNNING
int jobs_submitted;       // number of submitted jobs, i.e. number of bins set to SUBMITTED
char reading_dbs;         // 1: still HMMs to read in a database;  0: finshed reading all HMMs, no db left
const char DEBUG_THREADS=0; // Debugging flag

HMM* q = new HMM;         // Create query HMM with maximum of par.maxres match states
HMM* t[MAXBINS];          // Each bin has a template HMM allocated that was read from the database file
Hit* hit[MAXBINS];        // Each bin has an object of type Hit allocated with a separate dynamic programming matrix (memory!!)
HitList hitlist;          // list of hits with one Hit object for each pairwise comparison done
Hit hit_cur;              // current Hit element read from hitlist
int* format;              // format[bin] = 0 if in HHsearch format => add pcs; format[bin] = 1 if in HMMER format => no pcs
int read_from_db;         // The value of this flag is returned from HMM::Read(); 0:end of file  1:ok  2:skip HMM
int N_searched;  // Number of HMMs searched

char line[LINELEN]="";    // input line
int bin;                  // bin index


struct Thread_args // data to hand to WorkerLoop thread
{
  int thread_id;          // id of thread (for debugging)
  void (*function)(int);  // pointer to function (=job) to execute by WorkerLoop once old job is done
};

#ifdef PTHREAD
struct Thread_args thread_data[MAXTHREADS]; // store a threads thread_id and function to call (AlignByWorker, RealignByWorker)
pthread_t pthread[MAXTHREADS]; // info on thread's structures (needed by system)
pthread_attr_t joinable;       // attribute set for describing threads
int rc;                        // return code for threading commands

// With this condition variable the main thread signals to the worker threads that it has submitted a new job
pthread_cond_t new_job = PTHREAD_COND_INITIALIZER;

// Mutex assures exclusive access to bin_status[], jobs_sumitted, jobs_running,  and new_job by threads
pthread_mutex_t bin_status_mutex = PTHREAD_MUTEX_INITIALIZER;

// Mutex assures exclusive access to hitlist
pthread_mutex_t hitlist_mutex = PTHREAD_MUTEX_INITIALIZER;

// With this condition variable a worker thread signals to the main thread that it has finished a job
pthread_cond_t finished_job = PTHREAD_COND_INITIALIZER;
#endif


///////////////////////////////////////////////////////////////////////////////////////
//// For multi-threading: return a bin with the desired status, return -1 if no such bin found
//////////////////////////////////////////////////////////////////////////////////////
inline int PickBin(char status)
{
  for (int b=0; b<bins; b++) {if (bin_status[b]==status) return b;}
 return -1;
}

// Include hhworker.C here, because it needs some of the above variables
#include "hhworker.C"      // functions: AlignByWorker, RealignByWorker, WorkerLoop

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help(char all=0)
{
  printf("\n");
  printf("HHsearch %s\n",VERSION_AND_DATE);
  printf("Search a database of HMMs with a query alignment or query HMM\n");
  printf("%s",COPYRIGHT);
  printf("%s",HHSEARCH_REFERENCE);
  printf("\n");
  printf("Usage: %s -i query -d database [options]                       \n",program_name);
  printf(" -i <file>      input/query multiple sequence alignment (a2m, a3m, FASTA) or HMM\n");
  printf(" -d <file>      HMM database of concatenated HMMs in hhm, HMMER, or a3m format,\n");
  printf("                OR, if file has extension pal, list of HMM file names, one per\n");
  printf("                line. Multiple dbs, HMMs, or pal files with -d '<db1> <db2>...'\n");
  if (all) {
  printf("\n");
  printf("<file> may be 'stdin' or 'stdout' throughout.\n");
  }
  printf("\n");
  printf("Output options:                                                              \n");
  printf(" -o <file>      write results in standard format to file (default=<infile.hhr>)\n");
  if (all) {
  printf(" -Ofas <file>   write pairwise alignments of significant matches in FASTA format\n");
  printf("                Analogous for output in a3m, a2m, and psi format (e.g. -Oa3m)\n");
  printf(" -oa3m <file>   write MSA of significant matches in a3m format\n");
  printf("                Analogous for output in a2m, psi, and hhm format (e.g. -ohhm)\n");
  printf(" -e [0,1]       E-value cutoff for inclusion in multiple alignment (def=%G)    \n",par.e);
  printf(" -seq <int>     max. number of query/template sequences displayed (def=%i) \n",par.nseqdis);
  printf("                Beware of overflows! All these sequences are stored in memory.\n");
  printf(" -cons          show consensus sequence as master sequence of query MSA \n");
  }
  printf(" -nocons        don't show consensus sequence in alignments (default=show)     \n");
  printf(" -nopred        don't show predicted 2ndary structure in alignments (default=show)\n");
  printf(" -nodssp        don't show DSSP 2ndary structure in alignments (default=show)  \n");
  printf(" -ssconf        show confidences for predicted 2ndary structure in alignments\n");
  printf(" -p <float>     minimum probability in summary and alignment list (def=%G)   \n",par.p);
  printf(" -E <float>     maximum E-value in summary and alignment list (def=%G)       \n",par.E);
  printf(" -Z <int>       maximum number of lines in summary hit list (def=%i)         \n",par.Z);
  printf(" -z <int>       minimum number of lines in summary hit list (def=%i)         \n",par.z);
  printf(" -B <int>       maximum number of alignments in alignment list (def=%i)      \n",par.B);
  printf(" -b <int>       minimum number of alignments in alignment list (def=%i)      \n",par.b);
  if (all) {
  printf(" -aliw [40,..[  number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -dbstrlen      max length of database string to be printed in hhr file\n");
  }
  printf("\n");
  printf("Filter query multiple sequence alignment                                     \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[  filter MSA by selecting most diverse set of sequences, keeping \n");
  printf("                at least this many seqs in each MSA block of length 50 (def=%i) \n",par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",par.qsc);
  printf(" -neff [1,inf]  target diversity of alignment (default=off)\n");
  printf("\n");
  printf("Input alignment format:                                                       \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  if (all) {
  printf(" -tags          do NOT neutralize His-, C-myc-, FLAG-tags, and trypsin \n");
  printf("                recognition sequence to background distribution    \n");
  }
  printf("\n");
  printf("HMM-HMM alignment options:                                                    \n");
  printf(" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)   \n");
  printf(" -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                ness at alignment ends: 0:global  >0.1:local (default=%.2f)       \n",par.mact);
  printf(" -macins [0,1[  controls the cost of internal gap positions in the MAC algorithm.\n");
  printf("                0:dense alignments  1:gappy alignments (default=%.2f)\n",par.macins);
  
  printf(" -glob/-loc     use global/local alignment mode for searching/ranking (def=local)\n");
//   printf(" -vit          use Viterbi algorithm for searching/ranking (default)          \n");
//   printf(" -mac          use Maximum Accuracy MAC algorithm for searching/ranking\n");
//   printf(" -forward      use Forward probability for searching                       \n");
  printf(" -alt <int>     show up to this many significant alternative alignments(def=%i)\n",par.altali);
  if (all) {
  printf(" -vit           use Viterbi algorithm for searching/ranking (default)       \n");
  printf(" -mac           use Maximum Accuracy (MAC) algorithm for searching/ranking\n");
  printf(" -forward       use Forward probability for searching                       \n");
  printf(" -excl <range>  exclude query positions from the alignment, e.g. '1-33,97-168' \n");
  printf(" -shift [-1,1]  score offset (def=%-.2f)                                       \n",par.shift);
  printf(" -corr [0,1]    weight of term for pair correlations (def=%.2f)                \n",par.corr);
  printf(" -sc   <int>    amino acid score         (tja: template HMM at column j) (def=%i)\n",par.columnscore);
  printf("        0       = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
  printf("        1       = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
  printf("        2       = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
  printf("        3       = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
  printf("        5       local amino acid composition correction                     \n");
  }
  printf(" -ssm {0,..,4}  0:   no ss scoring                                             \n");
  printf("                1,2: ss scoring after or during alignment  [default=%1i]       \n",par.ssm);
  printf("                3,4: ss scoring after or during alignment, predicted vs. predicted \n");
  if (all) {
  printf(" -ssw  [0,1]    weight of ss score compared to column score (def=%-.2f)     \n",par.ssw);
  printf(" -ssa  [0,1]    SS substitution matrix = (1-ssa)*I + ssa*full-SS-substition-matrix [def=%-.2f)\n",par.ssa);
  printf("\n");
  printf("Gap cost options:                                                                      \n");
  printf(" -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                           \n",par.gapb);
  printf(" -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)          \n",par.gapd);
  printf(" -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)            \n",par.gape);
  printf(" -gapf ]0,inf]  factor to increase/reduce the gap open penalty for deletes (def=%-.2f) \n",par.gapf);
  printf(" -gapg ]0,inf]  factor to increase/reduce the gap open penalty for inserts (def=%-.2f) \n",par.gapg);
  printf(" -gaph ]0,inf]  factor to increase/reduce the gap extend penalty for deletes(def=%-.2f)\n",par.gaph);
  printf(" -gapi ]0,inf]  factor to increase/reduce the gap extend penalty for inserts(def=%-.2f)\n",par.gapi);
  printf(" -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f)      \n",par.egq);
  printf(" -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)   \n",par.egt);
  printf("\n");

  printf("Pseudocount (pc) options:                                                        \n");
  printf(" Context specific hhm pseudocounts:\n");
  printf("  -pc_hhm_contxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc_hhm_context_engine.admix);
  printf("               0: no pseudo counts:    tau = 0                                  \n");
  printf("               1: constant             tau = a                                  \n");
  printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
  printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
  printf("  -pc_hhm_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc_hhm_context_engine.pca);
  printf("  -pc_hhm_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",par.pc_hhm_context_engine.pcb);
  printf("  -pc_hhm_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",par.pc_hhm_context_engine.pcc);

  printf(" Context independent hhm pseudocounts (used for templates; used for query if contxt file is not available):\n");
  printf("  -pc_hhm_nocontxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc_hhm_nocontext_mode);
  printf("               0: no pseudo counts:    tau = 0                                  \n");
  printf("               1: constant             tau = a                                  \n");
  printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
//  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
  printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
  printf("  -pc_hhm_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc_hhm_nocontext_a);
  printf("  -pc_hhm_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",par.pc_hhm_nocontext_b);
  printf("  -pc_hhm_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",par.pc_hhm_nocontext_c);
  }

  printf("\n");
  printf("Context-specific pseudo-counts:                                                  \n");
  printf(" -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
  printf(" -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",par.clusterfile);

  if (all) {
  printf(" -csw  [0,inf]  weight of central position in cs pseudocount mode (def=%.1f)\n", par.csw);
  printf(" -csb  [0,1]    weight decay parameter for positions in cs pc mode (def=%.1f)\n", par.csb);
  }
  printf("\n");
  printf("Other options: \n");
  printf(" -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=1)\n");
#ifndef PTHREAD
  printf("(The -cpu option is inactive since POSIX threads ae not supported on your platform)\n");
#endif
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose   \n");
  if (all) {
  printf(" -maxres <int>  max number of HMM columns (def=%5i)             \n",par.maxres);
  printf(" -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n",par.maxmem);
  printf(" -scores <file> write scores for all pairwise comparisions to file         \n");
  printf(" -calm {0,..,3} empirical score calibration of 0:query 1:template 2:both   \n");
  printf("                default 3: neural network-based estimation of EVD params   \n");
  // printf(" -opt  <file>   parameter optimization mode (def=off): return sum of ranks \n");
  // printf("                of true positives (same superfamily) for minimization      \n");
  // printf("                and write result into file                                 \n");
  printf("\n");
  } else {
  printf("An extended list of options can be obtained by calling 'hhblits -help'\n");
  }
  printf("\n");
  printf("Example: %s -i a.1.1.1.a3m -d scop70_1.71.hhm \n",program_name);
  cout<<endl;

//   printf("More help:                                                         \n");
//   printf(" -h out        options for formatting ouput                        \n");
//   printf(" -h hmm        options for building HMM from multiple alignment    \n");
//   printf(" -h gap        options for setting gap penalties                   \n");
//   printf(" -h ali        options for HMM-HMM alignment                       \n");
//   printf(" -h other      other options \n");
//   printf(" -h all        all options \n");
 }




/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(int argc, char** argv)
{
  //Processing command line input
  for (int i=1; i<argc; i++)
    {
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-i"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -i\n"; exit(4);}
          else strcpy(par.infile,argv[i]);
        }
      else if (!strcmp(argv[i],"-d"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database file following -d\n"; exit(4);}
          else
            {
              par.dbfiles = new(char[strlen(argv[i])+1]);
              strcpy(par.dbfiles,argv[i]);
            }
        }
      else if (!strcmp(argv[i],"-o"))
        {
          par.append=0;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.outfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-Ofas"))
        {
          par.append=0;
          par.outformat=1;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.pairwisealisfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-Oa2m"))
        {
          par.append=0;
          par.outformat=2;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.pairwisealisfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-Oa3m"))
        {
          par.append=0;
          par.outformat=3;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.pairwisealisfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-oa3m"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Oa3m\n"; exit(4);}
          else strcpy(par.alnfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-ohhm"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -ohhm\n"; exit(4);}
          else strcpy(par.hhmfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-opsi"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Opsi\n"; exit(4);}
          else strcpy(par.psifile,argv[i]);
        }
      else if (!strcmp(argv[i],"-scores"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file following -scores\n"; exit(4);}
          else {strcpy(par.scorefile,argv[i]);}
        }
     else if (!strcmp(argv[i],"-atab") || !strcmp(argv[i],"-Aliout"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no query file following -atab\n"; exit(4);}
	  else strmcpy(par.alitabfile,argv[i],NAMELEN-1);
	}
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"-help"))
        { help(1); exit(0); }
      else if (!strcmp(argv[i],"-excl"))
        {
          if (++i>=argc) {help(); exit(4);}
          par.exclstr = new(char[strlen(argv[i])+1]);
          strcpy(par.exclstr,argv[i]);
        }
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-p") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-P") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-E") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-b") && (i<argc-1)) par.b = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-B") && (i<argc-1)) par.B = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-z") && (i<argc-1)) par.z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-Z") && (i<argc-1)) par.Z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-e") && (i<argc-1)) par.e = atof(argv[++i]);
      else if (!strncmp(argv[i],"-nocons",7)) par.showcons=0;
      else if (!strncmp(argv[i],"-nopred",7)) par.showpred=0;
      else if (!strncmp(argv[i],"-nodssp",7)) par.showdssp=0;
      else if (!strncmp(argv[i],"-ssconf",7)) par.showconf=1;
      else if (!strncmp(argv[i],"-cons",5)) par.cons=1;
      else if (!strncmp(argv[i],"-mark",5)) par.mark=1;
      else if (!strcmp(argv[i],"-seq") && (i<argc-1))  par.nseqdis=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-aliw") && (i<argc-1)) par.aliwidth=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-neff") && (i<argc-1)) par.Neff=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-Neff") && (i<argc-1)) par.Neff=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-Gonnet")) par.matrix=0;
      else if (!strncmp(argv[i],"-BLOSUM",7) || !strncmp(argv[i],"-Blosum",7))
        {
          if (!strcmp(argv[i]+7,"30")) par.matrix=30;
          else if (!strcmp(argv[i]+7,"40")) par.matrix=40;
          else if (!strcmp(argv[i]+7,"50")) par.matrix=50;
          else if (!strcmp(argv[i]+7,"65")) par.matrix=65;
          else if (!strcmp(argv[i]+7,"80")) par.matrix=80;
          else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
        }
      else if (!strcmp(argv[i],"-wg")) {par.wg=1;}
      else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pc_hhm_context_engine.admix=(Pseudocounts::Admix)atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pc_hhm_context_engine.pca=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pc_hhm_context_engine.pcb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pc_hhm_context_engine.pcc=atof(argv[++i]);
      //else if (!strcmp(argv[i],"-pcw") && (i<argc-1)) par.pcw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;}
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egq") && (i<argc-1)) par.egq=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egt") && (i<argc-1)) par.egt=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssgap")) par.ssgap=1;
      else if (!strcmp(argv[i],"-ssgapd") && (i<argc-1)) par.ssgapd=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssgape") && (i<argc-1)) par.ssgape=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssgapi") && (i<argc-1)) par.ssgapi=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-ssm") && (i<argc-1)) par.ssm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-ssw") && (i<argc-1)) par.ssw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssw_mac") && (i<argc-1)) par.ssw_realign=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssa") && (i<argc-1)) par.ssa=atof(argv[++i]);
      else if (!strcmp(argv[i],"-realign")) par.realign=1;
      else if (!strcmp(argv[i],"-norealign")) par.realign=0;
      else if (!strcmp(argv[i],"-mac") || !strcmp(argv[i],"-MAC")) par.forward=2;
      else if (!strcmp(argv[i],"-map") || !strcmp(argv[i],"-MAP")) par.forward=2;
      else if (!strcmp(argv[i],"-vit")) par.forward=0;
      else if (!strncmp(argv[i],"-glo",3)) {par.loc=0; if (par.mact>0.35 && par.mact<0.3502) {par.mact=0;} }
      else if (!strncmp(argv[i],"-loc",4)) par.loc=1;
      else if (!strncmp(argv[i],"-alt",4) && (i<argc-1)) par.altali=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-M") && (i<argc-1))
        if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1;
        else if(!strcmp(argv[i],"first"))  par.M=3;
        else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
        else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strcmp(argv[i],"-cal")) {par.calibrate=1; par.calm=0;}
      else if (!strcmp(argv[i],"-calm") && (i<argc-1)) par.calm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-shift") && (i<argc-1)) par.shift=atof(argv[++i]);
      else if ((!strcmp(argv[i],"-mact") || !strcmp(argv[i],"-mapt")) && (i<argc-1)) par.mact=atof(argv[++i]);
      else if (!strcmp(argv[i],"-macins") && (i<argc-1)) par.macins=atof(argv[++i]);
      else if (!strcmp(argv[i],"-sc") && (i<argc-1)) par.columnscore=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-scwin") && (i<argc-1)) {par.columnscore=5; par.half_window_size_local_aa_bg_freqs = imax(1,atoi(argv[++i]));}
//      TODO: is this necessary?
//      else if (!strcmp(argv[i],"-def")) ;
      else if (!strcmp(argv[i],"-maxres") && (i<argc-1)) {
	par.maxres=atoi(argv[++i]);
	par.maxcol=2*par.maxres;
      }
      else if (!strncmp(argv[i],"-cpu",4) && (i<argc-1)) threads=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-maxmem") && (i<argc-1)) {par.maxmem=atof(argv[++i]);}
      else if (!strcmp(argv[i],"-corr") && (i<argc-1)) par.corr=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ovlp") && (i<argc-1)) par.min_overlap=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-dbstrlen") && (i<argc-1)) par.maxdbstrlen=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-tags")) par.notags=0;
      else if (!strcmp(argv[i],"-notags")) par.notags=1;
      else if (!strncmp(argv[i],"-idummy",7) && (i<argc-1)) par.idummy=atoi(argv[++i]);
      else if (!strncmp(argv[i],"-premerge",9) && (i<argc-1)) par.premerge=atoi(argv[++i]);
      else if (!strncmp(argv[i],"-fdummy",7) && (i<argc-1)) par.fdummy=atof(argv[++i]);
      else if (!strcmp(argv[i],"-nocontxt")) par.nocontxt=1;
      else if (!strcmp(argv[i],"-csb") && (i<argc-1)) par.csb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-csw") && (i<argc-1)) par.csw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-contxt") || !strcmp(argv[i],"-cs"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -contxt\n"; exit(4);}
          else strcpy(par.clusterfile,argv[i]);
        }
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Realign hits with MAC algorithm
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void perform_realign(char *dbfiles[], int ndb)
{
  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
  int nhits=0;
  int N_aligned=0;
  
  // Longest allowable length of database HMM (backtrace: 5 chars, fwd, bwd: 1 double
  long int Lmaxmem=(par.maxmem*1024*1024*1024)/sizeof(double)/q->L/bins;
  long int Lmax=0;      // length of longest HMM to be realigned
    
  // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
  // This list can be sorted by ftellpos to access one template after the other efficiently during realignment
  Hash< List<Realign_hitpos>* >* phash_plist_realignhitpos;
  phash_plist_realignhitpos = new Hash< List<Realign_hitpos>* > (2311,NULL);
  
  // Some templates have several (suboptimal) alignments in hitlist. For realignment, we need to efficiently 
  // access all hit objects in hitlist belonging to one template (because we don't want to read templates twice)
  // We therefore need for each template (identified by its index between 0 and N_searched-1) a list of elements 
  // in hitlist that store the alignments with the template of that index. 
  // This list is pointed to by array_plist_phits[index].
  List<void*>** array_plist_phits; 
  array_plist_phits = new List<void*>*[N_searched];
  for (int index=0; index<N_searched; index++) array_plist_phits[index] = NULL; // initialize 
  
  // Store all dbfiles and ftell positions of templates to be displayed and realigned
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      if (nhits >= par.realign_max && nhits>=imax(par.B,par.Z)) break;
      if (hit_cur.Eval > par.e)
	{
	  if (nhits>=imax(par.B,par.Z)) continue;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) continue;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
	}

      if (hit_cur.L>Lmax) Lmax=hit_cur.L;
      if (hit_cur.L>Lmaxmem) {nhits++; continue;} // skip HMMs that require too much memory to be realigned

//    fprintf(stderr,"hit.name=%-15.15s  hit.index=%-5i hit.ftellpos=%-8i  hit.dbfile=%s\n",hit_cur.name,hit_cur.index,(unsigned int)hit_cur.ftellpos,hit_cur.dbfile);
      if (nhits>=par.premerge || hit_cur.irep>1) // realign the first premerge hits consecutively to query profile
	{
	  if (hit_cur.irep==1) 
	    {
	      // For each template (therefore irep==1), store template index and position on disk in a list
	      Realign_hitpos realign_hitpos;
	      realign_hitpos.ftellpos = hit_cur.ftellpos;    // stores position on disk of template for current hit
	      realign_hitpos.index = hit_cur.index;          // stores index of template of current hit
	      if (! phash_plist_realignhitpos->Contains(hit_cur.dbfile))
		{
		  List<Realign_hitpos>* newlist = new List<Realign_hitpos>;
		  phash_plist_realignhitpos->Add(hit_cur.dbfile,newlist);
		}
	      // Add template index and ftellpos to list which belongs to key dbfile in hash
	      phash_plist_realignhitpos->Show(hit_cur.dbfile)->Push(realign_hitpos);
	    }
	  if (! array_plist_phits[hit_cur.index]) // pointer at index is still NULL  
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
  if (Lmax>Lmaxmem)
    {
      Lmax=Lmaxmem;
      if (v>=1) 
	{
	  cerr<<"WARNING: Realigning sequences only up to length "<<Lmaxmem<<"."<<endl;
	  cerr<<"This is genarally unproboblematic but may lead to slightly sub-optimal alignments for these sequences."<<endl;
 	  cerr<<"You can increase available memory using the -maxmem <GB> option (currently "<<par.maxmem<<" GB)."<<endl; 
	  cerr<<"The maximum length realignable is approximately maxmem/query_length/(cpus+1)/8B."<<endl;
	}
    }
  
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs=1;   // needs to be set to 1 before threads are created
  for (bin=0; bin<bins; bin++)
    if (!hit[bin]->forward_allocated)  hit[bin]->AllocateForwardMatrix(q->L+2,Lmax+1);
    
  if (v>=2) printf("Realigning %i database HMMs using HMM-HMM Maximum Accuracy algorithm\n",nhits);
  v1 = v;
  if (v>0 && v<=3) v=1; else v-=2;  // Supress verbose output during iterative realignment and realignment
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Start premerge: align the first par.premerge templates?
  if (par.premerge>0)
    {
      
      // Read query alignment into Qali
      Alignment Qali;  // output A3M generated by merging A3M alignments for significant hits to the query alignment
      char qa3mfile[NAMELEN];
      char ta3mfile[NAMELEN];
      RemoveExtension(qa3mfile,par.infile); // directory??
      strcat(qa3mfile,".a3m");
      FILE* qa3mf=fopen(qa3mfile,"r");
      if (!qa3mf) OpenFileError(qa3mfile);
      Qali.Read(qa3mf,qa3mfile);
      fclose(qa3mf);
      delete[] Qali.longname;
      Qali.longname = new(char[strlen(q->longname)+1]);
      strcpy(Qali.longname,q->longname);
      strcpy(Qali.name,q->name);
      strcpy(Qali.fam,q->fam);
      RemovePathAndExtension(Qali.file,par.hhmfile);
      
      // If par.append==1 do not print query alignment
      if (par.append) Qali.MarkSeqsAsNonPrintable();
      
      if (v>=2) printf("Merging best hits to query alignment %s ...\n",qa3mfile);
      
      bin=0;
      nhits=0;
      hitlist.Reset();
      while (!hitlist.End() && nhits<par.premerge)
	{
	  hit_cur = hitlist.ReadNext();
	  if (nhits>=imax(par.B,par.Z)) break;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
	  
	  // if (hit_cur.irep>1) continue; // Align only the best hit of the first par.premerge templates // JS 13 Feb 13: commented out since this could lead to problems with hits that are then not realigned at all and missing posterior probs => remove entirely?
	 
	  if (hit_cur.L>Lmaxmem) {nhits++; continue;} //Don't align to long sequences due to memory limit
	  
	  // Open HMM database file dbfiles[idb]
	  FILE* dbf=fopen(hit_cur.dbfile,"rb");
	  if (!dbf) OpenFileError(hit_cur.dbfile);
	  read_from_db=1;
	  
	  // Forward stream position to start of next database HMM to be realigned
	  hit[bin]->index = hit_cur.index; // give hit a unique index for HMM
	  hit[bin]->ftellpos = hit_cur.ftellpos;
	  fseek(dbf,hit_cur.ftellpos,SEEK_SET);
	  hit[bin]->dbfile = new(char[strlen(hit_cur.dbfile)+1]);
	  strcpy(hit[bin]->dbfile,hit_cur.dbfile); // record db file name from which next HMM is read
	  hit[bin]->irep = 1;  // Needed for min_overlap calculation in InitializeForAlignment in hhhit.C
	  
	  char path[NAMELEN];
	  Pathname(path,hit_cur.dbfile);
	  
	  ///////////////////////////////////////////////////
	  // Read next HMM from database file
	  if (!fgetline(line,LINELEN,dbf)) {fprintf(stderr,"Error in %s: end of file %s reached prematurely!\n",par.argv[0],hit_cur.dbfile); exit(1);}
	  while (strscn(line)==NULL && fgetline(line,LINELEN,dbf)) {} // skip lines that contain only white space

	  if (!strncmp(line,"HMMER3",6))      // read HMMER3 format
	    {
	      format[bin] = 1;
	      read_from_db = t[bin]->ReadHMMer3(dbf,hit_cur.dbfile);
	      par.hmmer_used = true;
	    }
	  else if (!strncmp(line,"HMMER",5))      // read HMMER format
	    {
	      format[bin] = 1;
	      read_from_db = t[bin]->ReadHMMer(dbf,hit_cur.dbfile);
	      par.hmmer_used = true;
	    }
	  else if (!strncmp(line,"HH",2))     // read HHM format
	    {
	      format[bin] = 0;
	      read_from_db = t[bin]->Read(dbf,path);
	    }
	  else if (!strncmp(line,"NAME",4))  // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
	    {
	      format[bin] = 0;
	      fseek(dbf,hit_cur.ftellpos,SEEK_SET); // rewind to beginning of line
	      read_from_db = t[bin]->Read(dbf,path);
	    }
	  else if (line[0]=='#')             // read a3m alignment
	    {
	      Alignment tali;
	      tali.Read(dbf,hit_cur.dbfile,line);
	      tali.Compress(hit_cur.dbfile);
//            qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
	      tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
	      t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
	      tali.FrequenciesAndTransitions(t[bin]);
	      format[bin] = 0;
	    }
	  else {
	    cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<hit_cur.dbfile<<"\'. \n";
	    cerr<<"Context:\n'"<<line<<"\n";
	    fgetline(line,LINELEN,dbf); cerr<<line<<"\n";
	    fgetline(line,LINELEN,dbf); cerr<<line<<"'\n";
	    exit(1);
	  }
	  fclose(dbf);
	  
	  if (read_from_db!=1)
	    {
	      fprintf(stderr,"Error in %s: illegal format in %s while reading HMM %s\n",par.argv[0],hit_cur.dbfile,hit_cur.name);
	      continue;
	    }
	  if (v>=2) fprintf(stderr,"Realigning with %s ***** \n",t[bin]->name);
	  ///////////////////////////////////////////////////
	  
	  N_aligned++;
	  if (v1>=2 && !(N_aligned%10))
	    {
	      cout<<".";
	      if (!(N_aligned%500)) printf(" %-4i HMMs aligned\n",N_aligned);
	      cout.flush();
	    }
	  
	  // Prepare MAC comparison(s)
	  PrepareTemplateHMM(q,t[bin],format[bin]);
	  t[bin]->Log2LinTransitionProbs(1.0);
	  
	  // Align q to template in *hit[bin]
	  hit[bin]->Forward(q,t[bin]);
	  hit[bin]->Backward(q,t[bin]);
	  hit[bin]->MACAlignment(q,t[bin]);
	  hit[bin]->BacktraceMAC(q,t[bin]);
	  
	  // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
	  hit[bin]->score      = hit_cur.score;
	  hit[bin]->score_aass = hit_cur.score_aass;
	  hit[bin]->score_ss   = hit_cur.score_ss;
	  hit[bin]->Pval       = hit_cur.Pval;
	  hit[bin]->Pvalt      = hit_cur.Pvalt;
	  hit[bin]->logPval    = hit_cur.logPval;
	  hit[bin]->logPvalt   = hit_cur.logPvalt;
	  hit[bin]->Eval       = hit_cur.Eval;
	  hit[bin]->logEval    = hit_cur.logEval;
	  hit[bin]->Probab     = hit_cur.Probab;
	  hit[bin]->irep       = hit_cur.irep;
	  
	  // Replace original hit in hitlist with realigned hit
	  //hitlist.ReadCurrent().Delete();
	  hitlist.Delete().Delete();               // delete the list record and hit object
	  hitlist.Insert(*hit[bin]);
	  
	  // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
	  strcpy(ta3mfile,hit[bin]->file); // copy filename including path but without extension
	  strcat(ta3mfile,".a3m");
	  Alignment Tali;
	  FILE* ta3mf=fopen(ta3mfile,"r");
	  if (!ta3mf) OpenFileError(ta3mfile);
	  Tali.Read(ta3mf,ta3mfile); // Read template alignment into Tali
	  fclose(ta3mf);
	  Tali.Compress(ta3mfile); // Filter database alignment
	  Qali.MergeMasterSlave(*hit[bin],Tali,ta3mfile);
	  
	  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
	  Qali.Compress("merged A3M file");
	  
	  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
	  Qali.N_filtered = Qali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
	  
	  // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
	  Qali.FrequenciesAndTransitions(q);
	  
	  // Comput substitution matrix pseudocounts
	  if (par.nocontxt) { 
	    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	    q->PreparePseudocounts();
	    // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	    q->AddAminoAcidPseudocounts();
	  } else {
	    // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
	    // q->PrepareContextSpecificPseudocounts(); //OLD
	    q->AddContextSpecificPseudocounts(pc_hhm_context_engine, pc_hhm_context_mode);
	  }
	  
	  // Transform transition freqs to lin space if not already done
	  q->AddTransitionPseudocounts();
	  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
	  
	  q->CalculateAminoAcidBackground();

	  nhits++;
	}
    }
  // End premerge: align the first par.premerge templates?
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
#ifdef PTHREAD
  // Start threads for realignment
  for (int j=0; j<threads; j++)
    {
      thread_data[j].thread_id = j+1;
      thread_data[j].function  = &RealignByWorker;
      if (DEBUG_THREADS) fprintf(stderr,"Creating worker thread %i ...",j+1);
      pthread_create(&pthread[j], &joinable, WorkerLoop, (void*)&thread_data[j]);
      if (DEBUG_THREADS) fprintf(stderr," created!\n");
    }
#endif
  
  // Read all HMMs whose position is given in phash_plist_realignhitpos
  for (int idb=0; idb<ndb; idb++)
    {
      
      // Can we skip dbfiles[idb] because it contains no template to be realigned?
      if (! phash_plist_realignhitpos->Contains(dbfiles[idb])) continue;
      
      // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
      // This list is now sorted by ftellpos in ascending order to access one template after the other efficiently
      phash_plist_realignhitpos->Show(dbfiles[idb])->SortList();
      
      // Open HMM database file dbfiles[idb]
      FILE* dbf=fopen(dbfiles[idb],"rb");
      if (!dbf) OpenFileError(dbfiles[ndb]);
      read_from_db=1;
      
      ///////////////////////////////////////////////////////////////////////////////////////
      // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
      phash_plist_realignhitpos->Show(dbfiles[idb])->Reset();
      while (! phash_plist_realignhitpos->Show(dbfiles[idb])->End())
	{
	  // Submit jobs until no bin is free anymore
	  while (! phash_plist_realignhitpos->Show(dbfiles[idb])->End() && jobs_submitted+jobs_running<bins)
	    {
	      
	      // Allocate free bin
	      bin = PickBin(FREE);
	      if (bin<0) {
		fprintf(stderr,"Error during realignment: found no free bin! jobs running: %i  jobs_submitted:%i  threads:%i\n",jobs_running,jobs_submitted,threads);
		for (bin=0; bin<bins; bin++) fprintf(stderr,"bin_status[%i]=%i\n",bin,bin_status[bin]);
		exit(6);
	      }
	      
	      // Forward stream position to start of next database HMM to be realigned
	      Realign_hitpos hitpos_curr = phash_plist_realignhitpos->Show(dbfiles[idb])->ReadNext();
	      hit[bin]->index = hitpos_curr.index;      // give hit[bin] a unique index for HMM
	      fseek(dbf,hitpos_curr.ftellpos,SEEK_SET); // start to read at ftellpos for template
	      
	      // Give hit[bin] the pointer to the list of pointers to hitlist elements of same template (for realignment)
	      hit[bin]->plist_phits = array_plist_phits[hitpos_curr.index];
	      
	      // fprintf(stderr,"dbfile=%-40.40s  index=%-5i  ftellpos=%l\n",dbfiles[idb],hitpos_curr.index,hitpos_curr.ftellpos);
	      
	      char path[NAMELEN];
	      Pathname(path,dbfiles[idb]);
	      
	      ///////////////////////////////////////////////////
	      // Read next HMM from database file
	      if (!fgetline(line,LINELEN,dbf)) 
		{fprintf(stderr,"Error in %s: end of file %s reached prematurely!\n",par.argv[0],dbfiles[idb]); exit(1);}
	      while (strscn(line)==NULL && fgetline(line,LINELEN,dbf)) {} // skip lines that contain only white space

	      if (!strncmp(line,"HMMER3",6))      // read HMMER3 format
		{
		  format[bin] = 1;
		  read_from_db = t[bin]->ReadHMMer3(dbf,dbfiles[idb]);
		  par.hmmer_used = true;
		}
	      else if (!strncmp(line,"HMMER",5))      // read HMMER format
		{
		  format[bin] = 1;
		  read_from_db = t[bin]->ReadHMMer(dbf,dbfiles[idb]);
		  par.hmmer_used = true;
		}
	      else if (!strncmp(line,"HH",2))     // read HHM format
		{
		  format[bin] = 0;
		  read_from_db = t[bin]->Read(dbf,path);
		}
	      else if (!strncmp(line,"NAME",4))  // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
		{
		  format[bin] = 0;
		  fseek(dbf,hitpos_curr.ftellpos,SEEK_SET); // rewind to beginning of line
		  read_from_db = t[bin]->Read(dbf,path);
		}
	      else if (line[0]=='#')                 // read a3m alignment
		{
		  Alignment tali;
		  tali.Read(dbf,dbfiles[idb],line);
		  tali.Compress(dbfiles[idb]);
		  // qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
		  tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
		  t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
		  tali.FrequenciesAndTransitions(t[bin]);
		  format[bin] = 0;
		}
	      else {
		cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<dbfiles[idb]<<"\'. \n";
		cerr<<"Context:\n'"<<line<<"\n";
		fgetline(line,LINELEN,dbf); cerr<<line<<"\n";
		fgetline(line,LINELEN,dbf); cerr<<line<<"'\n";
		exit(1);

	      }
	      
	      
	      if (read_from_db==2) continue;  // skip current HMM or reached end of database
	      if (read_from_db==0) break;     // finished reading HMMs
	      if (v>=2) fprintf(stderr,"Realigning with %s\n",t[bin]->name);
	      ///////////////////////////////////////////////////
	      
	      hit[bin]->dbfile = new(char[strlen(dbfiles[idb])+1]);
	      strcpy(hit[bin]->dbfile,dbfiles[idb]); // record db file name from which next HMM is read
	      
	      N_aligned++;
	      if (v1>=2 && !(N_aligned%10))
		{
		  cout<<".";
		  if (!(N_aligned%500)) printf(" %-4i HMMs aligned\n",N_aligned);
		  cout.flush();
		}
#ifdef PTHREAD
	      // Lock access to bin_status
	      if (threads>0) rc = pthread_mutex_lock(&bin_status_mutex);
#endif
	      // Send the job in bin to a thread
	      bin_status[bin] = SUBMITTED;
	      jobs_submitted++;
	      
	      if (threads==0) // if no multi-threading mode, main thread executes job itself
		{
		  RealignByWorker(bin);
		  bin_status[bin] = FREE;
		  jobs_submitted--;
		  break;
		}
	      
#ifdef PTHREAD
	      // Restart threads waiting for a signal
	      rc = pthread_cond_signal(&new_job);
	      
	      // Unlock access to bin_status
	      rc = pthread_mutex_unlock(&bin_status_mutex);
	      
	      if (DEBUG_THREADS)
		fprintf(stderr,"Main: put job into bin %i  name=%-7.7s  Jobs running: %i  jobs_submitted:%i \n",bin,t[bin]->name,jobs_running,jobs_submitted);
#endif
	    }
	  
#ifdef PTHREAD
	  if (threads>0)
	    {
	      // Lock mutex
	      rc = pthread_mutex_lock(&bin_status_mutex);
	      
	      // Wait until job finishes and a bin becomes free
	      if (jobs_submitted+jobs_running>=bins)
		{
		  if (DEBUG_THREADS) fprintf(stderr,"Master thread is waiting for jobs to be finished...\n");
#ifdef WINDOWS
		  rc = pthread_cond_wait(&finished_job, &bin_status_mutex);
#else
#ifdef HH_MAC
		  
		  // If no submitted jobs are in the queue we have to wait for a new job ...
		  struct timespec ts;
		  struct timeval tv;
		  gettimeofday(&tv, NULL);
		  ts.tv_sec = tv.tv_sec + 1;
		  rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex,&ts);
#else
		  
		  // If no submitted jobs are in the queue we have to wait for a new job ...
		  struct timespec ts;
		  clock_gettime(CLOCK_REALTIME,&ts);
		  ts.tv_sec += 1;
		  rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex,&ts);
#endif
#endif
		}
	      // Unlock mutex
	      rc = pthread_mutex_unlock(&bin_status_mutex);
	    }
#endif
	  
	}
      // End while(1)
      ///////////////////////////////////////////////////////////////////////////////////////
      
      fclose(dbf);
    }
  reading_dbs=0;
  
#ifdef PTHREAD
  if (threads>0)
    {
      // No more HMMs in database => wait until all threads have finished
      if (DEBUG_THREADS)
	fprintf(stderr,"No more jobs read from database         Jobs running:%i  jobs_submitted:%i \n",jobs_running,jobs_submitted);
      
      // Free all threads still waiting for jobs
      rc = pthread_mutex_lock(&bin_status_mutex);
      rc = pthread_cond_broadcast(&new_job);
      rc = pthread_mutex_unlock(&bin_status_mutex); // Unlock mutex
      
      // Wait for all jobs to finish => join all jobs to main
      for (int j=0; j<threads; j++)
	{
	  int status;
	  pthread_join(pthread[j], (void **)&status);
	  if (DEBUG_THREADS) fprintf(stderr,"Thread %i finished its work\n",j+1);
	}
    }
#endif
  if (v1>=2) cout<<"\n";
  v=v1;

  
  // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
  if (*par.scorefile) {
    if (v>=3) printf("Printing scores file ...\n");
    hitlist.PrintScoreFile(q);
  }
  
  
  // Delete all hitlist entries with too short alignments
  nhits=0;
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      //        printf("Deleting alignment of %s with length %i? nhits=%-2i  par.B=%-3i  par.Z=%-3i par.e=%.2g par.b=%-3i  par.z=%-3i par.p=%.2g\n",hit_cur.name,hit_cur.matched_cols,nhits,par.B,par.Z,par.e,par.b,par.z,par.p);
      if (nhits>=imax(par.B,par.Z)) break;
      if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
      if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
      if (hit_cur.matched_cols < MINCOLS_REALIGN)
	{
	  if (v>=3) printf("Deleting alignment of %s with length %i\n",hit_cur.name,hit_cur.matched_cols);
	  hitlist.Delete().Delete();               // delete the list record and hit object
	  // // Make sure only realigned alignments get displayed! JS: Why? better unrealigned than none.
	  // if (par.B>par.Z) par.B--; else if (par.B==par.Z) {par.B--; par.Z--;} else par.Z--;
	  // if (par.b>par.z) par.b--; else if (par.b==par.z) {par.b--; par.z--;} else par.z--;
	}
      nhits++;
    }
  
  // Delete hash phash_plist_realignhitpos with lists
  phash_plist_realignhitpos->Reset();
  while (!phash_plist_realignhitpos->End())
    delete(phash_plist_realignhitpos->ReadNext()); // delete list to which phash_plist_realignhitpos->ReadNext() points
  delete(phash_plist_realignhitpos);
  
  // Delete array_plist_phits with lists
  for (int index=0; index<N_searched; index++) 
    if (array_plist_phits[index]) delete(array_plist_phits[index]); // delete list to which array[index] points
  delete[](array_plist_phits); 
}



/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  char* argv_conf[MAXOPT];       // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf=0;               // Number of arguments in argv_conf
  char inext[IDLEN]="";          // Extension of input file ich 5 Leute, die an Transkription arbeiten. (hhm or a3m)
  const char print_elapsed=0;

#ifdef PTHREAD
  pthread_attr_init(&joinable);  // initialize attribute set with default values
  if (pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)!=0) // set attribute 'joinable'
    cerr<<"Error "<<pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)<<": could not set detach state for thread attibute.\n";
#endif

  // Make command line input globally available
  par.argv=argv;
  par.argc=argc;
  RemovePathAndExtension(program_name,argv[0]);
  Pathname(program_path, argv[0]);

  // Enable changing verbose mode before defaults file and command line are processed
  for (int i=1; i<argc; i++)
    {
      if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1;
      else if (argc>1 && !strcmp(argv[i],"-v0")) v=0;
      else if (argc>1 && !strcmp(argv[i],"-v1")) v=1;
      else if (argc>2 && !strcmp(argv[i],"-v")) v=atoi(argv[i+1]);
    }

  par.SetDefaultPaths(program_path);

  // Read .hhdefaults file?
  if (par.readdefaultsfile)
    {
      // Process default otpions from .hhconfig file
      ReadDefaultsFile(argc_conf,argv_conf,program_path);
      ProcessArguments(argc_conf,argv_conf);
    }

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc,argv);

  // Check command line input and default values
  if (!*par.infile || !strcmp(par.infile,"")) // infile not given
    {help(); cerr<<endl<<"Error in "<<program_name<<": input file missing!\n"; exit(4);}
  if (!par.dbfiles) // pointer never set?
    {help(); cerr<<endl<<"Error in "<<program_name<<": no HMM database file given (-d file)\n"; exit(4);}
  if (threads<=1) threads=0;
  else if (threads>MAXTHREADS)
    {
      threads=MAXTHREADS;
      if (v>=1) fprintf(stderr,"WARNING: number of CPUs set to maximum value of %i\n",MAXTHREADS);
    }
  RemoveExtension(q->file,par.infile); // get rootname of infile (no directory path, no extension)
  Extension(inext,par.infile);        // get extension of infile
  if (!*par.outfile)      // outfile not given? Name it basename.hhm
    {
      RemoveExtension(par.outfile,par.infile);
      strcat(par.outfile,".hhr");
      if (v>=2) cout<<"Search results will be written to "<<par.outfile<<"\n";
    }

  // Check option compatibilities
  if (par.nseqdis>MAXSEQDIS-3-par.showcons) par.nseqdis=MAXSEQDIS-3-par.showcons; //3 reserved for secondary structure
  if (par.aliwidth<20) par.aliwidth=20;
  if (par.pc_hhm_context_engine.pca<0.001) par.pc_hhm_context_engine.pca=0.001; // to avoid log(0)
  if (par.b>par.B) par.B=par.b;
  if (par.z>par.Z) par.Z=par.z;
  if (par.maxmem<1.0) {cerr<<"Warning: setting -maxmem to its minimum allowed value of 1.0\n"; par.maxmem=1.0;}
  if (par.mact>=1.0) par.mact=0.999; else if (par.mact<0) par.mact=0.0;
  if (par.macins>=1.0) par.macins=0.999; else if (par.macins<0) par.macins=0.0;

  // Input parameters
  if (v>=3)
    {
      cout<<"Input file :   "<<par.infile<<"\n";
      cout<<"Database file: "<<par.dbfiles<<"\n";
      cout<<"Output file:   "<<par.outfile<<"\n";
    }

  // Prepare CS pseudocounts lib
  if (!par.nocontxt && *par.clusterfile) {
    InitializePseudocountsEngine();
  }

  // Set secondary structure substitution matrix
  if (par.ssm) SetSecStrucSubstitutionMatrix();

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix();

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  char input_format=0;
  ReadQueryFile(par.infile,input_format,q); 
  PrepareQueryHMM(input_format,q);

//  // Rescale matrix according to query aa composition? (Two iterations are sufficient)
//   if (par.pcm==4)
//     {
//       q->RescaleMatrix();
//       q->PreparePseudocounts();
//       q->AddAminoAcidPseudocounts();
//       SetSubstitutionMatrix();
//       q->RescaleMatrix();
//       q->PreparePseudocounts();
//       q->AddAminoAcidPseudocounts();
//     }

  // Reset lamda?
  if (par.calibrate>0) {q->lamda=LAMDA; q->mu=0.0;}

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags) q->NeutralizeTags();

  if (par.forward>=1)
    {
      if (v>=2 &&par.forward==2) printf("Using maximum accuracy (MAC) alignment algorithm ...\n");
    }
  else if (v>=3) printf("Using Viterbi algorithm ...\n");

  // Prepare multi-threading - reserve memory for threads, intialize, etc.
  if (threads==0) bins=1; else bins=iround(threads*1.2+0.5);
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs=1;   // needs to be set to 1 before threads are created
  for (bin=0; bin<bins; bin++)
    {
      t[bin]=new HMM;   // Each bin has a template HMM allocated that was read from the database file
      hit[bin]=new Hit; // Each bin has an object of type Hit allocated ...
      hit[bin]->AllocateBacktraceMatrix(q->L+2,par.maxres); // ...with a separate dynamic programming matrix (memory!!)
      if (par.forward>=1)
        hit[bin]->AllocateForwardMatrix(q->L+2,par.maxres);

      bin_status[bin] = FREE;
    }

#ifdef PTHREAD
  // Start threads for database search
  for (int j=0; j<threads; j++)
    {
      thread_data[j].thread_id = j+1;
      thread_data[j].function  = &AlignByWorker;
      if (DEBUG_THREADS) fprintf(stderr,"Creating worker thread %i ...",j+1);
      pthread_create(&pthread[j], &joinable, WorkerLoop, (void*)&thread_data[j]);
      if (DEBUG_THREADS) fprintf(stderr," created!\n");
    }
#endif
  format = new(int[bins]);


  // Cut par.dbfiles into separate database strings and read pal files

  int ndb=0;
  Hash<char>* doubled;
  doubled = new(Hash<char>);
  doubled->New(16381,0);
  char** dbfiles = new char*[par.maxnumdb+1];
  char* dbfile_cur=strscn_(par.dbfiles); // current file name
  char* dbfile_next; // next file name after current
  while (*dbfile_cur && ndb<par.maxnumdb) {

      // Cut off everything after next white space and store beginning of next database name in dbfile_next
      dbfile_next=strcut_(dbfile_cur);

      // Has the dbfiles[ndb] name not been seen before?
      if (! doubled->Contains(dbfile_cur)) {
          doubled->Add(dbfile_cur);
          if (!strcmp(dbfile_cur+strlen(dbfile_cur)-4,".pal"))
            {
              char dbfile[NAMELEN]="";    // input line
              FILE* palf = NULL;
              if (!strcmp(dbfile,"stdin")) palf=stdin;
              else
                {
                  palf = fopen(dbfile_cur,"rb");
                  if (!palf) OpenFileError(dbfile_cur);
                }
              while(fgetline(dbfile,NAMELEN,palf) && ndb<par.maxnumdb) // read HMM files in pal file
                {
                  if (! doubled->Contains(dbfile))
                    {
                      doubled->Add(dbfile);
                      dbfiles[ndb]=new(char[strlen(dbfile)+1]);
                      strcpy(dbfiles[ndb],dbfile);
                      if (ndb<5 && ndb>0 && access(dbfiles[ndb],R_OK)) OpenFileError(dbfiles[ndb]); // file not readable?
//	   	       printf("dbfiles[%i]='%s'\n",ndb,dbfiles[ndb]);
                      ndb++;
                    }
                  else if (v>=1) 
                    fprintf(stderr,"WARNING: skipping doubled datbase file %s\n",dbfile);
                }
              fclose(palf);
            }
          else
            {
              dbfiles[ndb]=new(char[strlen(dbfile_cur)+1]);
              strcpy(dbfiles[ndb],dbfile_cur);
              if (ndb<5 && ndb>0 && access(dbfiles[ndb],R_OK)) OpenFileError(dbfiles[ndb]); // file not readable?
//            printf("dbfiles[%i]='%s'\n",ndb,dbfiles[ndb]);
              ndb++;
            }
        }
      else if (v>=1)
        fprintf(stderr,"WARNING: skipping doubled datbase file %s\n",dbfile_cur);

     dbfile_cur=dbfile_next;
    }

  //fprintf (stderr,"ndb: %i  par.maxnumdb: %i\n",ndb,par.maxnumdb);

  if (v>=1 && ndb>=par.maxnumdb && *dbfiles[ndb-1])
    fprintf (stderr,"WARNING: maximum of %i allowed databases surpassed. Skipping the rest.\n",par.maxnumdb);



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Search databases

  // Initialize
  N_searched=0;
  v1=v;
  if (v>0 && v<=3) v=1; else v-=2;
  if (print_elapsed) ElapsedTimeSinceLastCall("(preparing for search)");

  // For all the databases given in -d '...' option ...
  for (int idb=0; idb<ndb; idb++)
    {

      // Open HMM database
//    cerr<<"\nReading db file "<<idb<<" dbfiles[idb]="<<dbfiles[idb]<<"\n";
      FILE* dbf=fopen(dbfiles[idb],"rb");
      if (!dbf) OpenFileError(dbfiles[idb]);
      read_from_db=1;

      ///////////////////////////////////////////////////////////////////////////////////////
      // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
      while(read_from_db!=0)
        {

          // Submit jobs until no bin is free anymore
          while (jobs_submitted+jobs_running<bins)
            {

              // Allocate free bin (no need to lock, since slave processes cannot change FREE to other status)
              bin = PickBin(FREE);
              if (bin<0) {
                fprintf(stderr,"Error during search: found no free bin! jobs running: %i  jobs_submitted:%i  threads:%i\n",jobs_running,jobs_submitted,threads);
                for (bin=0; bin<bins; bin++) fprintf(stderr,"bin_status[%i]=%i\n",bin,bin_status[bin]);
                exit(6);
              }
              hit[bin]->index = N_searched;          // give hit a unique index for HMM
              hit[bin]->ftellpos = ftell(dbf);       // record position in dbfile of next HMM to be read
              char path[NAMELEN];
              Pathname(path,dbfiles[idb]);

              ///////////////////////////////////////////////////
              // Read next HMM from database file
              if (!fgetline(line,LINELEN,dbf)) {read_from_db=0; break;}
	      while (strscn(line)==NULL && fgetline(line,LINELEN,dbf)) {} // skip lines that contain only white space
	      if (strscn(line)==NULL) break;
       
	      if (!strncmp(line,"HMMER3",6))      // read HMMER3 format
                {
                  format[bin] = 1;
                  read_from_db = t[bin]->ReadHMMer3(dbf,dbfiles[idb]);
		  par.hmmer_used = true;
                }
              else if (!strncmp(line,"HMMER",5))      // read HMMER format
                {
                  format[bin] = 1;
                  read_from_db = t[bin]->ReadHMMer(dbf,dbfiles[idb]);
		  par.hmmer_used = true;
                }
              else if (!strncmp(line,"HH",2))    // read HHM format
                {
                  format[bin] = 0;
                  read_from_db = t[bin]->Read(dbf,path);
                }
              else if (!strncmp(line,"NAME",4))  // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
                {
                  fseek(dbf,hit[bin]->ftellpos,SEEK_SET); // rewind to beginning of line
                  format[bin] = 0;
                  read_from_db = t[bin]->Read(dbf,path);
                }
              else if (line[0]=='#')             // read a3m alignment
                {
                  Alignment tali;
                  tali.Read(dbf,dbfiles[idb],line);
                  tali.Compress(dbfiles[idb]);
                  //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
                  tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
                  char wg=par.wg; par.wg=1; // use global weights
                  t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
                  tali.FrequenciesAndTransitions(t[bin]);
                  par.wg=wg; //reset global weights
                  format[bin] = 0;
                }
              else
                {
		  cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<dbfiles[idb]<<"\'. \n";
		  cerr<<"Context:\n'"<<line<<"\n";
		  fgetline(line,LINELEN,dbf); cerr<<line<<"\n";
		  fgetline(line,LINELEN,dbf); cerr<<line<<"'\n";
		  exit(1);
                }
	      
              if (read_from_db==2) continue;  // skip current HMM or reached end of database
              if (read_from_db==0) break;     // finished reading HMMs
              if (v>=4) printf("Aligning with %s\n",t[bin]->name);  /////////////////////v>=4
              ///////////////////////////////////////////////////

              hit[bin]->dbfile = new(char[strlen(dbfiles[idb])+1]);
              strcpy(hit[bin]->dbfile,dbfiles[idb]); // record db file name from which next HMM is read

              ++N_searched;
              if (v1>=2 && !(N_searched%20))
                {
                  cout<<".";
                  if (!(N_searched%1000)) printf(" %-4i HMMs searched\n",N_searched);
                  cout.flush();
                }

#ifdef PTHREAD
              // Lock access to bin_status
              if (threads>0) rc = pthread_mutex_lock(&bin_status_mutex);
#endif
              // Send the job in bin to a thread
              bin_status[bin] = SUBMITTED;
              jobs_submitted++;

              if (threads==0) // if no multi-threading mode, main thread executes job itself
                {
                  AlignByWorker(bin);
                  bin_status[bin] = FREE;
                  jobs_submitted--;
                  break;
                }

#ifdef PTHREAD
              // Restart threads waiting for a signal
              rc = pthread_cond_signal(&new_job);

              // Unlock access to bin_status
              rc = pthread_mutex_unlock(&bin_status_mutex);

              if (DEBUG_THREADS)
                fprintf(stderr,"Main: put job into bin %i  name=%-7.7s  Jobs running: %i  jobs_submitted:%i \n",bin,t[bin]->name,jobs_running,jobs_submitted);
#endif
            }

#ifdef PTHREAD
          if (threads>0)
            {
              // Lock mutex
              rc = pthread_mutex_lock(&bin_status_mutex);

              // Wait until job finishes and a bin becomes free
              if (jobs_submitted+jobs_running>=bins)
                {
                  if (DEBUG_THREADS) fprintf(stderr,"Master thread is waiting for jobs to be finished...\n");
#ifdef WINDOWS
                  rc = pthread_cond_wait(&finished_job, &bin_status_mutex);
#else
#ifdef HH_MAC

		  // If no submitted jobs are in the queue we have to wait for a new job ...
		  struct timespec ts;
		  struct timeval tv;
		  gettimeofday(&tv, NULL);
		  ts.tv_sec = tv.tv_sec + 1;
		  rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex,&ts);
#else
		  
		  // If no submitted jobs are in the queue we have to wait for a new job ...
		  struct timespec ts;
		  clock_gettime(CLOCK_REALTIME,&ts);
		  ts.tv_sec += 1;
		  rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex,&ts);
#endif
#endif
                }
              // Unlock mutex
              rc = pthread_mutex_unlock(&bin_status_mutex);
            }
#endif

        } // end while(read_from_db!=0)
      ///////////////////////////////////////////////////////////////////////////////////////

      fclose(dbf);
    }
  reading_dbs=0;


#ifdef PTHREAD
  if (threads>0)
    {
      // No more HMMs in database => wait until all threads have finished
      if (DEBUG_THREADS)
        fprintf(stderr,"No more jobs read from database         Jobs running:%i  jobs_submitted:%i \n",jobs_running,jobs_submitted);

      // Free all threads still waiting for jobs
      rc = pthread_mutex_lock(&bin_status_mutex);
      rc = pthread_cond_broadcast(&new_job);
      rc = pthread_mutex_unlock(&bin_status_mutex); // Unlock mutex

      // Wait for all jobs to finish => join all jobs to main
      for (int j=0; j<threads; j++)
        {
          int status;
          pthread_join(pthread[j], (void **)&status);
          if (DEBUG_THREADS) fprintf(stderr,"Thread %i finished its work\n",j+1);
        }
    }
#endif
  // End search databases
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  hitlist.N_searched=N_searched; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  if (v1>=2) cout<<"\n";
  v=v1;

  if (print_elapsed) ElapsedTimeSinceLastCall("(search through database)");

  // Sort list according to sortscore
  if (v>=3) printf("Sorting hit list ...\n");
  hitlist.SortList();

  // Fit EVD (with lamda, mu) to score distribution?
  if (par.calm==3)
    {
      hitlist.CalculatePvalues(q);  // Use NN prediction of lamda and mu
    }
  else if ((par.calm!=1 && q->lamda==0) || par.calibrate>0)
    {
      if (v>=2 && par.loc) printf("Fitting scores with EVD (first round) ...\n");
      hitlist.MaxLikelihoodEVD(q,3); // first ML fit: exclude 3 best superfamilies from fit

      if (v>=3) printf("Number of families present in database: %i\n",hitlist.fams); // DEBUG
      if (hitlist.fams>=100)
        {
          if (par.loc)
            {
              if (v>=2) printf("Fitting scores with EVD (second round) ...\n");
              hitlist.MaxLikelihoodEVD(q,0); // second ML fit: exclude superfamilies with E-value<MINEVALEXCL
              hitlist.ResortList();
            }
          else
            {
              if (v>=2)
                fprintf(stderr,"WARNING: E-values for global alignment option may be unreliable.\n");
              hitlist.ResortList();
            }
        }
      else
        {
          if (v)
            {
              fprintf(stderr,"\nWARNING: no E-values could be calculated.\n");
              fprintf(stderr,"To calculate E-values you have two options:\n");
              fprintf(stderr,"1. Calibrate your query profile HMM with a calibration database:\n");
              fprintf(stderr,"   > hhsearch -i yourHMM.hhm -d cal.hhm -cal\n");
              fprintf(stderr,"   This will insert a line in yourHMM.hhm with lamda and mu of the score distribution.\n");
              fprintf(stderr,"   cal.hhm contains 1220 HMMs from different SCOP superfamilies and is supplied with HHsearch.\n");
              fprintf(stderr,"   Instead of cal.hhm you may also use any SCOP database file, e.g. scop70_1.69\n");
              fprintf(stderr,"   Note that your HMM needs to be recalibrated when changing HMM-HMM alignment options.\n");
              fprintf(stderr,"2. Append cal.hhm to your own database:\n");
              fprintf(stderr,"   > cat cal.hhm >> yourDB.hhm\n");
              fprintf(stderr,"   But note that HMMs contained in cal.hmm will pop up among your hits.\n");
            }
        }
      if (par.calm==2) hitlist.GetPvalsFromCalibration(q);
    }
  else
    hitlist.GetPvalsFromCalibration(q);

  // Optimization mode?
  if (par.opt) hitlist.Optimize(q);

  // Set new ss weight for realign
  par.ssw = par.ssw_realign;

  
  // Realign hits with MAC algorithm
  if (par.realign && par.forward!=2) 
    {
      perform_realign(dbfiles,ndb); 
    } 
  else {
    // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
    if (*par.scorefile) {
      if (v>=3) printf("Printing scores file ...\n");
    }
  }
  
  // // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
  // if (*par.scorefile) {
  //   if (v>=3) printf("Printing scores file ...\n");
  //   hitlist.PrintScoreFile(q);
  // }
  
  // Print FASTA or A2M alignments?
  if (*par.pairwisealisfile) {
    if (v>=2) cout<<"Printing alignments in "<<(par.outformat==1? "FASTA" : par.outformat==2?"A2M" :"A3M")<<" format to "<<par.pairwisealisfile<<"\n";
    hitlist.PrintAlignments(q,par.pairwisealisfile,par.outformat);
  }
  
  // Warn, if HMMER files were used
  if (par.hmmer_used && v>=1)
    fprintf(stderr,"WARNING: Using HMMER files results in a drastically reduced sensitivity (>10%%).\nWe recommend to use HHMs build by hhmake.\n");

  // Print summary listing of hits
  if (v>=3) printf("Printing hit list ...\n");
  hitlist.PrintHitList(q,par.outfile);

  // Write only hit list to screen?
  if (v==2 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,1009); // write only hit list to screen

  // Print alignments of query sequences against hit sequences
  hitlist.PrintAlignments(q,par.outfile);

  // Write whole output file to screen? (max 10000 lines)
  if (v>=3 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,10009);

  // Write HMM to output file without pseudocounts
  if (par.calibrate)
    q->InsertCalibration(par.infile);


  // Generate result alignment or HMM file?
  if (*par.alnfile || *par.psifile || *par.hhmfile)
    {
      Alignment Qali;   // output A3M generated by merging A3M alignments for significant hits to the query alignment
      Hit hit;
      int nhits=0;

      // Read query alignment into Qali
      char qa3mfile[NAMELEN];
      char ta3mfile[NAMELEN];
      RemoveExtension(qa3mfile,par.infile); // directory??
      strcat(qa3mfile,".a3m");
      FILE* qa3mf=fopen(qa3mfile,"r");
      if (!qa3mf) OpenFileError(qa3mfile);
      Qali.Read(qa3mf,qa3mfile);
      fclose(qa3mf);

      // If par.append==1 do not print query alignment
      if (par.append) Qali.MarkSeqsAsNonPrintable();

      if (v>=2) printf("Merging hits to query alignment %s ...\n",qa3mfile);
      // If query consists of only one sequence:
      //     realign templates to query HMM enforced by sequences from up to the 10 best templates
      v1=v--;

      // For each template below threshold
      hitlist.Reset();
      while (!hitlist.End())
        {
          hit = hitlist.ReadNext();
          if (hit.Eval > 100.0*par.e) break; // E-value much too large
          if (hit.Eval > par.e) continue; // E-value too large

          // Read a3m alignment of hit from <file>.a3m file and merge into Qali alignment
          strcpy(ta3mfile,hit.file); // copy filename including path but without extension
          strcat(ta3mfile,".a3m");
 	  Alignment Tali;
	  FILE* ta3mf=fopen(ta3mfile,"r");
	  if (!ta3mf) OpenFileError(ta3mfile);
	  Tali.Read(ta3mf,ta3mfile); // Read template alignment into Tali
	  fclose(ta3mf);
	  Tali.Compress(ta3mfile); // Filter database alignment
	  Qali.MergeMasterSlave(hit,Tali,ta3mfile);
	  nhits++;
        }

      // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
      Qali.Compress("merged A3M file");

      // Sort out the nseqdis most dissimilar sequences for display in the result alignments
      Qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);

      v=v1;

      // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
      float const COV_ABS = 25;     // min. number of aligned residues
      int cov_tot = imax(imin((int)(COV_ABS / Qali.L * 100 + 0.5), 70), par.coverage);
      if (v>2) printf("Filter new alignment with cov %3i%%\n", cov_tot);
      Qali.N_filtered = Qali.Filter(par.max_seqid,cov_tot,par.qid,par.qsc,par.Ndiff);

      // Calculate (and write) output HMM?
      if (*par.hhmfile)
        {
          strcpy(Qali.longname,q->longname);
          strcpy(Qali.name,q->name);
          strcpy(Qali.fam,q->fam);
          RemovePathAndExtension(Qali.file,par.hhmfile);

          // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
	  HMM* Q = new HMM; // output HMM: generated from Qali
          Qali.FrequenciesAndTransitions(Q);

          // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
          Q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
          Q->CalculateAminoAcidBackground();

          // Write HMM to output file in HHsearch format?
          Q->WriteToFile(par.hhmfile);
	  delete Q;
        }

      // Write output A3M alignment?
      if (*par.alnfile) Qali.WriteToFile(par.alnfile,"a3m");

      // Write output PSI-BLAST-formatted alignment?
      if (*par.psifile) Qali.WriteToFile(par.psifile,"psi");

   }


  // Write alignments with posteriors etc. to alitabfile?
  if (*par.alitabfile) 
    {
      FILE* alitabf=NULL;
      if (strcmp(par.alitabfile,"stdout")) alitabf = fopen(par.alitabfile, "w"); else alitabf = stdout;
      if (!alitabf) OpenFileError(par.alitabfile);
      
      // Store all dbfiles and ftell positions of templates to be displayed and realigned
      int nhits=0;
      hitlist.Reset();
      while (!hitlist.End())
        {
          hit_cur = hitlist.ReadNext();
          if (nhits>=imax(par.B,par.Z)) break;
          if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
          if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;

	  fprintf(alitabf, ">%s\n", hit_cur.longname);
	  WriteToAlifile(alitabf,&hit_cur);
     	  nhits++;
	}
      fclose(alitabf);
    }
  
  // Delete memory for dynamic programming matrix
  for (bin=0; bin<bins; bin++)
    {
      hit[bin]->DeleteBacktraceMatrix(q->L+2);
      if (par.forward>=1 || par.realign)
        hit[bin]->DeleteForwardMatrix(q->L+2);
      delete hit[bin];
      delete t[bin];
     }
  if (par.dbfiles) delete[] par.dbfiles;
  if (dbfiles) {
    for (int idb=0; idb<ndb; idb++) delete[](dbfiles[idb]);
    delete[] dbfiles;
  }
  if (format) delete[](format);
 if (par.exclstr) delete[] par.exclstr;
  for (int n = 1; n < argc_conf; n++)
    delete[] argv_conf[n];
  delete doubled;

  DeletePseudocountsEngine();

  // Delete content of hits in hitlist
  hitlist.Reset();
  while (!hitlist.End())
    {
      //hitlist.ReadCurrent().Delete();
      hitlist.Delete().Delete(); // Delete list record and hit object
    }

  delete q;

#ifdef PTHREAD
  pthread_attr_destroy(&joinable);
  pthread_mutex_destroy(&bin_status_mutex);
  pthread_cond_destroy(&new_job);
  pthread_cond_destroy(&finished_job);
#endif


  if (print_elapsed) ElapsedTimeSinceLastCall("(sorting and formatting)");

  // Print 'Done!'
  FILE* outf=NULL;
  if (!strcmp(par.outfile,"stdout")) printf("Done!\n");
  else
    {
      if (*par.outfile)
        {
          outf=fopen(par.outfile,"a"); //open for append
          fprintf(outf,"Done!\n");
          fclose(outf);
        }
      if (v>=2) printf("Done\n");
    }

  exit(0);
} //end main

//////////////////////////////////////////////////////////////////////////////////////////////////////
// END OF MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////



