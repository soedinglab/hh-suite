// hhblits.C:
// Iterative search for a multiple alignment in a profile HMM database
//
// Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: command line error  6: internal logic error  7: internal numeric error

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
//     Nat. Methods 9:173-175 (2011); epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

//////////////////////////////////////////////////////////////////////////////////////////
//     This program contains, in file hhprefilter.C, code adapted from Michael Farrar
//     (http://sites.google.com/site/farrarmichael/smith-waterman). His code is marked
//     in the file hhprefilter.C.
//     The copy right of his code is shown below:

//     Copyright 2006, by Michael Farrar.  All rights reserved. The SWSSE2
//     program and documentation may not be sold or incorporated into a
//     commercial product, in whole or in part, without written consent of
//     Michael Farrar.
//
//     For further information regarding permission for use or reproduction, 
//     please contact Michael Farrar at:
//
//         farrar.michael@gmail.com
//
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
//     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
//     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
//     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
//     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
//     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//     Note by J. Soeding: Michael Farrar died unexpectedly in December 2010. 
//     Many thanks posthumously for your great code!

#include <map>

#define PTHREAD
#define MAIN
#define HHBLITS
// #define HHDEBUG
// #define DEBUG_THREADS
// #define WINDOWS

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <sstream>
#include <vector>
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror(), strerror(errno)
#include <cassert>
#include <stdexcept>

#ifdef PTHREAD
#include <pthread.h>  // POSIX pthread functions and data structures
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/time.h>

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
using std::string;
using std::stringstream;
using std::vector;
using std::pair;

extern "C" {
#include <ffindex.h>     // fast index-based database reading
}

#include "cs.h"          // context-specific pseudocounts
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"
#include "abstract_state_matrix.h"
cs::ContextLibrary<cs::AA> *cs_lib;

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure
#include "hhdecl.C"      // Constants, global variables, struct Parameters

std::map<std::string, unsigned char*> columnStateSequences;
ColumnStateScoring* columnStateScoring;

#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM
#include "hhhmm.h"       // class HMM
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

char line[LINELEN] = "";         // input line
int bin;                       // bin index
const char print_elapsed = 0;    // debug output for runtimes
char tmp_file[] = "/tmp/hhblitsXXXXXX"; // for runtime secondary structure prediction (only with -addss option)

// HHblits variables
const char HHBLITS_REFERENCE[] =
    "Remmert M., Biegert A., Hauser A., and Soding J.\nHHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.\nNat. Methods 9:173-175 (2011)\n";

int v1 = v;                               // verbose mode
int num_rounds = 2;                   // number of iterations
bool last_round = false;                // set to true in last iteration
bool already_seen_filter = true;       // Perform filtering of already seen HHMs
bool realign_old_hits = false; // Realign old hits in last round or use previous alignments

char input_format = 0; // Set to 1 if input in HMMER format (has already pseudocounts)

float neffmax = 10;                     // Break if Neff > Neffmax

char config_file[NAMELEN];
char infile[NAMELEN];
char alis_basename[NAMELEN];
char query_hhmfile[NAMELEN];             // -qhmm output file
bool alitab_scop = false;            // Write only SCOP alignments in alitabfile
char db_ext[NAMELEN];
int omp_threads = 2;                       // number of OpenMP threads to start

// Needed for fast index reading
size_t data_size;
FILE *dba3m_data_file;
FILE *dba3m_index_file;
FILE *dbhhm_data_file;
FILE *dbhhm_index_file;

char* dba3m_data;
char* dbhhm_data;
ffindex_index_t* dbhhm_index = NULL;
ffindex_index_t* dba3m_index = NULL;

//database filenames
char db_base[NAMELEN];                   // database basename
char dbcs_base[NAMELEN];
char dbcs_index_filename[NAMELEN];
char dbcs_data_filename[NAMELEN];

char dbhhm_base[NAMELEN];
char dbhhm_index_filename[NAMELEN];
char dbhhm_data_filename[NAMELEN];

char dba3m_base[NAMELEN];
char dba3m_index_filename[NAMELEN];
char dba3m_data_filename[NAMELEN];

char** dbfiles_new;
char** dbfiles_old;
int ndb_new = 0;
int ndb_old = 0;
Hash<Hit>* previous_hits;
Hash<char>* premerged_hits;

// HHsearch variables
const int MAXTHREADS = 256; // maximum number of threads (i.e. CPUs) for parallel computation
const int MAXBINS = 384;   // maximum number of bins (positions in thread queue)
enum bin_states {
  FREE = 0, SUBMITTED = 1, RUNNING = 2
};
int threads = 2; // number of compute pthreads during Viterbi and Realign (apart from the main thread which reads from db file) and # OpenMP threads; 0:no multithreading
int bins; // number of bins; jobs gets allocated to a FREE bin were they are waiting for execution by a thread
char bin_status[MAXBINS]; // The status for each bin is FREE, SUBMITTED, or RUNNING
int jobs_running;   // number of active jobs, i.e. number of bins set to RUNNING
int jobs_submitted; // number of submitted jobs, i.e. number of bins set to SUBMITTED
char reading_dbs; // 1: still HMMs to read in a database;  0: finshed reading all HMMs, no db left
const char DEBUG_THREADS = 0; // Debugging flag

HMM* q;              // Create query HMM with maximum of par.maxres match states
HMM* q_tmp; // Create query HMM with maximum of par.maxres match states (needed for prefiltering)
HMM* t[MAXBINS]; // Each bin has a template HMM allocated that was read from the database file
Hit* hit[MAXBINS]; // Each bin has an object of type Hit allocated with a separate dynamic programming matrix (memory!!)
Hit hit_cur;              // Current hit when going through hitlist
HitList hitlist; // list of hits with one Hit object for each pairwise comparison done
int* format; // format[bin] = 0 if in HHsearch format => add pcs; format[bin] = 1 if in HMMER format => no pcs
int read_from_db; // The value of this flag is returned from HMM::Read(); 0:end of file  1:ok  2:skip HMM
int N_searched;           // Number of HMMs searched
Alignment Qali; // output A3M generated by merging A3M alignments for significant hits to the query alignment
Alignment Qali_allseqs; // output A3M alignment with no sequence filtered out (only active with -all option)

struct Thread_args // data to hand to WorkerLoop thread
{
    int thread_id;          // id of thread (for debugging)
    void (*function)(int); // pointer to function (=job) to execute by WorkerLoop once old job is done
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

inline int PickBin(char status);

// Include hhworker.C and hhprefilter.C here, because it needs some of the above variables
#include "hhworker.C"      // functions: AlignByWorker, RealignByWorker, WorkerLoop
#include "hhprefilter.C"   // some prefilter functions

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help(char all = 0) {
  //      ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+---8-----+----9----+----0
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
      db_base);
  printf(
      " -n     [1,8]   number of iterations (default=%i)                              \n",
      num_rounds);
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
        " -noprefilt     disable all filter steps                                        \n");
    printf(
        " -noaddfilter   disable all filter steps (except for fast prefiltering)         \n");
    printf(
        " -nodbfilter    disable additional filtering of prefiltered HMMs                \n");
    printf(
        " -noblockfilter search complete matrix in Viterbi                               \n");
    printf(
        " -maxfilt       max number of hits allowed to pass 2nd prefilter (default=%i)  \n",
        par.maxnumdb);
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
      " -usecs         use column states of the templates in the database for scoring   \n");
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

    printf(" Context specific prefilter pseudocounts:\n");
    printf("  -pc_prefilter_contxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc_prefilter_context_engine.admix);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_prefilter_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc_prefilter_context_engine.pca);
    printf("  -pc_prefilter_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",par.pc_prefilter_context_engine.pcb);
    printf("  -pc_prefilter_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",par.pc_prefilter_context_engine.pcc);

    printf(" Context independent prefilter pseudocounts (used if context file is not available):\n");
    printf("  -pc_prefilter_nocontxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc_prefilter_nocontext_mode);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
  //  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_prefilter_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc_prefilter_nocontext_a);
    printf("  -pc_prefilter_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",par.pc_prefilter_nocontext_b);
    printf("  -pc_prefilter_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",par.pc_prefilter_nocontext_c);

    printf("\n");
    printf(" Context-specific pseudo-counts:                                                  \n");
    printf("  -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
    printf("  -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",par.clusterfile);
    //should not be in the section of pseudocounts ... associated to prefiltering ... and also to usecs (by markus)
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
      neffmax);
#ifdef PTHREAD
  printf(
      " -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=%i)      \n",
      threads);
#endif
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
#ifndef PTHREAD
  printf("(The -cpu option is inactive since HHblits was not compiled with POSIX thread support)\n");
#endif
  printf("\n");
  if (!all) {
    printf(
        "An extended list of options can be obtained by calling 'hhblits -help'\n");
  }
  printf("\n");
  printf("Example: %s -i query.fas -oa3m query.a3m -n 1  \n", program_name);
  cout << endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(int argc, char** argv) {
  //Processing command line input
  for (int i = 1; i < argc; i++) {
    if (v >= 4)
      cout << i << "  " << argv[i] << endl; //PRINT
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no query file following -i\n";
        exit(4);
      }
      else
        strcpy(par.infile, argv[i]);
    }
    else if (!strcmp(argv[i], "-d")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no database basename following -d\n";
        exit(4);
      }
      else
        strcpy(db_base, argv[i]);
    }
    else if (!strcmp(argv[i], "-contxt") || !strcmp(argv[i], "-context_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no lib following -contxt\n";
        exit(4);
      }
      else
        strcpy(par.clusterfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-cslib") || !strcmp(argv[i], "-cs_lib")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no lib following -cslib\n";
        exit(4);
      }
      else
        strcpy(par.cs_library, argv[i]);
    }
    else if (!strcmp(argv[i], "-psipred")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no directory following -psipred\n";
        exit(4);
      }
      else
        strcpy(par.psipred, argv[i]);
    }
    else if (!strcmp(argv[i], "-psipred_data")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no database directory following -psipred_data\n";
        exit(4);
      }
      else
        strcpy(par.psipred_data, argv[i]);
    }
    else if (!strcmp(argv[i], "-o")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa3m")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no output file following -oa3m\n";
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-ohhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no output file following -ohhm\n";
        exit(4);
      }
      else
        strcpy(par.hhmfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-opsi")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no output file following -opsi\n";
        exit(4);
      }
      else
        strcpy(par.psifile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oalis")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no file basename following -oalis\n";
        exit(4);
      }
      else
        strcpy(alis_basename, argv[i]);
    }
    else if (!strcmp(argv[i], "-Ofas")) {
      par.append = 0;
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help();
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
        help();
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
        help();
        cerr << endl << "Error in " << program_name
            << ": no output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-qhhm")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no filename following -qhhm\n";
        exit(4);
      }
      else
        strcpy(query_hhmfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-scores")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no file following -scores\n";
        exit(4);
      }
      else {
        strcpy(par.scorefile, argv[i]);
      }
    }
    else if (!strcmp(argv[i], "-db_ext")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no extension following -db_ext\n";
        exit(4);
      }
      else {
        strcpy(db_ext, argv[i]);
      }
    }
    else if (!strcmp(argv[i], "-atab")) {
      if (++i >= argc || argv[i][0] == '-') {
        help();
        cerr << endl << "Error in " << program_name
            << ": no file following -atab\n";
        exit(4);
      }
      else {
        strcpy(par.alitabfile, argv[i]);
      }
    }
    else if (!strcmp(argv[i], "-atab_scop"))
      alitab_scop = true;
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
      help(1);
      exit(0);
    }
    else if (!strcmp(argv[i], "-v") && (i < argc - 1) && argv[i + 1][0] != '-')
      v = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-v"))
      v = 2;
    else if (!strcmp(argv[i], "-v0"))
      v = 0;
    else if (!strcmp(argv[i], "-v1"))
      v = 1;
    else if (!strcmp(argv[i], "-n") && (i < argc - 1))
      num_rounds = atoi(argv[++i]);
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
      neffmax = atof(argv[++i]);
    else if ((!strcmp(argv[i], "-neff") || !strcmp(argv[i], "-Neff"))
        && (i < argc - 1))
      par.Neff = atof(argv[++i]);
    //pc hhm context variables
    else if (!strcmp(argv[i],"-pc_hhm_contxt_mode") && (i<argc-1)) par.pc_hhm_context_engine.admix=(Pseudocounts::Admix)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-pc_hhm_contxt_a") && (i<argc-1)) par.pc_hhm_context_engine.pca=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pc_hhm_contxt_b") && (i<argc-1)) par.pc_hhm_context_engine.pcb=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pc_hhm_contxt_c") && (i<argc-1)) par.pc_hhm_context_engine.pcc=atof(argv[++i]);
    //pc prefilter context variables
    else if (!strcmp(argv[i],"-pc_prefilter_contxt_mode") && (i<argc-1)) par.pc_prefilter_context_engine.admix=(Pseudocounts::Admix)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-pc_prefilter_contxt_a") && (i<argc-1)) par.pc_prefilter_context_engine.pca=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pc_prefilter_contxt_b") && (i<argc-1)) par.pc_prefilter_context_engine.pcb=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pc_prefilter_contxt_c") && (i<argc-1)) par.pc_prefilter_context_engine.pcc=atof(argv[++i]);
    //pc hhm nocontext variables
    else if (!strcmp(argv[i],"-pc_hhm_nocontxt_mode") && (i<argc-1)) par.pc_hhm_nocontext_mode=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-pc_hhm_nocontxt_a") && (i<argc-1)) par.pc_hhm_nocontext_a=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pc_hhm_nocontxt_b") && (i<argc-1)) par.pc_hhm_nocontext_b=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pc_hhm_nocontxt_c") && (i<argc-1)) par.pc_hhm_nocontext_c=atof(argv[++i]);
    //pc prefilter nocontext variables
    else if (!strcmp(argv[i],"-pc_prefilter_nocontxt_mode") && (i<argc-1)) par.pc_prefilter_nocontext_mode = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-pre_pca") && (i<argc-1)) par.pc_hhm_nocontext_a=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pre_pcb") && (i<argc-1)) par.pc_hhm_nocontext_b=atof(argv[++i]);
    else if (!strcmp(argv[i],"-pre_pcc") && (i<argc-1)) par.pc_hhm_nocontext_c=atof(argv[++i]);

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
    else if (!strcmp(argv[i], "-filterlen") && (i < argc - 1))
      early_stopping->length = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-filtercut") && (i < argc - 1))
      early_stopping->thresh = (double) atof(argv[++i]);
    else if (!strcmp(argv[i], "-noprefilt") || !strcmp(argv[i], "-nofilter")) {
      par.prefilter = false;
      already_seen_filter = false;
      par.early_stopping_filter = false;
      early_stopping->thresh = 0;
    }
    else if (!strcmp(argv[i], "-noaddfilter")) {
      already_seen_filter = false;
      par.early_stopping_filter = false;
      early_stopping->thresh = 0;
    }
    else if (!strcmp(argv[i], "-nodbfilter")) {
      early_stopping->thresh = 0;
    }
    else if (!strcmp(argv[i], "-noearlystoppingfilter")) {
      par.early_stopping_filter = false;
    }
    else if (!strcmp(argv[i], "-maxfilt") && (i < argc - 1))
      par.maxnumdb = par.maxnumdb_no_prefilter = atoi(argv[++i]);
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
      realign_old_hits = true;
    else if (!strcmp(argv[i], "-realign"))
      par.realign = 1;
    else if (!strcmp(argv[i], "-norealign"))
      par.realign = 0;
    else if (!strcmp(argv[i], "-ssm") && (i < argc - 1))
      par.ssm = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ssw") && (i < argc - 1))
      par.ssw = atof(argv[++i]);
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
      threads = atoi(argv[++i]);
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
    else if (!strcmp(argv[i], "-usecs")) {
      par.useCSScoring = true;
    }
    else
      cerr << endl << "WARNING: Ignoring unknown option " << argv[i]
          << " ...\n";
    if (v >= 4)
      cout << i << "  " << argv[i] << endl; //PRINT
  } // end of for-loop for command line input
}

///////////////////////////////////////////////////////////////////////////////////////
//// For multi-threading: return a bin with the desired status, return -1 if no such bin found
//////////////////////////////////////////////////////////////////////////////////////
inline int PickBin(char status) {
  for (int b = 0; b < bins; b++) {
    if (bin_status[b] == status)
      return b;
  }
  return -1;
}

///////////////////////////////////////////////////////////////////////////////////////
//// Do the pairwise comparison of q and t[bin] for the database search
//// Combination of RealignByWorker and AlignByWorker: 
//// Picks hits found in previous iterations and recalculates Viterbi scores using 
//// query profile from last iteration while KEEPING original (MAC) alignment.
//////////////////////////////////////////////////////////////////////////////////////
void PerformViterbiByWorker(int bin) {
  // Prepare q ant t and compare
  PrepareTemplateHMM(q, t[bin], format[bin]);

  // Do HMM-HMM comparison
  for (hit[bin]->irep = 1; hit[bin]->irep <= par.altali; hit[bin]->irep++) {
    // Break, if no previous_hit with irep is found
    hit[bin]->Viterbi(q, t[bin]);
    if (hit[bin]->irep > 1 && hit[bin]->score <= SMIN)
      break;
    hit[bin]->Backtrace(q, t[bin]);

    hit[bin]->score_sort = hit[bin]->score_aass;
    //printf("PerformViterbiByWorker:   %-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit[bin]->name,hit[bin]->fam,hit[bin]->irep,hit[bin]->score);

#ifdef PTHREAD
    pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif

    stringstream ss_tmp;
    ss_tmp << hit[bin]->file << "__" << hit[bin]->irep;

    if (previous_hits->Contains((char*) ss_tmp.str().c_str())) {
      //printf("Previous hits contains %s!\n",(char*)ss_tmp.str().c_str());
      hit_cur = previous_hits->Remove((char*) ss_tmp.str().c_str());
      previous_hits->Add((char*) ss_tmp.str().c_str(), *(hit[bin]));

      // Overwrite *hit[bin] with alignment, etc. of hit_cur
      hit_cur.score = hit[bin]->score;
      hit_cur.score_aass = hit[bin]->score_aass;
      hit_cur.score_ss = hit[bin]->score_ss;
      hit_cur.Pval = hit[bin]->Pval;
      hit_cur.Pvalt = hit[bin]->Pvalt;
      hit_cur.logPval = hit[bin]->logPval;
      hit_cur.logPvalt = hit[bin]->logPvalt;
      hit_cur.Eval = hit[bin]->Eval;
      hit_cur.logEval = hit[bin]->logEval;
      hit_cur.Probab = hit[bin]->Probab;

      hitlist.Push(hit_cur); // insert hit at beginning of list (last repeats first!)

    }
    else {
      // don't save alignments which where not found in previous rounds

      //printf("Don't save %s!\n",(char*)ss_tmp.str().c_str());
      //hitlist.Push(*(hit[bin]));          // insert hit at beginning of list (last repeats first!)
    }

#ifdef PTHREAD
    pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif

    if (hit[bin]->score <= SMIN)
      break;  // break if score for first hit is already worse than SMIN
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// This function is called if an HMM (<query>.hhm) instead of an MSA is given as query.
// It reads the query MSA either from <query>.a3m if it exists, or from the
//  representative sequences in the HHM/HMM file
/////////////////////////////////////////////////////////////////////////////////////
void ReadQueryA3MFile() {

  // Open query a3m MSA
  char qa3mfile[NAMELEN];
  RemoveExtension(qa3mfile, par.infile);
  strcat(qa3mfile, ".a3m");
  FILE* qa3mf = fopen(qa3mfile, "r");
  
  if (!qa3mf) {
    // <query>.a3m does not exist => extract query MSA from representative sequences in HHM
    if (v >= 1 && input_format == 0)  // HHM format
      printf(
          "Extracting representative sequences from %s to merge later with matched database sequences\n",
          par.infile);
    if (v >= 1 && input_format == 1) // HMMER format
      printf(
          "Extracting consensus sequence from %s to merge later with matched database sequences\n",
          par.infile);
    Qali.GetSeqsFromHMM(q);
    Qali.Compress(par.infile);
    if (Qali.L != q->L) {
      cerr << "Error in " << par.argv[0] << ": " << par.infile << " has "
          << q->L << " match states, while its representative sequences have "
          << Qali.L << "!\n";
      exit(1);
    }

    // if (num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile || *alis_basename)
    //  {
    //    if (input_format==0 && v>=1)   // HHM format
    //      cerr<<"WARNING: No alignment file found! Using only representative sequences in HHM file as starting MSA.\n";
    //    else if (input_format==1 && v>=1)
    //      cerr<<"WARNING: No alignment file found! Using only consensus sequence from HMMER file as starting MSA.\n";
    //  }
  }
  else {
    // <query>.a3m does exist => read it into Qali
    if (v >= 1)
      printf(
          "Reading query MSA from %s to merge later with matched database sequences\n",
          qa3mfile);
    Qali.Read(qa3mf, qa3mfile);
    Qali.Compress(qa3mfile);

    // Copy names from query HMM (don't use name of master sequence of MSA)
    delete[] Qali.longname;
    Qali.longname = new (char[strlen(q->longname) + 1]);
    strcpy(Qali.longname, q->longname);
    strcpy(Qali.name, q->name);
    strcpy(Qali.fam, q->fam);
    if (Qali.L != q->L) {
      printf(
          "Error in %s: query hhm %s has %i match states whereas query a3m %s has %i!\n",
          par.argv[0], qa3mfile, Qali.L, par.infile, q->L);
      exit(4);
    }
    fclose(qa3mf);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Perform Viterbi HMM-HMM search on all db HMM names in dbfiles
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DoViterbiSearch(char *dbfiles[], int ndb, bool alignByWorker = true) {
  // Search databases
  for (bin = 0; bin < bins; bin++) {
    hit[bin]->realign_around_viterbi = false;
  }

  double filter_cutoff = early_stopping->length * early_stopping->thresh;
  early_stopping->sum = early_stopping->length;
  early_stopping->counter = 0;
  for (int a = 0; a < early_stopping->length; a++) {
    early_stopping->evals[a] = 1.0;
  }

  // For all the databases comming through prefilter
  for (int idb = 0; idb < ndb; idb++) {
    // Check early stopping filter
    if (par.early_stopping_filter && early_stopping->sum < filter_cutoff) {
      if (v >= 4)
        printf(
            "Stop after DB-HHM %i from %i (filter_sum: %8.4f   cutoff: %8.4f)\n",
            idb, ndb, early_stopping->sum, filter_cutoff);
      printf("\n");
      break;
    }

    // Open HMM database
    FILE* dbf;
    dbf = ffindex_fopen_by_name(dbhhm_data, dbhhm_index, dbfiles[idb]);
    if (dbf == NULL) {
      if (dba3m_index_file == NULL) {
        cerr << endl << "Error opening " << dbfiles[idb]
            << ": A3M database missing\n";
        exit(4);
      }
      dbf = ffindex_fopen_by_name(dba3m_data, dba3m_index, dbfiles[idb]);
      if (dbf == NULL) {
        if (dbhhm_index_file == NULL) {
          cerr << endl << "Error opening " << dbfiles[idb]
              << ": HHM database missing\n";
          exit(4);
        }
        dbf = ffindex_fopen_by_name(dbhhm_data, dbhhm_index, dbfiles[idb]);
        if (dbf == NULL)
          OpenFileError(dbfiles[idb]);
      }
    }

    // Submit jobs if bin is free
    if (jobs_submitted + jobs_running < bins) {

      // Allocate free bin (no need to lock, since slave processes cannot change FREE to other status)
      bin = PickBin(FREE);
      if (bin < 0) {
        fprintf(stderr,
            "Error during search: found no free bin! jobs running: %i  jobs_submitted:%i  threads:%i\n",
            jobs_running, jobs_submitted, threads);

        for (bin = 0; bin < bins; bin++)
          fprintf(stderr, "bin_status[%i]=%i\n", bin, bin_status[bin]);
        exit(6);
      }
      hit[bin]->index = N_searched;          // give hit a unique index for HMM
      hit[bin]->ftellpos = ftell(dbf); // record position in dbfile of next HMM to be read
      //      fprintf(stderr,"dbfile=%-40.40s  index=%-5i  ftellpos=%i\n",dbfiles[idb],N_searched,(int) hit[bin]->ftellpos);

      char path[NAMELEN];
      Pathname(path, dbfiles[idb]);

      ///////////////////////////////////////////////////
      // Read next HMM from database file
      if (!fgetline(line, LINELEN, dbf)) {
        continue;
      }
      while (strscn(line) == NULL)
        fgetline(line, LINELEN, dbf); // skip lines that contain only white space

      if (!strncmp(line, "HMMER3", 6))      // read HMMER3 format
          {
        format[bin] = 1;
        t[bin]->ReadHMMer3(dbf, dbfiles[idb]);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HMMER", 5))      // read HMMER format
          {
        format[bin] = 1;
        t[bin]->ReadHMMer(dbf, dbfiles[idb]);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HH", 2))    // read HHM format
          {
        format[bin] = 0;
        t[bin]->Read(dbf, path);
      }
      else if (!strncmp(line, "NAME", 4)) // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
          {
        fseek(dbf, hit[bin]->ftellpos, SEEK_SET); // rewind to beginning of line
        format[bin] = 0;
        t[bin]->Read(dbf, path);
      }
      else if (line[0] == '#' || line[0] == '>')           // read a3m alignment
          {
        Alignment tali;
        tali.Read(dbf, dbfiles[idb], line);
        tali.Compress(dbfiles[idb]);
        //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
        tali.N_filtered = tali.Filter(par.max_seqid_db, par.coverage_db,
            par.qid_db, par.qsc_db, par.Ndiff_db);
        char wg = par.wg;
        par.wg = 1; // use global weights
        t[bin]->name[0] = t[bin]->longname[0] = t[bin]->fam[0] = '\0';
        tali.FrequenciesAndTransitions(t[bin]);
        par.wg = wg; //reset global weights
        format[bin] = 0;
      }
      else {
        cerr << endl << "Error in " << program_name
            << ": unrecognized HMM file format in \'" << dbfiles[idb]
            << "\'. \n";
        cerr << "Context:\n'" << line << "\n";
        fgetline(line, LINELEN, dbf);
        cerr << line << "\n";
        fgetline(line, LINELEN, dbf);
        cerr << line << "'\n";
        exit(1);
      }
      if (v >= 4)
        printf("Aligning with %s\n", t[bin]->name);  /////////////////////v>=4
      ///////////////////////////////////////////////////

      hit[bin]->dbfile = new (char[strlen(dbfiles[idb]) + 1]);
      strcpy(hit[bin]->dbfile, dbfiles[idb]); // record db file name from which next HMM is read

      ++N_searched;
      if (v1 >= 2 && !((idb + 1) % 20)) {
        cout << ".";
        if (!((idb + 1) % 1000))
          printf(" %-4i HMMs searched\n", (idb + 1));
        cout.flush();
      }

#ifdef PTHREAD
      // Lock access to bin_status
      if (threads > 0)
        rc = pthread_mutex_lock(&bin_status_mutex);
#endif
      // Send the job in bin to a thread
      bin_status[bin] = SUBMITTED;
      jobs_submitted++;

      if (threads == 0) // if no multi-threading mode, main thread executes job itself
          {
        if (alignByWorker)
          AlignByWorker(bin);
        else
          PerformViterbiByWorker(bin);
        bin_status[bin] = FREE;
        jobs_submitted--;
      }

#ifdef PTHREAD
      // Restart threads waiting for a signal
      rc = pthread_cond_signal(&new_job);

      // Unlock access to bin_status
      rc = pthread_mutex_unlock(&bin_status_mutex);

      if (DEBUG_THREADS)
        fprintf(stderr,
            "Main: put job into bin %i  name=%-7.7s  Jobs running: %i  jobs_submitted:%i \n",
            bin, t[bin]->name, jobs_running, jobs_submitted);
#endif
    }

#ifdef PTHREAD
    if (threads > 0) {
      // Lock mutex
      rc = pthread_mutex_lock(&bin_status_mutex);
      
      // Wait until job finishes and a bin becomes free
      if (jobs_submitted + jobs_running >= bins) {
        if (DEBUG_THREADS)
          fprintf(stderr,
              "Master thread is waiting for jobs to be finished...\n");
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
        clock_gettime(CLOCK_REALTIME, &ts);
        ts.tv_sec += 1;
        rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex, &ts);
#endif
#endif
      }
      // Unlock mutex
      rc = pthread_mutex_unlock(&bin_status_mutex);
    }
#endif

    fclose(dbf);
  }

  // Finished searching all database HHMs
  reading_dbs = 0;

#ifdef PTHREAD
  if (threads > 0) {
    // No more HMMs in database => wait until all threads have finished
    if (DEBUG_THREADS)
      fprintf(stderr,
          "No more jobs read from database         Jobs running:%i  jobs_submitted:%i \n",
          jobs_running, jobs_submitted);

    // Free all threads still waiting for jobs
    rc = pthread_mutex_lock(&bin_status_mutex);
    rc = pthread_cond_broadcast(&new_job);
    rc = pthread_mutex_unlock(&bin_status_mutex); // Unlock mutex

    // Wait for all jobs to finish => join all jobs to main
    for (int j = 0; j < threads; j++) {
      int status;
      pthread_join(pthread[j], (void **) &status);
      if (DEBUG_THREADS)
        fprintf(stderr, "Thread %i finished its work\n", j + 1);
    }
  }
#endif
}
// End DoViterbiSearch() of Viterbi HMM-HMM seach of database
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper around default Viterbi HMM-HMM search (function DoViterbiSearch)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ViterbiSearch(char *dbfiles[], int ndb, int db_size) {
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs = 1;   // needs to be set to 1 before threads are created
  for (bin = 0; bin < bins; bin++)
    bin_status[bin] = FREE;

#ifdef PTHREAD
  // Start threads for database search
  for (int j = 0; j < threads; j++) {
    thread_data[j].thread_id = j + 1;
    thread_data[j].function = &AlignByWorker;
    if (DEBUG_THREADS)
      fprintf(stderr, "Creating worker thread %i ...", j + 1);
    pthread_create(&pthread[j], &joinable, WorkerLoop, (void*) &thread_data[j]);
    if (DEBUG_THREADS)
      fprintf(stderr, " created!\n");
  }
#endif

  // Initialize
  v1 = v;
  if (v > 0 && v <= 3)
    v = 1;
  else
    v -= 2;
  if (print_elapsed)
    ElapsedTimeSinceLastCall("(preparing for search)");

  hitlist.N_searched = db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  DoViterbiSearch(dbfiles, ndb);

  if (v1 >= 2)
    cout << "\n";
  v = v1;

  // Sort list according to sortscore
  if (v >= 3)
    printf("Sorting hit list ...\n");
  hitlist.SortList();

  // Use NN prediction of lamda and mu
  hitlist.CalculatePvalues(q);

  // Calculate E-values as combination of P-value for Viterbi HMM-HMM comparison and prefilter E-value: E = Ndb P (Epre/Ndb)^alpha
  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(q);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant of ViterbiSearch() function for rescoring previously found HMMs with Viterbi algorithm.
// Perform Viterbi search on each hit object in global hash previous_hits, but keep old alignment
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RescoreWithViterbiKeepAlignment(int db_size) {
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs = 1;   // needs to be set to 1 before threads are created
  for (bin = 0; bin < bins; bin++)
    bin_status[bin] = FREE;

#ifdef PTHREAD
  // Start threads for database search
  for (int j = 0; j < threads; j++) {
    thread_data[j].thread_id = j + 1;
    thread_data[j].function = &PerformViterbiByWorker;
    if (DEBUG_THREADS)
      fprintf(stderr, "Creating worker thread %i ...", j + 1);
    pthread_create(&pthread[j], &joinable, WorkerLoop, (void*) &thread_data[j]);
    if (DEBUG_THREADS)
      fprintf(stderr, " created!\n");
  }
#endif

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
    hit_cur = previous_hits->ReadNext();
    if (hit_cur.irep == 1) {
      dbfiles[ndb] = new (char[strlen(hit_cur.dbfile) + 1]);
      strcpy(dbfiles[ndb], hit_cur.dbfile);
      ++ndb;
    }
  }
  
  hitlist.N_searched = db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  DoViterbiSearch(dbfiles, ndb, false);

  if (v1 >= 2)
    cout << "\n";
  v = v1;

  for (int n = 0; n < ndb; ++n)
    delete[] (dbfiles[n]);

  // Sort list according to sortscore
  if (v >= 3)
    printf("Sorting hit list ...\n");
  hitlist.SortList();

  hitlist.CalculatePvalues(q);  // Use NN prediction of lamda and mu

  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(q);
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Realign hits with MAC algorithm
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void perform_realign(char *dbfiles[], int ndb) {
  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
  int nhits = 0;
  int N_aligned = 0;

  // Longest allowable length of database HMM (backtrace: 5 chars, fwd, bwd: 1 double
  long int Lmaxmem = (par.maxmem * 1024 * 1024 * 1024) / sizeof(double) / q->L
      / bins;
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
    hit_cur = hitlist.ReadNext();
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
  
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs = 1;   // needs to be set to 1 before threads are created

  if (print_elapsed)
    ElapsedTimeSinceLastCall(
        "(prepare realign without reallocating/reseting forward/backwad matrices)");

  // (Re)allocate memory for forward matrix
  for (bin = 0; bin < bins; bin++) {
    // Free previously allocated memory (delete and reallocate, since Lmax may have increased)
    if (hit[bin]->forward_allocated)
      hit[bin]->DeleteForwardMatrix(q->L + 2);

    // Allocate memory for matrix and set to 0
    hit[bin]->AllocateForwardMatrix(q->L + 2, Lmax + 1);

    bin_status[bin] = FREE;
  }
  // // If the above reallocation and resetting to 0 is time-critical, it could be avoided by 
  // // replacing the above block with this block: (JS)
  // int Lmaxprev = 0;
  // if (hit[bin]->forward_allocated)  
  //   Lmaxprev = sizeof(hit[bin]->P_MM[0]) / sizeof(hit[bin]->P_MM[0][0]);
  // for (bin=0; bin<bins; bin++)
  //   {
  //     if (!hit[bin]->forward_allocated)  hit[bin]->AllocateForwardMatrix(q->L+2,Lmax+1);
  //     else if (Lmaxprev<Lmax) 
  //    {
  //      hit[bin]->DeleteForwardMatrix(q->L+2);
  //      hit[bin]->AllocateForwardMatrix(q->L+2,Lmax+1);
  //    }
  //     bin_status[bin] = FREE;
  //   }
  // // This is just an idea. This code would need to be tested!!! 
  // // First, it would need to be tested whether reseting P_MM to 0 during the allocation is needed at all.
  
  if (print_elapsed)
    ElapsedTimeSinceLastCall("(reallocate/reset forward/backwad matrices)");
  // if (print_elapsed) ElapsedTimeSinceLastCall("(prepare realign)");

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

    bin = 0;
    nhits = 0;
    hitlist.Reset();

    while (!hitlist.End() && nhits < par.premerge) {
      hit_cur = hitlist.ReadNext();
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

      //TODO
      // Open HMM database file dbfiles[idb]
      FILE* dbf;
      dbf = ffindex_fopen_by_name(dbhhm_data, dbhhm_index, hit_cur.dbfile);
      if (dbf == NULL) {
        if (dba3m_index_file != NULL) {
          dbf = ffindex_fopen_by_name(dba3m_data, dba3m_index, hit_cur.file);
        }
        else {
          cerr << endl << "Error opening " << hit_cur.dbfile
              << ": A3M database missing\n";
          exit(4);
        }
      }

      if (dbf == NULL)
        OpenFileError(hit_cur.dbfile);

      read_from_db = 1;

      // Forward stream position to start of next database HMM to be realigned
      hit[bin]->index = hit_cur.index; // give hit a unique index for HMM
      hit[bin]->ftellpos = hit_cur.ftellpos;
      fseek(dbf, hit_cur.ftellpos, SEEK_SET);
      hit[bin]->dbfile = new (char[strlen(hit_cur.dbfile) + 1]);
      strcpy(hit[bin]->dbfile, hit_cur.dbfile); // record db file name from which next HMM is read
      hit[bin]->irep = 1; // Needed for min_overlap calculation in InitializeForAlignment in hhhit.C

      char path[NAMELEN];
      Pathname(path, hit_cur.dbfile);

      ///////////////////////////////////////////////////
      // Read next HMM from database file
      if (!fgetline(line, LINELEN, dbf)) {
        fprintf(stderr, "Error in %s: end of file %s reached prematurely!\n",
            par.argv[0], hit_cur.dbfile);
        exit(1);
      }
      while (strscn(line) == NULL && fgetline(line, LINELEN, dbf)) {
      } // skip lines that contain only white space

      if (!strncmp(line, "HMMER3", 5))      // read HMMER3 format
          {
        format[bin] = 1;
        read_from_db = t[bin]->ReadHMMer3(dbf, hit_cur.dbfile);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HMMER", 5))      // read HMMER format
          {
        format[bin] = 1;
        read_from_db = t[bin]->ReadHMMer(dbf, hit_cur.dbfile);
        par.hmmer_used = true;
      }
      else if (!strncmp(line, "HH", 2))     // read HHM format
          {
        format[bin] = 0;
        read_from_db = t[bin]->Read(dbf, path);
      }
      else if (!strncmp(line, "NAME", 4)) // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
          {
        format[bin] = 0;
        fseek(dbf, hit_cur.ftellpos, SEEK_SET); // rewind to beginning of line
        read_from_db = t[bin]->Read(dbf, path);
      }
      else if (line[0] == '#' || line[0] == '>')           // read a3m alignment
          {
        Alignment tali;
        tali.Read(dbf, hit_cur.dbfile, line);
        tali.Compress(hit_cur.dbfile);
        //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
        tali.N_filtered = tali.Filter(par.max_seqid_db, par.coverage_db,
            par.qid_db, par.qsc_db, par.Ndiff_db);
        t[bin]->name[0] = t[bin]->longname[0] = t[bin]->fam[0] = '\0';
        tali.FrequenciesAndTransitions(t[bin]);
        format[bin] = 0;
      }
      else {
        cerr << endl << "Error in " << program_name
            << ": unrecognized HMM file format in \'" << hit_cur.dbfile
            << "\'. \n";
        cerr << "Context:\n'" << line << "\n";
        fgetline(line, LINELEN, dbf);
        cerr << line << "\n";
        fgetline(line, LINELEN, dbf);
        cerr << line << "'\n";
        exit(1);
      }
      fclose(dbf);

      if (read_from_db != 1) {
        cerr << "Error in " << par.argv[0] << ": wrong format while reading \'"
            << hit_cur.dbfile << ". Reached end of file while reading HMM "
            << hit_cur.name << " \n";
        exit(1);
      }

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
      PrepareTemplateHMM(q, t[bin], format[bin]);
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
      hit[bin]->Forward(q, t[bin]);
      hit[bin]->Backward(q, t[bin]);
      hit[bin]->MACAlignment(q, t[bin]);
      hit[bin]->BacktraceMAC(q, t[bin]);

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
      hitlist.Delete().Delete();        // delete the list record and hit object
      hitlist.Insert(*hit[bin]);

      // merge only when hit length > MINCOLS_REALIGN (don't merge 1 column matches)
      if (hit[bin]->matched_cols < MINCOLS_REALIGN)
        continue;

      //TODO
      // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
      // Reading in next db MSA and merging it onto Qali
      Alignment Tali;
      FILE* ta3mf = ffindex_fopen_by_name(dba3m_data, dba3m_index, hit[bin]->dbfile);
      if (ta3mf == NULL)
        OpenFileError(hit[bin]->dbfile);
      Tali.Read(ta3mf, hit[bin]->dbfile); // Read template alignment into Tali
      fclose(ta3mf);
      Tali.Compress(hit[bin]->dbfile); // Filter database alignment
      if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
        Qali_allseqs.MergeMasterSlave(*hit[bin], Tali, hit[bin]->dbfile);
      Tali.N_filtered = Tali.Filter(par.max_seqid_db, par.coverage_db,
          par.qid_db, par.qsc_db, par.Ndiff_db);
      Qali.MergeMasterSlave(*hit[bin], Tali, hit[bin]->dbfile);

      // //?????????????????????????????
      // // JS: Reading the db MSA twice for par.all seems pretty inefficient!!!!!!!!!!!
      // FILE* ta3mf;
      // ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
      // if (ta3mf == NULL) OpenFileError(ta3mfile);
      // Qali.MergeMasterSlave(*hit[bin],ta3mfile, ta3mf);
      // fclose(ta3mf);
      // if (par.allseqs)
      //   {
      //     ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
      //     Qali_allseqs.MergeMasterSlave(*hit[bin],ta3mfile, ta3mf, false); // filter db MSA = false
      //     fclose(ta3mf);
      //   }
      // //?????????????????????????????

      // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
      Qali.Compress("merged A3M file");

      // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
      Qali.N_filtered = Qali.Filter(par.max_seqid, par.coverage, par.qid,
          par.qsc, par.Ndiff);

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali.FrequenciesAndTransitions(q);

      stringstream ss_tmp;
      ss_tmp << hit[bin]->file << "__" << hit[bin]->irep;
      premerged_hits->Add((char*) ss_tmp.str().c_str());

      if (par.notags)
        q->NeutralizeTags();

      // Compute substitution matrix pseudocounts?
      if (par.nocontxt) {
        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        q->PreparePseudocounts();
        // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        q->AddAminoAcidPseudocounts();
      }
      else {
        // Add full context specific pseudocounts to query
        q->AddContextSpecificPseudocounts(pc_hhm_context_engine, pc_hhm_context_mode);
      }

      q->CalculateAminoAcidBackground();
      if (par.columnscore == 5 && !q->divided_by_local_bg_freqs)
        q->DivideBySqrtOfLocalBackgroundFreqs(
            par.half_window_size_local_aa_bg_freqs);

      // Transform transition freqs to lin space if not already done
      q->AddTransitionPseudocounts();
      q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done

    }
  }
  
  if (print_elapsed)
    ElapsedTimeSinceLastCall("(premerge)");
  // end premerge
  //////////////////////////////////////////////////////////////////////////////////

#ifdef PTHREAD
  // Start threads for realignment
  for (int j = 0; j < threads; j++) {
    thread_data[j].thread_id = j + 1;
    thread_data[j].function = &RealignByWorker;
    if (DEBUG_THREADS)
      fprintf(stderr, "Creating worker thread %i ...", j + 1);
    pthread_create(&pthread[j], &joinable, WorkerLoop, (void*) &thread_data[j]);
    if (DEBUG_THREADS)
      fprintf(stderr, " created!\n");
  }
#endif

  // Read all HMMs whose position is given in phash_plist_realignhitpos
  for (int idb = 0; idb < ndb; idb++) {

    // Can we skip dbfiles[idb] because it contains no template to be realigned?
    if (!phash_plist_realignhitpos->Contains(dbfiles[idb]))
      continue;

    // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
    // This list is now sorted by ftellpos in ascending order to access one template after the other efficiently
    phash_plist_realignhitpos->Show(dbfiles[idb])->SortList();

    //TODO
    // Open HMM database file dbfiles[idb]
    FILE* dbf;
    dbf = ffindex_fopen_by_name(dbhhm_data, dbhhm_index, dbfiles[idb]);
    if (dbf == NULL) {
      if (dba3m_index_file != NULL) {
        dbf = ffindex_fopen_by_name(dba3m_data, dba3m_index, dbfiles[idb]);
      }
      else {
        cerr << endl << "Error opening " << dbfiles[idb]
            << ": A3M database missing\n";
        exit(4);
      }
    }
    if (dbf == NULL)
      OpenFileError(dbfiles[idb]);

    read_from_db = 1;

    ///////////////////////////////////////////////////////////////////////////////////////
    // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
    phash_plist_realignhitpos->Show(dbfiles[idb])->Reset();
    while (!phash_plist_realignhitpos->Show(dbfiles[idb])->End()) {
      // Submit jobs until no bin is free anymore
      while (!phash_plist_realignhitpos->Show(dbfiles[idb])->End()
          && jobs_submitted + jobs_running < bins) {

        // Allocate free bin
        bin = PickBin(FREE);
        if (bin < 0) {
          fprintf(stderr,
              "Error during realignment: found no free bin! jobs running: %i  jobs_submitted:%i  threads:%i\n",
              jobs_running, jobs_submitted, threads);
          for (bin = 0; bin < bins; bin++)
            fprintf(stderr, "bin_status[%i]=%i\n", bin, bin_status[bin]);
          exit(6);
        }

        // Forward stream position to start of next database HMM to be realigned
        Realign_hitpos hitpos_curr = phash_plist_realignhitpos->Show(
            dbfiles[idb])->ReadNext();
        hit[bin]->index = hitpos_curr.index; // give hit[bin] a unique index for HMM
        fseek(dbf, hitpos_curr.ftellpos, SEEK_SET); // start to read at ftellpos for template

        // Give hit[bin] the pointer to the list of pointers to hitlist elements of same template (for realignment)
        hit[bin]->plist_phits = array_plist_phits[hitpos_curr.index];

        // fprintf(stderr,"dbfile=%-40.40s  index=%-5i  ftellpos=%l\n",dbfiles[idb],hitpos_curr.index,hitpos_curr.ftellpos);

        char path[NAMELEN];
        Pathname(path, dbfiles[idb]);

        ///////////////////////////////////////////////////
        // Read next HMM from database file
        if (!fgetline(line, LINELEN, dbf)) {
          fprintf(stderr, "Error in %s: end of file %s reached prematurely!\n",
              par.argv[0], dbfiles[idb]);
          exit(1);
        }
        while (strscn(line) == NULL && fgetline(line, LINELEN, dbf)) {
        } // skip lines that contain only white space

        if (!strncmp(line, "HMMER3", 5))      // read HMMER3 format
            {
          format[bin] = 1;
          read_from_db = t[bin]->ReadHMMer3(dbf, dbfiles[idb]);
          par.hmmer_used = true;
        }
        else if (!strncmp(line, "HMMER", 5))      // read HMMER format
            {
          format[bin] = 1;
          read_from_db = t[bin]->ReadHMMer(dbf, dbfiles[idb]);
          par.hmmer_used = true;
        }
        else if (!strncmp(line, "HH", 2))     // read HHM format
            {
          format[bin] = 0;
          read_from_db = t[bin]->Read(dbf, path);
        }
        else if (!strncmp(line, "NAME", 4)) // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
            {
          format[bin] = 0;
          fseek(dbf, hitpos_curr.ftellpos, SEEK_SET); // rewind to beginning of line
          read_from_db = t[bin]->Read(dbf, path);
        }
        else if (line[0] == '#' || line[0] == '>')         // read a3m alignment
            {
          Alignment tali;
          tali.Read(dbf, dbfiles[idb], line);
          tali.Compress(dbfiles[idb]);
          // qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
          tali.N_filtered = tali.Filter(par.max_seqid_db, par.coverage_db,
              par.qid_db, par.qsc_db, par.Ndiff_db);
          t[bin]->name[0] = t[bin]->longname[0] = t[bin]->fam[0] = '\0';
          tali.FrequenciesAndTransitions(t[bin]);
          format[bin] = 0;
        }
        else {
          cerr << endl << "Error in " << program_name
              << ": unrecognized HMM file format in \'" << dbfiles[idb]
              << "\'. \n";
          cerr << "Context:\n'" << line << "\n";
          fgetline(line, LINELEN, dbf);
          cerr << line << "\n";
          fgetline(line, LINELEN, dbf);
          cerr << line << "'\n";
          exit(1);
        }

        if (read_from_db == 2)
          continue;  // skip current HMM or reached end of database
        if (read_from_db == 0)
          break;     // finished reading HMMs
        if (v >= 2)
          fprintf(stderr, "Realigning with %s\n", t[bin]->name);
        ///////////////////////////////////////////////////

        hit[bin]->dbfile = new (char[strlen(dbfiles[idb]) + 1]);
        strcpy(hit[bin]->dbfile, dbfiles[idb]); // record db file name from which next HMM is read

        N_aligned++;
        if (v1 >= 2 && !(N_aligned % 10)) {
          cout << ".";
          if (!(N_aligned % 500))
            printf(" %-4i HMMs aligned\n", N_aligned);
          cout.flush();
        }
#ifdef PTHREAD
        // Lock access to bin_status
        if (threads > 0)
          rc = pthread_mutex_lock(&bin_status_mutex);
#endif
        // Send the job in bin to a thread
        bin_status[bin] = SUBMITTED;
        jobs_submitted++;

        if (threads == 0) // if no multi-threading mode, main thread executes job itself
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
          fprintf(stderr,
              "Main: put job into bin %i  name=%-7.7s  Jobs running: %i  jobs_submitted:%i \n",
              bin, t[bin]->name, jobs_running, jobs_submitted);
#endif
      }

#ifdef PTHREAD
      if (threads > 0) {
        // Lock mutex
        rc = pthread_mutex_lock(&bin_status_mutex);

        // Wait until job finishes and a bin becomes free
        if (jobs_submitted + jobs_running >= bins) {
          if (DEBUG_THREADS)
            fprintf(stderr,
                "Master thread is waiting for jobs to be finished...\n");
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
          clock_gettime(CLOCK_REALTIME, &ts);
          ts.tv_sec += 1;
          rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex, &ts);
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
  reading_dbs = 0;
  
#ifdef PTHREAD
  if (threads > 0) {
    // No more HMMs in database => wait until all threads have finished
    if (DEBUG_THREADS)
      fprintf(stderr,
          "No more jobs read from database         Jobs running:%i  jobs_submitted:%i \n",
          jobs_running, jobs_submitted);

    // Free all threads still waiting for jobs
    rc = pthread_mutex_lock(&bin_status_mutex);
    rc = pthread_cond_broadcast(&new_job);
    rc = pthread_mutex_unlock(&bin_status_mutex); // Unlock mutex

    // Wait for all jobs to finish => join all jobs to main
    for (int j = 0; j < threads; j++) {
      int status;
      pthread_join(pthread[j], (void **) &status);
      if (DEBUG_THREADS)
        fprintf(stderr, "Thread %i finished its work\n", j + 1);
    }
  }
#endif
  if (v1 >= 2)
    cout << "\n";
  v = v1;
  
  // Delete all hitlist entries with too short alignments
  nhits = 0;
  hitlist.Reset();
  while (!hitlist.End()) {
    hit_cur = hitlist.ReadNext();

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
      // printf("Deleting alignment of %s with length %i? irep=%i nhits=%-2i  par.B=%-3i  par.Z=%-3i par.e=%.2g par.b=%-3i  par.z=%-3i par.p=%.2g\n",hit_cur.name,hit_cur.matched_cols,hit_cur.irep,nhits,par.B,par.Z,par.e,par.b,par.z,par.p);

      if (v >= 3)
        printf("Deleting alignment of %s with length %i\n", hit_cur.name,
            hit_cur.matched_cols);
      hitlist.Delete().Delete();        // delete the list record and hit object
      // // Make sure only realigned alignments get displayed! JS: Why? better unrealigned than none.
      // if (last_round)
      // if (par.B>par.Z) par.B--; else if (par.B==par.Z) {par.B--; par.Z--;} else par.Z--;
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

	cuticle_init();

  int cluster_found = 0;
  int seqs_found = 0;
  char* argv_conf[MAXOPT]; // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf = 0;               // Number of arguments in argv_conf
  FILE* fin;

#ifdef PTHREAD
  pthread_attr_init(&joinable);  // initialize attribute set with default values
  if (pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE) != 0) // set attribute 'joinable'
    cerr << "Error "
        << pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)
        << ": could not set detach state for thread attibute.\n";
#endif

  par.premerge = 3;
  par.Ndiff = 1000;
  par.prefilter = true;
  par.early_stopping_filter = true;
  early_stopping = new Early_Stopping;
  early_stopping->length = 200;
  early_stopping->thresh = 0.01;
  strcpy(par.outfile, "");
  strcpy(db_ext, "hhm");
  N_searched = 0;
  previous_hits = new Hash<Hit>(1631, hit_cur);
  premerged_hits = new Hash<char>(1631);

  // Make command line input globally available
  par.argv = argv;
  par.argc = argc;
  RemovePathAndExtension(program_name, argv[0]);
  Pathname(program_path, argv[0]);

  // Enable changing verbose mode before command line are processed
  for (int i = 1; i < argc; i++) {
    if (argc > 1 && !strcmp(argv[i], "-v0"))
      v = 0;
    else if (argc > 1 && !strcmp(argv[i], "-v1"))
      v = 1;
    else if (argc > 2 && !strcmp(argv[i], "-v"))
      v = atoi(argv[i + 1]);
  }

  par.SetDefaultPaths(program_path);

  // Process default otpions from .hhdefaults file
  ReadDefaultsFile(argc_conf, argv_conf, program_path);
  ProcessArguments(argc_conf, argv_conf);
  
  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc, argv);

  // Check needed files
  if (!*par.infile || !strcmp(par.infile, "")) {
    help();
    cerr << endl << "Error in " << program_name << ": input file missing!\n";
    exit(4);
  }
  if (!*db_base) {
    help();
    cerr << endl << "Error in " << program_name
        << ": database missing (see -d)\n";
    exit(4);
  }
  if (par.addss == 1 && (!*par.psipred || !*par.psipred_data)) {
    help();
    cerr << endl << "Error in " << program_name
        << ": missing PSIPRED directory (see -psipred and -psipred_data).\nIf you don't need the predicted secondary structure, don't use the -addss option!\n";
    exit(4);
  }
  if (!par.nocontxt) {
    if (!strcmp(par.clusterfile, "")) {
      help();
      cerr << endl << "Error in " << program_name
          << ": context-specific library missing (see -contxt)\n";
      exit(4);
    }
    if (!strcmp(par.cs_library, "")) {
      help();
      cerr << endl << "Error in " << program_name
          << ": column state library (see -cslib)\n";
      exit(4);
    }
  }
  if (par.loc == 0 && num_rounds >= 2 && v >= 1)
    cerr
        << "WARNING: using -global alignment for iterative searches is deprecated since non-homologous sequence segments can easily enter the MSA and corrupt it.\n";
  if (num_rounds < 1)
    num_rounds = 1;
  else if (num_rounds > 8) {
    if (v >= 1)
      cerr << "WARNING: Number of iterations (" << num_rounds
          << ") to large => Set to 8 iterations\n";
    num_rounds = 8;
  }

  // Premerging can be very time-consuming on large database a3ms, such as from pdb70. 
  // Hence it is only done when iteratively searching against uniprot20 or nr20 with their much smaller MSAs:
  if (!(num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile
      || *alis_basename))
    par.premerge = 0;
  
  // No outfile given? Name it basename.hhm
  if (!*par.outfile) {     // outfile not given? Name it basename.hhm
    RemoveExtension(par.outfile, par.infile);
    strcat(par.outfile, ".hhr");
    if (v >= 2)
      cout << "Search results will be written to " << par.outfile << "\n";
  }

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

  fin = fopen(dba3m_data_filename, "r");
  if (fin) {
    fclose(fin);
  }
  else {
    if (errno == EOVERFLOW) {
      cerr << endl;
      cerr << "Error in " << program_name << ": A3M database  "
          << dba3m_data_filename << " too big (>2GB on 32bit system?):" << endl;
      exit(errno);
    }

    if (num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile
        || *alis_basename) {
      cerr << endl << "Error in " << program_name
          << ": Could not open A3M database " << dba3m_data_filename << ", "
          << strerror(errno) << " (needed to construct result MSA)" << endl;
      exit(4);
    }
    dba3m_data_filename[0] = 0;
  }

  q = new HMM;
  q_tmp = new HMM;

  dbfiles_new = new char*[par.maxnumdb_no_prefilter + 1];
  dbfiles_old = new char*[par.maxnumdb + 1];

  early_stopping->evals = new double[early_stopping->length];

  // Prepare index-based databases
  dbhhm_data_file = fopen(dbhhm_data_filename, "r");
  if (!dbhhm_data_file)
    OpenFileError(dbhhm_data_filename);

  dbhhm_index_file = fopen(dbhhm_index_filename, "r");
  if (!dbhhm_index_file)
    OpenFileError(dbhhm_index_filename);

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
      OpenFileError(dba3m_data_filename);

    filesize = CountLinesInFile(dba3m_index_filename);

    dba3m_index_file = fopen(dba3m_index_filename, "r");
    if (!dba3m_index_file)
      OpenFileError(dba3m_index_filename);

    dba3m_index = ffindex_index_parse(dba3m_index_file, filesize);
    if (dba3m_index == NULL) {
      cerr << "Error in " << par.argv[0] << ": could not read index file"
          << dba3m_index_filename << ". Is the file empty or corrupted?\n";
      exit(1);
    }
    dba3m_data = ffindex_mmap_data(dba3m_data_file, &data_size);
  }

  // Check for threads
  if (threads <= 1) {
    threads = 0;
    omp_threads = 1;
  }
  else if (threads > MAXTHREADS) {
    omp_threads = threads = MAXTHREADS;
    if (v >= 1)
      fprintf(stderr, "WARNING: number of CPUs set to maximum value of %i\n",
          MAXTHREADS);
  }
  else
    omp_threads = threads;

  // Set OpenMP threads
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  
  // Check option compatibilities
  if (par.nseqdis>MAXSEQDIS-3-par.showcons) par.nseqdis=MAXSEQDIS-3-par.showcons; //3 reserved for secondary structure
  if (par.aliwidth<20) par.aliwidth=20;
  if (par.pc_hhm_context_engine.pca<0.001) par.pc_hhm_context_engine.pca=0.001; // to avoid log(0)
  if (par.pc_prefilter_context_engine.pca<0.001) par.pc_prefilter_context_engine.pca=0.001; // to avoid log(0)
  if (par.b>par.B) par.B=par.b;
  if (par.z>par.Z) par.Z=par.z;
  if (par.maxmem<1.0) {cerr<<"WARNING: setting -maxmem to its minimum allowed value of 1.0\n"; par.maxmem=1.0;}
  if (par.mact>=1.0) par.mact=0.999; else if (par.mact<0) par.mact=0.0;
  if (par.macins>=1.0) par.macins=0.999; else if (par.macins<0) par.macins=0.0;

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix();

  // Set secondary structure substitution matrix
  if (par.ssm)
    SetSecStrucSubstitutionMatrix();

  // Prepare pseudocounts
  if (!par.nocontxt && *par.clusterfile) {
    char ext[100];
    Extension(ext, par.clusterfile);
    InitializePseudocountsEngine();
  }

  // Prepare column state lib (context size =1 )
  fin = fopen(par.cs_library, "r");
  if (!fin)
    OpenFileError(par.cs_library);
  cs_lib = new cs::ContextLibrary<cs::AA>(fin);
  fclose(fin);
  cs::TransformToLin(*cs_lib);
  
  if (print_elapsed)
    ElapsedTimeSinceLastCall("(prepare CS pseudocounts)");

  v1 = v;
  if (v > 0 && v <= 2)
    v = 1;
  else
    v--;

  // Read query input file (HHM, HMMER, or alignment format) without adding pseudocounts
  Qali.N_in = 0;
  ReadQueryFile(par.infile, input_format, q, &Qali);

  // If input file was not a sequence file (and hence a HHM or HMMER file) AND result MSA/HMM nee to be written, 
  // read in query sequences from a3m file or from representative seqs in HHM file or consensus seq in HMMER file 
  if (Qali.N_in == 0)
    if (num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile
        || *alis_basename || par.premerge > 0)
      ReadQueryA3MFile();
  if (Qali.N_in - Qali.N_ss > 1)
    par.premerge = 0;
  par.M = 1; // all database MSAs must be in A3M format

  if (par.allseqs) {
    Qali_allseqs = Qali; // make a *deep* copy of Qali!
    for (int k = 0; k < Qali_allseqs.N_in; ++k)
      Qali_allseqs.keep[k] = 1; // keep *all* sequences (reset filtering in Qali)
  }

  v = v1;

  if (print_elapsed)
    ElapsedTimeSinceLastCall("(initialize)");

  if (par.prefilter) {
    // Initialize Prefiltering (Get DBsize)
    init_prefilter();
  }
  // Set all HMMs in database as new_dbs
  else {
    init_no_prefiltering();
  }

  if (print_elapsed)
    ElapsedTimeSinceLastCall("(init prefilter)");

  // Input parameters
  if (v >= 3) {
    cout << "Input file       :   " << par.infile << "\n";
    cout << "Output file      :   " << par.outfile << "\n";
    cout << "Prefilter DB     :   " << dbcs_data_filename << " "
        << dbcs_index_filename << "\n";
    cout << "HHM DB           :   " << dbhhm_data_filename << " "
        << dbhhm_index_filename << "\n";
  }

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags)
    q->NeutralizeTags();

  // Prepare multi-threading - reserve memory for threads, intialize, etc.
  if (threads == 0)
    bins = 1;
  else
    bins = iround(threads * 1.2 + 0.5);
  for (bin = 0; bin < bins; bin++) {
    t[bin] = new HMM; // Each bin has a template HMM allocated that was read from the database file
    hit[bin] = new Hit; // Each bin has an object of type Hit allocated ...
    hit[bin]->AllocateBacktraceMatrix(q->L + 2, par.maxres); // ...with a separate dynamic programming matrix (memory!!)
  }
  format = new (int[bins]);

  if (print_elapsed)
    ElapsedTimeSinceLastCall("(finished init)");

  if (par.useCSScoring) {
    columnStateScoring = new ColumnStateScoring();
    columnStateScoring->number_column_states = cs::AS219::kSize;
    columnStateScoring->query_length = q->L;
    columnStateScoring->substitutionScores = new float*[q->L + 1];
    for (int i = 0; i <= q->L; i++) {
      columnStateScoring->substitutionScores[i] = new float[cs::AS219::kSize];
    }
  }
  else {
    columnStateScoring = NULL;
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  for (int round = 1; round <= num_rounds; round++) {

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
    if (*query_hhmfile) {
      v1 = v;
      if (v > 0 && v <= 3)
        v = 1;
      else
        v -= 2;

      // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
      q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
      q->CalculateAminoAcidBackground();

      q->WriteToFile(query_hhmfile);

      v = v1;
    }
    
    PrepareQueryHMM(input_format, q);

    if (print_elapsed)
      ElapsedTimeSinceLastCall("(before prefiltering (pseudocounts))");

    ////////////////////////////////////////////
    // Prefiltering
    ////////////////////////////////////////////

    if (par.prefilter) {
      if (v >= 2)
        printf("Prefiltering database\n");
      prefilter_db();  // in hhprefilter.C
    }
    
    if (print_elapsed)
      ElapsedTimeSinceLastCall("(prefiltering)");

    if (ndb_new == 0) {
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

    //precalculate column scores when using -usecs
    if (par.useCSScoring) {
      for (int i = 1; i <= q->L + 1; ++i) {
        for (int k = 0; k < columnStateScoring->number_column_states; ++k) {
          float sum = 0;
          for (int a = 0; a < 20; ++a) {
            sum += (q->p[i - 1][a] * cs_lib->operator [](k).probs[0][a]
                / q->pav[a]);
          }
          sum = flog2(sum);

          //Fitting of cs scores to non-heuristic scores
          sum = 0.6862403 * sum + 0.0342321 * pow(sum, 2)
              + 0.0002257 * pow(sum, 3) + 0.0006802;

          columnStateScoring->substitutionScores[i - 1][k] = fpow2(sum);
        }
      }
      //precalculateScores(q, columnStateScoring.number_column_states, columnStateScoring.substitutionScores);
    }

    // Main Viterbi HMM-HMM search
    // Starts with empty hitlist (hits of previous iterations were deleted) and creates a hitlist with the hits of this iteration
    ViterbiSearch(dbfiles_new, ndb_new, (ndb_new + ndb_old));

    if (print_elapsed)
      ElapsedTimeSinceLastCall("(Viterbi search)");

    // check for new hits or end with iteration
    int new_hits = 0;
    hitlist.Reset();
    while (!hitlist.End()) {
      hit_cur = hitlist.ReadNext();
      if (hit_cur.Eval > 100.0 * par.e)
        break; // E-value much too large
      if (hit_cur.Eval > par.e)
        continue; // E-value too large
      new_hits++;
    }

    if (new_hits == 0 || round == num_rounds) {
      last_round = true;
      if (round < num_rounds && v >= 2)
        printf("No new hits found in iteration %i => Stop searching\n", round);

      if (ndb_old > 0 && realign_old_hits) {
        if (v > 0) {
          printf("Rescoring previously found HMMs with Viterbi algorithm\n");
        }
        ViterbiSearch(dbfiles_old, ndb_old, (ndb_new + ndb_old));
        // Add dbfiles_old to dbfiles_new for realign
        for (int a = 0; a < ndb_old; a++) {
          dbfiles_new[ndb_new] = new (char[strlen(dbfiles_old[a]) + 1]);
          strcpy(dbfiles_new[ndb_new], dbfiles_old[a]);
          ndb_new++;
        }
      }
      else if (!realign_old_hits && previous_hits->Size() > 0) {
        if (v > 0) {
          printf("Rescoring previously found HMMs with Viterbi algorithm\n");
        }
        RescoreWithViterbiKeepAlignment(ndb_new + previous_hits->Size());

        if (print_elapsed)
          ElapsedTimeSinceLastCall("(Rescoring with Viterbi)");
      }
    }

    // Realign hits with MAC algorithm
    if (par.realign)
      perform_realign(dbfiles_new, ndb_new);

    if (print_elapsed)
      ElapsedTimeSinceLastCall("(realigning with MAC)");

    // Generate alignment for next iteration
    if (round < num_rounds || *par.alnfile || *par.psifile || *par.hhmfile
        || *alis_basename) {
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
          hit_cur = hitlist.ReadNext();
          if (hit_cur.Eval > 100.0 * par.e)
            break; // E-value much too large
          if (hit_cur.Eval > par.e)
            continue; // E-value too large
          if (hit_cur.matched_cols < MINCOLS_REALIGN)
            continue; // leave out too short alignments
          stringstream ss_tmp;
          ss_tmp << hit_cur.file << "__" << hit_cur.irep;
          if (previous_hits->Contains((char*) ss_tmp.str().c_str()))
            continue;  // Already in alignment

          // Add number of sequences in this cluster to total found
          seqs_found += SequencesInCluster(hit_cur.name); // read number after second '|'
          cluster_found++;

          // Skip merging this hit if hit alignment was already merged during premerging
          if (premerged_hits->Contains((char*) ss_tmp.str().c_str()))
            continue;

          //TODO
          // Read a3m alignment of hit from <file>.a3m file
          // Reading in next db MSA and merging it onto Qali
          FILE* ta3mf;
          Alignment Tali;
          ta3mf = ffindex_fopen_by_name(dba3m_data, dba3m_index, hit_cur.dbfile);
          if (ta3mf == NULL)
            OpenFileError(hit_cur.dbfile);
          Tali.Read(ta3mf, hit_cur.dbfile); // Read template alignment into Tali
          fclose(ta3mf);
          Tali.Compress(hit_cur.dbfile); // Filter database alignment
          if (par.allseqs) // need to keep *all* sequences in Qali_allseqs? => merge before filtering
            Qali_allseqs.MergeMasterSlave(hit_cur, Tali, hit_cur.dbfile);
          Tali.N_filtered = Tali.Filter(par.max_seqid_db, par.coverage_db,
              par.qid_db, par.qsc_db, par.Ndiff_db);
          Qali.MergeMasterSlave(hit_cur, Tali, hit_cur.dbfile);

          // //?????????????????????????????
          // // JS: Reading the db MSA twice for par.all seems pretty inefficient!!!!!!!!!!!
          // FILE* ta3mf;
          // ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
          // if (ta3mf == NULL) OpenFileError(ta3mfile);
          // Qali.MergeMasterSlave(hit_cur,ta3mfile, ta3mf);
          // fclose(ta3mf);
          // if (par.allseqs)
          //   {
          //     ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
          //     Qali_allseqs.MergeMasterSlave(hit_cur,ta3mfile, ta3mf, false); // filter db MSA = false
          //     fclose(ta3mf);
          //   }
          // //?????????????????????????????

          if (Qali.N_in >= MAXSEQ)
            break; // Maximum number of sequences reached
        }

        // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
        Qali.Compress("merged A3M file");

        // Sort out the nseqdis most dissimilacd r sequences for display in the result alignments
        Qali.FilterForDisplay(par.max_seqid, par.coverage, par.qid, par.qsc,
            par.nseqdis);

        // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
        float const COV_ABS = 25;     // min. number of aligned residues
        int cov_tot = imax(imin((int) (COV_ABS / Qali.L * 100 + 0.5), 70),
            par.coverage);
        if (v > 2)
          printf("Filter new alignment with cov %3i%%\n", cov_tot);
        Qali.N_filtered = Qali.Filter(par.max_seqid, cov_tot, par.qid, par.qsc,
            par.Ndiff);

        if (print_elapsed)
          ElapsedTimeSinceLastCall("(merge hits to Qali)");
      }

      // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
      Qali.FrequenciesAndTransitions(q, NULL, true);

      if (par.notags)
        q->NeutralizeTags();

      // Calculate SSpred if we need to print out alis after each iteration or if last iteration
      if (par.addss
          && (*alis_basename || round == num_rounds || new_hits == 0)) {
        char ss_pred[par.maxres];
        char ss_conf[par.maxres];

        CalculateSS(q, ss_pred, ss_conf);

        Qali.AddSSPrediction(ss_pred, ss_conf);

        if (print_elapsed)
          ElapsedTimeSinceLastCall("(calculate SS_Pred)");
      }

      if (print_elapsed)
        ElapsedTimeSinceLastCall("(Calculate AA frequencies and transitions)");

      if (*alis_basename) {
        stringstream ss_tmp;
        ss_tmp << alis_basename << "_" << round << ".a3m";
        if (par.allseqs)
          Qali_allseqs.WriteToFile(ss_tmp.str().c_str(), "a3m");
        else
          Qali.WriteToFile(ss_tmp.str().c_str(), "a3m");
      }

      v = v1;

    }
    else if (round == num_rounds) // Update counts for log
        {
      hitlist.Reset();
      while (!hitlist.End()) {
        hit_cur = hitlist.ReadNext();
        if (hit_cur.Eval > 100.0 * par.e)
          break; // E-value much too large
        if (hit_cur.Eval > par.e)
          continue; // E-value too large
        stringstream ss_tmp;
        ss_tmp << hit_cur.file << "__" << hit_cur.irep;
        if (previous_hits->Contains((char*) ss_tmp.str().c_str()))
          continue;  // Already in alignment

        // Add number of sequences in this cluster to total found
        seqs_found += SequencesInCluster(hit_cur.name); // read number after second '|'
        cluster_found++;
      }
    }

    if (v >= 2)
      printf(
          "%i sequences belonging to %i database HMMs found with an E-value < %-12.4g\n",
          seqs_found, cluster_found, par.e);

    if (v >= 2
        && (round < num_rounds || *par.alnfile || *par.psifile || *par.hhmfile
            || *alis_basename))
      printf(
          "Number of effective sequences of resulting query HMM: Neff = %4.2f\n",
          q->Neff_HMM);

    if (q->Neff_HMM > neffmax && round < num_rounds) {
      printf(
          "Diversity is above threshold (%4.2f). Stop searching! (Change threshold using -neffmax <float>.)\n",
          neffmax);
    }

    if (Qali.N_in >= MAXSEQ)
      printf(
          "Maximun number of sequences in query alignment reached (%i). Stop searching!\n",
          MAXSEQ);

    if (new_hits == 0 || round == num_rounds || q->Neff_HMM > neffmax
        || Qali.N_in >= MAXSEQ)
      break;

    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End()) {
      hit_cur = hitlist.ReadNext();
      char strtmp[NAMELEN + 6];
      sprintf(strtmp, "%s__%i%c", hit_cur.file, hit_cur.irep, '\0');
      if (!already_seen_filter || hit_cur.Eval > par.e
          || previous_hits->Contains(strtmp))
        hit_cur.Delete(); // Delete hit object (deep delete with Hit::Delete())
      else
        previous_hits->Add(strtmp, hit_cur);

      // // Old version by Michael => Delete
      // hit_cur = hitlist.ReadNext();
      // stringstream ss_tmp;
      // ss_tmp << hit_cur.file << "__" << hit_cur.irep;
      // if (!already_seen_filter || hit_cur.Eval > par.e || previous_hits->Contains((char*)ss_tmp.str().c_str()))
      //   hit_cur.Delete(); // Delete hit object (deep delete with Hit::Delete())
      // else
      //   previous_hits->Add((char*)ss_tmp.str().c_str(), hit_cur);

      hitlist.Delete(); // Delete list record (flat delete)

    }

    if (print_elapsed)
      ElapsedTimeSinceLastCall("(end of this round)");

  }   // end main loop overs search iterations
  //////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////
  // Result section
  //////////////////////////////////////////////////////////////////////////////////

  // Warn, if HMMER files were used
  if (par.hmmer_used && v >= 1)
    fprintf(stderr,
        "WARNING: Using HMMER files results in a drastically reduced sensitivity (>10%%).\nWe recommend to use HHMs build by hhmake.\n");

  // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
  if (*par.scorefile) {
    if (v >= 3)
      printf("Printing scores file ...\n");
    hitlist.PrintScoreFile(q);
  }

  // Print FASTA or A2M alignments?
  if (*par.pairwisealisfile) {
    if (v >= 2)
      cout << "Printing alignments in "
          << (par.outformat == 1 ? "FASTA" : par.outformat == 2 ? "A2M" : "A3M")
          << " format to " << par.pairwisealisfile << "\n";
    hitlist.PrintAlignments(q, par.pairwisealisfile, par.outformat);
  }

  // Write alignments in tabular layout to alitabfile
  if (*par.alitabfile)
    hitlist.WriteToAlifile(q, alitab_scop);

  // Print summary listing of hits
  if (v >= 3)
    printf("Printing hit list ...\n");
  hitlist.PrintHitList(q_tmp, par.outfile);

  // Write only hit list to screen?
  if (v == 2 && strcmp(par.outfile, "stdout"))
    WriteToScreen(par.outfile, 109); // write only hit list to screen

  // Print alignments of query sequences against hit sequences
  hitlist.PrintAlignments(q_tmp, par.outfile);

  // Write whole output file to screen? (max 10000 lines)
  if (v >= 3 && strcmp(par.outfile, "stdout"))
    WriteToScreen(par.outfile, 10009);

  // Generate result alignment or HMM file?
  if (*par.alnfile || *par.psifile || *par.hhmfile) {
    // Write output PSI-BLAST-formatted alignment?
    if (*par.psifile) {
      if (par.allseqs)
        Qali_allseqs.WriteToFile(par.psifile, "psi");
      else
        Qali.WriteToFile(par.psifile, "psi");
    }

    // Write output HHM file?
    if (*par.hhmfile) {
      // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
      q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
      q->CalculateAminoAcidBackground();

      q->WriteToFile(par.hhmfile);
    }

    // Write output A3M alignment?
    if (*par.alnfile) {
      if (par.allseqs)
        Qali_allseqs.WriteToFile(par.alnfile, "a3m");
      else
        Qali.WriteToFile(par.alnfile, "a3m");
    }

  }

  ////////////////////////////////////////////////////
  // Clean up 
  ////////////////////////////////////////////////////

  fclose(dbhhm_data_file);
  fclose(dbhhm_index_file);
  if (dba3m_index_file != NULL) {
    fclose(dba3m_data_file);
    fclose(dba3m_index_file);
  }
  free(dbhhm_index);
  free(dba3m_index);
  
  // Delete memory for dynamic programming matrix
  for (bin = 0; bin < bins; bin++) {
    hit[bin]->DeleteBacktraceMatrix(q->L + 2);
    if (hit[bin]->forward_allocated)
      hit[bin]->DeleteForwardMatrix(q->L + 2);
    delete hit[bin];
    delete t[bin];
  }
  delete q;
  delete q_tmp;
  if (format)
    delete[] (format);
  if (par.exclstr)
    delete[] par.exclstr;
  delete[] early_stopping->evals;
  delete early_stopping;
  for (int n = 1; n < argc_conf; n++)
    delete[] argv_conf[n];
  if (par.dbfiles)
    delete[] par.dbfiles;
  for (int idb = 0; idb < ndb_new; idb++)
    delete[] (dbfiles_new[idb]);
  for (int idb = 0; idb < ndb_old; idb++)
    delete[] (dbfiles_old[idb]);
  delete[] (dbfiles_new);
  delete[] (dbfiles_old);

  if (par.prefilter) {
    free(length);
    free(first);
    for (size_t n = 0; n < num_dbs; n++)
      delete[] (dbnames[n]);
    free(dbnames);
    fclose(db_data_file);
  }

  if (par.useCSScoring) {
    std::map<std::string, unsigned char*>::iterator it;
    for(it = columnStateSequences.begin(); it != columnStateSequences.end(); it++) {
      delete [] (*it).second;
    }
    columnStateSequences.clear();

    for (int i = 0; i <= columnStateScoring->query_length; i++) {
      delete [] columnStateScoring->substitutionScores[i];
    }
    delete [] columnStateScoring->substitutionScores;
  }

  DeletePseudocountsEngine();

  delete cs_lib;

  // Delete content of hits in hitlist
  hitlist.Reset();
  while (!hitlist.End())
    hitlist.Delete().Delete(); // Delete list record and hit object

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete(); // Delete hit object
  delete previous_hits;

  delete premerged_hits;

#ifdef PTHREAD
  pthread_attr_destroy(&joinable);
  pthread_mutex_destroy(&bin_status_mutex);
  pthread_cond_destroy(&new_job);
  pthread_cond_destroy(&finished_job);
#endif

  if (print_elapsed)
    ElapsedTimeSinceLastCall("(sorting and formatting)");

  // Print 'Done!'
  FILE* outf = NULL;
  if (!strcmp(par.outfile, "stdout"))
    printf("Done!\n");
  else {
    if (*par.outfile) {
      outf = fopen(par.outfile, "a"); //open for append
      fprintf(outf, "Done!\n");
      fclose(outf);
    }
    if (v >= 2)
      printf("Done\n");
  }

  exit(0);
} //end main

//////////////////////////////////////////////////////////////////////////////////////////////////////
// END OF MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////
