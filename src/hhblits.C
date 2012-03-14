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
#include <errno.h>    // perror()
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
#include "abstract_state_matrix.h"
cs::ContextLibrary<cs::AA> *cs_lib;

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure
#include "hhdecl.C"      // Constants, global variables, struct Parameters
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

char line[LINELEN]="";         // input line
int bin;                       // bin index
const char print_elapsed=0;    // debug output for runtimes
char tmp_file[]="/tmp/hhblitsXXXXXX";  // for runtime secondary structure prediction (only with -addss option)

// HHblits variables
const char HHBLITS_REFERENCE[] = "Remmert M., Biegert A., Hauser A., and Soding J.\nHHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.\nNat. Methods 9:173-175 (2011)\n";

int v1=v;                               // verbose mode
int num_rounds   = 2;                   // number of iterations
bool last_round = false;                // set to true in last iteration
bool already_seen_filter = true;        // Perform filtering of already seen HHMs
bool block_filter = true;               // Perform viterbi and forward algorithm only on unshaded tube given by prefiltering
bool realign_old_hits = false;          // Realign old hits in last round or use previous alignments

char input_format = 0;                  // Set to 1 if input in HMMER format (has already pseudocounts)

float neffmax = 10;                     // Break if Neff > Neffmax

int cpu = 2;                            // default: use 2 cores

char config_file[NAMELEN];
char infile[NAMELEN];
char alis_basename[NAMELEN];
char base_filename[NAMELEN];
char query_hhmfile[NAMELEN];

bool alitab_scop = false;                // Write only SCOP alignments in alitabfile

char db_ext[NAMELEN];

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

char db_base[NAMELEN];                   // database basename
char db[NAMELEN];                        // database with context-state sequences
char dba3m[NAMELEN];                     // database with A3M-files
char dbhhm[NAMELEN];                     // database with HHM-files

char** dbfiles_new;
char** dbfiles_old;
int ndb_new=0;
int ndb_old=0;
Hash<Hit>* previous_hits;
Hash<char>* premerged_hits;
 
// HHsearch variables
const int MAXTHREADS=256; // maximum number of threads (i.e. CPUs) for parallel computation
const int MAXBINS=384;    // maximum number of bins (positions in thread queue)
enum bin_states {FREE=0, SUBMITTED=1, RUNNING=2};
int threads=0;            // number of threads (apart from the main thread which reads from the databases file) 0:no multithreading
int bins;                 // number of bins; jobs gets allocated to a FREE bin were they are waiting for execution by a thread
char bin_status[MAXBINS]; // The status for each bin is FREE, SUBMITTED, or RUNNING
int jobs_running;         // number of active jobs, i.e. number of bins set to RUNNING
int jobs_submitted;       // number of submitted jobs, i.e. number of bins set to SUBMITTED
char reading_dbs;         // 1: still HMMs to read in a database;  0: finshed reading all HMMs, no db left
const char DEBUG_THREADS=0; // Debugging flag

HMM* q;                   // Create query HMM with maximum of par.maxres match states
HMM* q_tmp;               // Create query HMM with maximum of par.maxres match states (needed for prefiltering)
HMM* t[MAXBINS];          // Each bin has a template HMM allocated that was read from the database file
Hit* hit[MAXBINS];        // Each bin has an object of type Hit allocated with a separate dynamic programming matrix (memory!!)
Hit hit_cur;              // Current hit when going through hitlist
HitList hitlist;          // list of hits with one Hit object for each pairwise comparison done
int* format;              // format[bin] = 0 if in HHsearch format => add pcs; format[bin] = 1 if in HMMER format => no pcs
int read_from_db;         // The value of this flag is returned from HMM::Read(); 0:end of file  1:ok  2:skip HMM
int N_searched;           // Number of HMMs searched
Alignment Qali;           // output A3M generated by merging A3M alignments for significant hits to the query alignment
Alignment Qali_nodiff;    // output A3M alignment with no sequence filtered out (only active with -nodiff option)

struct Thread_args // data to hand to WorkerLoop thread
{
  int thread_id;          // id of thread (for debugging)
  void (*function)(int);  // pointer to function (=job) to execute by WorkerLoop once old job is done
};

#ifdef PTHREAD
// With this condition variable the main thread signals to the worker threads that it has submitted a new job
pthread_cond_t new_job = PTHREAD_COND_INITIALIZER;

// Mutex assures exclusive access to bin_status[], jobs_sumitted, jobs_running,  and new_job by threads
pthread_mutex_t bin_status_mutex = PTHREAD_MUTEX_INITIALIZER;

// Mutex assures exclusive access to hitlist
pthread_mutex_t hitlist_mutex = PTHREAD_MUTEX_INITIALIZER;

// With this condition variable a worker thread signals to the main thread that it has finished a job
pthread_cond_t finished_job = PTHREAD_COND_INITIALIZER;

struct Thread_args thread_data[MAXTHREADS]; // store a threads thread_id and function to call (AlignByWorker, RealignByWorker)
pthread_t pthread[MAXTHREADS]; // info on thread's structures (needed by system)
pthread_attr_t joinable;       // attribute set for describing threads
int rc;                        // return code for threading commands
#endif

inline int PickBin(char status);


// Include hhworker.C and hhprefilter.C here, because it needs some of the above variables
#include "hhworker.C"      // functions: AlignByWorker, RealignByWorker, WorkerLoop
#include "hhprefilter.C"   // some prefilter functions


/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help(char all=0)
{
  //      ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+---8-----+----9----+----0
  printf("\n");
  printf("HHblits %s:\nHMM-HMM-based lightning-fast iterative sequence search\n",VERSION_AND_DATE);
  printf("HHblits is a sensitive, general-purpose, iterative sequence search tool that represents\n");
  printf("both query and database sequences by HMMs. You can search HHblits databases starting\n");
  printf("with a single query sequence, a multiple sequence alignment (MSA), or an HMM. HHblits\n");
  printf("prints out a ranked list of database HMMs/MSAs and can also generate an MSA by merging\n");
  printf("the significant database HMMs/MSAs onto the query MSA.\n");
  printf("\n");
  printf("%s",HHBLITS_REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query [options] \n",program_name);
  printf(" -i <file>      input query (single FASTA-sequence, A3M- or FASTA-alignment, HMM-file)\n");
  printf("\n");
  printf("Options:                                                                       \n");
  printf(" -d    <base>   database basename (default=%s)                                 \n",db_base);
  printf(" -n     [1,8]   number of iterations (default=%i)                              \n",num_rounds); 
  printf(" -e     [0,1]   E-value cutoff for inclusion in result alignment (def=%G)      \n",par.e);
  printf("\n");
  printf("Input alignment format:                                                       \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               ' -' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  if (all) { 
  printf("Directory paths \n");
  printf(" -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",par.clusterfile);
  printf(" -cslib  <file> column state file for fast database prefiltering (default=%s)\n",par.cs_library);
  printf(" -psipred      <dir>  directory with PSIPRED executables (default=%s)  \n",par.psipred);
  printf(" -psipred_data <dir>  directory with PSIPRED data (default=%s) \n",par.psipred_data);
  printf("\n");
  }
  printf("\n");
  printf("Output options: \n");
  printf(" -o <file>      write results in standard format to file (default=<infile.hhr>)\n");
  printf(" -oa3m <file>   write multiple alignment of significant matches in a3m format\n");
  if (!all) {
  printf("                Analogous for a2m, fas, psi, hhm format (e.g. -ohhm, -ofas)\n");
  }
  if (all) {
  printf(" -ofas <file>   write MSA of significant matches in FASTA format\n");
  printf(" -opsi <file>   write MSA of significant matches in PSI format\n");
  printf(" -ohhm <file>   write HHM file for MSA of significant matches\n");
  }
  printf(" -oalis <name>  write MSAs in A3M format after each iteration\n");
  if (all) {
  printf(" -Ofas <file>   write pairwise alignments of significant matches in FASTA format\n");
  printf(" -qhhm <file>   write query input HHM file of last iteration (default=off)      \n");
  printf(" -seq <int>     max. number of query/template sequences displayed (default=%i)  \n",par.nseqdis);
  printf(" -addss         add predicted 2ndary structure in output alignments             \n");
  printf(" -aliw <int>    number of columns per line in alignment list (default=%i)       \n",par.aliwidth);
  printf(" -p [0,100]     minimum probability in summary and alignment list (default=%G)  \n",par.p);
  printf(" -E [0,inf[     maximum E-value in summary and alignment list (default=%G)      \n",par.E);
  printf(" -Z <int>       maximum number of lines in summary hit list (default=%i)        \n",par.Z);
  printf(" -z <int>       minimum number of lines in summary hit list (default=%i)        \n",par.z);
  printf(" -B <int>       maximum number of alignments in alignment list (default=%i)     \n",par.B);
  printf(" -b <int>       minimum number of alignments in alignment list (def=ault%i)     \n",par.b);
  printf("\n");
  printf("Prefilter options                                                               \n");
  printf(" -nofilter      disable all filter steps                                        \n");
  printf(" -noaddfilter   disable all filter steps (except for fast prefiltering)         \n");
  printf(" -nodbfilter    disable additional filtering of prefiltered HMMs                \n");
  printf(" -noblockfilter search complete matrix in Viterbi                               \n");
  printf(" -maxfilt       max number of hits allowed to pass 2nd prefilter (default=%i)  \n",par.maxnumdb);
  printf("\n");
  }
  printf("Filter result alignment (options can be combined):                              \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)                \n",par.max_seqid);
  printf(" -diff [0,inf[  filter query and db MSAs by selecting most diverse set of sequences,\n");
  printf("                keeping at least this many seq's in each MSA block of length 50 (def=%i)\n",par.Ndiff);
  printf(" -nodiff        do not filter sequences in output alignment (def=off)           \n");
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i)                       \n",par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i)              \n",par.qid);
  printf(" -neff [1,inf]  target diversity of alignment (default=off)                     \n");
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)                 \n",par.qsc);
  printf("\n");
  printf("HMM-HMM alignment options:                                                       \n");
  printf(" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)   \n");
  printf(" -mact [0,1[    posterior probability threshold for MAC re-alignment (def=%.3f)  \n",par.mact);
  printf("                Parameter controls alignment greediness: 0:global >0.1:local     \n");
  printf(" -glob/-loc     use global/local alignment mode for searching/ranking (def=local)\n");
  if (all) {
  printf(" -realign_max <int>  realign max. <int> hits (default=%i)                        \n",par.realign_max);  
  printf(" -alt <int>     show up to this many significant alternative alignments(def=%i)  \n",par.altali);
  printf(" -premerge <int> merge <int> hits to query MSA before aligning remaining hits (def=%i)\n",par.premerge);
  printf(" -shift [-1,1]  profile-profile score offset (def=%-.2f)                         \n",par.shift);
  printf(" -ssm  0-4      0:   no ss scoring                                               \n");
  printf("                1,2: ss scoring after or during alignment  [default=%1i]         \n",par.ssm);
  printf("                3,4: ss scoring after or during alignment, predicted vs. predicted\n");
  printf(" -ssw [0,1]     weight of ss score  (def=%-.2f)                                  \n",par.ssw);
  printf("\n");
  printf("Pseudocount options:                                                             \n");
  printf(" -pcm  0-2      Pseudocount mode (default=%-i)                                   \n",par.pcm);
  printf("                tau = substitution matrix pseudocount admixture                  \n");
  printf("                0: no pseudo counts:     tau = 0                                 \n");
  printf("                1: constant              tau = a                                 \n");
  printf("                2: divergence-dependent: tau = a/(1 + ((Neff-1)/b)^c)            \n");
  printf("                   Neff=( (Neff_q^d+Neff_t^d)/2 )^(1/d)                          \n");
  printf("                   Neff_q = av number of different AAs per column in query       \n");
  printf("                3: constant divergence pseudocounts                              \n");
  printf(" -pca  [0,1]    overall pseudocount admixture (def=%-.1f)                        \n",par.pca);
  printf(" -pcb  [1,inf[  threshold for Neff (def=%-.1f)                                   \n",par.pcb);
  printf(" -pcc  [0,3]    extinction exponent for tau(Neff)  (def=%-.1f)                   \n",par.pcc);
  printf(" -pcw  [0,3]    weight of pos-specificity for pcs  (def=%-.1f)                   \n",par.pcw);
  printf(" -pre_pca [0,1]   PREFILTER pseudocount admixture (def=%-.1f)                    \n",par.pre_pca);
  printf(" -pre_pcb [1,inf[ PREFILTER threshold for Neff (def=%-.1f)                       \n",par.pre_pcb);
  printf("\n");
  printf("Gap cost options:                                                                \n");
  printf(" -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                     \n",par.gapb);
  printf(" -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)    \n",par.gapd);
  printf(" -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)      \n",par.gape);
  printf(" -gapf ]0,inf]  factor to increase/reduce gap open penalty for deletes (def=%-.2f) \n",par.gapf);
  printf(" -gapg ]0,inf]  factor to increase/reduce gap open penalty for inserts (def=%-.2f) \n",par.gapg);
  printf(" -gaph ]0,inf]  factor to increase/reduce gap extend penalty for deletes(def=%-.2f)\n",par.gaph);
  printf(" -gapi ]0,inf]  factor to increase/reduce gap extend penalty for inserts(def=%-.2f)\n",par.gapi);
  printf(" -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f) \n",par.egq);
  printf(" -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)\n",par.egt);
  printf("\n");
  }
  printf("Other options:                                                                   \n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose (def=%i)\n",v);
#ifdef PTHREAD
  printf(" -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=%i)      \n",cpu);
#endif
  if (all) {
  printf(" -neffmax ]1,20] stop iterative search when number of effective sequences Neff in \n");
  printf("                evolving query MSA becomes larger than neffmax (default=%.1f) \n",neffmax); 
  printf(" -scores <file> write scores for all pairwise comparisions to file               \n");
  printf(" -atab   <file> write all alignments in tabular layout to file                   \n");
  printf(" -maxres <int>  max number of HMM columns (def=%5i)             \n",par.maxres);
  printf(" -maxmem [1,inf[ max available memory in GB (def=%.1f)          \n",par.maxmem);
  } 
#ifndef PTHREAD
  printf("(The -cpu option is inactive since HHblits was not compiled with POSIX thread support)\n");
#endif
  printf("\n");
  if (!all) {
  printf("An extended list of options can be obtained by using '-help all' as parameter    \n");
  }
  printf("\n");
  printf("Example: %s -i query.fas -oa3m query.a3m -n 2  \n",program_name);
  cout<<endl;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
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
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database basename following -d\n"; exit(4);}
          else
	    strcpy(db_base,argv[i]);
        }
      else if (!strcmp(argv[i],"-contxt") || !strcmp(argv[i],"-context_data"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no lib following -contxt\n"; exit(4);}
          else
	    strcpy(par.clusterfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-cslib") || !strcmp(argv[i],"-cs_lib"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no lib following -cslib\n"; exit(4);}
          else strcpy(par.cs_library,argv[i]);
        }
      else if (!strcmp(argv[i],"-psipred"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no directory following -psipred\n"; exit(4);}
          else
	    strcpy(par.psipred,argv[i]);
        }
      else if (!strcmp(argv[i],"-psipred_data"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database directory following -psipred_data\n"; exit(4);}
          else
	    strcpy(par.psipred_data,argv[i]);
        }
      else if (!strcmp(argv[i],"-o"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.outfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-oa3m"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -oa3m\n"; exit(4);}
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
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -opsi\n"; exit(4);}
          else strcpy(par.psifile,argv[i]);
        }
      else if (!strcmp(argv[i],"-oalis"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file basename following -oalis\n"; exit(4);}
          else strcpy(alis_basename,argv[i]);
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
      else if (!strcmp(argv[i],"-qhhm"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no filename following -qhhm\n"; exit(4);}
          else strcpy(query_hhmfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-scores"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file following -scores\n"; exit(4);}
          else {strcpy(par.scorefile,argv[i]);}
        }
      else if (!strcmp(argv[i],"-db_ext"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no extension following -db_ext\n"; exit(4);}
          else {strcpy(db_ext,argv[i]);}
        }
      else if (!strcmp(argv[i],"-atab"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file following -atab\n"; exit(4);}
          else {strcpy(par.alitabfile,argv[i]);}
        }
      else if (!strcmp(argv[i],"-atab_scop")) alitab_scop=true;
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"-help"))
        {
          if (++i>=argc || argv[i][0]=='-') {help(); exit(0);}
          if (!strcmp(argv[i],"all")) {help(1); exit(0);}
          else {help(); exit(0);}
        }
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-n") && (i<argc-1)) num_rounds = atoi(argv[++i]); 
      else if (!strncmp(argv[i],"-BLOSUM",7) || !strncmp(argv[i],"-Blosum",7))
        {
          if (!strcmp(argv[i]+7,"30")) par.matrix=30;
          else if (!strcmp(argv[i]+7,"40")) par.matrix=40;
          else if (!strcmp(argv[i]+7,"50")) par.matrix=50;
          else if (!strcmp(argv[i]+7,"62")) par.matrix=62;
          else if (!strcmp(argv[i]+7,"65")) par.matrix=65;
          else if (!strcmp(argv[i]+7,"80")) par.matrix=80;
          else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
        }
      else if (!strcmp(argv[i],"-M") && (i<argc-1))
        if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1;
        else if(!strcmp(argv[i],"first"))  par.M=3;
        else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
        else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strcmp(argv[i],"-p") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-P") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-E") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-b") && (i<argc-1)) par.b = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-B") && (i<argc-1)) par.B = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-z") && (i<argc-1)) par.z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-Z") && (i<argc-1)) par.Z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-realign_max") && (i<argc-1)) par.realign_max = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-e") && (i<argc-1)) par.e = atof(argv[++i]);
      else if (!strncmp(argv[i],"-nopred",7) || !strncmp(argv[i],"-noss",5)) par.showpred=0;
      else if (!strncmp(argv[i],"-addss",6)) par.addss=1;
      else if (!strcmp(argv[i],"-seq") && (i<argc-1))  par.nseqdis=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-aliw") && (i<argc-1)) par.aliwidth=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-nodiff")) par.nodiff=true;
      else if (!strcmp(argv[i],"-neffmax") && (i<argc-1)) neffmax=atof(argv[++i]); 
      else if ((!strcmp(argv[i],"-neff") || !strcmp(argv[i],"-Neff")) && (i<argc-1)) par.Neff=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pcm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pca=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pcb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pcc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcw") && (i<argc-1)) par.pcw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pre_pca") && (i<argc-1)) par.pre_pca=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pre_pcb") && (i<argc-1)) par.pre_pcb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;}
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egq") && (i<argc-1)) par.egq=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egt") && (i<argc-1)) par.egt=atof(argv[++i]);
      else if (!strcmp(argv[i],"-alphaa") && (i<argc-1)) par.alphaa=atof(argv[++i]);
      else if (!strcmp(argv[i],"-alphab") && (i<argc-1)) par.alphab=atof(argv[++i]);
      else if (!strcmp(argv[i],"-alphac") && (i<argc-1)) par.alphac=atof(argv[++i]);
      else if (!strcmp(argv[i],"-filterlen") && (i<argc-1)) 
	{
	  par.filter_length=atoi(argv[++i]);
	  delete par.filter_evals;
	  par.filter_evals=new double[par.filter_length];
	}
      else if (!strcmp(argv[i],"-filtercut") && (i<argc-1)) par.filter_thresh=(double)atof(argv[++i]);
      else if (!strcmp(argv[i],"-nofilter")) {par.prefilter=false; already_seen_filter=false; block_filter=false; par.early_stopping_filter=false; par.filter_thresh=0;}
      else if (!strcmp(argv[i],"-noaddfilter")) {already_seen_filter=false; block_filter=false; par.early_stopping_filter=false; par.filter_thresh=0;}
      else if (!strcmp(argv[i],"-nodbfilter")) {par.filter_thresh=0;}
      else if (!strcmp(argv[i],"-noblockfilter")) {block_filter=false;}
      else if (!strcmp(argv[i],"-noearlystoppingfilter")) {par.early_stopping_filter=false;}
      else if (!strcmp(argv[i],"-block_len") && (i<argc-1)) par.block_shading_space = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-maxfilt") && (i<argc-1)) par.maxnumdb = par.maxnumdb_no_prefilter = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-shading_mode") && (i<argc-1)) strcpy(par.block_shading_mode,argv[++i]);
      else if (!strcmp(argv[i],"-prepre_smax_thresh") && (i<argc-1)) par.preprefilter_smax_thresh = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pre_evalue_thresh") && (i<argc-1)) par.prefilter_evalue_thresh = atof(argv[++i]);
      else if (!strcmp(argv[i],"-pre_bitfactor") && (i<argc-1)) par.prefilter_bit_factor = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pre_gap_open") && (i<argc-1)) par.prefilter_gap_open = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pre_gap_extend") && (i<argc-1)) par.prefilter_gap_extend = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pre_score_offset") && (i<argc-1)) par.prefilter_score_offset = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-realignoldhits")) realign_old_hits=true;
      else if (!strcmp(argv[i],"-realign")) par.realign=1;
      else if (!strcmp(argv[i],"-norealign")) par.realign=0;
      else if (!strcmp(argv[i],"-ssm") && (i<argc-1)) par.ssm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-ssw") && (i<argc-1)) par.ssw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-maxres") && (i<argc-1)) {
	par.maxres=atoi(argv[++i]);
	par.maxcol=2*par.maxres;
      }
      else if (!strncmp(argv[i],"-glo",3)) {par.loc=0; if (par.mact>0.35 && par.mact<0.351) {par.mact=0;} }
      else if (!strncmp(argv[i],"-loc",4)) par.loc=1;
      else if (!strncmp(argv[i],"-alt",4) && (i<argc-1)) par.altali=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-shift") && (i<argc-1)) par.shift=atof(argv[++i]);
      else if ((!strcmp(argv[i],"-mact") || !strcmp(argv[i],"-mapt")) && (i<argc-1)) par.mact=atof(argv[++i]);
      else if (!strcmp(argv[i],"-scwin") && (i<argc-1)) {par.columnscore=5; par.half_window_size_local_aa_bg_freqs = imax(1,atoi(argv[++i]));}
      else if (!strncmp(argv[i],"-cpu",4) && (i<argc-1)) { threads=atoi(argv[++i]); cpu = threads;}
      else if (!strcmp(argv[i],"-maxmem") && (i<argc-1)) {par.maxmem=atof(argv[++i]);}
      else if (!strncmp(argv[i],"-premerge",9) && (i<argc-1)) par.premerge=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-csb") && (i<argc-1)) par.csb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-csw") && (i<argc-1)) par.csw=atof(argv[++i]);
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}

///////////////////////////////////////////////////////////////////////////////////////
//// For multi-threading: return a bin with the desired status, return -1 if no such bin found
//////////////////////////////////////////////////////////////////////////////////////
inline int PickBin(char status)
{
  for (int b=0; b<bins; b++) {if (bin_status[b]==status) return b;}
  return -1;
}

///////////////////////////////////////////////////////////////////////////////////////
//// Do the pairwise comparison of q and *(t[bin]) for the database search
//// Combination of RealignByWorker and AlignByWorker: 
//// Picks hits found in previous iterations and recalculates Viterbi scores using 
//// query profile from last iteration while KEEPING original (MAC) alignment.
//////////////////////////////////////////////////////////////////////////////////////
void PerformViterbiByWorker(int bin)
{
  // Prepare q ant t and compare
  PrepareTemplate(*q,*(t[bin]),format[bin]);

  // Do HMM-HMM comparison
  for (hit[bin]->irep=1; hit[bin]->irep<=par.altali; hit[bin]->irep++)
    {
      // Break, if no previous_hit with irep is found
      hit[bin]->Viterbi(*q,*(t[bin]));
      if (hit[bin]->irep>1 && hit[bin]->score <= SMIN) break;
      hit[bin]->Backtrace(*q,*(t[bin]));
      
      hit[bin]->score_sort = hit[bin]->score_aass;
      //printf("PerformViterbiByWorker:   %-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit[bin]->name,hit[bin]->fam,hit[bin]->irep,hit[bin]->score);

#ifdef PTHREAD
      pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif
     
      stringstream ss_tmp;
      ss_tmp << hit[bin]->name << "__" << hit[bin]->irep;

      if (previous_hits->Contains((char*)ss_tmp.str().c_str()))
	{
	  //printf("Previous hits contains %s!\n",(char*)ss_tmp.str().c_str());
	  hit_cur = previous_hits->Remove((char*)ss_tmp.str().c_str());
	  previous_hits->Add((char*)ss_tmp.str().c_str(), *(hit[bin]));
	  
	  // Overwrite *hit[bin] with alignment, etc. of hit_cur
	  hit_cur.score      = hit[bin]->score;
	  hit_cur.score_aass = hit[bin]->score_aass;
	  hit_cur.score_ss   = hit[bin]->score_ss;
	  hit_cur.Pval       = hit[bin]->Pval;
	  hit_cur.Pvalt      = hit[bin]->Pvalt;
	  hit_cur.logPval    = hit[bin]->logPval;
	  hit_cur.logPvalt   = hit[bin]->logPvalt;
	  hit_cur.Eval       = hit[bin]->Eval;
	  hit_cur.logEval    = hit[bin]->logEval;
	  hit_cur.Probab     = hit[bin]->Probab;

	  hitlist.Push(hit_cur);            // insert hit at beginning of list (last repeats first!)
	  
	}
      else
	{
	  // don't save alignments which where not found in previous rounds

	  //printf("Don't save %s!\n",(char*)ss_tmp.str().c_str());
	  //hitlist.Push(*(hit[bin]));          // insert hit at beginning of list (last repeats first!)
	}
	  

#ifdef PTHREAD
      pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif

      if (hit[bin]->score <= SMIN) break;  // break if score for first hit is already worse than SMIN
    }

  return;
}

//????????????????????????????????????????????????????????????????????????????????????????????????????????????
// Why is the framed code below needed?? Filetype will be checked in ReadInput()!
// Suggestion: 
// Rename this function as PrepareQa3mFile() and rename ReadInput() as ReadQueryFile(). 
// Throw out the framed code below, including the ReadInput() call,
// and replace the old call to ReadInputFile() in hhblits.C by { ReadQueryFile(); PrepareQa3mFile();} .
// Comment Michael Remmert: 
// "In hhblits.C gab es glaube ich 2 Gründe, warum ich den Filetype schon vorher checke:
// Zum einen habe ich da abgefangen, ob der Input im HMMER-Format ist (input_format == 1 && par.hmmer_used = true) 
// und das dann entsprechend ausgegeben.
// Zum anderen habe ich, wenn der Input aus nur einer Sequenz besteht, den Parameter par.M auf 3 gesetzt, 
// was hat den Vorteil, dass es bei einer einzelnen Sequenz, die nur aus Kleinbuchstaben besteht, 
// HHblits einfach alle Buchstaben der Sequenz als Match-States annimmt:
// if (num_seqs == 1 && par.M == 1) par.M=3; // if only single sequence in input file,  
//                                            //use par.M=3 (match states by first seq)"

/////////////////////////////////////////////////////////////////////////////////////
// Read input file
/////////////////////////////////////////////////////////////////////////////////////
void ReadInputFile()
{
  int num_seqs = 0;


  //????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // Remove framed code?

  FILE* qf=fopen(par.infile,"rb");
  if (!qf) OpenFileError(par.infile);

  char qa3mfile[NAMELEN];
  RemoveExtension(qa3mfile,par.infile);
  strcat(qa3mfile,".a3m");

  if (!fgetline(line,LINELEN,qf)) {help(); cerr<<endl<<"Error in "<<program_name<<": cannot read input file!\n"; exit(4);}
  if (!strncmp(line,"HMMER3",6) || !strncmp(line,"HMMER",5))  // HMMER/HMMER3 format
    {
      input_format = 1;
      par.hmmer_used = true;
      cerr<<endl<<"Error in "<<program_name<<": HMMER format as input not supported!\n";
      exit(1);
    }
  else if (!strncmp(line,"HH",2))     // HHM format
    {
      input_format = 0;
    }
  else if (line[0]=='#' || line[0]=='>')             // read sequence/alignment 
    {
      input_format = 0;
      strcpy(qa3mfile, par.infile);
      FILE* inf=fopen(par.infile,"r");
      while (fgetline(line,LINELEN,inf))
      	  if (line[0] == '>')
      	      num_seqs++;
      if (num_seqs == 1 && par.M == 1) par.M=3; // if only single sequence in input file, use par.M=3 (match states by first seq)
      fclose(inf);
    }
  else
    {
      cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<par.infile<<"\'. \n";
      cerr<<"Context:\n'"<<line<<"\n";
      fgetline(line,LINELEN,qf); cerr<<line<<"\n";
      fgetline(line,LINELEN,qf); cerr<<line<<"'\n";
      exit(1);
    }
  fclose(qf);


  v1=v;
  if (v>0 && v<=3) v=1; else v-=2;

  // Read input file (HMM or alignment format) without adding pseudocounts
  ReadInput(par.infile, *q);

  //????????????????????????????????????????????????????????????????????????????????????????????????????????????

  // Read in query alignment
  FILE* qa3mf=fopen(qa3mfile,"r");
  if (!qa3mf) 
    {
      // Read query alignment from HHM	
      Qali.GetSeqsFromHMM(*q);
      Qali.Compress("compress Qali");

      if (num_rounds > 1 || *par.alnfile || *par.psifile || *par.hhmfile || *alis_basename)
	{
	  if (input_format == 0 && v>=1)   // HHM format
	    cerr<<"WARNING: No alignment-file found, use only representative seqs from HHM-file as base alignment for evolving alignment!\n";
	  else if (input_format == 1 && v>=1)
	    cerr<<"WARNING: No alignment-file found, use only consensus sequence from HMMER-file as base alignment for evolving alignment!\n";
	}
    } 
  else 
    {
      if (num_seqs != 1) 
      	par.premerge=0;
      
      Qali.Read(qa3mf,qa3mfile);

      // Warn, if there are gaps in a single sequence
      if (num_seqs == 1 && par.M != 2) {
	int num_gaps = strtr(Qali.seq[0], "-", "-");
	if (num_gaps > 1 && v>=1) {  // 1 gap is always given at array pos 0
	  fprintf(stderr, "WARNING: Your input sequence contains gaps. These gaps will be ignored in this search!\nIf you wan't to make HHblits treat these as match states, you could start HHblits with the '-M 100' option.\n");
	}
      }

      Qali.Compress("compress Qali");
      fclose(qa3mf);

      delete[] Qali.longname;
      Qali.longname = new(char[strlen(q->longname)+1]);
      strcpy(Qali.longname,q->longname);
      strcpy(Qali.name,q->name);
      strcpy(Qali.fam,q->fam);
    }

  // Get basename
  RemoveExtension(base_filename,par.infile);
  if (!*par.outfile)      // outfile not given? Name it basename.hhm
    {
      strcpy(par.outfile,base_filename);
      strcat(par.outfile,".hhr");
      if (v>=2) cout<<"Search results will be written to "<<par.outfile<<"\n";
    }

  par.M=1;
  v=v1;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Perform Viterbi HMM-HMM search on all db HMM names in dbfiles
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void search_loop(char *dbfiles[], int ndb, bool alignByWorker=true)
{
  // Search databases
  for (bin=0; bin<bins; bin++) {hit[bin]->realign_around_viterbi=false;}

  double filter_cutoff = par.filter_length*par.filter_thresh;
  par.filter_sum=par.filter_length;
  par.filter_counter=0;
  for (int a=0; a<par.filter_length; a++) {par.filter_evals[a]=1;}

  // For all the databases comming through prefilter
  for (int idb=0; idb<ndb; idb++)
    {
      // Check early stopping filter
      if (par.early_stopping_filter && par.filter_sum < filter_cutoff)
	{
	  if (v>=4)
	    printf("Stop after DB-HHM %i from %i (filter_sum: %8.4f   cutoff: %8.4f)\n",idb,ndb,par.filter_sum,filter_cutoff);
	  printf("\n");
	  break;
	}
      
      // Open HMM database
      //cerr<<"\nReading db file "<<idb<<" dbfiles[idb]="<<dbfiles[idb]<<"\n";
      //FILE* dbf=fopen(dbfiles[idb],"rb");
      FILE* dbf;
      dbf = ffindex_fopen(dbhhm_data, dbhhm_index, dbfiles[idb]);
      if (dbf == NULL) {
	char filename[NAMELEN];
	RemoveExtension(filename, dbfiles[idb]);
	strcat(filename,".a3m");
	if(dba3m_index_file!=NULL) {
	  dbf = ffindex_fopen(dba3m_data, dba3m_index, filename);
	} else {
	  cerr<<endl<<"Error opening "<<dbfiles[idb]<<": A3M database missing\n"; exit(4);
	}	
	if (dbf == NULL) 
	  {
	    RemoveExtension(filename, dbfiles[idb]);
	    strcat(filename,".hmm");
	    if(dbhhm_index_file!=NULL) {
	      dbf = ffindex_fopen(dbhhm_data, dbhhm_index, filename);
	    } else {
	      cerr<<endl<<"Error opening "<<dbfiles[idb]<<": HHM database missing\n"; exit(4);
	    }
	    dbf = ffindex_fopen(dbhhm_data, dbhhm_index, dbfiles[idb]);
	  }
      }
      if (dbf == NULL) OpenFileError(dbfiles[idb]);
      
      // Submit jobs if bin is free
      if (jobs_submitted+jobs_running<bins)
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
	  //	  fprintf(stderr,"dbfile=%-40.40s  index=%-5i  ftellpos=%i\n",dbfiles[idb],N_searched,(int) hit[bin]->ftellpos);

	  char path[NAMELEN];
	  Pathname(path,dbfiles[idb]);

	  ///////////////////////////////////////////////////
	  // Read next HMM from database file
	  if (!fgetline(line,LINELEN,dbf)) {continue;}
	  while (strscn(line)==NULL && fgetline(line,LINELEN,dbf)) {} // skip lines that contain only white space

	  if (!strncmp(line,"HMMER3",6))      // read HMMER3 format
	    {
	      format[bin] = 1;
	      t[bin]->ReadHMMer3(dbf,dbfiles[idb]);
	      par.hmmer_used = true;
	    }
	  else if (!strncmp(line,"HMMER",5))      // read HMMER format
	    {
	      format[bin] = 1;
	      t[bin]->ReadHMMer(dbf,dbfiles[idb]);
	      par.hmmer_used = true;
	    }
	  else if (!strncmp(line,"HH",2))    // read HHM format
	    {
	      format[bin] = 0;
	      t[bin]->Read(dbf,path);
	    }
	  else if (!strncmp(line,"NAME",4))  // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
	    {
	      fseek(dbf,hit[bin]->ftellpos,SEEK_SET); // rewind to beginning of line
	      format[bin] = 0;
	      t[bin]->Read(dbf,path);
	    }
	  else if (line[0]=='#' || line[0]=='>')             // read a3m alignment
	    {
	      Alignment tali;
	      tali.Read(dbf,dbfiles[idb],line);
	      tali.Compress(dbfiles[idb]);
	      //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
	      tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
	      char wg=par.wg; par.wg=1; // use global weights
	      t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
	      tali.FrequenciesAndTransitions(*(t[bin]));
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
	  if (v>=4) printf("Aligning with %s\n",t[bin]->name);  /////////////////////v>=4
	  ///////////////////////////////////////////////////
	  
	  hit[bin]->dbfile = new(char[strlen(dbfiles[idb])+1]);
	  strcpy(hit[bin]->dbfile,dbfiles[idb]); // record db file name from which next HMM is read
	  
	  ++N_searched;
	  if (v1>=2 && !((idb+1)%20))
	    {
	      cout<<".";
	      if (!((idb+1)%1000)) printf(" %-4i HMMs searched\n",(idb+1));
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
      
      fclose(dbf);
    }

  // Finished searching all database HHMs
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
}
// End search_loop() of Viterbi HMM-HMM seach of database
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper around default Viterbi HMM-HMM search (function search_loop)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void search_database(char *dbfiles[], int ndb, int db_size)
{
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs=1;   // needs to be set to 1 before threads are created
  for (bin=0; bin<bins; bin++)
    bin_status[bin] = FREE;

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

  // Initialize
  v1=v;
  if (v>0 && v<=3) v=1; else v-=2;
  if (print_elapsed) ElapsedTimeSinceLastCall("(preparing for search)");

  hitlist.N_searched=db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  search_loop(dbfiles,ndb);

  if (v1>=2) cout<<"\n";
  v=v1;

  if (print_elapsed) ElapsedTimeSinceLastCall("(search through database)");

  // Sort list according to sortscore
  if (v>=3) printf("Sorting hit list ...\n");
  hitlist.SortList();

  // Use NN prediction of lamda and mu
  hitlist.CalculatePvalues(*q);  

  // Calculate E-values as combination of P-value for Viterbi HMM-HMM comparison and prefilter E-value: E = Ndb P (Epre/Ndb)^alpha
  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(*q);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variant of search_database() function:
// Perform Viterbi search on each hit object in global hash previous_hits, but keep old alignment
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void perform_viterbi_search(int db_size)
{
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs=1;   // needs to be set to 1 before threads are created
  for (bin=0; bin<bins; bin++)
    bin_status[bin] = FREE;

#ifdef PTHREAD
  // Start threads for database search
  for (int j=0; j<threads; j++)
    {
      thread_data[j].thread_id = j+1;
      thread_data[j].function  = &PerformViterbiByWorker;
      if (DEBUG_THREADS) fprintf(stderr,"Creating worker thread %i ...",j+1);
      pthread_create(&pthread[j], &joinable, WorkerLoop, (void*)&thread_data[j]);
      if (DEBUG_THREADS) fprintf(stderr," created!\n");
    }
#endif

  // Initialize
  v1=v;
  if (v>0 && v<=3) v=1; else v-=2;

  char *dbfiles[db_size+1];
  int ndb = 0;

  par.block_shading->Reset();
  while (!par.block_shading->End())
    delete[] (par.block_shading->ReadNext()); 
  par.block_shading->New(16381,NULL);

  // Get dbfiles of previous hits
  previous_hits->Reset();
  while (!previous_hits->End())
    {
      hit_cur = previous_hits->ReadNext();
      if (hit_cur.irep==1)  
	{
	  dbfiles[ndb]=new(char[strlen(hit_cur.dbfile)+1]);
	  strcpy(dbfiles[ndb],hit_cur.dbfile);
	  ++ndb;
	}

      // Seach only around previous HMM-HMM alignment (not the prefilter alignment as in search_database())
      if (block_filter)
	{
	  // Hash par.block_shading contains pointer to int array which is here called block.
	  //        This array contains start and stop positions of alignment for shading
	  // Hash par.block_shading_counter contains currently next free position in par.block_shading int array

	  //printf("Viterbi hit %s   q: %i-%i   t: %i-%i\n",hit_cur.name, hit_cur.i1, hit_cur.i2, hit_cur.j1, hit_cur.j2);
	  int* block;
	  int counter;
	  if (par.block_shading->Contains(hit_cur.name))
	    {
	      block = par.block_shading->Show(hit_cur.name); 
	      counter = par.block_shading_counter->Remove(hit_cur.name);
	    }
	  else
	    {
	      block = new(int[400]);
	      counter = 0;
	    }
	  block[counter++] = hit_cur.i1;
	  block[counter++] = hit_cur.i2;
	  block[counter++] = hit_cur.j1;
	  block[counter++] = hit_cur.j2;
	  par.block_shading_counter->Add(hit_cur.name,counter);
	  if (!par.block_shading->Contains(hit_cur.name))
	    par.block_shading->Add(hit_cur.name,block);
	  // printf("Add to block shading   key: %s    data:",hit_cur.name);
	  // for (int i = 0; i < counter; i++)
	  //   printf(" %i,",block[i]);
	  // printf("\n");
	}
    }
  
  hitlist.N_searched=db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  //////////////////////////////////////////////////////////
  // Start Viterbi search through db HMMs listed in dbfiles
  search_loop(dbfiles,ndb,false);

  if (v1>=2) cout<<"\n";
  v=v1;

  for (int n=0;n<ndb;++n) delete[](dbfiles[n]);

  if (print_elapsed) ElapsedTimeSinceLastCall("(search through database)");

  // Sort list according to sortscore
  if (v>=3) printf("Sorting hit list ...\n");
  hitlist.SortList();

  hitlist.CalculatePvalues(*q);  // Use NN prediction of lamda and mu
  
  if (par.prefilter)
    hitlist.CalculateHHblitsEvalues(*q);
  
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Realign hits with MAC algorithm
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void perform_realign(char *dbfiles[], int ndb)
{
  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
  int nhits=0;
  int N_aligned=0;

  // Longest allowable length of database HMM (backtrace: 5 chars, fwd: 1 double, bwd: 1 double 
  long int Lmaxmem=((par.maxmem-0.5)*1024*1024*1024)/(2*sizeof(double)+8)/q->L/bins;
  long int Lmax=0;      // length of longest HMM to be realigned
    
  par.block_shading->Reset();
  while (!par.block_shading->End())
    delete[] (par.block_shading->ReadNext()); 
  par.block_shading->New(16381,NULL);
  par.block_shading_counter->New(16381,0);

  // phash_plist_realignhitpos->Show(dbfile) is pointer to list with template indices and their ftell positions.
  // This list can be sorted by ftellpos to access one template after the other efficiently during realignment
  Hash< List<Realign_hitpos>* >* phash_plist_realignhitpos;
  phash_plist_realignhitpos = new Hash< List<Realign_hitpos>* > (30031,NULL);
  
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
      if (hit_cur.L<=Lmaxmem)
	{
	  // Seach only around viterbi hit
	  if (block_filter)
	    {
	      // Hash par.block_shading contains pointer to int array which is here called block.
	      //        This array contains start and stop positions of alignment for shading
	      // Hash par.block_shading_counter contains currently next free position in par.block_shading int array
	      
	      //printf("Viterbi hit %s   q: %i-%i   t: %i-%i\n",hit_cur.name, hit_cur.i1, hit_cur.i2, hit_cur.j1, hit_cur.j2);
	      int* block;
	      int counter;
	      if (par.block_shading->Contains(hit_cur.name))
		{
		  block = par.block_shading->Show(hit_cur.name);
		  counter = par.block_shading_counter->Remove(hit_cur.name);
		}
	      else
		{
		  block = new(int[400]);
		  counter = 0;
		}
	      block[counter++] = hit_cur.i1;
	      block[counter++] = hit_cur.i2;
	      block[counter++] = hit_cur.j1;
	      block[counter++] = hit_cur.j2;
	      par.block_shading_counter->Add(hit_cur.name,counter);
	      if (!par.block_shading->Contains(hit_cur.name))
		par.block_shading->Add(hit_cur.name,block);
	      // printf("Add to block shading in realign   key: %s    data:",hit_cur.name);
	      // for (int i = 0; i < counter; i++)
	      // 	printf(" %i,",block[i]);
	      // printf("\n");
	    }
	  
	      //fprintf(stderr,"hit.name=%-15.15s  hit.index=%-5i hit.ftellpos=%-8i  hit.dbfile=%s\n",hit_cur.name,hit_cur.index,(unsigned int)hit_cur.ftellpos,hit_cur.dbfile);

	  if (nhits>=par.premerge) // realign the first premerge hits consecutively to query profile
	    {
	      if (hit_cur.irep==1) 
		{
		  // For each template (therefore irep==1), store template index and position on disk in a list
		  Realign_hitpos realign_hitpos;
		  realign_hitpos.ftellpos = hit_cur.ftellpos;    // stores position on disk of template for current hit
		  realign_hitpos.index = hit_cur.index;          // stores index of template of current hit
		  if (!phash_plist_realignhitpos->Contains(hit_cur.dbfile))
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
	  
	}
      nhits++;
    }
  if (Lmax>Lmaxmem)
    {
      Lmax=Lmaxmem;
      if (v>=1) 
	{
	  cerr<<"WARNING: Realigning sequences only up to length "<<Lmaxmem<<"."<<endl;
	  cerr<<"This is genarally unproboblematic but may lead to slightly sub-optimal alignments for longer sequences."<<endl;
 	  cerr<<"You can increase available memory using the -maxmem <GB> option (currently "<<par.maxmem<<" GB)."<<endl; // still to be implemented
	  cerr<<"The maximum length realignable is approximately (maxmem-0.5GB)/query_length/(cpus+1)/24B."<<endl;
	}
    }
  


  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs=1;   // needs to be set to 1 before threads are created
  for (bin=0; bin<bins; bin++)
    {
      // Free previously allocated memory
      if (hit[bin]->forward_allocated)
	hit[bin]->DeleteForwardMatrix(q->L+2);
      if (hit[bin]->backward_allocated)
	hit[bin]->DeleteBackwardMatrix(q->L+2);
      
      hit[bin]->AllocateForwardMatrix(q->L+2,Lmax+1);
      hit[bin]->AllocateBackwardMatrix(q->L+2,Lmax+1);

      bin_status[bin] = FREE;
    }
 
  if (print_elapsed) ElapsedTimeSinceLastCall("(prepare realign)");

  if (v>=2)
      printf("Realigning %i HMMs using HMM-HMM Maximum Accuracy algorithm\n",phash_plist_realignhitpos->Size());

  v1=v;
  if (v>0 && v<=3) v=1; else v-=2;  // Supress verbose output during iterative realignment and realignment

  // Align the first par.premerge templates?
  if (par.premerge>0)
    {
      if (v>=2) printf("Merging %i best hits to query alignment ...\n",par.premerge);
      
      bin=0;
      nhits=0;
      hitlist.Reset();
      while (!hitlist.End() && nhits<par.premerge)
	{
	  hit_cur = hitlist.ReadNext();
	  if (nhits>=imax(par.B,par.Z)) break;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
	  
	  nhits++;

	  if (hit_cur.L>Lmaxmem) continue;  // Don't align to long sequences due to memory limit
	  if (hit_cur.Eval > par.e) continue; // Don't align hits with an E-value below the inclusion threshold

	  // Open HMM database file dbfiles[idb]
	  FILE* dbf;
	  dbf = ffindex_fopen(dbhhm_data, dbhhm_index, hit_cur.dbfile);
	  if (dbf == NULL) {
	    char filename[NAMELEN];
	    RemoveExtension(filename, hit_cur.file);
	    strcat(filename,".a3m");
	    if(dba3m_index_file!=NULL) {
	      dbf = ffindex_fopen(dba3m_data, dba3m_index, filename);
	    } else {
	      cerr<<endl<<"Error opening "<<hit_cur.dbfile<<": A3M database missing\n"; exit(4);
	    }
	  }
	  if (dbf == NULL) OpenFileError(hit_cur.dbfile);

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

	  if (!strncmp(line,"HMMER3",5))      // read HMMER3 format
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
	  else if (line[0]=='#' || line[0]=='>')             // read a3m alignment
	    {
	      Alignment tali;
	      tali.Read(dbf,hit_cur.dbfile,line);
	      tali.Compress(hit_cur.dbfile);
	      //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
	      tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
	      t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
	      tali.FrequenciesAndTransitions(*(t[bin]));
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
	      cerr<<"Error in "<<par.argv[0]<<": wrong format while reading \'"<<hit_cur.dbfile<<". Reached end of file while reading HMM "<<hit_cur.name<<" \n";
	      exit(1);
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
	  PrepareTemplate(*q,*(t[bin]),format[bin]);
	  t[bin]->Log2LinTransitionProbs(1.0);

	  // Realign only around previous Viterbi hit
	  hit[bin]->i1 = hit_cur.i1;
	  hit[bin]->i2 = hit_cur.i2;
	  hit[bin]->j1 = hit_cur.j1;
	  hit[bin]->j2 = hit_cur.j2;
	  hit[bin]->nsteps = hit_cur.nsteps;
	  hit[bin]->i = hit_cur.i;
	  hit[bin]->j = hit_cur.j;
	  hit[bin]->realign_around_viterbi=true;

	  // Align q to template in *hit[bin]
	  hit[bin]->Forward(*q,*(t[bin]));
	  hit[bin]->Backward(*q,*(t[bin]));
	  hit[bin]->MACAlignment(*q,*(t[bin]));
	  hit[bin]->BacktraceMAC(*q,*(t[bin]));
	  
	  // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
	  hit[bin]->score      = hit_cur.score;
	  hit[bin]->score_ss   = hit_cur.score_ss;
	  hit[bin]->score_aass = hit_cur.score_aass;
	  hit[bin]->score_sort = hit_cur.score_sort;
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
	  
	  // merge only when hit length > MINCOLS_REALIGN (don't merge 1 column matches)
	  if (hit[bin]->matched_cols < MINCOLS_REALIGN) continue;

	  // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
	  char ta3mfile[NAMELEN];
	  //strcpy(ta3mfile,hit[bin]->file); // copy filename including path but without extension
	  RemoveExtension(ta3mfile,hit[bin]->dbfile);
	  strcat(ta3mfile,".a3m");
	  FILE* ta3mf;
	  ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
	  if (ta3mf == NULL) OpenFileError(ta3mfile);

	  Qali.MergeMasterSlave(*hit[bin],ta3mfile, ta3mf);
	  if (par.nodiff)
	    {
	      fclose(ta3mf);
	      ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
	      Qali_nodiff.MergeMasterSlave(*hit[bin],ta3mfile, ta3mf, false);
	    }

	  fclose(ta3mf);
	  
	  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
	  Qali.Compress("merged A3M file");
	  
	  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
	  Qali.N_filtered = Qali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);

	  // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
	  Qali.FrequenciesAndTransitions(*q);

	  stringstream ss_tmp;
	  ss_tmp << hit[bin]->name << "__" << hit[bin]->irep;
	  premerged_hits->Add((char*)ss_tmp.str().c_str());

	  if (par.notags) q->NeutralizeTags();
	  
	  if (!*par.clusterfile) { //compute context-specific pseudocounts?
	    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	    q->PreparePseudocounts();
	    // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	    q->AddAminoAcidPseudocounts();
	  } else {
	    // Add full context specific pseudocounts to query
	    q->AddContextSpecificPseudocounts();
	  }

	  q->CalculateAminoAcidBackground();
	  if (par.columnscore == 5 && !q->divided_by_local_bg_freqs) q->DivideBySqrtOfLocalBackgroundFreqs(par.half_window_size_local_aa_bg_freqs);

	  // Transform transition freqs to lin space if not already done
	  q->AddTransitionPseudocounts();
	  q->Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done

	}
    }
  
  if (print_elapsed) ElapsedTimeSinceLastCall("(premerge)");

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
      FILE* dbf;
      dbf = ffindex_fopen(dbhhm_data, dbhhm_index, dbfiles[idb]);
      if (dbf == NULL) {
	char filename[NAMELEN];
	RemoveExtension(filename, dbfiles[idb]);
	strcat(filename,".a3m");
	if(dba3m_index_file!=NULL) {
	  dbf = ffindex_fopen(dba3m_data, dba3m_index, filename);
	} else {
	  cerr<<endl<<"Error opening "<<dbfiles[idb]<<": A3M database missing\n"; exit(4);
	}
      }
      if (dbf == NULL) OpenFileError(dbfiles[idb]);

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
	      if (!fgetline(line,LINELEN,dbf)) {fprintf(stderr,"Error in %s: end of file %s reached prematurely!\n",par.argv[0],dbfiles[idb]); exit(1);}
	      while (strscn(line)==NULL && fgetline(line,LINELEN,dbf)) {} // skip lines that contain only white space

	      if (!strncmp(line,"HMMER3",5))      // read HMMER3 format
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
	      else if (line[0]=='#' || line[0]=='>')                 // read a3m alignment
		{
		  Alignment tali;
		  tali.Read(dbf,dbfiles[idb],line);
		  tali.Compress(dbfiles[idb]);
		  // qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
		  tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
		  t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
		  tali.FrequenciesAndTransitions(*(t[bin]));
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
  
  // Delete all hitlist entries with too short alignments
  nhits=0;
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      //printf("Deleting alignment of %s with length %i? irep=%i nhits=%-2i  par.B=%-3i  par.Z=%-3i par.e=%.2g par.b=%-3i  par.z=%-3i par.p=%.2g\n",hit_cur.name,hit_cur.matched_cols,hit_cur.irep,nhits,par.B,par.Z,par.e,par.b,par.z,par.p);

      if (nhits > par.realign_max && nhits>=imax(par.B,par.Z)) break;
      if (hit_cur.Eval > par.e)
	{
	  if (nhits>=imax(par.B,par.Z)) continue;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) continue;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
	}

      if (hit_cur.matched_cols < MINCOLS_REALIGN)
	{
	  if (v>=3) printf("Deleting alignment of %s with length %i\n",hit_cur.name,hit_cur.matched_cols);
	  //hitlist.ReadCurrent().Delete();
	  hitlist.Delete().Delete();               // delete the list record and hit object
	  // Make sure only realigned alignments get displayed!
	  if (last_round) {
	    if (par.B>par.Z) par.B--; else if (par.B==par.Z) {par.B--; par.Z--;} else par.Z--;
	    //if (par.b>par.z) par.b--; else if (par.b==par.z) {par.b--; par.z--;} else par.z--;
	  }
	}
      else nhits++;
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  int cluster_found = 0;
  int seqs_found = 0;
  char* argv_conf[MAXOPT];       // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf=0;               // Number of arguments in argv_conf
  FILE* fin;

#ifdef PTHREAD
  pthread_attr_init(&joinable);  // initialize attribute set with default values
  if (pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)!=0) // set attribute 'joinable'
    cerr<<"Error "<<pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)<<": could not set detach state for thread attibute.\n";
#endif

  par.premerge = 3;
  par.Ndiff = 1000;
  par.prefilter=true;
  par.early_stopping_filter=true;
  par.filter_thresh=0.01;
  par.filter_evals=new double[par.filter_length];
  strcpy(par.outfile,"");
  strcpy(db_ext,"hhm");
  N_searched=0;
  previous_hits = new Hash<Hit>(1631,hit_cur);
  premerged_hits = new Hash<char>(1631);
  
  // Make command line input globally available
  par.argv=argv;
  par.argc=argc;
  RemovePathAndExtension(program_name,argv[0]);
  Pathname(program_path, argv[0]);

  // Enable changing verbose mode before command line are processed
  for (int i=1; i<argc; i++)
    {
      if (argc>1 && !strcmp(argv[i],"-v0")) v=0;
      else if (argc>1 && !strcmp(argv[i],"-v1")) v=1;
      else if (argc>2 && !strcmp(argv[i],"-v")) v=atoi(argv[i+1]);
    }

  par.SetDefaultPaths(program_path);

  // Process default otpions from .hhdefaults file
  ReadDefaultsFile(argc_conf,argv_conf,program_path);
  ProcessArguments(argc_conf,argv_conf);
  
  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc,argv);

  // Check needed files
  if (!*par.infile || !strcmp(par.infile,"") || !strcmp(par.infile,"stdin")) // infile not given
    {help(); cerr<<endl<<"Error in "<<program_name<<": input file missing!\n"; exit(4);}
  if (!*db_base)
    {help(); cerr<<endl<<"Error in "<<program_name<<": database missing (see -d)\n"; exit(4);}
  if (par.addss==1 && (!*par.psipred || !*par.psipred_data))
    {help(); cerr<<endl<<"Error in "<<program_name<<": missing PSIPRED directory (see -psipred and -psipred_data).\nIf you don't need the predicted secondary structure, don't use the -addss option!\n"; exit(4);}
  if (!strcmp(par.clusterfile,""))
    {help(); cerr<<endl<<"Error in "<<program_name<<": context-specific library missing (see -contxt)\n"; exit(4);}
  if (!strcmp(par.cs_library,""))
    {help(); cerr<<endl<<"Error in "<<program_name<<": column state library (see -cslib)\n"; exit(4);}
  if (par.loc==0 && num_rounds>=2 && v>=1) cerr<<"WARNING: using -global alignment for iterative searches is deprecated since non-homologous sequence segments can easily enter the MSA and corrupt it.\n";
    if (num_rounds < 1) num_rounds=1; 
    else if (num_rounds > 8) 
      {
	if (v>=1) cerr<<"WARNING: Number of iterations ("<<num_rounds<<") to large => Set to 8 iterations\n";
	num_rounds=8; 
      }


  // Set databases
  strcpy(db,db_base);
  strcat(db,".cs219");

  strcpy(dbhhm,db_base);
  strcat(dbhhm,"_hhm_db");

  strcpy(dba3m,db_base);
  strcat(dba3m,"_a3m_db");

  fin = fopen(dba3m, "r");
  if (fin) {
    fclose(fin);
  } else {
    if(errno == EOVERFLOW)
    {
      cerr << endl;
      cerr <<"ERROR in "<< program_name <<": A3M database too big (>2GB on 32bit system?):";
      cerr << endl;
      cerr << dba3m;
      cerr << endl;
      exit(errno);
    }
    if (num_rounds > 1)
      {help(); cerr<<endl<<"Error in "<<program_name<<": A3M database missing (needed for more than 1 search iteration)\n"; exit(4);}
    if (*par.alnfile || *par.psifile || *par.hhmfile || *alis_basename)
      {help(); cerr<<endl<<"Error in "<<program_name<<": A3M database missing (needed for output alignment)\n"; exit(4);}
    dba3m[0] = 0;
  }

  q = new HMM;
  q_tmp = new HMM;

  par.block_shading = new Hash<int*>;
  par.block_shading_counter = new Hash<int>;
  par.block_shading->New(16381,NULL);
  par.block_shading_counter->New(16381,0);

  dbfiles_new = new char*[par.maxnumdb_no_prefilter+1];
  dbfiles_old = new char*[par.maxnumdb+1];

  // Prepare index-based databases
  char filename[NAMELEN];
  dbhhm_data_file = fopen(dbhhm, "r");
  if (!dbhhm_data_file) OpenFileError(dbhhm);
  strcpy(filename, dbhhm);
  strcat(filename, ".index");

  int filesize;
  filesize = CountLinesInFile(filename); 

  dbhhm_index_file = fopen(filename, "r");
  if (!dbhhm_index_file) OpenFileError(filename);

  dbhhm_index = ffindex_index_parse(dbhhm_index_file, filesize);
  if (dbhhm_index==NULL) {
    cerr<<"Error in "<<par.argv[0]<<": could not read index file"<<filename<<". Is the file empty or corrupted?\n";
    exit(1);
  }
  dbhhm_data = ffindex_mmap_data(dbhhm_data_file, &data_size);

  if (!*dba3m) {
    dba3m_data_file = dba3m_index_file = NULL;
    dba3m_index = NULL;
    // set premerge = 0 (no a3m database)
    par.premerge = 0;
  } else {
    dba3m_data_file = fopen(dba3m, "r");
    if (!dba3m_data_file) OpenFileError(dba3m);
    strcpy(filename, dba3m);
    strcat(filename, ".index");

    filesize = CountLinesInFile(filename);

    dba3m_index_file = fopen(filename, "r");
    if (!dba3m_index_file) OpenFileError(filename);

    dba3m_index = ffindex_index_parse(dba3m_index_file,filesize);
    if (dba3m_index==NULL) {
      cerr<<"Error in "<<par.argv[0]<<": could not read index file"<<filename<<". Is the file empty or corrupted?\n";
      exit(1);
    }
    dba3m_data = ffindex_mmap_data(dba3m_data_file, &data_size);
  }

  // Check for threads
  if (threads<=1) threads=0;
  else if (threads>MAXTHREADS)
    {
      threads=MAXTHREADS;
      if (v>=1) fprintf(stderr,"WARNING: number of CPUs set to maximum value of %i\n",MAXTHREADS);
    }

  // Set OpenMP threads
#ifdef _OPENMP
  cpu = imin(cpu,omp_get_max_threads());
  omp_set_num_threads(cpu);
#endif
  
  // Check option compatibilities
  if (par.nseqdis>MAXSEQDIS-3-par.showcons) par.nseqdis=MAXSEQDIS-3-par.showcons; //3 reserved for secondary structure
  if (par.aliwidth<20) par.aliwidth=20;
  if (par.pca<0.001) par.pca=0.001; // to avoid log(0)
  if (par.pre_pca<0.001) par.pre_pca=0.001; // to avoid log(0)
  if (par.b>par.B) par.B=par.b;
  if (par.z>par.Z) par.Z=par.z;
  if (par.maxmem<1.0) {cerr<<"WARNING: setting -maxmem to its minimum allowed value of 1.0\n"; par.maxmem=1.0;}

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix();

// Set secondary structure substitution matrix
  if (par.ssm) SetSecStrucSubstitutionMatrix();

  // Prepare context state pseudocounts lib
  if (*par.clusterfile) {
    fin = fopen(par.clusterfile, "r");
    if (!fin) OpenFileError(par.clusterfile);
    context_lib = new cs::ContextLibrary<cs::AA>(fin);
    fclose(fin);
    cs::TransformToLog(*context_lib);
    
    lib_pc = new cs::LibraryPseudocounts<cs::AA>(*context_lib, par.csw, par.csb);
  }

  // Prepare column state lib (context size =1 )
  fin = fopen(par.cs_library, "r");
  if (!fin) OpenFileError(par.cs_library);
  cs_lib = new cs::ContextLibrary<cs::AA>(fin);
  fclose(fin);
  cs::TransformToLin(*cs_lib);
  
  if (print_elapsed) ElapsedTimeSinceLastCall("(prepare CS pseudocounts)"); 

  // Read input file
  ReadInputFile();

  if (print_elapsed) ElapsedTimeSinceLastCall("(initialize)");
  
  if (par.prefilter)
    {
      // Initialize Prefiltering (Get DBsize)
      if (v>=2) printf("Reading in column state sequences for prefiltering\n");
      init_prefilter();
    }
  else // Set all HMMs in database as new_dbs
    {
      init_no_prefiltering();
    }

  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefilter)"); 

  // Input parameters
  if (v>=3)
    {
      cout<<"Input file       :   "<<par.infile<<"\n";
      cout<<"Output file      :   "<<par.outfile<<"\n";
      cout<<"Prefilter DB     :   "<<db<<"\n";
      cout<<"HHM DB           :   "<<dbhhm<<"\n";
    }

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags) q->NeutralizeTags();

  // Prepare multi-threading - reserve memory for threads, intialize, etc.
  if (threads==0) bins=1; else bins=iround(threads*1.2+0.5);
  for (bin=0; bin<bins; bin++)
    {
      t[bin]=new HMM;   // Each bin has a template HMM allocated that was read from the database file
      hit[bin]=new Hit; // Each bin has an object of type Hit allocated ...
      hit[bin]->AllocateBacktraceMatrix(q->L+2,par.maxres); // ...with a separate dynamic programming matrix (memory!!)
    }
  format = new(int[bins]);

  if (print_elapsed) ElapsedTimeSinceLastCall("(finished init)");

  if (par.nodiff)
    Qali_nodiff = Qali;

  //////////////////////////////////////////////////////////
  // Main loop overs search iterations
  //////////////////////////////////////////////////////////

  for (int round = 1; round <= num_rounds; round++) {

    if (v>=2) printf("\nIteration %i\n",round);

    // Settings for different rounds
    if (par.premerge > 0 && round > 1 && previous_hits->Size() >= par.premerge)
      {
	if (v>3) printf("Set premerge to 0! (premerge: %i   iteration: %i   hits.Size: %i)\n",par.premerge,round,previous_hits->Size());
	par.premerge = 0;
      }
    else 
      {
	par.premerge -= previous_hits->Size();
      }

    // Save HMM without pseudocounts for prefilter query-profile
    *q_tmp = *q;

    // Write query HHM file?
    if (*query_hhmfile) 
      {
	v1=v;
	if (v>0 && v<=3) v=1; else v-=2;
	
	// Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
	q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
	q->CalculateAminoAcidBackground();

	q->WriteToFile(query_hhmfile);

	v=v1;
      }
    
    // Add Pseudocounts, if no HMMER input
    if (input_format == 0)
      {
	// Transform transition freqs to lin space if not already done
	q->AddTransitionPseudocounts();
	
	if (!*par.clusterfile) { //compute context-specific pseudocounts?
	  // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	  q->PreparePseudocounts();
	  // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	  q->AddAminoAcidPseudocounts();
	} else {
	  // Add full context specific pseudocounts to query
	  q->AddContextSpecificPseudocounts();
	}
      } 
    else 
      {
	q->AddAminoAcidPseudocounts(0);
      }
    
    q->CalculateAminoAcidBackground();

    if (par.columnscore == 5 && !q->divided_by_local_bg_freqs) q->DivideBySqrtOfLocalBackgroundFreqs(par.half_window_size_local_aa_bg_freqs);
    
    if (print_elapsed) ElapsedTimeSinceLastCall("(before prefiltering (pseudocounts))");

    ///////////////////////////////////////////////////////////////////////////////
    // Prefiltering
    ///////////////////////////////////////////////////////////////////////////////

    if (par.prefilter)
      if (v>=2) printf("Prefiltering database\n");
      prefilter_db();  // in hhprefilter.C
    
    if (print_elapsed) ElapsedTimeSinceLastCall("(prefiltering)"); 

    if (ndb_new == 0) {
      printf("No HMMs pass prefilter => Stop searching!\n");
      break;
    }

    // Search datbases
    if (v>=2) {
      printf("HMMs passed 2nd prefilter (gapped profile-profile alignment)   : %6i\n", (ndb_new+ndb_old));
      printf("HMMs passed 2nd prefilter and not found in previous iterations : %6i\n", ndb_new);
      printf("Scoring %i HMMs using HMM-HMM Viterbi alignment\n", ndb_new);
    }

    // Main Viterbi HMM-HMM search
    // Starts with empty hitlist (hits of previous iterations were deleted) and creates a hitlist with the hits of this iteration
    search_database(dbfiles_new,ndb_new,(ndb_new + ndb_old));

    // check for new hits or end with iteration
    int new_hits = 0;
    hitlist.Reset();
    while (!hitlist.End())
      {
	hit_cur = hitlist.ReadNext();
	if (hit_cur.Eval > 100.0*par.e) break; // E-value much too large
	if (hit_cur.Eval > par.e) continue; // E-value too large
	new_hits++;
      }

    if (new_hits == 0 || round == num_rounds) 
      {
	last_round = true;
	if (round < num_rounds && v>=2)
	  printf("No new hits found in iteration %i => Stop searching\n",round);

	if (ndb_old > 0 && realign_old_hits)
	  {
	    printf("Rescoring previously found HMMs with Viterbi algorithm\n");
	    search_database(dbfiles_old,ndb_old,(ndb_new + ndb_old));
	    // Add dbfiles_old to dbfiles_new for realign
	    for (int a = 0; a < ndb_old; a++) 
	      {
		dbfiles_new[ndb_new]=new(char[strlen(dbfiles_old[a])+1]);
		strcpy(dbfiles_new[ndb_new],dbfiles_old[a]);
		ndb_new++;
	      }
	  }
	else if (!realign_old_hits && previous_hits->Size() > 0)
	  {
	    printf("Rescoring previously found HMMs with Viterbi algorithm\n");
	    perform_viterbi_search(ndb_new+previous_hits->Size());
	  }
      }

    // Realign hits with MAC algorithm
    if (par.realign) perform_realign(dbfiles_new,ndb_new);

    // Generate alignment for next iteration
    if (round < num_rounds || *par.alnfile || *par.psifile || *par.hhmfile || *alis_basename)
      {
	v1=v;
	if (v>0 && v<=3) v=1; else v-=2;
	
	// If new hits found, merge hits to query alignment
	if (new_hits != 0)
	  {
	    char ta3mfile[NAMELEN];
	    
	    if (v>=1) printf("Merging hits to query profile ...\n");
	    
	    // For each template below threshold
	    hitlist.Reset();
	    while (!hitlist.End())
	      {
		hit_cur = hitlist.ReadNext();
		if (hit_cur.Eval > 100.0*par.e) break; // E-value much too large
		if (hit_cur.Eval > par.e) continue; // E-value too large
		if (hit_cur.matched_cols < MINCOLS_REALIGN) continue; // leave out too short alignments
		stringstream ss_tmp;
		ss_tmp << hit_cur.name << "__" << hit_cur.irep;
		if (previous_hits->Contains((char*)ss_tmp.str().c_str())) continue;  // Already in alignment

		// Update counts
		cluster_found++;
		if (!strncmp(hit_cur.name,"cl|",3) || !strncmp(hit_cur.name,"UP20|",5) || !strncmp(hit_cur.name,"NR20|",5))   // kClust formatted database (NR20, ...)
		  {
		    char *ptr = hit_cur.name;
		    seqs_found += strint(ptr);
		  }
		else
		  seqs_found++;

		// Skip mering this hit if hit alignment was already merged during premerging
		if (premerged_hits->Contains((char*)ss_tmp.str().c_str())) continue;

		// Read a3m alignment of hit from <file>.a3m file and merge into Qali alignment
		RemoveExtension(ta3mfile,hit_cur.dbfile);
		strcat(ta3mfile,".a3m");
		FILE* ta3mf;
		ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
		if (ta3mf == NULL) OpenFileError(ta3mfile);
		Qali.MergeMasterSlave(hit_cur,ta3mfile, ta3mf);
		fclose(ta3mf);

		if (par.nodiff)
		  {
		    ta3mf = ffindex_fopen(dba3m_data, dba3m_index, ta3mfile);
		    Qali_nodiff.MergeMasterSlave(hit_cur,ta3mfile, ta3mf, false);
		    fclose(ta3mf);
		  }



		if (Qali.N_in>=MAXSEQ) break; // Maximum number of sequences reached 
	      }

	    // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
	    Qali.Compress("merged A3M file");

	    // Sort out the nseqdis most dissimilar sequences for display in the output alignments
	    Qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
	    
	    // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
	    float const COV_ABS = 25;     // min. number of aligned residues
	    int cov_tot = imax(imin((int)(COV_ABS / Qali.L * 100 + 0.5), 70), par.coverage);
	    if (v>2) printf("Filter new alignment with cov %3i%%\n", cov_tot);
	    Qali.N_filtered = Qali.Filter(par.max_seqid,cov_tot,par.qid,par.qsc,par.Ndiff);

	    if (print_elapsed) ElapsedTimeSinceLastCall("(merge hits to Qali)");
	  }
	  
	// Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
	Qali.FrequenciesAndTransitions(*q,NULL,true);

	if (par.notags) q->NeutralizeTags();

	// if needed, calculate SSpred
	if (par.addss && Qali.L>25 && (*alis_basename || round == num_rounds || new_hits == 0))
	  {
	    char ss_pred[par.maxres];
	    char ss_conf[par.maxres];
	    
	    CalculateSS(*q, ss_pred, ss_conf);
	    
	    Qali.AddSSPrediction(ss_pred, ss_conf);
	    
	    if (print_elapsed) ElapsedTimeSinceLastCall("(calculate SS_Pred)");
	  }
	
	if (print_elapsed) ElapsedTimeSinceLastCall("(Calculate AA frequencies and transitions)");

	if (*alis_basename)
	  {
	    stringstream ss_tmp;
	    ss_tmp << alis_basename << "_" << round << ".a3m";
	    if (par.nodiff)
	      Qali_nodiff.WriteToFile(ss_tmp.str().c_str(),"a3m");
	    else
	      Qali.WriteToFile(ss_tmp.str().c_str(),"a3m");
	  }

	v=v1;
	
      }
    else if (round == num_rounds) // Update counts for log
      {
	hitlist.Reset();
	while (!hitlist.End())
	  {
	    hit_cur = hitlist.ReadNext();
	    if (hit_cur.Eval > 100.0*par.e) break; // E-value much too large
	    if (hit_cur.Eval > par.e) continue; // E-value too large
	    stringstream ss_tmp;
	    ss_tmp << hit_cur.name << "__" << hit_cur.irep;
	    if (previous_hits->Contains((char*)ss_tmp.str().c_str())) continue;  // Already in alignment

	    // Update counts
	    cluster_found++;
	    if (!strncmp(hit_cur.name,"cl|",3) || !strncmp(hit_cur.name,"UP20|",5) || !strncmp(hit_cur.name,"NR20|",5))   // kClust formatted database (NR20, ...)
	      {
		char *ptr = hit_cur.name;
		seqs_found += strint(ptr);
	      }
	    else
	      seqs_found++;

	  }
      }
    
    if (v>=2)
      printf("%i sequences belonging to %i database HMMs found with an E-value < %-6.4g\n",seqs_found,cluster_found, par.e);

    if (v>=2 && (round < num_rounds || *par.alnfile || *par.psifile || *par.hhmfile || *alis_basename))
      printf("Number of effective sequences of resulting query HMM: Neff = %4.2f\n", q->Neff_HMM);

    if (q->Neff_HMM > neffmax && round < num_rounds) {
      printf("Diversity is above threshold (%4.2f). Stop searching! (Change threshold using -neffmax <float>.)\n", neffmax);
    }

    if (Qali.N_in>=MAXSEQ)
      printf("Maximun number of sequences in query alignment reached (%i). Stop searching!\n", MAXSEQ);

    if (new_hits == 0 || round == num_rounds || q->Neff_HMM > neffmax || Qali.N_in>=MAXSEQ) 
      break;

    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End())
      {
	hit_cur = hitlist.ReadNext();
	char strtmp[NAMELEN+6];
	sprintf(strtmp,"%s__%i",hit_cur.name,hit_cur.irep);
	if (!already_seen_filter || hit_cur.Eval > par.e || previous_hits->Contains(strtmp))
	  hit_cur.Delete(); // Delete hit object (deep delete with Hit::Delete())
	else
	  previous_hits->Add(strtmp,hit_cur);
	
	// Old version by Michael => Delete
	// hit_cur = hitlist.ReadNext();
	// stringstream ss_tmp;
	// ss_tmp << hit_cur.name << "__" << hit_cur.irep;
	// if (!already_seen_filter || hit_cur.Eval > par.e || previous_hits->Contains((char*)ss_tmp.str().c_str()))
	//   hit_cur.Delete(); // Delete hit object (deep delete with Hit::Delete())
	// else
	//   previous_hits->Add((char*)ss_tmp.str().c_str(), hit_cur);

	hitlist.Delete(); // Delete list record (flat delete)

      }

    if (print_elapsed) ElapsedTimeSinceLastCall("(end of this round)");

  } // end for-loop rounds
  

  //////////////////////////////////////////////////////////
  // Result section
  //////////////////////////////////////////////////////////

  // Warn, if HMMER files were used
  if (par.hmmer_used && v>=1)
    fprintf(stderr,"WARNING: Using HMMER files results in a drastically reduced sensitivity (>10%%).\nWe recommend to use HHMs build by hhmake.\n");

  // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
  if (*par.scorefile) {
    if (v>=3) printf("Printing scores file ...\n");
    hitlist.PrintScoreFile(*q);
  }

  // Print FASTA or A2M alignments?
  if (*par.pairwisealisfile) {
    if (v>=2) cout<<"Printing alignments in "<<(par.outformat==1? "FASTA" : par.outformat==2?"A2M" :"A3M")<<" format to "<<par.pairwisealisfile<<"\n";
    hitlist.PrintAlignments(*q,par.pairwisealisfile,par.outformat);
  }

  // Write alignments in tabular layout to alitabfile
  if (*par.alitabfile) 
    hitlist.WriteToAlifile(*q,alitab_scop);

  // Print summary listing of hits
  if (v>=3) printf("Printing hit list ...\n");
  hitlist.PrintHitList(*q_tmp,par.outfile);

  // Write only hit list to screen?
  if (v==2 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,109); // write only hit list to screen

  // Print alignments of query sequences against hit sequences
  hitlist.PrintAlignments(*q_tmp,par.outfile);

  // Write whole output file to screen? (max 10000 lines)
  if (v>=3 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,10009);

  // Generate output alignment or HMM file?
  if (*par.alnfile || *par.psifile || *par.hhmfile)
    {
      // Write output PSI-BLAST-formatted alignment?
      if (*par.psifile) 
	{
	  if (par.nodiff)
	    Qali_nodiff.WriteToFile(par.psifile,"psi");
	  else
	    Qali.WriteToFile(par.psifile,"psi");
	}

      // Write output HHM file?
      if (*par.hhmfile) 
	{
	  // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
	  q->AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
	  q->CalculateAminoAcidBackground();

	  q->WriteToFile(par.hhmfile);
	}

      // Write output A3M alignment?
      if (*par.alnfile) 
	{
	  if (par.nodiff)
	    Qali_nodiff.WriteToFile(par.alnfile,"a3m");
	  else
	    Qali.WriteToFile(par.alnfile,"a3m");
	}
      
    }

  ////////////////////////////////////////////////////
  // Clean up 
  ////////////////////////////////////////////////////

  fclose(dbhhm_data_file);
  fclose(dbhhm_index_file);
  if (dba3m_index_file!=NULL) {
    fclose(dba3m_data_file);
    fclose(dba3m_index_file);
  }
  free(dbhhm_index);
  free(dba3m_index);
  
  // Delete memory for dynamic programming matrix
  for (bin=0; bin<bins; bin++)
    {
      hit[bin]->DeleteBacktraceMatrix(q->L+2);
      if (hit[bin]->forward_allocated)
	hit[bin]->DeleteForwardMatrix(q->L+2);
      if (hit[bin]->backward_allocated)
	hit[bin]->DeleteBackwardMatrix(q->L+2);

      delete hit[bin];
      delete t[bin];
    }
  delete q;
  delete q_tmp;
  if (format) delete[](format);
  if (par.exclstr) delete[] par.exclstr;
  delete[] par.filter_evals;
  for (int n = 1; n < argc_conf; n++)
    delete[] argv_conf[n];
  if (par.dbfiles) delete[] par.dbfiles;
  for (int idb=0; idb<ndb_new; idb++) delete[](dbfiles_new[idb]);
  for (int idb=0; idb<ndb_old; idb++) delete[](dbfiles_old[idb]);
  delete[](dbfiles_new);
  delete[](dbfiles_old);
  par.block_shading->Reset();
  while (!par.block_shading->End())
    delete[] (par.block_shading->ReadNext()); 
  delete par.block_shading;
  delete par.block_shading_counter;
  
  if (par.prefilter)
    {
      free(X);
      free(length);
      free(first);
      for (int n = 0; n < num_dbs; n++)
	delete[](dbnames[n]);
      delete[](dbnames);
    }

  if (*par.clusterfile) {
    delete context_lib;
    delete lib_pc;
  }

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
