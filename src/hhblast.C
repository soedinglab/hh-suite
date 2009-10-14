// hhblast.C:
// Iterative search for a multiple alignment in a profile HMM database
// Compile:              g++ hhblast.C -o hhblast -O3 -lpthread -lrt -fno-strict-aliasing
// Compile for Valgrind: g++ hhblast.C -o hhblast_valgrind -lpthread -lrt -O -g
//
// Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: command line error  6: internal logic error  7: internal numeric error

////#define WINDOWS
#define PTHREAD
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <sstream>
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

#include <sys/time.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure
#include "hhdecl.C"      // Constants, global variables, struct Parameters
#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM

// includes needed for context specific pseudocounts
#include "amino_acid.cpp"
#include "sequence.cpp"
#include "profile.cpp"
#include "cluster.cpp"
#include "simple_cluster.cpp"
#include "matrix.cpp"
#include "cs_counts.cpp"

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
const char print_elapsed=0;
string command;
char tmp_file[]="/tmp/hhblastXXXXXX";
char dummydb [NAMELEN];
char* ptr;                // pointer for string manipulation

// HHblast variables

const char HHBLAST_VERSION[]="version 1.4.4 (September 2009)";
const char HHBLAST_REFERENCE[]="to be published.\n";
const char HHBLAST_COPYRIGHT[]="(C) Michael Remmert and Johannes Soeding\n";

int num_rounds   = 2;                  // number of iterations
float e_psi      = 100;                // E-value cutoff for prefiltering
const int N_PSI  = 10000;              // number of max. PSI-BLAST hits to take as HMM-DB
bool nodiff = false;                   // if true, do not filter in last round
bool filter = true;                    // Perform filtering of already seen HHMs
bool block_filter = true;              // Perform viterbi and forward algorithm only on block given by prefiltering
bool realign_old_hits = false;         // Realign old hits in last round or use previous alignments

int cpu = 1;

char config_file[NAMELEN];
char a3m_infile[NAMELEN];
char hhm_infile[NAMELEN];
char psi_infile[NAMELEN];
char alis_basename[NAMELEN];
char base_filename[NAMELEN];
char query_hhmfile[NAMELEN];

string tmp_psifile="";
string tmp_infile="";
// Read from config-file:
char db[NAMELEN];                    // BLAST formatted database with consensus sequences
char dbhhm[NAMELEN];                 // directory with database HMMs
char hh[NAMELEN];                    // directory with hhblast
char pre_mode[NAMELEN];              // prefilter mode (csblast or blast)
char blast[NAMELEN];                 // BLAST binaries (not needed with csBLAST)
char csblast[NAMELEN];               // csBLAST binaries (not needed with BLAST)
char csblast_db[NAMELEN];            // csBLAST database (not needed with BLAST)
char psipred[NAMELEN];               // PsiPred binaries
char psipred_data[NAMELEN];          // PsiPred data

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
class Posindex // stores position in file stream and unique index of template
{
public:
  long ftellpos;
  int index;
  int operator<(const Posindex& posindex) {return ftellpos<posindex.ftellpos;}
};

HMM q;                    // Create query  HMM with maximum of MAXRES match states
HMM qrev;                 // Create query  HMM with maximum of MAXRES match states
HMM* t[MAXBINS];          // Each bin has a template HMM allocated that was read from the database file
Hit* hit[MAXBINS];        // Each bin has an object of type Hit allocated with a separate dynamic programming matrix (memory!!)
HitList hitlist;          // list of hits with one Hit object for each pairwise comparison done
int* format;              // format[bin] = 0 if in HHsearch format => add pcs; format[bin] = 1 if in HMMER format => no pcs
int read_from_db;         // The value of this flag is returned from HMM::Read(); 0:end of file  1:ok  2:skip HMM
int N_searched;           // Number of HMMs searched
Alignment Qali;           // output A3M generated by merging A3M alignments for significant hits to the query alignment

const int MAXNUMDB=2*N_PSI;
Hash<char>* doubled;
int ndb_new=0;
char* dbfiles_new[MAXNUMDB+1];
int ndb_old=0;
char* dbfiles_old[MAXNUMDB+1];
Hash<Hit>* previous_hits;
 
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

///////////////////////////////////////////////////////////////////////////////////////
//// Do the pairwise comparison of q and *(t[bin]) for the database search
//////////////////////////////////////////////////////////////////////////////////////
void PerformViterbiByWorker(int bin)
{
  Hit hit_cur;

  // Prepare q ant t and compare
  PrepareTemplate(q,*(t[bin]),format[bin]);

  // Do HMM-HMM comparison
  for (hit[bin]->irep=1; hit[bin]->irep<=par.altali; hit[bin]->irep++)
    {
      // Break, if no previous_hit with irep is found
      hit[bin]->Viterbi(q,*(t[bin]));
      if (hit[bin]->irep>1 && hit[bin]->score <= SMIN) break;
      hit[bin]->Backtrace(q,*(t[bin]));
      
      hit[bin]->score_sort = hit[bin]->score_aass;
      //printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit[bin]->name,hit[bin]->fam,hit[bin]->irep,hit[bin]->score);

#ifdef PTHREAD
      pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif
     
      stringstream ss_tmp;
      ss_tmp << hit[bin]->name << "__" << hit[bin]->irep;

      if (previous_hits->Contains((char*)ss_tmp.str().c_str()))
	{
	  hit_cur = previous_hits->Remove((char*)ss_tmp.str().c_str());   // Remove hit from hash -> add to hitlist
	  //previous_hits->Add((char*)ss_tmp.str().c_str(), *(new Hit));
	  previous_hits->Add((char*)ss_tmp.str().c_str(), *(hit[bin]));
	  
	  // Overwrite *hit[bin] with alignment, etc. of hit_cur
	  hit_cur.score      = hit[bin]->score;
	  hit_cur.score_aass = hit[bin]->score_aass;
	  hit_cur.score_ss   = hit[bin]->score_ss;
	  hit_cur.Pval       = hit[bin]->Pval;
	  hit_cur.Pvalt      = hit[bin]->Pvalt;
	  hit_cur.logPval    = hit[bin]->logPval;
	  hit_cur.logPvalt   = hit[bin]->logPvalt;
	  hit_cur.logP1val   = hit[bin]->logP1val;
	  hit_cur.Eval       = hit[bin]->Eval;
	  hit_cur.logEval    = hit[bin]->logEval;
	  hit_cur.E1val      = hit[bin]->E1val;
	  hit_cur.Probab     = hit[bin]->Probab;

	  hitlist.Push(hit_cur);            // insert hit at beginning of list (last repeats first!)
	  
	}
      else
	hitlist.Push(*(hit[bin]));          // insert hit at beginning of list (last repeats first!)
	  

#ifdef PTHREAD
      pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif

      if (hit[bin]->score <= SMIN) break;  // break if score for first hit is already worse than SMIN
    }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Execute system command
/////////////////////////////////////////////////////////////////////////////////////
void runSystem(string cmd)
{
  if (v>2)
    cout << "Command: " << cmd << "!\n";
  int res = system(cmd.c_str());
  if (res!=0) 
    {
      cerr << endl << "ERROR when executing: " << cmd << "!\n";
      exit(1);
    }
    
}

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHblast %s\n",HHBLAST_VERSION);
  printf("Fast homology detection method HHblast to iteratively search a filtered NR HMM database.\n");
  printf("%s",HHBLAST_REFERENCE);
  printf("%s",HHBLAST_COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query [options]                                                             \n",program_name);
  printf("\n");
  printf(" -i <file>      input query sequence (FASTA format, not needed with -psi or -a3m)        \n");
  printf("\n");
  printf("Options:                                                                                 \n");
  printf(" -conf  <file>  config file for databases and bioprogs path (default=%s)                 \n",config_file); 
  printf(" -db    <file>  BLAST formatted database with consensus sequences (default=%s)           \n",db);
  printf(" -dbhhm <dir>   directory with database HMMs (default=%s)                                \n",dbhhm);
  printf(" -n     [1,8]   number of rounds (default=%i)                                            \n",num_rounds); 
  printf(" -e_hh  [0,1]   E-value cutoff for inclusion in result alignment (def=%G)                \n",par.e); 
  printf(" -e_psi [1,inf[ E-value cutoff for prefiltering (default=%-.2f)                          \n",e_psi); 
  printf("\n");
  printf("Input options:                                                                           \n");
  printf(" -a3m <file>    A3M alignment to restart HHblast ('A3M' format)                          \n");
  printf(" -psi <file>    PSI-BLAST alignment file to restart HHblast ('PSI' format)               \n");
  printf(" -hhm <file>    HMM File to restart HHblast (in combination with -psi or -a3m)           \n");
  printf("\n");
  printf("Output options:                                                                          \n");
  printf(" -o <file>      write results in standard format to file (default=<infile.hhr>)          \n");
  printf(" -oa3m <file>   write pairwise alignments in A3M format (default=none)                   \n");
  printf(" -opsi <file>   write pairwise alignments in PSI format (default=none)                   \n");
  printf(" -ohhm <file>   write HHM file of the pairwise alignments (default=none)                 \n");
  printf(" -oalis <base>  write pairwise alignments in A3M format after each round (default=none)  \n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose (def=%i)    \n",v);
  printf("\n");
  printf("HMM-HMM alignment options:                                                               \n");
  printf(" -realign       realign displayed hits with max. accuracy (MAC) algorithm                \n");
  printf(" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)           \n");
  printf(" -mact [0,1[    posterior probability threshold for MAC re-alignment (def=%.3f)          \n",par.mact);
  printf("                Parameter controls alignment greediness: 0:global >0.1:local             \n");
  printf(" -glob/-loc     use global/local alignment mode for searching/ranking (def=local)        \n");
  printf("\n");
  printf("Other options:                                                                           \n");
  printf(" -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=1)               \n");
#ifndef PTHREAD
  printf("(The -cpu option is inactive since POSIX threads ae not supported on your platform)      \n");
#endif
  printf("\n");
  printf("An extended list of options can be obtained by using '--help all' as parameter           \n");
  printf("\n\n");
  printf("Example: %s -i query.fas -oa3m query.a3m -n 2                                            \n",program_name);
  cout<<endl;
}

void help_all()
{
  printf("\n");
  printf("HHblast %s\n",HHBLAST_VERSION);
  printf("Fast homology detection method HHblast to iteratively search a filtered NR HMM database.\n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n\n");
  printf("Usage: %s -i query [options]                                                             \n",program_name);
  printf(" -i <file>      input query sequence (FASTA format, not needed with -psi or -a3m)        \n");
  printf("\n");
  printf("Options:                                                                                 \n");
  printf(" -conf <file>   config file for databases and bioprogs path (default=%s)                 \n",config_file); 
  printf(" -db <file>     BLAST formatted database with consensus sequences (default=%s)           \n",db);
  printf(" -dbhhm <dir>   directory with database HMMs (default=%s)                                \n",dbhhm);
  printf(" -n     [1,5]   number of rounds (default=%i)                                            \n",num_rounds); 
  printf(" -e_hh  [0,1]   E-value cutoff for inclusion in result alignment (default=%G)            \n",par.e); 
  printf(" -e_psi [1,inf[ E-value cutoff for prefiltering (default=%-.0f)                          \n",e_psi); 
  printf("\n");
  printf("Input options:                                                                           \n");
  printf(" -a3m <file>    A3M alignment to restart HHblast ('A3M' format)                          \n");
  printf(" -psi <file>    PSI-BLAST alignment file to restart HHblast ('PSI' format)               \n");
  printf(" -hhm <file>    HMM File to restart HHblast (in combination with -psi or -a3m)           \n");
  printf("\n");
  printf("Output options:                                                                          \n");
  printf(" -o <file>      write results in standard format to file (default=<infile.hhr>)          \n");
  printf(" -oa3m <file>   write pairwise alignments in A3M format (default=none)                   \n");
  printf(" -opsi <file>   write pairwise alignments in PSI format (default=none)                   \n");
  printf(" -ohhm <file>   write HHM file of the pairwise alignments (default=none)                 \n");
  printf(" -oalis <base>  write pairwise alignments in A3M format after each round (default=none)  \n");
  printf(" -qhhm <file>   write query input HHM file of last round (default=none)                  \n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose (def=%i)    \n",v);
  printf(" -seq <int>     max. number of query/template sequences displayed (def=%i)               \n",par.nseqdis);
  printf(" -nopred        don't add predicted 2ndary structure in output alignments                \n");
  printf(" -aliw <int>    number of columns per line in alignment list (def=%i)                    \n",par.aliwidth);
  printf(" -p <float>     minimum probability in summary and alignment list (def=%G)               \n",par.p);
  printf(" -E <float>     maximum E-value in summary and alignment list (def=%G)                   \n",par.E);
  printf(" -Z <int>       maximum number of lines in summary hit list (def=%i)                     \n",par.Z);
  printf(" -z <int>       minimum number of lines in summary hit list (def=%i)                     \n",par.z);
  printf(" -B <int>       maximum number of alignments in alignment list (def=%i)                  \n",par.B);
  printf(" -b <int>       minimum number of alignments in alignment list (def=%i)                  \n",par.b);
  printf("\n");
  printf("Directories for needed programs                                                          \n");
  printf(" -pre_mode    <mode>  prefilter mode (blast or csblast) (default=%s)                     \n",pre_mode);
  printf(" -hh           <dir>  directory with HHBLAST executables and reformat.pl (default=%s)    \n",hh);
  printf(" -blast        <dir>  directory with BLAST executables (default=%s)                      \n",blast);
  printf(" -csblast      <dir>  directory with csBLAST executables (default=%s)                    \n",csblast);
  printf(" -csblast_db   <dir>  directory with csBLAST database (default=%s)                       \n",csblast_db);
  printf(" -psipred      <dir>  directory with PsiPred executables (default=%s)                    \n",psipred);
  printf(" -psipred_data <dir>  directory with PsiPred data (default=%s)                           \n",psipred_data);
  printf("\n");
  printf("Filter options                                                                           \n");
  printf(" -nofilter      disable all filter steps (except for PSI-BLAST prefiltering)             \n");
  printf(" -nodbfilter    disable additional filtering of prefiltered HMMs                         \n");
  printf(" -noblockfilter search complete matrix in Viterbi                                        \n");
  printf("\n");
  printf("Filter result alignment (options can be combined):                                       \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)                         \n",par.max_seqid);
  printf(" -diff [0,inf[  filter most diverse set of sequences, keeping at least this              \n");
  printf("                many sequences in each block of >50 columns (def=%i)                     \n",par.Ndiff);
  printf(" -nodiff        do not filter sequences in last iteration (def=off)                      \n");
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i)                                \n",par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i)                       \n",par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)                          \n",par.qsc);
  printf("\n");
  printf("HMM-HMM alignment options:                                                               \n");
  printf(" -realign       realign displayed hits with max. accuracy (MAC) algorithm                \n");
  printf(" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)           \n");
  printf(" -mact [0,1[    posterior probability threshold for MAC re-alignment (def=%.3f)          \n",par.mact);
  printf("                Parameter controls alignment greediness: 0:global >0.1:local             \n");
  printf(" -glob/-loc     use global/local alignment mode for searching/ranking (def=local)        \n");
  printf(" -alt <int>     show up to this many significant alternative alignments(def=%i)          \n",par.altali);
  printf(" -jdummy [0,20] ... (default=%i)                                                         \n",par.jdummy);       
  printf("\n");
  printf("Pseudocount options:                                                                     \n");
  printf(" -pcm  0-2      Pseudocount mode (default=%-i)                                           \n",par.pcm);
  printf("                tau = substitution matrix pseudocount admixture                          \n");
  printf("                0: no pseudo counts:     tau = 0                                         \n");
  printf("                1: constant              tau = a                                         \n");
  printf("                2: divergence-dependent: tau = a/(1 + ((Neff-1)/b)^c)                    \n");
  printf("                   Neff=( (Neff_q^d+Neff_t^d)/2 )^(1/d)                                  \n");
  printf("                   Neff_q = av number of different AAs per column in query               \n");
  printf("                3: constant divergence pseudocounts                                      \n");
  printf("                4: divergence-dependent, with composition-adjusted matrix                \n");
  printf(" -pca  [0,1]    overall pseudocount admixture (def=%-.1f)                                \n",par.pca);
  printf(" -pcb  [1,inf[  threshold for Neff) (def=%-.1f)                                          \n",par.pcb);
  printf(" -pcc  [0,3]    extinction exponent for tau(Neff)  (def=%-.1f)                           \n",par.pcc);
  printf(" -pcw  [0,3]    weight of pos-specificity for pcs  (def=%-.1f)                           \n",par.pcw);
  printf("\n");
  printf("Gap cost options:                                                                        \n");
  printf(" -gapb [0,inf[  transition pseudocount admixture (def=%-.2f)                             \n",par.gapb);
  printf(" -gapd [0,inf[  Transition pseudocount admixture for opening gap (default=%-.2f)         \n",par.gapd);
  printf(" -gape [0,1.5]  Transition pseudocount admixture for extending gap (def=%-.2f)           \n",par.gape);
  printf(" -gapf ]0,inf]  factor for inc./reducing the gap open penalty for deletes (def=%-.2f)    \n",par.gapf);
  printf(" -gapg ]0,inf]  factor for inc./reducing the gap open penalty for deletes (def=%-.2f)    \n",par.gapg);
  printf(" -gaph ]0,inf]  factor for inc./reducing the gap extension penalty for deletes(def=%-.2f)\n",par.gaph);
  printf(" -gapi ]0,inf]  factor for inc./reducing the gap extension penalty for inserts(def=%-.2f)\n",par.gapi);
  printf(" -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f)        \n",par.egq);
  printf(" -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)     \n",par.egt);
  printf("\n");
  printf("Other options:                                                                           \n");
  printf(" -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=1)               \n");
#ifndef PTHREAD
  printf("(The -cpu option is inactive since POSIX threads ae not supported on your platform)      \n");
#endif
  printf("\n");
  printf("An extended list of options can be obtained by using '--help all' as parameter           \n");
  printf("\n\n");
  printf("Example: %s -i query.fas -oa3m query.a3m -n 2                                            \n",program_name);
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
      else if (!strcmp(argv[i],"-a3m"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -a3m\n"; exit(4);}
          else strcpy(a3m_infile,argv[i]);
	  strcpy(par.infile,"");
        }
      else if (!strcmp(argv[i],"-hhm"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -hhm\n"; exit(4);}
          else strcpy(hhm_infile,argv[i]);
        }
      else if (!strcmp(argv[i],"-psi"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -psi\n"; exit(4);}
          else strcpy(psi_infile,argv[i]);
	  strcpy(par.infile,"");
        }
      else if (!strcmp(argv[i],"-db"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database file following -db\n"; exit(4);}
          else
              strcpy(db,argv[i]);
        }
      else if (!strcmp(argv[i],"-dbhhm"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database directory following -dbhhm\n"; exit(4);}
          else
              strcpy(dbhhm,argv[i]);
        }
      else if (!strcmp(argv[i],"-hh"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no directory following -hh\n"; exit(4);}
          else
              strcpy(hh,argv[i]);
        }
      else if (!strcmp(argv[i],"-pre_mode"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no mode following -pre_mode\n"; exit(4);}
          else
              strcpy(pre_mode,argv[i]);
        }
      else if (!strcmp(argv[i],"-blast"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no directory following -blast\n"; exit(4);}
          else
              strcpy(blast,argv[i]);
        }
      else if (!strcmp(argv[i],"-csblast"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no directory following -csblast\n"; exit(4);}
          else
              strcpy(csblast,argv[i]);
        }
      else if (!strcmp(argv[i],"-csblast_db"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no directory following -csblast_db\n"; exit(4);}
          else
              strcpy(csblast_db,argv[i]);
        }
      else if (!strcmp(argv[i],"-psipred"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no directory following -psipred\n"; exit(4);}
          else
              strcpy(psipred,argv[i]);
        }
      else if (!strcmp(argv[i],"-psipred_data"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database directory following -psipred_data\n"; exit(4);}
          else
              strcpy(psipred_data,argv[i]);
        }
      else if (!strcmp(argv[i],"-conf"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no config file following -conf\n"; exit(4);}
          else strcpy(config_file,argv[i]);
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
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help"))
        {
          if (++i>=argc || argv[i][0]=='-') {help(); exit(0);}
          if (!strcmp(argv[i],"all")) {help_all(); exit(0);}
          else {help(); exit(0);}
        }
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-n") && (i<argc-1)) 
	{ 
	  num_rounds = atoi(argv[++i]); 
	  if (num_rounds < 1) 
	    num_rounds=1; 
	  if (num_rounds > 8) 
	    {
	      cerr<<endl<<"WARNING! Number of rounds ("<<num_rounds<<") to large => Set to 8 rounds\n";
	      num_rounds=8; 
	    }
	}
      else if (!strcmp(argv[i],"-p") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-P") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-E") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-b") && (i<argc-1)) par.b = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-B") && (i<argc-1)) par.B = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-z") && (i<argc-1)) par.z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-Z") && (i<argc-1)) par.Z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-e_hh") && (i<argc-1)) par.e = atof(argv[++i]);
      else if (!strcmp(argv[i],"-e_psi") && (i<argc-1)) e_psi = atof(argv[++i]);
      else if (!strncmp(argv[i],"-nopred",7)) par.showpred=0;
      else if (!strcmp(argv[i],"-seq") && (i<argc-1))  par.nseqdis=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-aliw") && (i<argc-1)) par.aliwidth=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-nodiff")) nodiff=true;
      else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pcm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pca=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pcb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pcc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcw") && (i<argc-1)) par.pcw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;}
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egq") && (i<argc-1)) par.egq=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egt") && (i<argc-1)) par.egt=atof(argv[++i]);
      else if (!strcmp(argv[i],"-filterlen") && (i<argc-1)) 
	{
	  par.filter_length=atoi(argv[++i]);
	  delete par.filter_evals;
	  par.filter_evals=new double[par.filter_length];
	}
      else if (!strcmp(argv[i],"-filtercut") && (i<argc-1)) par.filter_thresh=(double)atof(argv[++i]);
      else if (!strcmp(argv[i],"-nofilter")) {filter=false; block_filter=false; par.filter_thresh=0;}
      else if (!strcmp(argv[i],"-nodbfilter")) {par.filter_thresh=0;}
      else if (!strcmp(argv[i],"-noblockfilter")) {block_filter=false;}
      else if (!strcmp(argv[i],"-block_len") && (i<argc-1)) par.block_shading_space = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-shading_mode") && (i<argc-1)) strcpy(par.block_shading_mode,argv[++i]);
      else if (!strcmp(argv[i],"-realignoldhits")) realign_old_hits=true;
      else if (!strcmp(argv[i],"-realign")) par.realign=1;
      else if (!strcmp(argv[i],"-norealign")) par.realign=0;
      else if (!strncmp(argv[i],"-glo",3)) {par.loc=0; if (par.mact>0.3 && par.mact<0.301) {par.mact=0;} }
      else if (!strncmp(argv[i],"-loc",4)) par.loc=1;
      else if (!strncmp(argv[i],"-alt",4) && (i<argc-1)) par.altali=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-mact") && (i<argc-1)) par.mact=atof(argv[++i]);
      else if (!strcmp(argv[i],"-mapt") && (i<argc-1)) par.mact=atof(argv[++i]);
      else if (!strncmp(argv[i],"-cpu",4) && (i<argc-1)) { threads=atoi(argv[++i]); cpu = threads; }
      else if (!strncmp(argv[i],"-jdummy",7) && (i<argc-1)) par.jdummy=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-csb") && (i<argc-1)) par.csb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-csw") && (i<argc-1)) par.csw=atof(argv[++i]);
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input

  if (*csblast_db) {
    //cerr<<"CSBLAST DB given!\n";
    strcpy(par.clusterfile,hh);
    strcat(par.clusterfile,"/nr30_neff2.5_1psi_N1000000_W13_K4000_wcenter1.6_wmax1.36_beta0.85.prf");
    //strcpy(par.clusterfile,"/cluster/bioprogs/hhblast/clusters.prf");
    // TODO: change back, when Andreas has rewritten hhhmm.C
    //strcpy(par.clusterfile,csblast_db);
  } else {
    cerr<<"CSBLAST DB not given!\n";
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Read config file
/////////////////////////////////////////////////////////////////////////////////////
void ReadConfigFile(char filename[])
{
  char* c;         //pointer to scan line read in for end of argument
  FILE* configf=NULL;

  // Open config file
  configf = fopen(filename,"r");
  if (!configf) return;

  // Read in options until end-of-file
  while (fgets(line,LINELEN,configf))
    {
      if (!strcmp(line,"\n") || line[0] == '#')
	continue;

      // Analyze line
      c=line;
      
      // find parameter name
      char param[NAMELEN];
      c = strwrd(param, c);
      // ignore '='
      c++; 
      // find parameter value
      char value[NAMELEN];
      c = strwrd(value, c);

      if (!*param || !*value) {
	 cerr<<endl<<"WARNING: Ignoring line '"<<line<<"' in config-file (WRONG FORMAT)!\n";
	 continue;
      }

      if (!strcmp(param, "db")) strcpy(db,value);
      else if (!strcmp(param, "dbhhm")) strcpy(dbhhm,value);
      else if (!strcmp(param, "hh")) strcpy(hh,value);
      else if (!strcmp(param, "pre_mode")) strcpy(pre_mode,value);
      else if (!strcmp(param, "blast")) strcpy(blast,value);
      else if (!strcmp(param, "csblast")) strcpy(csblast,value);
      else if (!strcmp(param, "csblast_db")) strcpy(csblast_db,value);
      else if (!strcmp(param, "psipred")) strcpy(psipred,value);
      else if (!strcmp(param, "psipred_data")) strcpy(psipred_data,value);
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<param<<" in config file...\n";

    }

  fclose(configf);
}

/////////////////////////////////////////////////////////////////////////////////////
// Check given input files
/////////////////////////////////////////////////////////////////////////////////////
void CheckInputFiles()
{
  if (!*par.infile || !strcmp(par.infile,"") || !strcmp(par.infile,"stdin")) // infile not given
    {
      if (!*a3m_infile && !*psi_infile)
	{help(); cerr<<endl<<"Error in "<<program_name<<": input file missing!\n"; exit(4);}
      
      if (*a3m_infile && strcmp(a3m_infile,"")) 
	{
	  strcpy(par.infile,a3m_infile);
	  if (*psi_infile && strcmp(psi_infile,"")) 
	    {
	      command = (string)hh + "/reformat.pl psi psi " + (string)psi_infile + " " + tmp_psifile + " -r -M first > /dev/null";
	    }
	  else
	    {
	      command = (string)hh + "/reformat.pl a3m psi " + (string)a3m_infile + " " + tmp_psifile + " -r > /dev/null";
	    }
	  runSystem(command);
	  strcpy(par.infile,a3m_infile);
	  RemoveExtension(base_filename,a3m_infile);
	}
      else
	{
	  if (*psi_infile && strcmp(psi_infile,"")) 
	    {
	      command = (string)hh + "/reformat.pl psi psi " + (string)psi_infile + " " + tmp_psifile + " -r -M first > /dev/null";
	      runSystem(command);
	      command = (string)hh + "/reformat.pl psi a3m " + (string)psi_infile + " " + (string)tmp_file + ".a3m -M first > /dev/null";
	      runSystem(command);
	      strcpy(par.infile,tmp_file);
	      strcat(par.infile,".a3m");
	      RemoveExtension(base_filename,psi_infile);
	    }
	  else
	    {
	      help(); 
	      cerr<<endl<<"Error in "<<program_name<<": input file missing!\n"; 
	      exit(4);
	    }
	}

      // Extract query sequence
      FILE* inf = NULL;
      FILE* outf = NULL;
      inf = fopen(par.infile,"rb");
      outf = fopen((tmp_infile).c_str(),"w");
      if (!inf) OpenFileError(par.infile);
      if (!outf) OpenFileError((tmp_infile).c_str());
      while(fgetline(line,LINELEN,inf))
	{
	  if (line[0] == '>' && strncmp(line,">ss",3) && strncmp(line,">sa",3))
	    break;
	}
      fprintf(outf,"%s\n",line);
      while(fgetline(line,LINELEN,inf))
	{
	  if (line[0] == '>')
	    break;
	  fprintf(outf,"%s\n",line);
	}
      fclose(inf);
      fclose(outf);
      
    }
  else 
    {
      // TODO! Check for single sequence in FASTA!!!

      // Copy infile to tmp_file.fas as input for the BLAST prefilter searches
      command = (string)hh + "/reformat.pl fas fas " + (string)par.infile + " " + tmp_infile + " > /dev/null";
      runSystem(command);
      RemoveExtension(base_filename,par.infile);
      strcpy(par.infile,tmp_infile.c_str());
      // Create simple PSI-file
      command = (string)hh + "/reformat.pl fas psi " + (string)par.infile + " " + tmp_psifile + " -r -M first > /dev/null";
      runSystem(command);
    }

  int v1=v;
  if (v<=3) v=1; else v-=2;

  // Read Query alignment
  FILE* qa3mf=fopen(par.infile,"r");
  if (!qa3mf) OpenFileError(par.infile);
  Qali.Read(qa3mf,par.infile);
  Qali.Compress("compress Qali");
  //Qali.N_filtered = Qali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
  fclose(qa3mf);

  v=v1;

  // Check if hhm file is given
  if (*hhm_infile && strcmp(hhm_infile,"")) 
    strcpy(par.infile,hhm_infile);

  if (!*par.outfile)      // outfile not given? Name it basename.hhm
    {
      strcpy(par.outfile,base_filename);
      strcat(par.outfile,".hhr");
      if (v>=2) cout<<"Search results will be written to "<<par.outfile<<"\n";
    }

}

// Calculate secondary structure prediction with PSIpred
void CalculateSS(char *ss_pred, char *ss_conf)
{
  // Initialize
  strcpy(ss_pred," ");
  strcpy(ss_conf," ");
  char rootname[NAMELEN];
  RemovePath(rootname,tmp_file);

  // Create dummy-DB if not exists
  if (!*dummydb)
    {
      strcpy(dummydb,tmp_file);
      strcat(dummydb,"_dummy_db");
      command = "cp " + tmp_infile + " " + (string)dummydb;
      runSystem(command);
      command = (string)blast + "/formatdb -i " + (string)dummydb + " -l /dev/null > /dev/null";
      runSystem(command);
    }

  // Create BLAST checkpoint file
  command = (string)blast + "/blastpgp -b 1 -j 1 -h 0.001 -d " + (string)dummydb + " -i " + tmp_infile + " -B " + tmp_psifile + " -C " + (string)tmp_file + ".chk 1> /dev/null 2> /dev/null";
  runSystem(command);
  command = "echo " + (string)rootname + ".chk > " + (string)tmp_file + ".pn";
  runSystem(command);
  command = "echo " + (string)rootname + ".fas > " + (string)tmp_file + ".sn";
  runSystem(command);
  command =  (string)blast + "/makemat -P " + (string)tmp_file;
  runSystem(command);

  // Run PSIpred
  command = (string)psipred + "/psipred " + (string)tmp_file + ".mtx " + (string)psipred_data + "/weights.dat " + (string)psipred_data + "/weights.dat2 " + (string)psipred_data + "/weights.dat3 " + (string)psipred_data + "/weights.dat4 > " + (string)tmp_file + ".ss";
  runSystem(command);
  command = (string)psipred + "/psipass2 " + (string)psipred_data + "/weights_p2.dat 1 0.98 1.09 " + (string)tmp_file + ".ss2 " + (string)tmp_file + ".ss > " + (string)tmp_file + ".horiz";
  runSystem(command);

  // Read results
  char filename[NAMELEN];
  strcpy(filename,tmp_file);
  strcat(filename,".horiz");
  FILE* horizf = fopen(filename,"r");
  if (!horizf) return;

  while (fgets(line,LINELEN,horizf))
    {
      char tmp_seq[NAMELEN]="";
      char* ptr=line;
      if (!strncmp(line,"Conf:",5))
	{
	  ptr+=5;
	  strwrd(tmp_seq,ptr);
	  strcat(ss_conf,tmp_seq);
	}
      if (!strncmp(line,"Pred:",5))
	{
	  ptr+=5;
	  strwrd(tmp_seq,ptr);
	  strcat(ss_pred,tmp_seq);
	}
    }
  fclose(horizf);

  // NEEDED???????
  //$ss_conf=~tr/0-9/0/c; # replace all non-numerical symbols with a 0

  if (v>3)
    {
      printf("SS-pred: %s\n",ss_pred);
      printf("SS-conf: %s\n",ss_conf);
    }

}

void search_loop(char *dbfiles[], int ndb, bool alignByWorker=true)
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Search databases
  double filter_cutoff = par.filter_length*par.filter_thresh;
  par.filter_sum=par.filter_length;
  par.filter_counter=0;
  for (int a=0; a<par.filter_length; a++) {par.filter_evals[a]=1;}

  // For all the databases given in -d '...' option ...
  for (int idb=0; idb<ndb; idb++)
    {

      // Check early stopping filter
      if (par.filter_sum < filter_cutoff)
	{
	  if (v>4)
	    printf("Stop after DB-HHM %i from %i\n",idb,ndb);
	  break;
	}
      
      // Open HMM database
      //cerr<<"\nReading db file "<<idb<<" dbfiles[idb]="<<dbfiles[idb]<<"\n";
      FILE* dbf=fopen(dbfiles[idb],"rb");
      if (!dbf) OpenFileError(dbfiles[idb]);

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
	  char path[NAMELEN];
	  Pathname(path,dbfiles[idb]);

	  ///////////////////////////////////////////////////
	  // Read next HMM from database file
	  if (!fgetline(line,LINELEN,dbf)) {continue;}
	  if (!strncmp(line,"HMMER",5))      // read HMMER format
	    {
	      format[bin] = 1;
	      t[bin]->ReadHMMer(dbf,dbfiles[idb]);
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
	  else if (line[0]=='#')             // read a3m alignment
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
	      cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<dbfiles[idb]<<"\'\n";
	      cerr<<"line = "<<line<<"\n";
	      exit(1);
	    }
	  if (v>=4) printf("Aligning with %s\n",t[bin]->name);  /////////////////////v>=4
	  ///////////////////////////////////////////////////
	  
	  hit[bin]->dbfile = new(char[strlen(dbfiles[idb])+1]);
	  strcpy(hit[bin]->dbfile,dbfiles[idb]); // record db file name from which next HMM is read
	  
	  ++N_searched;
	  if (v>=1 && !((idb+1)%20))
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
	      // If no submitted jobs are in the queue we have to wait for a new job ...
	      struct timespec ts;
	      clock_gettime(CLOCK_REALTIME,&ts);
	      ts.tv_sec += 1;
	      rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex,&ts);
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
  // End search databases
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

/////////////////////////////////////////////////////////////////////////////////////
// Perform Viterbi search on each hit object in previous_hits, but keep old alignment
/////////////////////////////////////////////////////////////////////////////////////
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
  int v1=v;
  if (v<=3) v=1; else v-=2;

  char *dbfiles[MAXNUMDB+1];
  int ndb = 0;

  par.block_shading->Reset();
  while (!par.block_shading->End())
      delete[] (par.block_shading->ReadNext()); 
  par.block_shading->New(16381,NULL);

  // Get dbfiles of previous hits
  previous_hits->Reset();
  while (!previous_hits->End())
    {
      Hit hit_cur = previous_hits->ReadNext();
      if (hit_cur.irep==1)  
	{
	  dbfiles[ndb]=new(char[strlen(hit_cur.dbfile)+1]);
	  strcpy(dbfiles[ndb],hit_cur.dbfile);
	  ++ndb;
	}
      // Seach only around viterbi hit
      if (block_filter)
	{
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

  search_loop(dbfiles,ndb,false);

  if (v1>=1) cout<<"\n";
  v=v1;

  if (print_elapsed) ElapsedTimeSinceLastCall("(search through database)");

  // Sort list according to sortscore
  if (v>=3) printf("Sorting hit list ...\n");
  hitlist.SortList();

  hitlist.CalculatePvalues(q);  // Use NN prediction of lamda and mu
  
  hitlist.CalculateHHblastEvalues(q);
  
}

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
  int v1=v;
  if (v<=3) v=1; else v-=2;
  if (print_elapsed) ElapsedTimeSinceLastCall("(preparing for search)");

  hitlist.N_searched=db_size; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  search_loop(dbfiles,ndb);

  if (v1>=1) cout<<"\n";
  v=v1;

  if (print_elapsed) ElapsedTimeSinceLastCall("(search through database)");

  // Sort list according to sortscore
  if (v>=3) printf("Sorting hit list ...\n");
  hitlist.SortList();

  hitlist.CalculatePvalues(q);  // Use NN prediction of lamda and mu

  hitlist.CalculateHHblastEvalues(q);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Realign hits with MAC algorithm

void perform_realign(char *dbfiles[], int ndb)
{
  q.Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done
      
  Hash< List<Posindex>* >* realign; // realign->Show(dbfile) is list with ftell positions for templates in dbfile to be realigned
  realign = new(Hash< List<Posindex>* >);
  realign->New(3601,NULL);
  par.block_shading->Reset();
  while (!par.block_shading->End())
      delete[] (par.block_shading->ReadNext()); 
  par.block_shading->New(16381,NULL);
  par.block_shading_counter->New(16381,NULL);
  Hit hit_cur;
  const float MEMSPACE_DYNPROG = 2.0*1024.0*1024.0*1024.0;
  int nhits=0;
  int Lmax=0;      // length of longest HMM to be realigned
  int Lmaxmem=(int)((float)MEMSPACE_DYNPROG/q.L/6.0/8.0/bins); // longest allowable length of database HMM
  int N_aligned=0;
  
  // Store all dbfiles and ftell positions of templates to be displayed and realigned
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      if (nhits>=imax(par.B,par.Z)) break;
      if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
      if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
      if (hit_cur.L>Lmax) Lmax=hit_cur.L;
      if (hit_cur.L<=Lmaxmem)
	{
	  // Seach only around viterbi hit
	  if (block_filter)
	    {
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
	  if (nhits>=par.jdummy || hit_cur.irep>1 || hit_cur.Eval > par.e) // realign the first jdummy hits consecutively to query profile
	    {
	      Posindex posindex;
	      posindex.ftellpos = hit_cur.ftellpos;
	      posindex.index = hit_cur.index;
	      if (realign->Contains(hit_cur.dbfile))
		  realign->Show(hit_cur.dbfile)->Push(posindex);
	      else
		{
		  List<Posindex>* newlist = new(List<Posindex>);
		  newlist->Push(posindex);
		  realign->Add(hit_cur.dbfile,newlist);
		}
	    }
	}
      nhits++;
    }
  if (Lmax>Lmaxmem)
    {
      Lmax=Lmaxmem;
      if (v>=1) 
	{
	  cerr<<"WARNING: Realigning sequences only up to length "<<Lmaxmem<<" due to limited memory."<<endl;
	  if (bins>1) cerr<<"Note: you can reduce memory requirement by lowering N in the -cpu N option."<<endl;
	}
    }
  
  // Initialize and allocate space for dynamic programming
  jobs_running = 0;
  jobs_submitted = 0;
  reading_dbs=1;   // needs to be set to 1 before threads are created
  for (bin=0; bin<bins; bin++)
    bin_status[bin] = FREE;
 
  if (print_elapsed) ElapsedTimeSinceLastCall("(prepare realign)");
 
  if (v>=1) printf("Realigning %i query-template alignments with maximum accuracy (MAC) algorithm ...\n",nhits);

  int v1=v;
  if (v<=3) v=1; else v-=1;  // Supress verbose output during iterative realignment and realignment
  
  // Align the first par.jdummy templates?
  if (par.jdummy>0)
    {
      if (v>=2) printf("Merging %i best hits to query alignment ...\n",par.jdummy);
      
      bin=0;
      nhits=0;
      hitlist.Reset();
      while (!hitlist.End() && nhits<par.jdummy)
	{
	  hit_cur = hitlist.ReadNext();
	  if (nhits>=imax(par.B,par.Z)) break;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
	  if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
	  nhits++;

	  if (hit_cur.irep>1) continue;  // Align only the best hit of the first par.jdummy templates
	  if (hit_cur.L>Lmaxmem) continue;  // Don't align to long sequences due to memory limit
	  if (hit_cur.Eval > par.e) continue; // Don't align hits with an E-value below the inclusion threshold

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
	  if (!fgetline(line,LINELEN,dbf)) {fprintf(stderr,"Error: end of file %s reached prematurely!\n",hit_cur.dbfile); exit(1);}
	  if (!strncmp(line,"HMMER",5))      // read HMMER format
	    {
	      format[bin] = 1;
	      read_from_db = t[bin]->ReadHMMer(dbf,hit_cur.dbfile);
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
	      //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
	      tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
	      t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
	      tali.FrequenciesAndTransitions(*(t[bin]));
	      format[bin] = 0;
	    }
	  else {
	    cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<hit_cur.dbfile<<"\'\n";
	    cerr<<"line = "<<line<<"\n";
	    exit(1);
	  }
	  fclose(dbf);
	  
	  if (read_from_db!=1)
	    {
	      fprintf(stderr,"Error: illegal format in %s while reading HMM %s\n",hit_cur.dbfile,hit_cur.name);
	      continue;
	    }

	  if (v>=2) fprintf(stderr,"Realigning with %s ***** \n",t[bin]->name);
	  ///////////////////////////////////////////////////
	  
	  N_aligned++;
	  if (v>=1 && !(N_aligned%10))
	    {
	      cout<<".";
	      if (!(N_aligned%500)) printf(" %-4i HMMs aligned\n",N_aligned);
	      cout.flush();
	    }
	  
	  // Prepare MAC comparison(s)
	  PrepareTemplate(q,*(t[bin]),format[bin]);
	  t[bin]->Log2LinTransitionProbs(1.0);
	  
	  // Align q to template in *hit[bin]
	  hit[bin]->Forward(q,*(t[bin]));
	  hit[bin]->Backward(q,*(t[bin]));
	  hit[bin]->MACAlignment(q,*(t[bin]));
	  hit[bin]->BacktraceMAC(q,*(t[bin]));
	  
	  // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
	  hit[bin]->score      = hit_cur.score;
	  hit[bin]->score_aass = hit_cur.score_aass;
	  hit[bin]->score_ss   = hit_cur.score_ss;
	  hit[bin]->Pval       = hit_cur.Pval;
	  hit[bin]->Pvalt      = hit_cur.Pvalt;
	  hit[bin]->logPval    = hit_cur.logPval;
	  hit[bin]->logPvalt   = hit_cur.logPvalt;
	  hit[bin]->logP1val   = hit_cur.logP1val;
	  hit[bin]->Eval       = hit_cur.Eval;
	  hit[bin]->logEval    = hit_cur.logEval;
	  hit[bin]->E1val      = hit_cur.E1val;
	  hit[bin]->Probab     = hit_cur.Probab;
	  hit[bin]->irep       = hit_cur.irep;
	  
	  // Replace original hit in hitlist with realigned hit
	  //hitlist.ReadCurrent().Delete();
	  hitlist.Delete().Delete();               // delete the list record and hit object
	  hitlist.Insert(*hit[bin]);
	  
	  // Read a3m alignment of hit and merge with Qali according to Q-T-alignment in hit[bin]
	  char ta3mfile[NAMELEN];
	  strcpy(ta3mfile,hit[bin]->file); // copy filename including path but without extension
	  strcat(ta3mfile,".a3m");
	  Qali.MergeMasterSlave(*hit[bin],ta3mfile);
	  
	  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
	  Qali.Compress("merged A3M file");
	  
	  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
	  Qali.N_filtered = Qali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
	  
	  // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
	  Qali.FrequenciesAndTransitions(q);
	  
	  if (!*par.clusterfile) { //compute context-specific pseudocounts?
	    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
	    q.PreparePseudocounts();
	  } else {
	    // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
	    q.PrepareContextSpecificPseudocounts();
	  }
	  
	  // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
	  q.AddAminoAcidPseudocounts();
	  q.CalculateAminoAcidBackground();
	  
	  // Transform transition freqs to lin space if not already done
	  q.AddTransitionPseudocounts();
	  q.Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done

	}
    }
  
  if (print_elapsed) ElapsedTimeSinceLastCall("(jdummy)");

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

  // Read all HMMs whose position is given in list realign_pos
  for (int idb=0; idb<ndb; idb++)
    {

      // Can we skip dbfiles[idb] because it contains no template to be realigned?
      if (! realign->Contains(dbfiles[idb])) continue;
      
      // Should depend on mode -Ohhm: no resort; else resort !!!!!!!!! ////////////////////////////////////////
      // Sort dbfile list of ftell positions in ascending order
      realign->Show(dbfiles[idb])->SortList();
      
      // Open HMM database file dbfiles[idb]
      FILE* dbf=fopen(dbfiles[idb],"rb");
      if (!dbf) OpenFileError(dbfiles[ndb]);
      read_from_db=1;
      int index_prev=-1;
      
      ///////////////////////////////////////////////////////////////////////////////////////
      // The loop (reads HMMs from the database file and) submits jobs into free bins as soon as they become available
      realign->Show(dbfiles[idb])->Reset();
      while (! realign->Show(dbfiles[idb])->End())
	{
	  // Submit jobs until no bin is free anymore
	  while (! realign->Show(dbfiles[idb])->End() && jobs_submitted+jobs_running<bins)
	    {
	      
	      // Allocate free bin
	      bin = PickBin(FREE);
	      if (bin<0) {
		fprintf(stderr,"Error during realignment: found no free bin! jobs running: %i  jobs_submitted:%i  threads:%i\n",jobs_running,jobs_submitted,threads);
		for (bin=0; bin<bins; bin++) fprintf(stderr,"bin_status[%i]=%i\n",bin,bin_status[bin]);
		exit(6);
	      }
	      
	      // Forward stream position to start of next database HMM to be realigned
	      hit[bin]->index = realign->Show(dbfiles[idb])->ReadNext().index;  // give hit a unique index for HMM
	      if (hit[bin]->index <= index_prev) continue;
	      index_prev = hit[bin]->index;
	      fseek(dbf,realign->Show(dbfiles[idb])->ReadCurrent().ftellpos,SEEK_SET);
	      
	      //                fprintf(stderr,"dbfile=%-40.40s  index=%-5i  ftellpos=%i\n",dbfiles[idb],realign->Show(dbfiles[idb])->ReadCurrent().index,(unsigned int) realign->Show(dbfiles[idb])->ReadCurrent().ftellpos);
	      
	      char path[NAMELEN];
	      Pathname(path,dbfiles[idb]);
	      
	      ///////////////////////////////////////////////////
	      // Read next HMM from database file
	      if (!fgetline(line,LINELEN,dbf)) {fprintf(stderr,"Error: end of file %s reached prematurely!\n",dbfiles[idb]); exit(1);}
	      if (!strncmp(line,"HMMER",5))      // read HMMER format
		{
		  format[bin] = 1;
		  read_from_db = t[bin]->ReadHMMer(dbf,dbfiles[idb]);
		}
	      else if (!strncmp(line,"HH",2))     // read HHM format
		{
		  format[bin] = 0;
		  read_from_db = t[bin]->Read(dbf,path);
		}
	      else if (!strncmp(line,"NAME",4))  // The following lines are for backward compatibility of HHM format version 1.2 with 1.1
		{
		  format[bin] = 0;
		  fseek(dbf,realign->Show(dbfiles[idb])->ReadCurrent().ftellpos,SEEK_SET); // rewind to beginning of line
		  read_from_db = t[bin]->Read(dbf,path);
		}
	      else if (line[0]=='#')                 // read a3m alignment
		{
		  Alignment tali;
		  tali.Read(dbf,dbfiles[idb],line);
		  tali.Compress(dbfiles[idb]);
		  //                  qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
		  tali.N_filtered = tali.Filter(par.max_seqid_db,par.coverage_db,par.qid_db,par.qsc_db,par.Ndiff_db);
		  t[bin]->name[0]=t[bin]->longname[0]=t[bin]->fam[0]='\0';
		  tali.FrequenciesAndTransitions(*(t[bin]));
		  format[bin] = 0;
		}
	      else {
		cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<dbfiles[idb]<<"\'\n";
		cerr<<"line = "<<line<<"\n";
		exit(1);
	      }
	      
	      
	      if (read_from_db==2) continue;  // skip current HMM or reached end of database
	      if (read_from_db==0) break;     // finished reading HMMs
	      if (v>=2) fprintf(stderr,"Realigning with %s\n",t[bin]->name);
	      ///////////////////////////////////////////////////
	      
	      hit[bin]->dbfile = new(char[strlen(dbfiles[idb])+1]);
	      strcpy(hit[bin]->dbfile,dbfiles[idb]); // record db file name from which next HMM is read
	      
	      N_aligned++;
	      if (v>=1 && !(N_aligned%10))
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
		  // If no submitted jobs are in the queue we have to wait for a new job, but max. 1 second ...
		  struct timespec ts;
		  clock_gettime(CLOCK_REALTIME,&ts);
		  ts.tv_sec += 1;
		  rc = pthread_cond_timedwait(&finished_job, &bin_status_mutex,&ts);
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
  if (v1>=1) cout<<"\n";
  v=v1;
  
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
	  //hitlist.ReadCurrent().Delete();
	  hitlist.Delete().Delete();               // delete the list record and hit object
	  // Make sure only realigned alignments get displayed!
	  if (par.B>par.Z) par.B--; else if (par.B==par.Z) {par.B--; par.Z--;} else par.Z--;
	  if (par.b>par.z) par.b--; else if (par.b==par.z) {par.b--; par.z--;} else par.z--;
	}
      else nhits++;
    }

  // Delete realign hash with lists
  realign->Reset();
  while (!realign->End())
    delete(realign->ReadNext()); // delete List<Posindex> to which realign->ReadNext() points
  delete(realign);

  if (print_elapsed) ElapsedTimeSinceLastCall("(realign)");
}


/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  Hit hit_cur;

#ifdef PTHREAD
  pthread_attr_init(&joinable);  // initialize attribute set with default values
  if (pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)!=0) // set attribute 'joinable'
    cerr<<"Error "<<pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)<<": could not set detach state for thread attibute.\n";
#endif

  SetDefaults();
  par.jdummy = 3;
  par.Ndiff = 1000;
  par.filter_thresh=0.05;
  strcpy(pre_mode,"csblast");
  strcpy(par.outfile,"");
  N_searched=0;
  previous_hits = new Hash<Hit>(1631,hit_cur);
  par.block_shading = new Hash<int*>;
  par.block_shading_counter = new Hash<int>;
  
  // Make command line input globally available
  par.argv=argv;
  par.argc=argc;
  RemovePathAndExtension(program_name,argv[0]);

  // Set default for config-file
  Pathname(config_file,argv[0]);
  strcat(config_file,"hhblast.config");

  // Enable changing verbose mode before command line are processed
  for (int i=1; i<argc; i++)
    {
      if (!strcmp(argv[i],"-conf"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no config file following -conf\n"; exit(4);}
          else strcpy(config_file,argv[i]);
        }
      else if (argc>1 && !strcmp(argv[i],"-v0")) v=0;
      else if (argc>1 && !strcmp(argv[i],"-v1")) v=1;
      else if (argc>2 && !strcmp(argv[i],"-v")) v=atoi(argv[i+1]);
    }

  // Read config file?
  if (is_regular_file(config_file))
    {
      // Process default otpions from .hhconfig file
      ReadConfigFile(config_file);
    }

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc,argv);

  // Check needed files
  if (!*db || !*dbhhm)
    {help(); cerr<<endl<<"Error in "<<program_name<<": database missing (see config-file or -db and -dbhhm)\n"; exit(4);}
  if (!strcmp(pre_mode,"csblast") && (!*csblast || !*csblast_db))
    {help(); cerr<<endl<<"Error in "<<program_name<<": missing csBLAST directory (see config-file or -csblast and -csblast_db)\n"; exit(4);}
  if (!*blast)
    {help(); cerr<<endl<<"Error in "<<program_name<<": missing BLAST directory (see config-file or -blast)\n"; exit(4);}
  if (par.showpred==1 && (!*psipred || !*psipred_data))
    {help(); cerr<<endl<<"Error in "<<program_name<<": missing PsiPred directory (see config-file or -psipred and -psipred_data)\n"; exit(4);}
 
  // Create tmp-directory
  if (mkstemp(tmp_file) == -1) {
    cerr << "ERROR! Could not create tmp-file!\n"; 
    exit(4);
  }
  tmp_infile = (string)tmp_file + ".fas";
  tmp_psifile = (string)tmp_file + ".psi";

  // Check input files
  CheckInputFiles();
  
  // Check for threads
  if (threads<=1) threads=0;
  else if (threads>MAXTHREADS)
    {
      threads=MAXTHREADS;
      if (v>=1) fprintf(stderr,"WARNING: number of CPUs set to maximum value of %i\n",MAXTHREADS);
    }
  
  // Check option compatibilities
  if (par.nseqdis>MAXSEQDIS-3-par.showcons) par.nseqdis=MAXSEQDIS-3-par.showcons; //3 reserved for secondary structure
  if (par.aliwidth<20) par.aliwidth=20;
  if (par.pca<0.001) par.pca=0.001; // to avoid log(0)
  if (par.b>par.B) par.B=par.b;
  if (par.z>par.Z) par.Z=par.z;

  // E-value shift of csBLAST
  if (*par.clusterfile && !strcmp(csblast,"/cluster/bioprogs/csblast-2.0.3-linux64/bin")) {
    cout << "CS-BLAST 2.0.3\n";
    e_psi = e_psi / 2;
  }

  // Get Prefilter Pvalue (Evalue / DBsize)
  float dbsize = 0;
  FILE *stream;
  command = (string)blast + "/fastacmd -d " + (string)db + " -I";
  stream = popen(command.c_str(), "r");
  while (fgets(line, LINELEN, stream))
    {
      char tmp_number[NAMELEN];
      char tmp_name[NAMELEN];
      ptr=strscn(line); 
      ptr=strwrd(tmp_number,ptr);
      ptr=strwrd(tmp_name,ptr);
      if (!strcmp(tmp_name,"sequences;")) {
	strtrd(tmp_number,",");
	dbsize = atof(tmp_number);
	break;
      }
    }
  pclose(stream);
  if (dbsize == 0)
    {cerr<<endl<<"Error in "<<program_name<<": Could not determine DB-size of prefilter db ("<<db<<")\n"; exit(4);}
  par.hhblast_prefilter_logpval=-log(e_psi / dbsize);

  // Input parameters
  if (v>=3)
    {
      cout<<"Input file       :   "<<par.infile<<"\n";
      cout<<"Output file      :   "<<par.outfile<<"\n";
      cout<<"Prefilter mode   :   "<<pre_mode<<"\n";
      cout<<"Prefilter DB     :   "<<db<<"\n";
      cout<<"HHM DB           :   "<<dbhhm<<"\n";
      cout<<"Prefilter Pval   :   "<<(e_psi / dbsize)<<"\n";
    }

  // Set secondary structure substitution matrix
  if (par.ssm) SetSecStrucSubstitutionMatrix();

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix();

  int v1=v;
  if (v<=3) v=1; else v-=2;

  // Read input file (HMM or alignment format) without adding pseudocounts
  ReadInput(par.infile, q);
  
  v=v1;

  // set Qali names
  Qali.longname = new(char[strlen(q.longname)+1]);
  strcpy(Qali.longname,q.longname);
  strcpy(Qali.name,q.name);
  strcpy(Qali.fam,q.fam);

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags) q.NeutralizeTags();

  // Prepare multi-threading - reserve memory for threads, intialize, etc.
  if (threads==0) bins=1; else bins=iround(threads*1.2+0.5);
  for (bin=0; bin<bins; bin++)
    {
      t[bin]=new HMM;   // Each bin has a template HMM allocated that was read from the database file
      hit[bin]=new Hit; // Each bin has an object of type Hit allocated ...
      hit[bin]->AllocateBacktraceMatrix(q.L+2,MAXRES); // ...with a separate dynamic programming matrix (memory!!)
      if (par.realign==1)
	{
	  hit[bin]->AllocateForwardMatrix(q.L+2,MAXRES);
	  hit[bin]->AllocateBackwardMatrix(q.L+2,MAXRES);
	}

    }
  format = new(int[bins]);

  if (print_elapsed) ElapsedTimeSinceLastCall("(initialize)");

  //////////////////////////////////////////////////////////
  // Main loop
  //////////////////////////////////////////////////////////

  if (v>=2) printf("\n******************************************************\n* Building alignment for query with %i rounds HHblast *\n******************************************************\n\n",num_rounds);

  for (int round = 1; round <= num_rounds; round++) {

    if (v>=2) printf("\nRound: %i\n",round);

    // Settings for different rounds
    if (par.jdummy > 0 && round > 1 && previous_hits->Size() > (par.jdummy-1))
      {
	if (v>=3) printf("Set jdummy to 0! (jdummy: %i   round: %i   hits.Size: %i)\n",par.jdummy,round,previous_hits->Size());
	par.jdummy = 0;
      }
    else 
      {
	par.jdummy -= previous_hits->Size();
      }

    if (round == num_rounds && nodiff)
      {
	if (v>=4) printf("Set Ndiff to 0!\n");
	par.Ndiff = 0;
      }

    // Write query HHM file?
    if (*query_hhmfile) 
      {
	// Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
	q.AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
	q.CalculateAminoAcidBackground();

	q.WriteToFile(query_hhmfile);
      }
    
    // Add Pseudocounts
    if (!*par.clusterfile) { //compute context-specific pseudocounts?
      // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
      q.PreparePseudocounts();
    } else {
      // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
      q.PrepareContextSpecificPseudocounts();
    }
    
    // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    q.AddAminoAcidPseudocounts();
    q.CalculateAminoAcidBackground();
    
    // Transform transition freqs to lin space if not already done
    q.AddTransitionPseudocounts();
    
    // Prefilter with csBLAST/PSI-BLAST
    if (v>=2) printf("Pre-filtering with %s ...\n",(strcmp(pre_mode,"csblast"))?"PSI-BLAST":"CS-BLAST");
    stringstream ss;
    if (round == 1 && !is_regular_file(tmp_psifile.c_str())) 
      {
    	if (strcmp(pre_mode,"csblast")) 
    	  ss << blast << "/blastpgp -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 2> /dev/null";
    	else
    	  ss << csblast << "/csblast -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -D " << csblast_db << " --blast-path " << blast << " --no-penalty 2> /dev/null";
      } 
    else
      {
    	if (strcmp(pre_mode,"csblast")) 
    	  ss << blast << "/blastpgp -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -B " << tmp_psifile << " 2> /dev/null";
    	else
    	  ss << csblast << "/csblast -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -B " << tmp_psifile << " -D " << csblast_db << " --blast-path " << blast << " --no-penalty 2> /dev/null";
      }
    
    command = ss.str();
    if (v>=4) cout << "Command: " << command << "\n";
    
    // and extract prefiltered HHMs
    for (int idb=0; idb<ndb_new; idb++) delete[](dbfiles_new[idb]);
    for (int idb=0; idb<ndb_old; idb++) delete[](dbfiles_old[idb]);
    if(doubled) delete doubled;
    doubled = new(Hash<char>);
    doubled->New(16381,0);
    par.block_shading->Reset();
    while (!par.block_shading->End())
	delete[] (par.block_shading->ReadNext()); 
    par.block_shading->New(16381,NULL);
    par.block_shading_counter->New(16381,NULL);
    ndb_new = ndb_old = 0;
    int pos;
    int count=0;
    char actual_hit[NAMELEN];
    int* block = new(int[400]);
    strcpy(actual_hit,"");

    stream = popen(command.c_str(), "r");
    while (fgets(line, LINELEN, stream))
      {
	ptr=strscn(line);
	ptr=strcut(ptr);
	if (ptr==NULL)
	  {cerr<<endl<<"Error in "<<program_name<<": Cannot parse result of pre-filtering!\n"; exit(6);}
	char tmp_name[NAMELEN];
	char db_name[NAMELEN];
	char tmp[NAMELEN];
	ptr=strwrd(tmp_name,ptr);
	
	if (block_filter)
	  {
	    // Write block for template
	    if (strcmp(actual_hit,"") && strcmp(actual_hit,tmp_name))   // New template
	      {
		//printf("Add to block shading   key: %s    data:",actual_hit);
		// for (int i = 0; i < count; i++)
		//   printf(" %i,",block[i]);
		// printf("\n");
		par.block_shading->Add(actual_hit,block);
		par.block_shading_counter->Add(actual_hit,count);
		block = new(int[400]);
		count = 0;
	      }
	    if (count >= 400) { continue; }
	    strcpy(actual_hit,tmp_name);
	    // Get block of HSP
	    ptr=strwrd(tmp,ptr); // sequence identity
	    pos=strint(ptr); // ali length
	    pos=strint(ptr); // mismatches
	    pos=strint(ptr); // gap openings
	    pos=strint(ptr); // query start
	    block[count++]=pos;
	    pos=strint(ptr); // query end
	    block[count++]=pos;
	    pos=strint(ptr); // subject start
	    block[count++]=pos;
	    pos=strint(ptr); // subject end
	    block[count++]=pos;	    
	  }
	
	if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	  {
	    substr(tmp,tmp_name,3,11);
	    substr(db_name,tmp,0,1);
	    strcat(db_name,"/");
	    strcat(db_name,tmp);
	    strcat(db_name,".db");
	  }
	else                              // other database
	  {
	    strcpy(db_name,tmp_name);
	    strtr(db_name,"|", "_");
	    strcat(db_name,".hhm");
	  }

	if (! doubled->Contains(db_name))
	  {
	    doubled->Add(db_name);
	    // check, if DB was searched in previous rounds 
	    strcat(tmp_name,"__1");  // irep=1
	    if (previous_hits->Contains(tmp_name))
	      {
		dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
		strcpy(dbfiles_old[ndb_old],dbhhm);
		strcat(dbfiles_old[ndb_old],"/");
		strcat(dbfiles_old[ndb_old],db_name);
		if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
		ndb_old++;
	      }
	    else 
	      {
		dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
		strcpy(dbfiles_new[ndb_new],dbhhm);
		strcat(dbfiles_new[ndb_new],"/");
		strcat(dbfiles_new[ndb_new],db_name);
		if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
		ndb_new++;
	      }
	  }
      }
    pclose(stream);
    if (block_filter && strcmp(actual_hit,""))   // New template
      {
	// printf("Add to block shading   key: %s    data:",actual_hit);
	// for (int i = 0; i < count; i++)
	//   printf(" %i,",block[i]);
	// printf("\n");
	par.block_shading->Add(actual_hit,block);
	par.block_shading_counter->Add(actual_hit,count);
      }

    if (v>=3) printf("Number of new extracted HMMs: %i\n",ndb_new);
    if (v>=3) printf("Number of extracted HMMs (previous searched): %i\n",ndb_old);

    if (print_elapsed) ElapsedTimeSinceLastCall("(prefiltering)"); 

    // Search datbases
    if (v>=2) printf("HMM search ...\n");

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
	if (round < num_rounds && v>=2)
	  printf("No new hits found in round %i => Stop searching\n",round);

	if (ndb_old > 0 && realign_old_hits)
	  {
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
	    perform_viterbi_search(ndb_new+previous_hits->Size());
	  }
      }

    // Realign all hits with MAC algorithm?
    if (par.realign)
      perform_realign(dbfiles_new,ndb_new);

    // Generate alignment for next iteration
    if (round < num_rounds || *par.alnfile || *par.psifile || *par.hhmfile || *alis_basename)
      {
	int v1=v;
	if (v<=3) v=1; else v-=2;
	
	// If new hits found, merge hits to query alignment
	if (new_hits != 0)
	  {
	    char ta3mfile[NAMELEN];
	    
	    if (v>=2) printf("Merging hits to query alignment ...\n");
	    
	    // For each template below threshold
	    hitlist.Reset();
	    while (!hitlist.End())
	      {
		hit_cur = hitlist.ReadNext();
		if (hit_cur.Eval > 100.0*par.e) break; // E-value much too large
		if (hit_cur.Eval > par.e) continue; // E-value too large
		stringstream ss_tmp;
		ss_tmp << hit_cur.name << "__" << hit_cur.irep;
		if (previous_hits->Contains((char*)ss_tmp.str().c_str())) continue;  // Already in alignment
		
		// Read a3m alignment of hit from <file>.a3m file and merge into Qali alignment
		strcpy(ta3mfile,hit_cur.file); // copy filename including path but without extension
		strcat(ta3mfile,".a3m");
		Qali.MergeMasterSlave(hit_cur,ta3mfile);
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
	    
	    // Write PSI-alignment for next round prefiltering
	    Qali.WriteToFile(tmp_psifile.c_str(),"psi");
	    
	    if (print_elapsed) ElapsedTimeSinceLastCall("(write PSI-alignment)");
	  }
	  
	// if needed, calculate SSpred
	if (par.showpred && Qali.L>25 && (*alis_basename || round == num_rounds || new_hits == 0))
	  {
	    char ss_pred[MAXRES];
	    char ss_conf[MAXRES];
	    
	    CalculateSS(ss_pred, ss_conf);
	    
	    Qali.AddSSPrediction(ss_pred, ss_conf);
	    
	    if (print_elapsed) ElapsedTimeSinceLastCall("(calculate SS_Pred)");
	  }
	
	// Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
    	Qali.FrequenciesAndTransitions(q,NULL,true);
    	
	if (print_elapsed) ElapsedTimeSinceLastCall("(Calculate AA frequencies and transitions)");

	if (*alis_basename)
	  {
	    stringstream ss_tmp;
	    ss_tmp << alis_basename << "_" << round << ".a3m";
	    Qali.WriteToFile(ss_tmp.str().c_str(),"a3m");
	  }

	v=v1;
	
      }
    
    if (new_hits == 0 || round == num_rounds) 
      break;

    if (v>=2 && filter) printf("%i new hits found in round %i\n",new_hits,round);

    // Write good hits to previous_hits hash and clear hitlist
    hitlist.Reset();
    while (!hitlist.End())
      {
	hit_cur = hitlist.ReadNext();
	stringstream ss_tmp;
	ss_tmp << hit_cur.name << "__" << hit_cur.irep;
	if (!filter || hit_cur.Eval > par.e || previous_hits->Contains((char*)ss_tmp.str().c_str()))
	  hit_cur.Delete(); // Delete hit object
	else
	  previous_hits->Add((char*)ss_tmp.str().c_str(), hit_cur);

	hitlist.Delete(); // Delete list record

      }

    if (print_elapsed) ElapsedTimeSinceLastCall("(end of this round)");

  } // end for-loop rounds
  

  //////////////////////////////////////////////////////////
  // Result section
  //////////////////////////////////////////////////////////

  // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
  if (*par.scorefile) {
    if (v>=3) printf("Printing scores file ...\n");
    hitlist.PrintScoreFile(q);
  }

  // Print summary listing of hits
  if (v>=3) printf("Printing hit list ...\n");
  hitlist.PrintHitList(q,par.outfile);

  // Write only hit list to screen?
  if (v==2 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,109); // write only hit list to screen

  // Print alignments of query sequences against hit sequences
  hitlist.PrintAlignments(q,par.outfile);

  // Write whole output file to screen? (max 10000 lines)
  if (v>=3 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,10009);

  // Generate output alignment or HMM file?
  if (*par.alnfile || *par.psifile || *par.hhmfile)
    {
      // Write output PSI-BLAST-formatted alignment?
      if (*par.psifile) Qali.WriteToFile(par.psifile,"psi");

      // Write output HHM file?
      if (*par.hhmfile) 
      {
	// Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
	q.AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
	q.CalculateAminoAcidBackground();

	q.WriteToFile(par.hhmfile);
      }

      // Write output A3M alignment?
      if (*par.alnfile) Qali.WriteToFile(par.alnfile,"a3m");
      
    }
  
  // Delete memory for dynamic programming matrix
  for (bin=0; bin<bins; bin++)
    {
      hit[bin]->DeleteBacktraceMatrix(q.L+2);
      if (par.realign)
	{
	  hit[bin]->DeleteForwardMatrix(q.L+2);
	  hit[bin]->DeleteBackwardMatrix(q.L+2);
	}
      delete hit[bin];
      delete t[bin];
    }
  if (par.dbfiles) delete[] par.dbfiles;
  for (int idb=0; idb<ndb_new; idb++) delete[](dbfiles_new[idb]);
  for (int idb=0; idb<ndb_old; idb++) delete[](dbfiles_old[idb]);
  if (format) delete[](format);
  if (par.blafile) delete[] par.blafile;
  if (par.exclstr) delete[] par.exclstr;
  delete doubled;
  par.block_shading->Reset();
  while (!par.block_shading->End())
    delete[] (par.block_shading->ReadNext()); 
  delete par.block_shading;
  previous_hits->Reset();
  while (!previous_hits->End())
    {
      previous_hits->ReadNext().Delete(); // Delete hit object
    }
  delete previous_hits;
  
  // Delete content of hits in hitlist
  hitlist.Reset();
  while (!hitlist.End())
    {
      hitlist.Delete().Delete(); // Delete list record and hit object
    }

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

  // Remove temp-files
  command = "rm " + (string)tmp_file + "*";
  runSystem(command);
  //cout << "Command: " << command << "!\n";

  exit(0);
} //end main

//////////////////////////////////////////////////////////////////////////////////////////////////////
// END OF MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////



