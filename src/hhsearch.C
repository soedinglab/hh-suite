// hhsearch.C:
// Search for a multiple alignment (transformed into HMM) in a profile HMM database
// Compile:              g++ hhsearch.C -o hhsearch -O3 -lpthread -lrt -fno-strict-aliasing
// Compile with efence:  g++ hhsearch.C -o hhsearch -lefence -lpthread -lrt -O -g
// Compile for Valgrind: g++ hhsearch.C -o hhsearch2 -lpthread -lrt -O -g
// With wnlib:           g++ hhsearch.C /home/soeding/programs/wnlib/acc/text.a  -o hhsearch -O3 -lpthread -lrt -fno-strict-aliasing -g -I/home/soeding/programs/wnlib/acc/h/ -L/home/soeding/programs/electric-fence-2.1.13/
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
//#include <new>
//#include "efence.h"
//#include "efence.c"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

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
int N_searched;  // Number of HMMs searched



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
#endif


/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHsearch %s\n",VERSION_AND_DATE);
  printf("Search a database of HMMs with a query alignment or query HMM\n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query -d database [options]                       \n",program_name);
  printf(" -i <file>     input query alignment (A2M, A3M, FASTA) or HMM\n");
  printf(" -d <file>     HMM database of concatenated HMMs in hhm, HMMER, or A3M format,\n");
  printf("               OR, if file has extension pal, list of HMM file names, one per\n");
  printf("               line. Multiple dbs, HMMs, or pal files with -d '<db1> <db2>...'\n");
  printf("\n");
  printf("Output options:                                                              \n");
  printf(" -cal          calibrate query HMM (write mu and lamda into hhm file)        \n");
  printf(" -o <file>     write results in standard format to file (default=<infile.hhr>)\n");
  printf(" -ofas <file>  write pairwise alignments in FASTA (-oa2m: A2M, -oa3m: A3M) format\n");
  printf(" -v <int>      verbose mode: 0:no screen output  1:only warings  2: verbose   \n");
  printf(" -seq <int>    max. number of query/template sequences displayed (def=%i) \n",par.nseqdis);
  printf(" -nocons       don't show consensus sequence in alignments (default=show)     \n");
  printf(" -nopred       don't show predicted 2ndary structure in alignments (default=show)\n");
  printf(" -nodssp       don't show DSSP 2ndary structure in alignments (default=show)  \n");
  printf(" -ssconf       show confidences for predicted 2ndary structure in alignments\n");
  printf(" -aliw <int>   number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -p <float>    minimum probability in summary and alignment list (def=%G)   \n",par.p);
  printf(" -E <float>    maximum E-value in summary and alignment list (def=%G)       \n",par.E);
  printf(" -Z <int>      maximum number of lines in summary hit list (def=%i)         \n",par.Z);
  printf(" -z <int>      minimum number of lines in summary hit list (def=%i)         \n",par.z);
  printf(" -B <int>      maximum number of alignments in alignment list (def=%i)      \n",par.B);
  printf(" -b <int>      minimum number of alignments in alignment list (def=%i)      \n",par.b);
  printf("               Remark: you may use 'stdin' and 'stdout' instead of file names\n");
  printf("\n");
  printf("Filter input alignment (options can be combined):                             \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
  printf("\n");
  printf("Input alignment format:                                                       \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");
  printf("HMM-HMM alignment options:                                                    \n");
  printf(" -realign      realign displayed hits with max. accuracy (MAC) algorithm \n");
  printf(" -norealign    do NOT realign displayed hits with MAC algorithm (def=realign)\n");
  printf(" -mact [0,1[   posterior probability threshold for MAC re-alignment (def=%.3f)\n",par.mact);
  printf("               Parameter controls alignment greediness: 0:global >0.1:local\n");
  printf(" -glob/-loc    use global/local alignment mode for searching/ranking (def=local)\n");
//   printf(" -vit          use Viterbi algorithm for searching/ranking (default)          \n");
//   printf(" -mac          use Maximum Accuracy MAC algorithm for searching/ranking\n");
//   printf(" -forward      use Forward probability for searching                       \n");
  printf(" -alt <int>    show up to this many significant alternative alignments(def=%i)\n",par.altali);
  printf(" -excl <range> exclude query positions from the alignment, e.g. '1-33,97-168' \n");
  printf(" -shift [-1,1] score offset (def=%-.2f)                                       \n",par.shift);
  printf(" -corr [0,1]   weight of term for pair correlations (def=%.2f)                \n",par.corr);
  printf(" -ssm  0-4     0:   no ss scoring                                             \n");
  printf("               1,2: ss scoring after or during alignment  [default=%1i]       \n",par.ssm);
  printf("               3,4: ss scoring after or during alignment, predicted vs. predicted \n");
  printf(" -ssw [0,1]    weight of ss score  (def=%-.2f)                                \n",par.ssw);
  printf("\n");
  printf("Other options:                                                                \n");
//  printf(" -w  <file>    calculate scores transitively with this database weight file  \n");
  printf(" -def          read default options from ./.hhdefaults or <home>/.hhdefault.  \n");
  printf("               Write 'hhsearch', 'hhmake' and/or 'hhfilter' etc. in one line, \n");
  printf("               followed by its list of options, one per line.\n");
  printf(" -cpu <int>    number of CPUs to use (for shared memory SMPs) (default=1)\n");
#ifndef PTHREAD
  printf("(The -cpu option is inactive since POSIX threads ae not supported on your platform)\n");
#endif
  printf("\n");
  printf("Example: %s -i a.1.1.1.a3m -d scop70_1.71.hhm             \n",program_name);
  cout<<endl;

//  printf("Network scoring options:                                                 \n");
//  printf(" -netn <int>   max number of intermediate HMMs on comparison paths (def=%i)\n",par.netn);
//  printf(" -neta [0,1]   strength of transitive vs direct comparison (def=%f)        \n",par.neta);
//  printf(" -netb [0,1]   strength of covariance weighting (0=no weighting) (def=%f)  \n",par.netb);

//   printf("More help:                                                         \n");
//   printf(" -h out        options for formatting ouput                        \n");
//   printf(" -h hmm        options for building HMM from multiple alignment    \n");
//   printf(" -h gap        options for setting gap penalties                   \n");
//   printf(" -h ali        options for HMM-HMM alignment                       \n");
//   printf(" -h other      other options \n");
//   printf(" -h all        all options \n");
 }

void help_out()
{
  printf("\n");
  printf("Output options:                                                           \n");
  printf(" -o <file>     write output alignment to file (default=%s)\n",par.outfile);
  printf("               Omit this option to write to standard output\n");
  printf(" -v            verbose mode (default: show only warnings)                 \n");
  printf(" -v 0          suppress all screen output                                 \n");
  printf(" -p <float>    minimum probability in summary and alignment list (def=%G) \n",par.p);
  printf(" -E <float>    maximum E-value in summary and alignment list (def=%G)     \n",par.E);
  printf(" -Z <int>      maximum number of lines in summary hit list (def=%i)       \n",par.Z);
  printf(" -z <int>      minimum number of lines in summary hit list (def=%i)       \n",par.z);
  printf(" -B <int>      maximum number of alignments in alignment list (def=%i)    \n",par.B);
  printf(" -b <int>      minimum number of alignments in alignment list (def=%i)    \n",par.b);
  printf(" -seq  [1,inf[ max. number of query/template sequences displayed  (def=%i)  \n",par.nseqdis);
  printf(" -nocons       don't show consensus sequence in alignments (default=show) \n");
  printf(" -nopred       don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(" -nodssp       don't show DSSP 2ndary structure in alignments (default=show) \n");
  printf(" -ssconf       show confidences for predicted 2ndary structure in alignments\n");
  printf(" -aliw [40,..[ number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -cal          calibrate query HMM (write mu and lamda into hmm file)     \n");
  printf(" -dbstrlen     max length of database string to be printed in hhr file\n");
}

void help_hmm()
{
  printf("\n");
  printf("Filter input alignment (options can be combined):                         \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
//   printf(" -csc  [0,100] minimum score per column with core alignment (def=%-.2f)\n",par.coresc);
//   printf(" -qscc [0,100] minimum score per column of core sequence with query (def=%-.2f)\n",par.qsc_core);
  printf("                                                                          \n");
  printf("HMM-building options:                                                     \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf(" -tags         do NOT neutralize His-, C-myc-, FLAG-tags, and \n");
  printf("               trypsin recognition sequence to background distribution    \n");
  printf("                                                                          \n");
  printf("Pseudocount options:                                                      \n");
  printf(" -Gonnet       use the Gonnet substitution matrix (default)               \n");
  printf(" -BlosumXX     use a Blosum substitution matrix (XX=30,40,50,65, or 80)    \n");
  printf(" -pcm  0-2     Pseudocount mode (default=%-i)                             \n",par.pcm);
  printf("               tau = substitution matrix pseudocount admixture            \n");
  printf("               0: no pseudo counts:     tau = 0                           \n");
  printf("               1: constant              tau = a                           \n");
  printf("               2: divergence-dependent: tau = a/(1 + ((Neff-1)/b)^c)      \n");
  printf("                  Neff=( (Neff_q^d+Neff_t^d)/2 )^(1/d)                    \n");
  printf("                  Neff_q = av number of different AAs per column in query \n");
  printf("               3: constant divergence pseudocounts \n");
  printf("               4: divergence-dependent, with composition-adjusted matrix  \n");
  printf(" -pca  [0,1]   overall pseudocount admixture (def=%-.1f)                      \n",par.pca);
  printf(" -pcb  [1,inf[ threshold for Neff) (def=%-.1f)                     \n",par.pcb);
  printf(" -pcc  [0,3]   extinction exponent for tau(Neff)  (def=%-.1f)     \n",par.pcc);
  printf(" -pcw  [0,3]   weight of pos-specificity for pcs  (def=%-.1f)      \n",par.pcw);
  printf(" -cs   <file>  use context-specific pseudocounts computed from given context library (def=off)\n");
  printf(" -csw  [0,inf] weight of central position in cs pseudocount mode (def=%.1f)\n", par.csw);
  printf(" -csb  [0,1]   weight decay parameter for positions in cs pc mode (def=%.1f)\n", par.csb);
}

void help_gap()
{
  printf("\n");
  printf("Gap cost options:                                                         \n");
  printf(" -gapb [0,inf[ transition pseudocount admixture (def=%-.2f)               \n",par.gapb);
  printf(" -gapd [0,inf[ Transition pseudocount admixture for opening gap (default=%-.2f)\n",par.gapd);
  printf(" -gape [0,1.5] Transition pseudocount admixture for extending gap (def=%-.2f)\n",par.gape);
  printf(" -gapf ]0,inf] factor for increasing/reducing the gap open penalty for deletes (def=%-.2f)\n",par.gapf);
  printf(" -gapg ]0,inf] factor for increasing/reducing the gap open penalty for deletes (def=%-.2f)\n",par.gapg);
  printf(" -gaph ]0,inf] factor for increasing/reducing the gap extension penalty for deletes(def=%-.2f)\n",par.gaph);
  printf(" -gapi ]0,inf] factor for increasing/reducing the gap extension penalty for inserts(def=%-.2f)\n",par.gapi);
  printf(" -egq  [0,inf[ penalty (bits) for end gaps aligned to query residues (def=%-.2f)\n",par.egq);
  printf(" -egt  [0,inf[ penalty (bits) for end gaps aligned to template residues (def=%-.2f)\n",par.egt);
}

void help_ali()
{
  printf("\n");
  printf("Alignment options:  \n");
  printf(" -realign      realign displayed hits with MAC algorithm \n");
  printf(" -norealign    do NOT realign displayed hits with MAC algorithm (def=realign)\n");
  printf(" -mact [0,1]   posterior prob. threshold in MAC (re-)alignment (def=%-.3f) \n",par.mact);
  printf(" -glob/-loc    use global/local alignment mode for searching/ranking (def=local)\n");
  printf(" -vit          use Viterbi algorithm for searching/ranking (default)       \n");
  printf(" -mac          use Maximum Accuracy (MAC) algorithm for searching/ranking\n");
  printf(" -forward      use Forward probability for searching                       \n");
  printf("               This controls alignment greediness: 0:global ~0.2-0.99:local\n");
  printf(" -alt <int>    show up to this number of alternative alignments (def=%i)  \n",par.altali);
  printf(" -excl <range> exclude query positions from the alignment, e.g. '1-33,97-168'\n");
  printf(" -sc   <int>   amino acid score         (tja: template HMM at column j) (def=%i)\n",par.columnscore);
  printf("        0      = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
  printf("        1      = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
  printf("        2      = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
  printf("        3      = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
  printf(" -corr [0,1]   weight of term for pair correlations (def=%.2f)             \n",par.corr);
  printf(" -shift [-1,1] score offset (def=%-.3f)                                    \n",par.shift);
  printf(" -r            repeat identification: multiple hits not treated as independent\n");
  printf(" -ssm  0-4     0:no ss scoring [default=%i]               \n",par.ssm);
  printf("               1:ss scoring after alignment                                \n");
  printf("               2:ss scoring during alignment (default)                     \n");
  printf("               3:ss scoring after alignment; use only psipred (not dssp)   \n");
  printf("               4:ss scoring during alignment use only psipred (not dssp)   \n");
  printf(" -ssw  [0,1]   weight of ss score compared to column score (def=%-.2f)     \n",par.ssw);
  printf(" -ssa  [0,1]   SS substitution matrix = (1-ssa)*I + ssa*full-SS-substition-matrix [def=%-.2f)\n",par.ssa);
  printf(" -ssgap        Gap opening within SS elements costs iX bits after ith residue,\n");
  printf("               where X is %f bit by default and can be changed with -ssgapd\n",par.ssgapd);
  printf(" -ssgapd       Controls additional penalty for opening gap within SS elements\n");
  printf(" -excl <range> exclude query positions from the alignment, e.g. '1-33,97-168'\n");
}

void help_other()
{
  printf("\n");
  printf("Other options: \n");
  printf(" -calm 0|1|2   use score calibration of  0:query  1:template  2:both (def=%i)\n",par.calm);
  printf(" -opt  <file>  parameter optimization mode (def=off): return sum of ranks \n");
  printf("               of true positives (same superfamily) for minimization      \n");
  printf("               and write result into file                                 \n");
  printf(" -scores <file> write scores for all pairwise comparisions to file         \n");
}

void help_all()
{
  help();
  help_out();
  help_hmm();
  help_gap();
  help_ali();
  help_other();
  printf("\n");
  printf("Default options can be specified in './.hhdefaults' or '~/.hhdefaults'\n");
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
      else if (!strcmp(argv[i],"-ofas"))
        {
          par.append=0;
          par.outformat=1;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.pairwisealisfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-oa2m"))
        {
          par.append=0;
          par.outformat=2;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.pairwisealisfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-oa3m"))
        {
          par.append=0;
          par.outformat=3;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.pairwisealisfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-Oa3m"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Oa3m\n"; exit(4);}
          else strcpy(par.alnfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-Ohhm"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Ohhm\n"; exit(4);}
          else strcpy(par.hhmfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-Opsi"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Opsi\n"; exit(4);}
          else strcpy(par.psifile,argv[i]);
        }
      else if (!strcmp(argv[i],"-w"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no weights file following -a\n"; exit(4);}
          else {strcpy(par.wfile,argv[i]); par.trans=1;}
        }
      else if (!strcmp(argv[i],"-w2"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no weights file following -a\n"; exit(4);}
          else {strcpy(par.wfile,argv[i]); par.trans=2;}
        }
      else if (!strcmp(argv[i],"-w3"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no weights file following -a\n"; exit(4);}
          else {strcpy(par.wfile,argv[i]); par.trans=3;}
        }
      else if (!strcmp(argv[i],"-w4"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no weights file following -a\n"; exit(4);}
          else {strcpy(par.wfile,argv[i]); par.trans=4;}
        }
      else if (!strcmp(argv[i],"-opt"))
        {
          par.opt=1;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file following -opt\n"; exit(4);}
          else strcpy(par.buffer,argv[i]);
        }
      else if (!strcmp(argv[i],"-scores"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file following -scores\n"; exit(4);}
          else {strcpy(par.scorefile,argv[i]);}
        }
      else if (!strncmp(argv[i],"-bla",4))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no file following -bla\n"; exit(4);}
          else {
              par.blafile = new(char[strlen(argv[i])+1]);
              strcpy(par.blafile,argv[i]);
          }
        }
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help"))
        {
          if (++i>=argc || argv[i][0]=='-') {help(); exit(0);}
          if (!strcmp(argv[i],"out")) {help_out(); exit(0);}
          if (!strcmp(argv[i],"hmm")) {help_hmm(); exit(0);}
          if (!strcmp(argv[i],"gap")) {help_gap(); exit(0);}
          if (!strcmp(argv[i],"ali")) {help_ali(); exit(0);}
          if (!strcmp(argv[i],"other")) {help_other(); exit(0);}
          if (!strcmp(argv[i],"all")) {help_all(); exit(0);}
          else {help(); exit(0);}
        }
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
      else if (!strcmp(argv[i],"-qscc") && (i<argc-1))    par.qsc_core=atof(argv[++i]);
      else if (!strcmp(argv[i],"-csc") && (i<argc-1))     par.coresc=atof(argv[++i]);
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
      else if (!strcmp(argv[i],"-neffa") && (i<argc-1)) par.neffa=atof(argv[++i]);
      else if (!strcmp(argv[i],"-neffb") && (i<argc-1)) par.neffb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-Etrans") && (i<argc-1)) par.Emax_trans=atof(argv[++i]);
      else if (!strcmp(argv[i],"-wtrans") && (i<argc-1)) par.wtrans=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssgap")) par.ssgap=1;
      else if (!strcmp(argv[i],"-ssgapd") && (i<argc-1)) par.ssgapd=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssgape") && (i<argc-1)) par.ssgape=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssgapi") && (i<argc-1)) par.ssgapi=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-ssm") && (i<argc-1)) par.ssm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-ssw") && (i<argc-1)) par.ssw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssa") && (i<argc-1)) par.ssa=atof(argv[++i]);
      else if (!strcmp(argv[i],"-realign")) par.realign=1;
      else if (!strcmp(argv[i],"-norealign")) par.realign=0;
      else if (!strcmp(argv[i],"-forward")) par.forward=1;
      else if (!strcmp(argv[i],"-mac") || !strcmp(argv[i],"-MAC")) par.forward=2;
      else if (!strcmp(argv[i],"-map") || !strcmp(argv[i],"-MAP")) par.forward=2;
      else if (!strcmp(argv[i],"-vit")) par.forward=0;
      else if (!strncmp(argv[i],"-glo",3)) {par.loc=0; if (par.mact>0.3 && par.mact<0.301) {par.mact=0;} }
      else if (!strncmp(argv[i],"-loc",4)) par.loc=1;
      else if (!strncmp(argv[i],"-alt",4) && (i<argc-1)) par.altali=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-r")) par.repmode=1;
      else if (!strcmp(argv[i],"-M") && (i<argc-1))
        if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1;
        else if(!strcmp(argv[i],"first"))  par.M=3;
        else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
        else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strcmp(argv[i],"-cal")) par.calibrate=1;
      else if (!strcmp(argv[i],"-calm") && (i<argc-1)) par.calm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-shift") && (i<argc-1)) par.shift=atof(argv[++i]);
      else if (!strcmp(argv[i],"-mact") && (i<argc-1)) par.mact=atof(argv[++i]);
      else if (!strcmp(argv[i],"-mapt") && (i<argc-1)) par.mact=atof(argv[++i]);
      else if (!strcmp(argv[i],"-sc") && (i<argc-1)) par.columnscore=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-def")) ;
      else if (!strncmp(argv[i],"-cpu",4) && (i<argc-1)) threads=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-corr") && (i<argc-1)) par.corr=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ovlp") && (i<argc-1)) par.min_overlap=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-dbstrlen") && (i<argc-1)) par.maxdbstrlen=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-tags")) par.notags=0;
      else if (!strcmp(argv[i],"-notags")) par.notags=1;
      else if (!strcmp(argv[i],"-mode") && (i<argc-1)) par.mode=atoi(argv[++i]);
      else if (!strncmp(argv[i],"-idummy",7) && (i<argc-1)) par.idummy=atoi(argv[++i]);
      else if (!strncmp(argv[i],"-jdummy",7) && (i<argc-1)) par.jdummy=atoi(argv[++i]);
      else if (!strncmp(argv[i],"-fdummy",7) && (i<argc-1)) par.fdummy=atof(argv[++i]);
      else if (!strncmp(argv[i],"-hhb_pval",8) && (i<argc-1)) par.hhblast_prefilter_pval=-log(atof(argv[++i]));
      else if (!strcmp(argv[i],"-csb") && (i<argc-1)) par.csb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-csw") && (i<argc-1)) par.csw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-cs"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -cs\n"; exit(4);}
          else strcpy(par.clusterfile,argv[i]);
        }
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
//////////////////////////////////////////////////////////////////////////////////////
void AlignByWorker(int bin)
{

  // Prepare q ant t and compare
  PrepareTemplate(q,*(t[bin]),format[bin]);

  // Do HMM-HMM comparison, store results if score>SMIN, and try next best alignment
  for (hit[bin]->irep=1; hit[bin]->irep<=par.altali; hit[bin]->irep++)
    {
      if (par.forward==0)
        {
          hit[bin]->Viterbi(q,*(t[bin]));
          if (hit[bin]->irep>1 && hit[bin]->score <= SMIN) break;
          hit[bin]->Backtrace(q,*(t[bin]));
        }
      else if (par.forward==1)
        {
          hit[bin]->Forward(q,*(t[bin]));
          hit[bin]->StochasticBacktrace(q,*(t[bin]),1); // the 1 selects maximization instead of stochastic backtracing
        }
      else if (par.forward==2)
        {
          hit[bin]->Forward(q,*(t[bin]));
          hit[bin]->Backward(q,*(t[bin]));
          hit[bin]->MACAlignment(q,*(t[bin]));
          hit[bin]->BacktraceMAC(q,*(t[bin]));
        }
      hit[bin]->score_sort = hit[bin]->score_aass;
      //printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit[bin]->name,hit[bin]->fam,hit[bin]->irep,hit[bin]->score);

#ifdef PTHREAD
      pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif
      hitlist.Push(*(hit[bin]));            // insert hit at beginning of list (last repeats first!)
#ifdef PTHREAD
      pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif

      if (par.forward>0) break; // find only best alignment for forward algorithm and stochastic sampling
      if (hit[bin]->score <= SMIN) break;  // break if score for first hit is already worse than SMIN
    }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////
//// Realign q and *(t[bin]) with the MAC algorithm (after the database search)
//////////////////////////////////////////////////////////////////////////////////////
void RealignByWorker(int bin)
{
  Hit hit_cur;
  int nhits=0;
  int pos;
  hit[bin]->irep=1;

  // Prepare MAC comparison(s)
  PrepareTemplate(q,*(t[bin]),format[bin]);
  t[bin]->Log2LinTransitionProbs(1.0);

#ifdef PTHREAD
          pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif

  // Search positions in hitlist with correct index of template
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      if (nhits>=imax(par.B,par.Z)) break;
      //    fprintf(stderr,"t->name=%s  hit_cur.name=%s  hit[bin]->irep=%i  nhits=%i  hit_cur.index=%i  hit[bin]->index=%i\n",t[bin]->name,hit_cur.name,hit_cur.irep,nhits,hit_cur.index,hit[bin]->index);
      if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
      if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
      if (hit_cur.index==hit[bin]->index) // found position with correct template
        {
          //      fprintf(stderr,"  t->name=%s   hit_cur.irep=%i  hit[bin]->irep=%i  nhits=%i\n",t[bin]->name,hit_cur.irep,hit[bin]->irep,nhits);
          pos = hitlist.GetPos();
#ifdef PTHREAD
          pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif
          // Align q to template in *hit[bin]
          hit[bin]->Forward(q,*(t[bin]));
          hit[bin]->Backward(q,*(t[bin]));
          hit[bin]->MACAlignment(q,*(t[bin]));
          hit[bin]->BacktraceMAC(q,*(t[bin]));
#ifdef PTHREAD
          pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif
          hit_cur = hitlist.Read(pos);

          // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
          hit[bin]->score      = hit_cur.score;
          hit[bin]->score_aass = hit_cur.score_aass;
          hit[bin]->score_ss   = hit_cur.score_ss; // comment out?? => Andrea
          hit[bin]->Pval       = hit_cur.Pval;
          hit[bin]->Pvalt      = hit_cur.Pvalt;
          hit[bin]->logPval    = hit_cur.logPval;
          hit[bin]->logPvalt   = hit_cur.logPvalt;
          hit[bin]->logP1val   = hit_cur.logP1val;
          hit[bin]->Eval       = hit_cur.Eval;
          hit[bin]->logEval    = hit_cur.logEval;
          hit[bin]->E1val      = hit_cur.E1val;
          hit[bin]->Probab     = hit_cur.Probab;

          // Replace original hit in hitlist with realigned hit
          //hitlist.ReadCurrent().Delete();
          hitlist.Delete().Delete();                // delete list record and hit object
          hitlist.Insert(*hit[bin]);
          hit[bin]->irep++;
        }
      nhits++;
    }
#ifdef PTHREAD
  pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif


  if (hit[bin]->irep==1)
    {
      fprintf(stderr,"*************************************************\n");
      fprintf(stderr,"\nError: could not find template %s in hit list (index:%i dbfile:%s ftell:%i\n\n",hit[bin]->name, hit[bin]->index,hit[bin]->dbfile,(unsigned int)hit[bin]->ftellpos);
      fprintf(stderr,"*************************************************\n");
    }

    return;
}



#ifdef PTHREAD

///////////////////////////////////////////////////////////////////////////////////////
//// This is the main thread loop that waits for new jobs (i.e. pairwise alignment) and executes them
//////////////////////////////////////////////////////////////////////////////////////
void* WorkerLoop(void* data)
{
  int thread_id = (*((Thread_args*)data)).thread_id;  // typecast 'data' from pointer-to-void to pointer-to-Thread_args
  void (*ExecuteJob)(int) = (*((Thread_args*)data)).function; // dito; data.function is a pointer to a function(int) that returns void
  int rc;                            // return code for threading commands
  int bin;                           // bin index

  // Lock access to bin_status
  rc = pthread_mutex_lock(&bin_status_mutex);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Take jobs from one of the SUBMITTED bins and execute them. If no submitted jobs found, wait for signal 'new_job'
  while (reading_dbs || jobs_submitted>0)
    {

     if (jobs_submitted>0 ) // submitted job in one of the bins?
        {
          bin = PickBin(SUBMITTED);
          if (DEBUG_THREADS)
            fprintf(stderr,"Thread %3i:   start job in bin %i       jobs running: %i  jobs_submitted:%i \n",thread_id,bin,jobs_running,jobs_submitted);
          jobs_running++;            // the order of the following three lines is important, since...
          jobs_submitted--;          // ... the main thread tries to find a free bin...
          bin_status[bin] = RUNNING; // ... if jobs_running+jobs_submitted<bins !

          // Execute job
          rc = pthread_mutex_unlock(&bin_status_mutex); // unlock access to bin_status
          if (DEBUG_THREADS)
            fprintf(stderr,"Thread %3i:   executing HMM %-10.10s jobs running: %i  jobs_submitted:%i \n",thread_id,t[bin]->name,jobs_running,jobs_submitted);
          ExecuteJob(bin);
          rc = pthread_mutex_lock(&bin_status_mutex); // lock access to bin_status

          bin_status[bin] = FREE;
          jobs_running--;

          // Signal completion of job
          rc = pthread_cond_signal(&finished_job);

          if (DEBUG_THREADS)
            fprintf(stderr,"Thread %3i:   finished job in bin %i    jobs running: %i  jobs_submitted:%i \n",thread_id,bin,jobs_running,jobs_submitted);
        }
     else
       {

         if (DEBUG_THREADS)
           fprintf(stderr,"Thread %3i:   waiting for new job ...   jobs running: %i  jobs_submitted:%i \n",thread_id,jobs_running,jobs_submitted);
         rc = pthread_cond_wait(&new_job, &bin_status_mutex);
       }

      // Unlock access to bin_status
      rc = pthread_mutex_unlock(&bin_status_mutex);
      // Lock access to bin_status
      rc = pthread_mutex_lock(&bin_status_mutex);
    }
  ///////////////////////////////////////////////////////////////////////////////////////

  // Unlock access to bin_status
  rc = pthread_mutex_unlock(&bin_status_mutex);

  // Exit thread automatically
  pthread_exit(NULL);
  return NULL;
}

#endif

/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PR_OGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  char line[LINELEN]="";         // input line
  char* argv_conf[MAXOPT];       // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf;                 // Number of arguments in argv_conf
  char inext[IDLEN]="";          // Extension of input file (hhm or a3m)
  int bin;                       // bin index

#ifdef PTHREAD
  struct Thread_args thread_data[MAXTHREADS]; // store a threads thread_id and function to call (AlignByWorker, RealignByWorker)
  pthread_t pthread[MAXTHREADS]; // info on thread's structures (needed by system)
  pthread_attr_t joinable;       // attribute set for describing threads
  int rc;                        // return code for threading commands
  pthread_attr_init(&joinable);  // initialize attribute set with default values
  if (pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)!=0) // set attribute 'joinable'
    cerr<<"Error "<<pthread_attr_setdetachstate(&joinable, PTHREAD_CREATE_JOINABLE)<<": could not set detach state for thread attibute.\n";
#endif

  const char print_elapsed=0;

  SetDefaults();
  par.neffa=1.0;
  par.neffb=10.0;

  // Make command line input globally available
  par.argv=argv;
  par.argc=argc;
  RemovePathAndExtension(program_name,argv[0]);

  // Enable changing verbose mode before defaults file and command line are processed
  for (int i=1; i<argc; i++)
    {
      if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1;
      else if (argc>1 && !strcmp(argv[i],"-v0")) v=0;
      else if (argc>1 && !strcmp(argv[i],"-v1")) v=1;
      else if (argc>2 && !strcmp(argv[i],"-v")) v=atoi(argv[i+1]);
    }

  // Read .hhdefaults file?
  if (par.readdefaultsfile)
    {
      // Process default otpions from .hhconfig file
      ReadDefaultsFile(argc_conf,argv_conf);
      ProcessArguments(argc_conf,argv_conf);
    }

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc,argv);

  // Check command line input and default values
  if (!*par.infile) // string empty?
    {help(); cerr<<endl<<"Error in "<<program_name<<": no query alignment file given (-i file)\n"; exit(4);}
  if (!par.dbfiles) // pointer never set?
    {help(); cerr<<endl<<"Error in "<<program_name<<": no HMM database file given (-d file)\n"; exit(4);}
  if (threads<=1) threads=0;
  else if (threads>MAXTHREADS)
    {
      threads=MAXTHREADS;
      if (v>=1) fprintf(stderr,"WARNING: number of CPUs set to maximum value of %i\n",MAXTHREADS);
    }
  RemoveExtension(q.file,par.infile); // get rootname of infile (no directory path, no extension)
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
  if (par.pca<0.001) par.pca=0.001; // to avoid log(0)
  if (par.b>par.B) par.B=par.b;
  if (par.z>par.Z) par.Z=par.z;

  // Input parameters
  if (v>=3)
    {
      cout<<"Input file :   "<<par.infile<<"\n";
      cout<<"Database file: "<<par.dbfiles<<"\n";
      cout<<"Output file:   "<<par.outfile<<"\n";
    }

  // Set secondary structure substitution matrix
  if (par.ssm) SetSecStrucSubstitutionMatrix();

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix();

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  ReadAndPrepare(par.infile, q);

//  // Rescale matrix according to query aa composition? (Two iterations are sufficient)
//   if (par.pcm==4)
//     {
//       q.RescaleMatrix();
//       q.PreparePseudocounts();
//       q.AddAminoAcidPseudocounts();
//       SetSubstitutionMatrix();
//       q.RescaleMatrix();
//       q.PreparePseudocounts();
//       q.AddAminoAcidPseudocounts();
//     }

  // Reset lamda?
  if (par.calibrate>0 || par.trans>0) {q.lamda=LAMDA; q.mu=0.0;}

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags) q.NeutralizeTags();

  if (par.forward>=1)
    {
      if (v>=2 && par.forward==1) printf("Using Forward algorithm ...\n");
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
      hit[bin]->AllocateBacktraceMatrix(q.L+2,MAXRES); // ...with a separate dynamic programming matrix (memory!!)
      if (par.forward>=1)
        hit[bin]->AllocateForwardMatrix(q.L+2,MAXRES);
      if (par.forward==2)
        hit[bin]->AllocateBackwardMatrix(q.L+2,MAXRES);

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
  const int MAXNUMDB=16383;
  int ndb=0;
  Hash<char>* doubled;
  doubled = new(Hash<char>);
  doubled->New(16381,0);
  char* dbfiles[MAXNUMDB+1];
  char* dbfile_cur=strscn_(par.dbfiles); // current file name
  char* dbfile_next; // next file name after current
  while (*dbfile_cur && ndb<MAXNUMDB)
    {
      // Cut off everything after next white space and store beginning of next database name in dbfile_next
      dbfile_next=strcut_(dbfile_cur);

      // Has the dbfiles[ndb] name not been seen before?
      if (! doubled->Contains(dbfile_cur))
        {
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
              while(fgetline(dbfile,NAMELEN,palf) && ndb<MAXNUMDB) // read HMM files in pal file
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
                  else
                    printf("WARNING: skipping doubled datbase file %s\n",dbfile);
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
      else
        printf("WARNING: skipping doubled datbase file %s\n",dbfile_cur);

     dbfile_cur=dbfile_next;
    }

  //fprintf (stderr,"ndb: %i  MAXNUMDB: %i\n",ndb,MAXNUMDB);

  if (v>=1 && ndb>=MAXNUMDB && *dbfiles[ndb-1])
    fprintf (stderr,"WARNING: maximum of %i allowed databases surpassed. Skipping the rest.\n",MAXNUMDB);



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Search databases

  // Initialize
  N_searched=0;
  int v1=v;
  if (v<=3) v=1; else v-=2;
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
              if (!strncmp(line,"HMMER",5))      // read HMMER format
                {
                  format[bin] = 1;
                  read_from_db = t[bin]->ReadHMMer(dbf,dbfiles[idb]);
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
                  tali.N_filtered = tali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
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
              if (read_from_db==2) continue;  // skip current HMM or reached end of database
              if (read_from_db==0) break;     // finished reading HMMs
              if (v>=4) printf("Aligning with %s\n",t[bin]->name);  /////////////////////v>=4
              ///////////////////////////////////////////////////

              hit[bin]->dbfile = new(char[strlen(dbfiles[idb])+1]);
              strcpy(hit[bin]->dbfile,dbfiles[idb]); // record db file name from which next HMM is read

              ++N_searched;
              if (v>=1 && !(N_searched%20))
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

  if (v1>=1) cout<<"\n";
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
  else if ((par.calm!=1 && q.lamda==0) || par.calibrate>0 || par.trans>0)
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
                fprintf(stderr,"\nWARNING: E-values for global alignment option may be unreliable.\n");
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // HHblast: read BLAST file with P-values?
//   if (par.blafile)
//     {
//       hitlist.ReadBlastFile(q);

//       // Calculate composite P-values/E-values for HHblast from HHsearch and PSI-BLAST P-values
//       hitlist.CalculateHHblastEvalues(q);
//     }

  // Calculate corrected E-values for HHblast by NN prediction of correlation of PSI-BLAST and HHsearch P-values
  if (par.hhblast_prefilter_pval)
    {
      hitlist.CalculateHHblastEvalues(q);
    }

  // Optimization mode?
  if (par.opt) {hitlist.Optimize(q,par.buffer);}


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Realign all hits with MAC algorithm?
  if (par.realign && par.forward!=2)
    {

      q.Log2LinTransitionProbs(1.0); // transform transition freqs to lin space if not already done

      Hash< List<Posindex>* >* realign; // realign->Show(dbfile) is list with ftell positions for templates in dbfile to be realigned
      realign = new(Hash< List<Posindex>* >);
      realign->New(3601,NULL);
      Hit hit_cur;
      const float MEMSPACE_DYNPROG = 512*1024*1024;
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
//            fprintf(stderr,"hit.name=%-15.15s  hit.index=%-5i hit.ftellpos=%-8i  hit.dbfile=%s\n",hit_cur.name,hit_cur.index,(unsigned int)hit_cur.ftellpos,hit_cur.dbfile);
              if (nhits>=par.jdummy || hit_cur.irep>1) // realign the first jdummy hits consecutively to query profile
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
          if (v>=1) {
            cerr<<"WARNING: Realigning sequences only up to length "<<Lmaxmem<<" due to limited memory."<<endl;
            if (bins>1) cerr<<"Note: you can reduce memory requirement by lowering N in the -cpu N option."<<endl;
          }
        }

      // Initialize and allocate space for dynamic programming
      jobs_running = 0;
      jobs_submitted = 0;
      reading_dbs=1;   // needs to be set to 1 before threads are created
      for (bin=0; bin<bins; bin++)
        {
          if (par.forward==0)
            hit[bin]->AllocateForwardMatrix(q.L+2,Lmax+1);
          if (par.forward<=1)
            hit[bin]->AllocateBackwardMatrix(q.L+2,Lmax+1);
          bin_status[bin] = FREE;
        }


      if (v>=1) printf("Realigning %i query-template alignments with maximum accuracy (MAC) algorithm ...\n",nhits);
      if (v<=3) v=1; else v-=1;  // Supress verbose output during iterative realignment and realignment

      // Align the first par.jdummy templates?
      if (par.jdummy>0)
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
          Qali.longname = new(char[strlen(q.longname)+1]);
          strcpy(Qali.longname,q.longname);
          strcpy(Qali.name,q.name);
          strcpy(Qali.fam,q.fam);
          RemovePathAndExtension(Qali.file,par.hhmfile);

          if (v>=2) printf("Merging best hits to query alignment %s ...\n",qa3mfile);

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
	      if (hit_cur.L>Lmaxmem) continue;  //Don't align to long sequences due to memory limit

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
                  tali.N_filtered = tali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
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
              nhits++;
            }
        }

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
                      tali.N_filtered = tali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
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
      // End Realign all hits with MAC algorithm?
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    } else {

    // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
    if (*par.scorefile) {
      if (v>=3) printf("Printing scores file ...\n");
      hitlist.PrintScoreFile(q);
    }
  }



  // Print for each HMM: n  score  -log2(Pval)  L  name  (n=5:same name 4:same fam 3:same sf...)
  if (*par.scorefile) {
    if (v>=3) printf("Printing scores file ...\n");
    hitlist.PrintScoreFile(q);
  }

  // Print FASTA or A2M alignments?
  if (*par.pairwisealisfile) {
    if (v>=2) cout<<"Printing alignments in "<<(par.outformat==1? "FASTA" : par.outformat==2?"A2M" :"A3M")<<" format to "<<par.pairwisealisfile<<"\n";
    hitlist.PrintAlignments(q,par.pairwisealisfile,par.outformat);
  }

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
    q.InsertCalibration(par.infile);


  // Generate output alignment or HMM file?
  if (*par.alnfile || *par.psifile || *par.hhmfile)
    {
      Alignment Qali;  // output A3M generated by merging A3M alignments for significant hits to the query alignment
      HMM Q;           // output HMM: generated from Qali
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
      if (v>=2) printf("Merging hits to query alignment %s ...\n",qa3mfile);
      // If query consists of only one sequence:
      //     realign templates to query HMM enforced by sequences from up to the 10 best templates
      int v1=v--;

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
          Qali.MergeMasterSlave(hit,ta3mfile);
          nhits++;
        }

      // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
      Qali.Compress("merged A3M file");

      // Sort out the nseqdis most dissimilar sequences for display in the output alignments
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
          strcpy(Qali.longname,q.longname);
          strcpy(Qali.name,q.name);
          strcpy(Qali.fam,q.fam);
          RemovePathAndExtension(Qali.file,par.hhmfile);

          // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
          Qali.FrequenciesAndTransitions(Q);

          // Add *no* amino acid pseudocounts to query. This is necessary to copy f[i][a] to p[i][a]
          Q.AddAminoAcidPseudocounts(0, 0.0, 0.0, 1.0);
          Q.CalculateAminoAcidBackground();

          // Write HMM to output file in HHsearch format?
          Q.WriteToFile(par.hhmfile);
        }

      // Write output A3M alignment?
      if (*par.alnfile) Qali.WriteToFile(par.alnfile,"a3m");

      // Write output PSI-BLAST-formatted alignment?
      if (*par.psifile) Qali.WriteToFile(par.psifile,"psi");
   }

  // Delete memory for dynamic programming matrix
  for (bin=0; bin<bins; bin++)
    {
      hit[bin]->DeleteBacktraceMatrix(q.L+2);
      if (par.forward>=1 || par.realign)
        hit[bin]->DeleteForwardMatrix(q.L+2);
      if (par.forward==2 || par.realign)
        hit[bin]->DeleteBackwardMatrix(q.L+2);
      delete hit[bin];
      delete t[bin];
     }
  if (par.dbfiles) delete[] par.dbfiles;
  for (int idb=0; idb<ndb; idb++) delete[](dbfiles[idb]);
  if (format) delete[](format);
  if (par.blafile) delete[] par.blafile;
  if (par.exclstr) delete[] par.exclstr;
  delete doubled;

  // Delete content of hits in hitlist
  Hit hit_cur;
  hitlist.Reset();
  while (!hitlist.End())
    {
      //hitlist.ReadCurrent().Delete();
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

  exit(0);
} //end main

//////////////////////////////////////////////////////////////////////////////////////////////////////
// END OF MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////



