/*
 * hhdecl.h
 *
 *  Created on: Mar 28, 2014
 *      Author: meiermark
 */

#ifndef HHDECL_H_
#define HHDECL_H_

#include <unistd.h>

#include "crf.h"
#include "aa.h"
#include "crf_pseudocounts-inl.h"
#include "library_pseudocounts-inl.h"
#include "log.h"
#include "util.h"

/////////////////////////////////////////////////////////////////////////////////////
//// Global variable declarations
/////////////////////////////////////////////////////////////////////////////////////

const char REFERENCE[] = "Steinegger M, Meier M, Mirdita M, Vöhringer H, Haunsberger S J, and Söding J (2019)\nHH-suite3 for fast remote homology detection and deep protein annotation.\nBMC Bioinformatics, doi:10.1186/s12859-019-3019-7\n";
const char COPYRIGHT[] = "(c) The HH-suite development team\n";

const int LINELEN=524288; //max length of line read in from input files; must be >= MAXCOL
const int MAXSEQDIS=10238;//max number of sequences stored in 'hit' objects and displayed in output alignment
const int IDLEN=255;     //max length of scop hierarchy id and pdb-id
const int DESCLEN=32765;//max length of sequence description (longname)
const int NAMELEN=(PATH_MAX>512? PATH_MAX:512); //max length of file names etc., defined in limits.h
const int NAA=20;       //number of amino acids (0-19)
const int NTRANS=7;     //number of transitions recorded in HMM (M2M,M2I,M2D,I2M,I2I,D2M,D2D)
const int NCOLMIN=10;   //min number of cols in subalignment for calculating pos-specific weights w[k][i]
const int ANY=20;       //number representing an X (any amino acid) internally
const int GAP=21;       //number representing a gap internally
const int FWD_BKW_PATHWITDH=40;       //cell off path width around viterbi alignment
const int ENDGAP=22;    //Important to distinguish because end gaps do not contribute to tansition counts
const int HMMSCALE=1000;//Scaling number for log2-values in HMMs
const int MAXPROF=32766;//Maximum number of HMM scores for fitting EVD
const float MAXENDGAPFRAC=0.1; //For weighting: include only columns into subalignment i that have a max fraction of seqs with endgap
const float LAMDA=0.388; //lamda in score EVD used for -local mode in length correction: S = S-log(Lq*Lt)/LAMDA)
const float LAMDA_GLOB=0.42; //lamda in score EVD used for -global mode
const int SELFEXCL=3;   // exclude self-alignments with j-i<SELFEXCL
const float PLTY_GAPOPEN=6.0f; // for -qsc option (filter for min similarity to query): 6 bits to open gap
const float PLTY_GAPEXTD=1.0f; // for -qsc option (filter for min similarity to query): 1 bit to extend gap
const int MINCOLS_REALIGN=6; // hits with MAC alignments with fewer matched columns will be deleted in hhsearch hitlist; must be at least 2 to avoid nonsense MAC alignments starting from the left/upper edge
const float LOG1000=log(1000.0);
const float POSTERIOR_PROBABILITY_THRESHOLD = 0.01;
const int VITERBI_PATH_WIDTH=40;

// Secondary structure
const int NDSSP=8;      //number of different ss states determined by dssp: 0-7 (0: no state available)
const int NSSPRED=4;    //number of different ss states predicted by psipred: 0-3 (0: no prediction availabe)
const int MAXCF=11;     //number of different confidence values: 0-10 (0: no prediction availabe)

// const char aa[]="ARNDCQEGHILKMFPSTWYVX-";
//Amino acids Sorted by alphabet     -> internal numbers a
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  X
const int s2a[]={ 0, 4, 3, 6,13, 7, 8, 9,11,10,12, 2,14, 5, 1,15,16,19,17,18,20};

//Internal numbers a for amino acids -> amino acids Sorted by alphabet:
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
const int a2s[]={ 0,14,11, 2, 1,13, 3, 5, 6, 7, 9, 8,10, 4,12,15,16,18,19,17,20};

enum transitions {M2M,M2I,M2D,I2M,I2I,D2M,D2D}; // index for transitions within a HMM

#ifdef __GNUC__
#define __DEPRECATED__ __attribute__((deprecated))
#elif defined(_MSC_VER)
#define __DEPRECATED__ __declspec(deprecated)
#else
#pragma message("WARNING: No DEPRECATED for this compiler")
#define __DEPRECATED__
#endif


enum pair_states {DEPRECATED_STOP=0,DEPRECATED_MM=2,DEPRECATED_GD=3,DEPRECATED_IM=4,DEPRECATED_DG=5,DEPRECATED_MI=6};
//__DEPRECATED__
const pair_states STOP = pair_states(DEPRECATED_STOP);
//__DEPRECATED__
const pair_states MM = pair_states(DEPRECATED_MM);
//__DEPRECATED__
const pair_states GD = pair_states(DEPRECATED_GD);
//__DEPRECATED__
const pair_states IM = pair_states(DEPRECATED_IM);
//__DEPRECATED__
const pair_states DG = pair_states(DEPRECATED_DG);
//__DEPRECATED__
const pair_states MI = pair_states(DEPRECATED_MI);


//states for the interim filter of the query msa during merging
enum InterimFilterStates {INTERIM_FILTER_NONE=0, INTERIM_FILTER_FULL=1};

// Pseudocounts
namespace Pseudocounts {
enum Admix {
  ConstantAdmix = 1,
  HHsearchAdmix = 2,
  CSBlastAdmix  = 3
};

struct Params {
  Params(
      Admix m     = ConstantAdmix,
      double a    = 1.0,
      double b    = 1.0,
      double c    = 1.0,
      double neff = 0.0)
    : admix(m), pca(a), pcb(b), pcc(c), target_neff(neff) {}

  cs::Admix* CreateAdmix() {
    switch (admix) {
      case ConstantAdmix:
        return new cs::ConstantAdmix(pca);
        break;
      case HHsearchAdmix:
        return new cs::HHsearchAdmix(pca, pcb, pcc);
        break;
      case CSBlastAdmix:
        return new cs::CSBlastAdmix(pca, pcb);
        break;
      default:
        return NULL;
    }
  }

  Admix admix;        // admixture mode
  double pca;         // admixture paramter a
  double pcb;         // admixture paramter b
  double pcc;         // admixture parameter c needed for HHsearchAdmix
  double target_neff; // target diversity adjusted by optimizing a
};

};

class Parameters          // Parameters for gap penalties and pseudocounts
{
public:
  Parameters(const int argc, const char** argv);

  const char** argv;            //command line parameters
  const char argc;              //dimension of argv

  LogLevel v;

  char infile[NAMELEN];   // input filename
  char outfile[NAMELEN];  // output filename
  char matrices_output_file[NAMELEN];
  bool filter_matrices;
  char pairwisealisfile[NAMELEN]; // output filename with pairwise alignments
  char alisbasename[NAMELEN];
  char alnfile[NAMELEN];  // name of output alignment file in A3M format (for iterative search)
  char hhmfile[NAMELEN];  // name of output HHM file for (iterative search)
  char psifile[NAMELEN];  // name of output alignmen file in PSI-BLAST format (iterative search)
  char scorefile[NAMELEN];// table of scores etc for all HMMs in searched database
  char m8file[NAMELEN];   // blast tab format for all HMMs in searched database
  char indexfile[NAMELEN];// optional file containing indeices of aligned residues in given alignment
  std::vector<std::string> tfiles;    // template filenames (in hhalign)
  char alitabfile[NAMELEN]; // where to write pairs of aligned residues (-atab option)
  char* exclstr;          // optional string containing list of excluded residues, e.g. '1-33,97-168'
  char* template_exclstr;
  int aliwidth;           // number of characters per line in output alignments for HMM search
  float p;                // minimum probability for inclusion in hit list and alignments
  double E;               // maximum E-value for inclusion in hit list and alignment list
  double e;               // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
  int Z;                  // max number of lines in hit list
  int z;                  // min number of lines in hit list
  int B;                  // max number of lines in alignment list
  int b;                  // min number of lines in alignment list
  char showcons;           // in query-template alignments  0: don't show consensus sequence   1:show
  char showdssp;           // in query-template alignments  0: don't show ss_dssp lines        1:show
  char showpred;           // in query-template alignments  0: don't show ss_pred and ss_conf lines  1:show
  char showconf;           // in query-template alignments  0: don't show ss_conf lines        1:show
  char cons;              // if set to 1, include consensus as first representative sequence of HMM
  int nseqdis;            // maximum number of query or template sequences in output alignments
  char mark;              // which sequences to mark for display in output alignments? 0: auto; 1:all
  char append;            // append to output file? (hhmake)
  char outformat;         // 0: hhr  1: FASTA  2:A2M   3:A3M
                          //0:MAC alignment, master-slave  1:MAC blending, master-slave  2:MAC alignment, combining


  int max_seqid;          // Maximum sequence identity with all other sequences in alignment
  int qid;                // Minimum sequence identity with query sequence (sequence 0)
  float qsc;              // Minimum score per column with query sequence (sequence 0)
  int coverage;           // Minimum coverage threshold
  int Ndiff;              // Pick Ndiff most different sequences that passed the other filter thresholds
  bool allseqs;           // if true, do not filter in output alignment; show all sequences

  int Mgaps;              // Maximum percentage of gaps for match states
  int M;                  // Match state assignment by  1:upper/lower case  2:percentage rule  3:marked sequence
  int M_template;
  char matrix;            // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50

  char wg;                // 0: use local sequence weights   1: use global ones

  Pseudocounts::Params pc_hhm_context_engine;       // Pseudocounts parameters for query hhm if context given
  Pseudocounts::Params pc_prefilter_context_engine; // Pseudocounts parameters for prefiltering if context given

  //pseudocount variables if no context is used
  int pc_hhm_nocontext_mode;              // Admixture method
  float pc_hhm_nocontext_a;               // Admixture parameter a
  float pc_hhm_nocontext_b;               // Admixture parameter b
  float pc_hhm_nocontext_c;               // Admixture parameter c

  //pseudocount variables for the prefilter if no context is used
  int pc_prefilter_nocontext_mode;           // Admixture method
  float pc_prefilter_nocontext_a;            // Admixture parameter a
  float pc_prefilter_nocontext_b;            // Admixture parameter b
  float pc_prefilter_nocontext_c;            // Admixture parameter c

  float gapb;             // Diversity threshold for adding pseudocounts to transitions from M state
  float gapd;             // Gap open penalty factor for deletions
  float gape;             // Gap extend penalty: factor to multiply hmmer values (def=1)
  float gapf;             // factor for increasing/reducing the gap opening penalty for deletes
  float gapg;             // factor for increasing/reducing the gap opening penalty for inserts
  float gaph;             // factor for increasing/reducing the gap extension penalty for deletes
  float gapi;             // factor for increasing/reducing the gap extension penalty for inserts

  float egq;              // penalty for end gaps when query not fully covered
  float egt;              // penalty for end gaps when template not fully covered

  float Neff;

  char ssm;               // SS comparison mode: 0:no ss scoring  1:ss scoring AFTER alignment  2:ss score in column score
  float ssw;              // SS weight as compared to column score
  float ssw_realign;      // SS weight as compared to column score for realign
  float ssa;              // SS state evolution matrix M1 = (1-ssa)*I + ssa*M0

  char loc;               // 0: local alignment (wrt. query), 1: global alignement
  char realign;           // realign database hits to be displayed with MAC algorithm
  int altali;             // find up to this many possibly overlapping alignments
  float smin;             //Minimum score of hit needed to search for another repeat of same profile: p=exp(-(4-mu)/lamda)=0.01
  int columnscore;        // 0: no aa comp corr  1: 1/2(qav+tav) 2: template av freqs 3: query av freqs 4:...
  int half_window_size_local_aa_bg_freqs; // half-window size to average local aa background frequencies
  float corr;             // Weight of correlations between scores with |i-j|<=4
  float shift;            // Score offset for match-match states
  double mact;            // Probability threshold (negative offset) in MAC alignment determining greediness at ends of alignment
  float mac_min_length;   // Minimum length of MAC hit in comparison with the original Viterbi hit (in terms of matched_cols)
  int realign_max;        // Realign max ... hits
  float maxmem;           // maximum available memory in GB for realignment (approximately)

  int min_overlap;        // all cells of dyn. programming matrix with L_T-j+i or L_Q-i+j < min_overlap will be ignored
  char notags;            // neutralize His-tags, FLAG tags, C-myc tags?

  //TODO: a const would be nicer
  unsigned int maxdbstrlen; // maximum length of database string to be printed in 'Command' line of hhr file

  int maxcol;             // max number of columns in sequence/MSA input files; must be <= LINELEN and >= maxres
  int maxres;             // max number of states in HMM; must be <= LINELEN
  int maxseq;             // max number of sequences in MSA
  int maxnumdb;           // max number of hits allowed past prefilter

  bool hmmer_used;        // True, if a HMMER database is used

  // parameters for context-specific pseudocounts
  float csb;
  float csw;
  std::string clusterfile;
  bool nocontxt;

  // HHblits
  int dbsize;           // number of clusters of input database

  // HHblits Evalue calculation  (alpha = a + b(Neff(T) - 1)(1 - c(Neff(Q) - 1)) )
  float alphaa;
  float alphab;
  float alphac;

  // For filtering database alignments in HHsearch and HHblits
  // JS: What are these used for? They are set to the values without _db anyway.
  int max_seqid_db;
  int qid_db;
  float qsc_db;
  int coverage_db;
  int Ndiff_db;

  // HHblits context state prefilter
  std::string cs_library;

  // HHblits prefilter
  bool prefilter;             // perform prefiltering in HHblits?

  //early stopping stuff
  bool early_stopping_filter; // Break HMM search, when the sum of the last N HMM-hit-Evalues is below threshold
  double filter_thresh;    // Threshold for early stopping

  // For HHblits prefiltering with SSE2
  short prefilter_gap_open;
  short prefilter_gap_extend;
  int prefilter_score_offset;
  int prefilter_bit_factor;
  double prefilter_evalue_thresh;
  double prefilter_evalue_coarse_thresh;
  int preprefilter_smax_thresh;

  int min_prefilter_hits;

  size_t max_number_matrices;

  //hhblits specific variables
  int num_rounds;
  std::vector<std::string> db_bases;
  // Perform filtering of already seen HHMs
  bool already_seen_filter;
  // Realign old hits in last round or use previous alignments
  bool realign_old_hits;
  float neffmax;
  int threads;

  InterimFilterStates interim_filter;
};


#endif /* HHDECL_H_ */
