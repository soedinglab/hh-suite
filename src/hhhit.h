// hhhit.h
#ifndef HHHIT_H_
#define HHHIT_H_

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
#include <vector>

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

#include "list.h"
#include "hhhmm.h"
#include "hhdecl.h"
#include "hhutil.h"
#include "util.h"
#include "HHDatabase.h"

#include "hhhit-inl.h"

/////////////////////////////////////////////////////////////////////////////////////
// // Describes an alignment of two profiles. Used as list element in Hits : List<Hit> 
/////////////////////////////////////////////////////////////////////////////////////
class Hit
{
 public:  
  char* longname;       // Name of HMM
  char* name;           // One-word name of HMM
  char* file;           // Basename (with path, without extension) of alignment file that was used to construct the HMM
                        // (path from db-file is prepended to FILE record in HMM file!)
  char fam[IDLEN];      // family ID (derived from name) (FAM field)
  char sfam[IDLEN];     // superfamily ID (derived from name) 
  char fold[IDLEN];     // fold ID (derived from name)
  char cl[IDLEN];       // class ID (derived from name)

  char* dbfile;         // full database file name from which HMM was read
  long ftellpos;        // start position of HMM in database file
  int index;            // index of HMM in order of reading in (first=0)
  HHDatabaseEntry* entry;
  List<void*>* plist_phits; // points to a list of pointers to hitlist elements of same template (for realignment)
  
  float score;          // Score of alignment (i.e. of Viterbi path)
  float score_sort;     // score to sort hits in output list (negative means first/best!)
  float score_aass;     // first: just hit.score, then hit.logPval-SSSCORE2NATLOG*hit.score_ss;(negative means best!)
  float score_ss;       // Part of score due to secondary structure
  double Pval;           // P-value for whole protein based on score distribution of query
  double Pvalt;          // P-value for whole protein based on score distribution of template
  double logPval;        // natural logarithm of Pval
  double logPvalt;       // natural logarithm of Pvalt
  double Eval;           // E-value for whole protein
  double logEval;        // natural logarithm of Eval
  float Probab;         // probability in % for a positive (depends only on score)
  double Pforward;      // scaled total forward probability : Pforward * Product_{i=1}^{Lq+1}(scale[i])
  
  int L;                // Number of match states in template
  short int irep;       // Index  of single-repeat hit (1: highest scoring repeat hit)
  char lastrep;         // is current hit last (sub)optimal alignment? 0: no  1: yes
  
  int n_display;        // number of sequences stored for display of alignment 
  char** sname;         // names of stored sequences 
  char** seq;           // residues of stored sequences (first at pos 1)
  int nss_dssp;         // index of dssp secondary structure sequence in seq[]
  int nsa_dssp;         // index of of dssp solvent accessibility in seq[]
  int nss_pred;         // index of dssp secondary structure sequence in seq[]
  int nss_conf;         // index of dssp secondary structure sequence in seq[]
  int nfirst;           // index of query sequence in seq[]
  int ncons;            // index of consensus sequence
  
  float** P_MM;        // Posterior probability matrix, filled in Forward and Backward algorithms
  int nsteps;           // index for last step in Viterbi path; (first=1)
  int* i;               // i[step] = query match state at step of Viterbi path
  int* j;               // j[step] = template match state at step of Viterbi path
  std::vector<int>* alt_i; // Path of alternative alignments (query positions)
  std::vector<int>* alt_j; // Path of alternative alignments (template positions)
  char* states;         // state at step of Viterbi path  0: Start  1: M(MM)  2: A(-D)  3: B(IM)  4: C(D-)  5 D(MI)
  float* S;             // S[step] = match-match score contribution at alignment step
  float* S_ss;          // S_ss[step] = secondary structure score contribution
  float* P_posterior;   // P_posterior[step] = posterior prob for MM states (otherwise zero)
  int i1;               // First aligned residue in query
  int i2;               // Last aligned residue in query
  int j1;               // First aligned residue in template 
  int j2;               // Last aligned residue in template
  int matched_cols;     // number of matched columns in alignment against query
  int ssm1;             // SS scoring AFTER  alignment? 0:no  1:yes; t->dssp q->psipred  2:yes; q->dssp t->psipred
  int ssm2;             // SS scoring DURING alignment? 0:no  1:yes; t->dssp q->psipred  2:yes; q->dssp t->psipred
  char self;            // 0: align two different HMMs  1: align HMM with itself
  int min_overlap;      // Minimum overlap between query and template
  float sum_of_probs;   // sum of probabilities for Maximum ACcuracy alignment (if dssp states defined, only aligned pairs with defined dssp state contribute to sum)
  float Neff_HMM;       // Diversity of underlying alignment

  bool realign_around_viterbi;
  bool forward_allocated;

  // Constructor (only set pointers to NULL)
  Hit();
  ~Hit(){};
  
  // Free all allocated memory (to delete list of hits)
  void Delete();

  // Allocate/delete memory for dynamic programming matrix
  void AllocateBacktraceMatrix(int Nq, int Nt);
  void DeleteBacktraceMatrix(int Nq);
  void AllocateForwardMatrix(int Nq, int Nt);
  void DeleteForwardMatrix(int Nq);
  
  void AllocateIndices(int len);
  void DeleteIndices();

  // Compare an HMM with overlapping subalignments
  void Viterbi(HMM* q, HMM* t, const char loc, const char ssm, const int maxres, const int par_min_overlap, const float shift, const float egt, const float egq, const float ssw, const char* exclstr, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

  // Compare two HMMs with each other in lin space
  void Forward(HMM* q, HMM* t, const char ssm, const int par_min_overlap, const char loc, const float shift, const float ssw, const char* exclstr, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

  // Compare two HMMs with each other in lin space
  void Backward(HMM* q, HMM* t, const char loc, const float shift, const float ssw, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

   // Find maximum accuracy alignment (after running Forward and Backward algorithms)
  void MACAlignment(HMM* q, HMM* t, const char loc, const double mact, const double macins);

  // Trace back alignment of two profiles based on matrices btr[][]
  void Backtrace(HMM* q, HMM* t, const float corr, const float ssw, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

  // Trace back MAC alignment of two profiles based on matrix btr[][]
  void BacktraceMAC(HMM* q, HMM* t, const float corr, const float ssw, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

  // Calculate score between columns i and j of two HMMs (query and template)
  inline float ProbFwd(float* qi, float* tj) {
    return ScalarProd20(qi, tj); //
  }

  //Calculate score between columns i and j of two HMMs (query and template)
  inline float Score(float* qi, float* tj) {
    return fast_log2(ProbFwd(qi, tj));
  }

  // Calculate secondary structure score between columns i and j of two HMMs (query and template)
  inline float ScoreSS(HMM* q, HMM* t, int i, int j, int ssm, const float ssw, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]) {
    switch (ssm) //SS scoring during alignment
    {
      case 0: // no SS scoring during alignment
        return 0.0;
      case 1: // t has dssp information, q has psipred information
        return ssw
            * S73[(int) t->ss_dssp[j]][(int) q->ss_pred[i]][(int) q->ss_conf[i]];
      case 2: // q has dssp information, t has psipred information
        return ssw
            * S73[(int) q->ss_dssp[i]][(int) t->ss_pred[j]][(int) t->ss_conf[j]];
      case 3: // q has dssp information, t has psipred information
        return ssw
            * S33[(int) q->ss_pred[i]][(int) q->ss_conf[i]][(int) t->ss_pred[j]][(int) t->ss_conf[j]];
        //     case 4: // q has dssp information, t has dssp information
        //       return par.ssw*S77[ (int)t->ss_dssp[j]][ (int)t->ss_conf[j]];
    }
    return 0.0;
  }

  // Calculate secondary structure score between columns i and j of two HMMs (query and template)
  inline float ScoreSS(HMM* q, HMM* t, int i, int j, const float ssw, const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]) {
    return ScoreSS(q, t, i, j, ssm2, ssw, S73, S33);
  }

  // Calculate score for a given alignment
  void ScoreAlignment(HMM* q, HMM* t, int steps);

  // Comparison (used to sort list of hits)
  int operator<(const Hit& hit2)  {return score_sort<hit2.score_sort;}

  static int compare_score_sort(const Hit& hit1, const Hit& hit2) {
    return hit1.score_sort < hit2.score_sort;
  }

  static int compare_sum_of_probs(const Hit& hit1, const Hit& hit2) {
    return hit1.sum_of_probs < hit2.sum_of_probs;
  }

  // Calculate Evalue, score_aass, Proba from logPval and score_ss
  // Calculate Evalue, score_aass, Proba from logPval and score_ss
  inline void CalcEvalScoreProbab(int N_searched, float lamda, const char loc, const char ssm, const float ssw) {
    Eval = exp(logPval + log(N_searched));
    logEval = logPval + log(N_searched);
    // P-value = 1 - exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
    score_aass = (logPval < -10.0 ? logPval : log(-log(1 - Pval))) / 0.45
        - fmin(lamda * score_ss, fmax(0.0, 0.2 * (score - 8.0))) / 0.45 - 3.0;
    score_sort = score_aass;
    Probab = CalcProbab(loc, ssm, ssw);
  }

  /* // Merge HMM with next aligned HMM   */
  /* void MergeHMM(HMM* Q, HMM* t, float wk[]); */
  
private:
  char state;          // 0: Start/STOP state  1: MM state  2: GD state (-D)  3: IM state  4: DG state (D-)  5 MI state
  char** btr;          // backtracing matrix for all  5 pair states in one bit representation: btr[i][j] = 0|MI|DG|IM|GD|MM = 0|1|1|1|1|111
  char** cell_off;     // cell_off[i][j]=1 means this cell will get score -infinity
  double* scale;       // 

  void InitializeBacktrace(HMM* q, HMM* t);
  void InitializeForAlignment(HMM* q, HMM* t, const int min_overlap, const char ssm, const char* exclstr, bool vit=true);

  // Calculate probability of true positive : p_TP(score)/( p_TP(score)+p_FP(score) )
  // TP: same superfamily OR MAXSUB score >=0.1
  inline double CalcProbab(const char loc, const char ssm, const float ssw) {
    double s = -score_aass;
    double t;
    if (s > 200)
      return 100.0;
    if (loc) {
      if (ssm && (ssm1 || ssm2) && ssw > 0) {
        // local with SS
        const double a = sqrt(6000.0);
        const double b = 2.0 * 2.5;
        const double c = sqrt(0.12);
        const double d = 2.0 * 32.0;
        t = a * exp(-s / b) + c * exp(-s / d);
      }
      else {
        // local no SS
        const double a = sqrt(4000.0);
        const double b = 2.0 * 2.5;
        const double c = sqrt(0.15);
        const double d = 2.0 * 34.0;
        t = a * exp(-s / b) + c * exp(-s / d);
      }
    }
    else {
      if (ssm > 0 && ssw > 0) {
        // global with SS
        const double a = sqrt(4000.0);
        const double b = 2.0 * 3.0;
        const double c = sqrt(0.13);
        const double d = 2.0 * 34.0;
        t = a * exp(-s / b) + c * exp(-s / d);
      }
      else {
        // global no SS
        const double a = sqrt(6000.0);
        const double b = 2.0 * 2.5;
        const double c = sqrt(0.10);
        const double d = 2.0 * 37.0;
        t = a * exp(-s / b) + c * exp(-s / d);
      }

    }

    return 100.0 / (1.0 + t * t); // ??? JS Jul'12
  }

};


double Pvalue(double x, double a[]);
double Pvalue(float x, float lamda, float mu);
double logPvalue(float x, float lamda, float mu);
double logPvalue(float x, double a[]);


#endif
