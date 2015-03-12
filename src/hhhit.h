// hhhit.h
#ifndef HHHIT_H_
#define HHHIT_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <ctype.h>
#include <vector>

class Hit;

#include "list.h"
#include "hhhmm.h"
#include "hhdecl.h"
#include "hhutil.h"
#include "util.h"
#include "hhdatabase.h"
#include "log.h"

#include "hhhit-inl.h"

struct Posterior_Triple {
	int query_pos;
	int template_pos;
	float posterior_probability;

	Posterior_Triple(int query_pos, int template_pos, float posterior_probability) {
		this->query_pos = query_pos;
		this->template_pos = template_pos;
		this->posterior_probability = posterior_probability;
	}
};

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

  HHEntry* entry;
  
  float* forward_profile;
  float* backward_profile;

  size_t backward_entries;
  float** backward_matrix;
  size_t forward_entries;
  float** forward_matrix;
  size_t posterior_entries;
  float** posterior_matrix;

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
  float sum_of_probs;   // sum of probabilities for Maximum ACcuracy alignment (if dssp states defined, only aligned pairs with defined dssp state contribute to sum)
  float Neff_HMM;       // Diversity of underlying alignment

  bool realign_around_viterbi;

  float qsc;

  // Constructor (only set pointers to NULL)
  Hit();
  ~Hit(){};
  
  // Free all allocated memory (to delete list of hits)
  void Delete();

  void AllocateIndices(int len);
  void DeleteIndices();

  // Comparison (used to sort list of hits)
  int operator<(const Hit& hit2)  {return score_sort<hit2.score_sort;}

  static int compare_score_sort(const Hit& hit1, const Hit& hit2) {
    return hit1.score_sort < hit2.score_sort;
  }

  static int compare_sum_of_probs(const Hit& hit1, const Hit& hit2) {
    return hit1.sum_of_probs > hit2.sum_of_probs;
  }

  static bool compare_evalue(const Hit& hit1, const Hit& hit2) {
    return hit1.Eval < hit2.Eval;
  }

  void initHitFromHMM(HMM * t, const int nseqdis);

  float calculateSimilarity(HMM* q, const float S[20][20]);

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
  char state;          // 0: Start/STOP state  1: MM state  2: GD state (-D)  3: IM state  4: DG state (D-)  5 MI state
  int min_overlap;
  
private:
  // Calculate probability of true positive : p_TP(score)/( p_TP(score)+p_FP(score) )
  // TP: same superfamily OR MAXSUB score >=0.1
  inline double CalcProbab(const char loc, const char ssm, const float ssw) {
    double s = -score_aass;
    double t = 0;
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

struct ViterbiScores {
  ViterbiScores(){
    score = 0.0;
    score_aass = 0.0;
    score_ss = 0.0;
    Pval = 0.0;
    Pvalt = 0.0;
    logPval = 0.0;
    logPvalt = 0.0;
    Eval = 0.0;
    logEval = 0.0;
    Probab = 0.0;
  };

  ViterbiScores(Hit& hit) {
    score = hit.score;
    score_aass = hit.score_aass;
    score_ss = hit.score_ss;
    Pval = hit.Pval;
    Pvalt = hit.Pvalt;
    logPval = hit.logPval;
    logPvalt = hit.logPvalt;
    Eval = hit.Eval;
    logEval = hit.logEval;
    Probab = hit.Probab;
  }

  float score;
  float score_aass;
  float score_ss;
  double Pval;
  double Pvalt;
  double logPval;
  double logPvalt;
  double Eval;
  double logEval;
  float Probab;
};


double Pvalue(double x, double a[]);
double Pvalue(float x, float lamda, float mu);
double logPvalue(float x, float lamda, float mu);
double logPvalue(float x, double a[]);
int compareHitLengths(const void * a, const void * b);

#endif
