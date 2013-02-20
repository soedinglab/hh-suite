// hhhit.h

#ifdef HH_SSE3
#ifdef __SUNPRO_C
#include <sunmedia_intrin.h>
#else
#include <emmintrin.h>
#include <pmmintrin.h>
#endif
#endif

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
  
  int nsteps;           // index for last step in Viterbi path; (first=1)
  int* i;               // i[step] = query match state at step of Viterbi path
  int* j;               // j[step] = template match state at step of Viterbi path
  List<int>* alt_i;      // Path of alternative alignments (query positions)
  List<int>* alt_j;      // Path of alternative alignments (template positions)
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
  bool backward_allocated;

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
  void AllocateBackwardMatrix(int Nq, int Nt);
  void DeleteBackwardMatrix(int Nq);
  
  void AllocateIndices(int len);
  void DeleteIndices();

  // Compare an HMM with overlapping subalignments
  void Viterbi(HMM* q, HMM* t, float** Sstruc=NULL);

  // Compare two HMMs with each other in lin space
  void Forward(HMM* q, HMM* t, float** Pstruc=NULL);

  // Compare two HMMs with each other in lin space
  void Backward(HMM* q, HMM* t);

   // Find maximum accuracy alignment (after running Forward and Backward algorithms)
  void MACAlignment(HMM* q, HMM* t);

  // Trace back alignment of two profiles based on matrices bXX[][]
  void Backtrace(HMM* q, HMM* t);

  // Trace back alignment of two profiles based on matrices bXX[][]
  void StochasticBacktrace(HMM* q, HMM* t, char maximize=0);

  // Trace back MAC alignment of two profiles based on matrix bMM[][]
  void BacktraceMAC(HMM* q, HMM* t);

  // Calculate secondary structure score between columns i and j of two HMMs (query and template)
  inline float ScoreSS(HMM* q, HMM* t, int i, int j, int ssm);

  // Calculate secondary structure score between columns i and j of two HMMs (query and template)
  inline float ScoreSS(HMM* q, HMM* t, int i, int j);

  // Calculate in log2 space the amino acid similarity score between columns i and j of two HMMs (query and template)
  inline float Score(float* qi, float* tj);

  // Calculate in lin space the amino acid similarity score between columns i and j of two HMMs (query and template)
  inline float ProbFwd(float* qi, float* tj);

  // Calculate score for a given alignment
  void ScoreAlignment(HMM* q, HMM* t, int steps);

  // Comparison (used to sort list of hits)
  int operator<(const Hit& hit2)  {return score_sort<hit2.score_sort;}

  // Calculate Evalue, score_aass, Proba from logPval and score_ss
  void CalcEvalScoreProbab(int N_searched, float lamda);

  /* // Merge HMM with next aligned HMM   */
  /* void MergeHMM(HMM* Q, HMM* t, float wk[]); */

  double** B_MM;        // Backward matrices
  
private:
  char state;          // 0: Start/STOP state  1: MM state  2: GD state (-D)  3: IM state  4: DG state (D-)  5 MI state
  char** bMM;          // (backtracing) bMM[i][j] = STOP:start of alignment  MM:prev was MM  GD:prev was GD etc
  char** bGD;          // (backtracing) bMM[i][j] = STOP:start of alignment  MM:prev was MM  SAME:prev was GD
  char** bDG;          // (backtracing)
  char** bIM;          // (backtracing)
  char** bMI;          // (backtracing)
  char** cell_off;     // cell_off[i][j]=1 means this cell will get score -infinity

  double** F_MM;        // Forward matrices 
  /*
  double** F_GD;        // F_XY[i][j] * Prod_1^i(scale[i]) 
  double** F_DG;        //   = Sum_x1..xl{ P(HMMs aligned up to Xi||Yj co-emmitted x1..xl ) / (Prod_k=1^l f(x_k)) }   
  double** F_IM;        // end gaps are not penalized!
  double** F_MI;        // 
  */
  double* scale;        // 

  void InitializeBacktrace(HMM* q, HMM* t);
  void InitializeForAlignment(HMM* q, HMM* t, bool vit=true);
  double CalcProbab();

};


double Pvalue(double x, double a[]);
double Pvalue(float x, float lamda, float mu);
double logPvalue(float x, float lamda, float mu);
double logPvalue(float x, double a[]);




