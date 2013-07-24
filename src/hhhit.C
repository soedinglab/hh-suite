// hhhit.C

#ifndef MAIN
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
#include "util.C"     // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.h"     // list data structure
#include "hash.h"     // hash data structure
#include "hhdecl.C"      // constants, class 
#include "hhutil.C"      // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhhmm.h"       // class HMM
#include "hhalignment.h" // class Alignment
#include "hhhitlist.h"   // class HitList
#endif

#define CALCULATE_MAX6(max, var1, var2, var3, var4, var5, var6, varb) \
  if (var1>var2) { max=var1; varb=STOP;}			          \
  else           { max=var2; varb=MM;};					\
  if (var3>max)  { max=var3; varb=GD;};					\
  if (var4>max)  { max=var4; varb=IM;};					\
  if (var5>max)  { max=var5; varb=DG;};					\
  if (var6>max)  { max=var6; varb=MI;}; 

#define CALCULATE_MAX4(max, var1, var2, var3, var4, varb)	\
  if (var1>var2) { max=var1; varb=STOP;}			\
  else           { max=var2; varb=MM;};				\
  if (var3>max)  { max=var3; varb=MI;};				\
  if (var4>max)  { max=var4; varb=IM;}; 

// Generate random number in [0,1[
#define frand() ((float) rand()/(RAND_MAX+1.0))

// Function declarations
inline float max2(const float& xMM, const float& xX, char& b, unsigned char bit); 
inline int pickprob2(const double& xMM, const double& xX, const int& state); 
inline int pickprob3_GD(const double& xMM, const double& xDG, const double& xGD); 
inline int pickprob3_IM(const double& xMM, const double& xMI, const double& xIM); 
inline int pickprob6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI); 
inline int pickmax2(const double& xMM, const double& xX, const int& state); 
inline int pickmax3_GD(const double& xMM, const double& xDG, const double& xGD); 
inline int pickmax3_IM(const double& xMM, const double& xMI, const double& xIM); 
inline int pickmax6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI); 
inline double Pvalue(double x, double a[]);
inline double Pvalue(float x, float lamda, float mu);
inline double logPvalue(float x, float lamda, float mu);
inline double logPvalue(float x, double a[]);
inline double Probab();
#ifdef HH_SSE2
inline __m128 _mm_flog2_ps(__m128 X); // Fast SSE2 log2 for four floats
#endif

/////////////////////////////////////////////////////////////////////////////////////
//// Constructor
/////////////////////////////////////////////////////////////////////////////////////
Hit::Hit()
{
  longname = name = file = dbfile = NULL;
  sname = NULL;
  seq = NULL;  
  btr = NULL;
  self = 0;
  i = j = NULL;
  // alt_i = new List<int>();
  // alt_j = new List<int>();
  alt_i = alt_j = NULL;
  states = NULL;
  S = S_ss = P_posterior = NULL;
  P_MM=NULL;
  cell_off = NULL;
  scale = NULL;
  sum_of_probs=0.0; 
  Neff_HMM=0.0;
  realign_around_viterbi=false;
  forward_allocated = false;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Free all allocated memory (to delete list of hits) 
//// Do NOT delete DP matrices
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Delete() {
  if (i)
    delete[] i;
  if (j)
    delete[] j;

  if (irep == 1) {
    if (alt_i)
      delete alt_i;
    if (alt_j)
      delete alt_j;
  }
  
  if (states)
    delete[] states;
  if (S)
    delete[] S;
  if (S_ss)
    delete[] S_ss;
  if (P_posterior)
    delete[] P_posterior;
  //  delete[] l;    
  i = j = NULL;
  states = NULL;
  S = S_ss = P_posterior = NULL;

  delete[] longname;
  delete[] name;
  delete[] file;
  delete[] dbfile;
  if (sname) {
    for (int k = 0; k < n_display; ++k)
      delete[] sname[k];
    delete[] sname;
  }
  if (seq) {
    for (int k = 0; k < n_display; ++k)
      delete[] seq[k];
    delete[] seq;
  }

  longname = name = file = NULL;
  dbfile = NULL;
  sname = NULL;
  seq = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateBacktraceMatrix(int Nq, int Nt)
{
  int i;
  btr=new(char*[Nq]);
  cell_off=new(char*[Nq]);
  for (i=0; i<Nq; ++i) 
    {
      btr[i]=new(char[Nt]);
      cell_off[i]=new(char[Nt]);
      if (!btr[i] || !cell_off[i]) 
    {
      fprintf(stderr,"Error in %s: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",par.argv[0],i+1,Nq);
      fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
      fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
      exit(3);
    }
    }
}


void Hit::DeleteBacktraceMatrix(int Nq)
{
  int i;
  for (i=0; i<Nq; ++i) 
    {
      delete[] btr[i];
      delete[] cell_off[i];
    }
  delete[] btr;
  delete[] cell_off;
  btr = NULL;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for Forward dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateForwardMatrix(int Nq, int Nt)
{
  P_MM=new(float*[Nq]);
  scale=new(double[Nq+1]); // need Nq+3?
  for (int i=0; i<Nq; ++i) 
    {
      P_MM[i] = new(float[Nt]);
      if (!P_MM[i]) 
    {
      fprintf(stderr,"Error in %s: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",par.argv[0],i+1,Nq);
      fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
      fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
      exit(3);
    }
      for (int j=0; j<Nt; ++j) 
    P_MM[i][j]=0.0; // This might be time-consuming! Is it necessary???? JS

    }
  forward_allocated = true;
}


void Hit::DeleteForwardMatrix(int Nq)
{
  for (int i=0; i<Nq; ++i) 
    delete[] P_MM[i];
  delete[] P_MM;
  delete[] scale;
  P_MM = NULL;
  forward_allocated = false;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for indices by given alignment
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateIndices(int len) {
  i = new (int[len]);
  j = new (int[len]);
}

void Hit::DeleteIndices() {
  delete[] i;
  delete[] j;
}

/////////////////////////////////////////////////////////////////////////////////////
// Compare HMMs with one another and look for sub-optimal alignments that share no pair with previous ones
// The function is called with q and t
// If q and t are equal (self==1), only the upper right part of the matrix is calculated: j>=i+3
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Viterbi(HMM* q, HMM* t) {
  // Linear topology of query (and template) HMM:
  // 1. The HMM HMM has L+2 columns. Columns 1 to L contain 
  //    a match state, a delete state and an insert state each.
  // 2. The Start state is M0, the virtual match state in column i=0 (j=0). (Therefore X[k][0]=ANY)
  //    This column has only a match state and it has only a transitions to the next match state.
  // 3. The End state is M(L+1), the virtual match state in column i=L+1.(j=L+1) (Therefore X[k][L+1]=ANY)
  //    Column L has no transitions to the delete state: tr[L][M2D]=tr[L][D2D]=0.
  // 4. Transitions I->D and D->I are ignored, since they do not appear in PsiBlast alignments 
  //    (as long as the gap opening penalty d is higher than the best match score S(a,b)). 
  
  // Pairwise alignment of two HMMs:
  // 1. Pair-states for the alignment of two HMMs are 
  //    MM (Q:Match T:Match) , GD (Q:Gap T:Delete), IM (Q:Insert T:Match),  DG (Q:Delelte, T:Match) , MI (Q:Match T:Insert) 
  // 2. Transitions are allowed only between the MM-state and each of the four other states.
  
  // Saving space:
  // The best score ending in pair state XY sXY[i][j] is calculated from left to right (j=1->t->L) 
  // and top to bottom (i=1->q->L). To save space, only the last row of scores calculated is kept in memory.
  // (The backtracing matrices are kept entirely in memory [O(t->L*q->L)]).
  // When the calculation has proceeded up to the point where the scores for cell (i,j) are caculated,
  //    sXY[i-1][j'] = sXY[j']   for j'>=j (A below)  
  //    sXY[i][j']   = sXY[j']   for j'<j  (B below)
  //    sXY[i-1][j-1]= sXY_i_1_j_1         (C below) 
  //    sXY[i][j]    = sXY_i_j             (D below)
  //                   j-1   
  //                     j
  // i-1:               CAAAAAAAAAAAAAAAAAA
  //  i :   BBBBBBBBBBBBBD
  //
  // The backtracing information is kept for all 5 pair states in a matrix btr[i][j] with a single byte per cell.
  // The last 3 bits store the previous state for the MM state, the successively higher bits are for GD, IM, DG, MI states:
  // btr[i][j] = 0|MI|DG|IM|GD|MM = 0|1|1|1|1|111

  unsigned char* csSeq = NULL;
  if (par.useCSScoring) {
    std::string id(t->name);
    csSeq = columnStateSequences[id];
  }
  
  // Variable declarations
  float __attribute__((aligned(16))) Si[par.maxres]; // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
  float sMM[par.maxres]; // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
  float sGD[par.maxres]; // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete)
  float sDG[par.maxres]; // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
  float sIM[par.maxres]; // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
  float sMI[par.maxres]; // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins)
  float smin = (par.loc ? 0 : -FLT_MAX); //used to distinguish between SW and NW algorithms in maximization
  int i, j;      //query and template match state indices
  float sMM_i_j = 0, sMI_i_j, sIM_i_j, sGD_i_j, sDG_i_j;
  float sMM_i_1_j_1, sMI_i_1_j_1, sIM_i_1_j_1, sGD_i_1_j_1, sDG_i_1_j_1;
  int jmin, jmax;

  // Reset crossed out cells?
  if (irep == 1)
    InitializeForAlignment(q, t);

  // Initialization of top row, i.e. cells (0,j)
  for (j = 0; j <= t->L; ++j) {
    sMM[j] = (self ? 0 : -j * par.egt);
    sIM[j] = sMI[j] = sDG[j] = sGD[j] = -FLT_MAX;
  }
  score = -INT_MAX;
  i2 = j2 = 0;
  btr[0][0] = STOP;

  // Viterbi algorithm

  // Loop through query positions i
  for (i = 1; i <= q->L; ++i) {
    //       if (v>=5) printf("\n");

    if (self) {
      // If q is compared to itself, ignore cells below diagonal+SELFEXCL
      jmin = i + SELFEXCL;
      jmax = t->L;
      if (jmin > jmax)
        continue;
    }
    else {
      // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
      jmin = imax(1, i + min_overlap - q->L); // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq}
      jmax = imin(t->L, i - min_overlap + t->L); // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt}
    }

    // Initialize cells
    if (jmin == 1) {
      sMM_i_1_j_1 = -(i - 1) * par.egq;  // initialize at (i-1,0)
      sMM[0] = -i * par.egq;           // initialize at (i,0)
      sIM_i_1_j_1 = sMI_i_1_j_1 = sDG_i_1_j_1 = sGD_i_1_j_1 = -FLT_MAX; // initialize at (i-1,jmin-1)
    }
    else {
      // Initialize at (i-1,jmin-1) if lower left triagonal is excluded due to min overlap
      sMM_i_1_j_1 = sMM[jmin - 1];     // initialize at (i-1,jmin-1)
      sIM_i_1_j_1 = sIM[jmin - 1];     // initialize at (i-1,jmin-1)
      sMI_i_1_j_1 = sMI[jmin - 1];     // initialize at (i-1,jmin-1)
      sDG_i_1_j_1 = sDG[jmin - 1];     // initialize at (i-1,jmin-1)
      sGD_i_1_j_1 = sGD[jmin - 1];     // initialize at (i-1,jmin-1)
      sMM[jmin - 1] = -FLT_MAX;        // initialize at (i,jmin-1)
    }
    if (jmax < t->L) // initialize at (i-1,jmmax) if upper right triagonal is excluded due to min overlap
      sMM[jmax] = sIM[jmax] = sMI[jmax] = sDG[jmax] = sGD[jmax] = -FLT_MAX;
    sIM[jmin - 1] = sMI[jmin - 1] = sDG[jmin - 1] = sGD[jmin - 1] = -FLT_MAX; // initialize at (i,jmin-1)

    // Precalculate amino acid profile-profile scores
    if (!par.useCSScoring || !csSeq) {
      #ifdef HH_SSE2
        for (j=jmin; j<=jmax; ++j)
        Si[j] = ProbFwd(q->p[i],t->p[j]);
        __m128* sij = (__m128*)(Si+(jmin/4)*4);
        for (j=jmin/4; j<=jmax/4; ++j, ++sij)
        _mm_store_ps((float*)(sij),_mm_flog2_ps(*sij));
      #else
        for (j = jmin; j <= jmax; ++j)
          Si[j] = Score(q->p[i], t->p[j]);
      #endif
    }
    else {
      for (j = jmin; j <= jmax; ++j) {
        Si[j] = flog2(columnStateScoring->substitutionScores[i][csSeq[j]]);
      }
    }


    // Loop through template positions j
    for (j=jmin; j<=jmax; ++j) {

      if (cell_off[i][j]) {
        sMM_i_1_j_1 = sMM[j]; // sMM_i_1_j_1 (for j->j+1) = sMM(i-1,(j+1)-1) = sMM[j]
        sGD_i_1_j_1 = sGD[j];
        sIM_i_1_j_1 = sIM[j];
        sDG_i_1_j_1 = sDG[j];
        sMI_i_1_j_1 = sMI[j];
        sMM[j]=sMI[j]=sIM[j]=sDG[j]=sGD[j]=-FLT_MAX; // sMM[j] = sMM(i,j) is cell_off
      }
      else {
        // Recursion relations
        //          printf("S[%i][%i]=%4.1f  ",i,j,Score(q->p[i],t->p[j])); // DEBUG!!

        CALCULATE_MAX6( sMM_i_j,
                smin,
                sMM_i_1_j_1 + q->tr[i-1][M2M] + t->tr[j-1][M2M],
                sGD_i_1_j_1 + q->tr[i-1][M2M] + t->tr[j-1][D2M],
                sIM_i_1_j_1 + q->tr[i-1][I2M] + t->tr[j-1][M2M],
                sDG_i_1_j_1 + q->tr[i-1][D2M] + t->tr[j-1][M2M],
                sMI_i_1_j_1 + q->tr[i-1][M2M] + t->tr[j-1][I2M],
                btr[i][j]
                );

        sMM_i_j += Si[j] + ScoreSS(q,t,i,j) + par.shift;


        sGD_i_j = max2
       (
       sMM[j-1] + t->tr[j-1][M2D], // MM->GD gap opening in query
       sGD[j-1] + t->tr[j-1][D2D], // GD->GD gap extension in query
       btr[i][j], 0x08
       );
        sIM_i_j = max2
      (
       sMM[j-1] + q->tr[i][M2I] + t->tr[j-1][M2M] ,
       sIM[j-1] + q->tr[i][I2I] + t->tr[j-1][M2M], // IM->IM gap extension in query
       btr[i][j], 0x10
       );
        sDG_i_j = max2
      (
       sMM[j] + q->tr[i-1][M2D],
       sDG[j] + q->tr[i-1][D2D], //gap extension (DD) in query
       btr[i][j], 0x20
       );
        sMI_i_j = max2
      (
       sMM[j] + q->tr[i-1][M2M] + t->tr[j][M2I], // MM->MI gap opening M2I in template
       sMI[j] + q->tr[i-1][M2M] + t->tr[j][I2I], // MI->MI gap extension I2I in template
       btr[i][j], 0x40
       );

        sMM_i_1_j_1 = sMM[j];
        sGD_i_1_j_1 = sGD[j];
        sIM_i_1_j_1 = sIM[j];
        sDG_i_1_j_1 = sDG[j];
        sMI_i_1_j_1 = sMI[j];
        sMM[j] = sMM_i_j;
        sGD[j] = sGD_i_j;
        sIM[j] = sIM_i_j;
        sDG[j] = sDG_i_j;
        sMI[j] = sMI_i_j;

        // Find maximum score; global alignment: maxize only over last row and last column
        if(sMM_i_j>score && (par.loc || i==q->L)) { i2=i; j2=j; score=sMM_i_j; }

      } // end if

  } //end for j

    // if global alignment: look for best cell in last column
    if (!par.loc && sMM_i_j>score) {
      i2=i;
      j2=jmax;
      score=sMM_i_j;
    }

  } // end for i

  state=MM; // state with maximum score is MM state

  //printf("Template=%-12.12s  i=%-4i j=%-4i score=%6.3f\n",t->name,i2,j2,score);  ///??? DEBUG

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Compare two HMMs with Forward Algorithm in lin-space (~ 2x faster than in log-space)
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Forward(HMM* q, HMM* t) {
  unsigned char* csSeq = NULL;
  if (par.useCSScoring) {
    std::string id(t->name);
    csSeq = columnStateSequences[id];
  }

  // Variable declarations
  int i, j;      // query and template match state indices
  double pmin = (par.loc ? 1.0 : 0.0); // used to distinguish between SW and NW algorithms in maximization
  double Cshift = pow(2.0, par.shift); // score offset transformed into factor in lin-space
  double Pmax_i;                        // maximum of F_MM in row i
  double scale_prod = 1.0;                // Prod_i=1^i (scale[i])
  int jmin;

  double F_MM_prev[t->L + 1];
  double F_GD_prev[t->L + 1];
  double F_DG_prev[t->L + 1];
  double F_IM_prev[t->L + 1];
  double F_MI_prev[t->L + 1];

  double F_MM_curr[t->L + 1];
  double F_GD_curr[t->L + 1];
  double F_DG_curr[t->L + 1];
  double F_IM_curr[t->L + 1];
  double F_MI_curr[t->L + 1];

  // First alignment of this pair of HMMs?
  if (irep == 1) {
    q->tr[0][M2D] = q->tr[0][M2I] = 0.0;
    q->tr[0][I2M] = q->tr[0][I2I] = 0.0;
    q->tr[0][D2M] = q->tr[0][D2D] = 0.0;
    t->tr[0][M2M] = 1.0;
    t->tr[0][M2D] = t->tr[0][M2I] = 0.0;
    t->tr[0][I2M] = t->tr[0][I2I] = 0.0;
    t->tr[0][D2M] = t->tr[0][D2D] = 0.0;
    q->tr[q->L][M2M] = 1.0;
    q->tr[q->L][M2D] = q->tr[q->L][M2I] = 0.0;
    q->tr[q->L][I2M] = q->tr[q->L][I2I] = 0.0;
    q->tr[q->L][D2M] = 1.0;
    q->tr[q->L][D2D] = 0.0;
    t->tr[t->L][M2M] = 1.0;
    t->tr[t->L][M2D] = t->tr[t->L][M2I] = 0.0;
    t->tr[t->L][I2M] = t->tr[t->L][I2I] = 0.0;
    t->tr[t->L][D2M] = 1.0;
    t->tr[t->L][D2D] = 0.0;
    InitializeForAlignment(q, t, false);
  }

  if (realign_around_viterbi) {
    if (irep > 1)
      InitializeForAlignment(q, t, false);

    int step;
    // fprintf(stderr,"\nViterbi-hit (Index: %i  Irep: %i) Query: %4i-%4i   Template %4i-%4i\n",index,irep,i1,i2,j1,j2);
    // fprintf(stderr,"Step Query Templ\n");
    // for (step=nsteps; step>=1; step--)
    // 	fprintf(stderr,"%4i  %4i %4i\n",step,this->i[step],this->j[step]);

    // Cross out regions
    for (i = 1; i <= q->L; ++i)
      for (j = 1; j <= t->L; ++j)
        if (!((i < i1 && j < j1) || (i > i2 && j > j2)))
          cell_off[i][j] = 1;

    // Clear Viterbi path
    for (step = nsteps; step >= 1; step--) {
      int path_width = 40;
      for (i = imax(1, this->i[step] - path_width);
          i <= imin(q->L, this->i[step] + path_width); ++i)
        cell_off[i][this->j[step]] = 0;
      for (j = imax(1, this->j[step] - path_width);
          j <= imin(t->L, this->j[step] + path_width); ++j)
        cell_off[this->i[step]][j] = 0;
    }

    // Mask previous found alternative alignments
    if (alt_i && alt_j) {
      alt_i->Reset();
      alt_j->Reset();
      while (!alt_i->End()) {
        i=alt_i->ReadNext();
        j=alt_j->ReadNext();

        for (int ii=imax(i-2,1); ii<=imin(i+2,q->L); ++ii)
          cell_off[ii][j]=1;
        for (int jj=imax(j-2,1); jj<=imin(j+2,t->L); ++jj)
          cell_off[i][jj]=1;
      }
    }

    // DEBUG!!!
    // fprintf(stdout,"\nCell_off Matrix:\n");
    // fprintf(stdout,"----+----|----+----|----+----|----+----|----+----50---+----|----+----|----+----|----+----|----+---100");
    // fprintf(stdout, "---+----|----+----|----+----|----+----|----+---150---+----|----+----|----+----|----+----|----+---200");
    // fprintf(stdout, "---+----|----+----|----+----|----+----|----+---250---+----|----+----|----+----|----+----|----+---300");
    // fprintf(stdout, "---+----|----+----|----+----|----+---340\n");
    // for (j=1; j<=t->L; ++j)
    //    {
    //      for (i=1; i<=q->L; ++i)
    //        fprintf(stderr,"%1i", cell_off[i][j]);
    //      fprintf(stderr,"\n");
    //    }
  }


  // Initialize F_XX_prev (representing i=1) and P_MM[1][j]
  F_MM_prev[0] = F_IM_prev[0] = F_GD_prev[0] = F_MI_prev[0] = F_DG_prev[0] =0.0;

  for (j=1; j<=t->L; ++j)
  {
    if (cell_off[1][j])
      F_MM_prev[j] = F_MI_prev[j] = F_DG_prev[j] = F_IM_prev[j] = F_GD_prev[j] = 0.0;
    else
    {
      float substitutionScore =
          (par.useCSScoring && csSeq) ?
              columnStateScoring->substitutionScores[1][csSeq[j]] :
              ProbFwd(q->p[1], t->p[j]);

      F_MM_prev[j] = P_MM[1][j] = substitutionScore * fpow2(ScoreSS(q,t,1,j)) * Cshift;
      F_MI_prev[j] = F_DG_prev[j] = 0.0;
      F_IM_prev[j] = F_MM_prev[j-1] * q->tr[1][M2I] * t->tr[j-1][M2M] + F_IM_prev[j-1] * q->tr[1][I2I] * t->tr[j-1][M2M];
      F_GD_prev[j] = F_MM_prev[j-1] * t->tr[j-1][M2D]                 + F_GD_prev[j-1] * t->tr[j-1][D2D];
    }
  }

  scale[0]=scale[1]=scale[2]=1.0;

  // Forward algorithm

  // Loop through query positions i
  for (i=2; i<=q->L; ++i) {
    if (self)
      jmin = imin(i+SELFEXCL+1,t->L);
    else
      jmin=1;

    if (scale_prod<DBL_MIN*100)
      scale_prod = 0.0;
    else
      scale_prod *= scale[i];

    // Initialize cells at (i,0)
    if (cell_off[i][jmin])
      F_MM_curr[jmin] = F_MI_curr[jmin] = F_DG_curr[jmin] = F_IM_curr[jmin] = F_GD_curr[jmin] = 0.0;
    else
    {
      float substitutionScore =
          (par.useCSScoring && csSeq) ?
              columnStateScoring->substitutionScores[i][csSeq[jmin]] :
              ProbFwd(q->p[i], t->p[jmin]);

      F_MM_curr[jmin] = scale_prod * substitutionScore * fpow2(ScoreSS(q,t,i,jmin)) * Cshift;
      F_IM_curr[jmin] = F_GD_curr[jmin] = 0.0;
      F_MI_curr[jmin] = scale[i] * (F_MM_prev[jmin] * q->tr[i-1][M2M] * t->tr[jmin][M2I] + F_MI_prev[jmin] * q->tr[i-1][M2M] * t->tr[jmin][I2I]);
      F_DG_curr[jmin] = scale[i] * (F_MM_prev[jmin] * q->tr[i-1][M2D]                   + F_DG_prev[jmin] * q->tr[i-1][D2D]);
    }

    /* copy back */
    P_MM[i][jmin] = F_MM_curr[jmin];

    Pmax_i=0;

    // Loop through template positions j
    for (j=jmin+1; j<=t->L; ++j) {

      // Recursion relations
      if (cell_off[i][j])
        F_MM_curr[j] = F_MI_curr[j] = F_DG_curr[j] = F_IM_curr[j] = F_GD_curr[j] = 0.0;
      else {
        float substitutionScore =
            (par.useCSScoring && csSeq) ?
                columnStateScoring->substitutionScores[i][csSeq[j]] :
                ProbFwd(q->p[i], t->p[j]);

        F_MM_curr[j] = substitutionScore * fpow2(ScoreSS(q,t,i,j)) * Cshift * scale[i] *
          ( pmin
            + F_MM_prev[j-1] * q->tr[i-1][M2M] * t->tr[j-1][M2M] // BB -> MM (BB = Begin/Begin, for local alignment)
            + F_GD_prev[j-1] * q->tr[i-1][M2M] * t->tr[j-1][D2M] // GD -> MM
            + F_IM_prev[j-1] * q->tr[i-1][I2M] * t->tr[j-1][M2M] // IM -> MM
            + F_DG_prev[j-1] * q->tr[i-1][D2M] * t->tr[j-1][M2M] // DG -> MM
            + F_MI_prev[j-1] * q->tr[i-1][M2M] * t->tr[j-1][I2M] // MI -> MM
          );
        F_GD_curr[j] =
          (   F_MM_curr[j-1] * t->tr[j-1][M2D]                    // GD -> MM
            + F_GD_curr[j-1] * t->tr[j-1][D2D]                    // GD -> GD
    //       + F_DG_curr[j-1] * t->tr[j-1][M2D] * q->tr[i][D2M] )  // DG -> GD (only when structure scores given)
          );
        F_IM_curr[j] =
          (   F_MM_curr[j-1] * q->tr[i][M2I] * t->tr[j-1][M2M]     // MM -> IM
            + F_IM_curr[j-1] * q->tr[i][I2I] * t->tr[j-1][M2M]     // IM -> IM
    //       + F_MI_curr[j-1] * q->tr[i][M2I] * t->tr[j-1][I2M] )   // MI -> IM (only when structure scores given)
          );
        F_DG_curr[j] = scale[i] *
          (   F_MM_prev[j] * q->tr[i-1][M2D]                    // DG -> MM
            + F_DG_prev[j] * q->tr[i-1][D2D]                    // DG -> DG
          ) ;
        F_MI_curr[j] = scale[i] *
          (   F_MM_prev[j] * q->tr[i-1][M2M] * t->tr[j][M2I]     // MI -> MM
            + F_MI_prev[j] * q->tr[i-1][M2M] * t->tr[j][I2I]     // MI -> MI
          );

        if(F_MM_curr[j]>Pmax_i)
          Pmax_i=F_MM_curr[j];
      } // end else

    } //end for j

    /* F_MM_prev = F_MM_curr */
    for (int jj = 0; jj <= t->L; jj++) {
      F_MM_prev[jj] = F_MM_curr[jj];
      F_MI_prev[jj] = F_MI_curr[jj];
      F_IM_prev[jj] = F_IM_curr[jj];
      F_DG_prev[jj] = F_DG_curr[jj];
      F_GD_prev[jj] = F_GD_curr[jj];

      // Fill posterior probability matrix with forward score
      P_MM[i][jj] = F_MM_curr[jj];
    }

    pmin *= scale[i];
    if (pmin<DBL_MIN*100)
      pmin = 0.0;

    scale[i+1] = 1.0/(Pmax_i+1.0);
//    scale[i+1] = 1.0;   // to debug scaling


  } // end for i

// Calculate P_forward * Product_{i=1}^{Lq+1}(scale[i])
  if (par.loc) {
    Pforward = 1.0; // alignment contains no residues (see Mueckstein, Stadler et al.)
    // Loop through query positions i
    for (i=1; i<=q->L; ++i) {
      if (self)
        jmin = imin(i+SELFEXCL+1,t->L);
      else
        jmin=1;

      for (j=jmin; j<=t->L; ++j) // Loop through template positions j
        Pforward += P_MM[i][j];

      Pforward *= scale[i+1];
    }
  }
  // global alignment
  else {
    Pforward = 0.0;
    for (i=1; i<q->L; ++i) Pforward = (Pforward + P_MM[i][t->L]) * scale[i+1];
    for (j=1; j<=t->L; ++j) Pforward += P_MM[q->L][j];
    Pforward *= scale[q->L+1];
  }

  // Calculate log2(P_forward)
  score = log2(Pforward)-10.0f;
  for (i=1; i<=q->L+1; ++i)
    score -= log2(scale[i]);

  if (par.loc) {
    if (self)
      score -= log(0.5*t->L*q->L)/LAMDA+14.; // +14.0 to get approx same mean as for -global
    else
      score -= log(t->L*q->L)/LAMDA+14.; // +14.0 to get approx same mean as for -global
  }

// Debugging output
/*
if (v>=4)
  {
    const int i0=0, i1=q->L;
    const int j0=0, j1=t->L;
    scale_prod=1;
    fprintf(stderr,"\nFwd      scale     ");
    for (j=j0; j<=j1; ++j) fprintf(stderr,"%3i     ",j);
    fprintf(stderr,"\n");
    for (i=i0; i<=i1; ++i)
  {
    scale_prod *= scale[i];
    fprintf(stderr,"%3i: %9.3G ",i,1/scale_prod);
    for (j=j0; j<=j1; ++j)
      fprintf(stderr,"%7.4f ",(P_MM[i][j]+F_MI[i][j]+F_IM[i][j]+F_DG[i][j]+F_GD[i][j]));
    fprintf(stderr,"\n");
    //      printf(" MM  %9.5f ",1/scale[i]);
    //      for (j=j0; j<=j1; ++j)
    //        printf("%7.4f ",P_MM[i][j]);
    //      printf("\n");
  }
    fprintf(stderr,"Template=%-12.12s  score=%6.3f i2=%i  j2=%i \n",t->name,score,i2,j2);
    fprintf(stderr,"\nForward total probability ratio: %8.3G\n",Pforward);
  }
  */
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Compare two HMMs with Backward Algorithm (in lin-space, 2x faster), for use in MAC alignment 
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Backward(HMM* q, HMM* t) {
  unsigned char* csSeq = NULL;
  if (par.useCSScoring) {
    std::string id(t->name);
    csSeq = columnStateSequences[id];
  }

  // Variable declarations
  int i,j;      // query and template match state indices
  double pmin;  // this is the scaled 1 in the SW algorithm that represents a starting alignment
  double Cshift = pow(2.0,par.shift);   // score offset transformed into factor in lin-space
  double scale_prod=scale[q->L+1];
  int jmin;
  
  double B_MM_prev[t->L + 1];
  double B_DG_prev[t->L + 1];
  double B_MI_prev[t->L + 1];

  double B_MM_curr[t->L + 1];
  double B_GD_curr[t->L + 1];
  double B_DG_curr[t->L + 1];
  double B_IM_curr[t->L + 1];
  double B_MI_curr[t->L + 1];

  // Initialization of top row, i.e. cells (0,j)
  for (int j=t->L; j>=1; j--) 
  {
    if (cell_off[q->L][j]) 
      P_MM[q->L][j] = B_MM_prev[j] = 0.0;
    else 
      {
    B_MM_prev[j] = scale[q->L+1];
    P_MM[q->L][j] *= scale[q->L+1]/Pforward;
      }
    B_MI_prev[j] = B_DG_prev[j] = 0.0;
  }
  if (par.loc) pmin = scale[q->L+1]; else pmin = 0.0; // transform pmin (for local alignment) to scale of present (i'th) row

  // Backward algorithm
  // Loop through query positions i
  for (i=q->L-1; i>=1; i--) {
    //       if (v>=5) printf("\n");

    if (self) jmin = imin(i+SELFEXCL,t->L); else jmin=1; // jmin = i+SELFEXCL and not (i+SELFEXCL+1) to set matrix element at boundary to zero

    // Initialize cells at (i,t->L+1)
    scale_prod *= scale[i+1];
    if (scale_prod<DBL_MIN*100)
      scale_prod = 0.0;

    if (cell_off[i][t->L])
      P_MM[i][t->L] = B_MM_curr[t->L] = 0.0;
    else {
      B_MM_curr[t->L] = scale_prod;
      P_MM[i][t->L] *= scale_prod/Pforward;
    }

    B_IM_curr[t->L] = B_MI_curr[t->L] = B_DG_curr[t->L] = B_GD_curr[t->L] = 0.0;

    pmin *= scale[i+1]; // transform pmin (for local alignment) to scale of present (i'th) row
    if (pmin<DBL_MIN*100)
      pmin = 0.0;

    // Loop through template positions j
    for (j=t->L-1; j>=jmin; j--) {
    // Recursion relations
    //          printf("S[%i][%i]=%4.1f  ",i,j,Score(q->p[i],t->p[j]));

      if (cell_off[i][j])
        B_MM_curr[j] = B_GD_curr[j] = B_IM_curr[j] = B_DG_curr[j] = B_MI_curr[j] = 0.0;
      else {
        float substitutionScore =
            (par.useCSScoring && csSeq) ?
                columnStateScoring->substitutionScores[i+1][csSeq[j+1]] :
                ProbFwd(q->p[i+1], t->p[j+1]);

        double pmatch = B_MM_prev[j+1] * substitutionScore * fpow2(ScoreSS(q,t,i+1,j+1)) * Cshift * scale[i+1];
        B_MM_curr[j] =
        (
         + pmin                                                    // MM -> EE (End/End, for local alignment)
         + pmatch       * q->tr[i][M2M]   * t->tr[j][M2M]              // MM -> MM
         + B_GD_curr[j+1]                 * t->tr[j][M2D]              // MM -> GD (q->tr[i][M2M] is already contained in GD->MM)
         + B_IM_curr[j+1] * q->tr[i][M2I] * t->tr[j][M2M]              // MM -> IM
         + B_DG_prev[j] * q->tr[i][M2D]                 * scale[i+1] // MM -> DG (t->tr[j][M2M] is already contained in DG->MM)
         + B_MI_prev[j] * q->tr[i][M2M]   * t->tr[j][M2I] * scale[i+1] // MM -> MI
         );

          B_GD_curr[j] =
        (
         + pmatch       * q->tr[i][M2M] * t->tr[j][D2M]              // GD -> MM
         + B_GD_curr[j+1]               * t->tr[j][D2D]              // DG -> DG
         );

          B_IM_curr[j] =
        (
         + pmatch         * q->tr[i][I2M] * t->tr[j][M2M]            // IM -> MM
         + B_IM_curr[j+1] * q->tr[i][I2I] * t->tr[j][M2M]            // IM -> IM
         );

          B_DG_curr[j] =
        (
         + pmatch       * q->tr[i][D2M] * t->tr[j][M2M]              // DG -> MM
         + B_DG_prev[j] * q->tr[i][D2D]                * scale[i+1]  // DG -> DG
  //             + B_GD[i][j+1] * q->tr[i][D2M] * t->tr[j][M2D]              // DG -> GD
         );

          B_MI_curr[j] =
        (
         + pmatch       * q->tr[i][M2M] * t->tr[j][I2M]              // MI -> MM
         + B_MI_prev[j] * q->tr[i][M2M] * t->tr[j][I2I] * scale[i+1] // MI -> MI
    //           + B_IM[i][j+1] * q->tr[i][M2I] * t->tr[j][I2M]              // MI -> IM
         );

      } // end else

      // Calculate posterior probability from Forward and Backward matrix elements
      P_MM[i][j] *= B_MM_curr[j]/Pforward;
    } //end for j

    for(int jj = 0; jj <= t->L; jj++) {
      B_MM_prev[jj] = B_MM_curr[jj];
      B_DG_prev[jj] = B_DG_curr[jj];
      B_MI_prev[jj] = B_MI_curr[jj];
    }
  } // end for i
  
  /*
  // Debugging output
  if (v>=6)
    {
      const int i0=0, i1=q->L;
      const int j0=0, j1=t->L;
      double scale_prod[q->L+2];
      scale_prod[q->L] = scale[q->L+1];
      for (i=q->L-1; i>=1; i--) scale_prod[i] = scale_prod[i+1] * scale[i+1];

      printf("\nBwd      scale     ");
      for (j=j0; j<=j1; ++j) printf("%3i     ",j);
      printf("\n");
      for (i=i0; i<=i1; ++i)
    {
      printf("%3i: %9.3G ",i,1/scale_prod[i]);
      for (j=j0; j<=j1; ++j)
        printf("%7.4f ",(B_MM[i][j]+B_MI[i][j]+B_IM[i][j]+B_DG[i][j]+B_GD[i][j]) * (ProbFwd(q->p[i],t->p[j])*fpow2(ScoreSS(q,t,i,j)) * Cshift));
      printf("\n");

      //      printf("MM   %9.5f ",1/scale[i]);
      //      for (j=j0; j<=j1; ++j)
      //        printf("%7.4f ",B_MM[i][j] * (ProbFwd(q->p[i],t->p[j])*fpow2(ScoreSS(q,t,i,j)) * Cshift));
      //      printf("\n");
    }
      printf("\nPost     scale     ");
      for (j=j0; j<=j1; ++j) printf("%3i     ",j);
      printf("\n");
      for (i=i0; i<=i1; ++i) 
    {
      printf("%3i: %9.3G ",i,1/scale_prod[i]);
      for (j=j0; j<=j1; ++j)
        printf("%7.4f ",B_MM[i][j]*P_MM[i][j]/Pforward);
      printf("\n");
    }
      printf("\n");
    }

  if (v>=4) printf("\nForward total probability ratio: %8.3G\n",Pforward);
  */

  return;
}



/////////////////////////////////////////////////////////////////////////////////////
// Maximum Accuracy alignment 
/////////////////////////////////////////////////////////////////////////////////////
void Hit::MACAlignment(HMM* q, HMM* t) {
  // Use Forward and Backward matrices to find that alignment which 
  // maximizes the expected number of correctly aligned pairs of residues (mact=0)
  // or, more generally, which maximizes the expectation value of the number of 
  // correctly aligned pairs minus (mact x number of aligned pairs)
  // "Correctly aligned" can be based on posterior probabilities calculated with
  // a local or a global version of the Forward-Backward algorithm.

  int i,j;           // query and template match state indices
  int jmin,jmax;     // range of dynamic programming for j
  double S_prev[t->L + 1];    // scores
  double S_curr[t->L + 1];    // scores
  double score_MAC;   // score of the best MAC alignment
  const double GAPPENALTY = 0.5*(1.0 - (double) par.macins)*par.mact;

  // float S[1000][2000]; //DEBUG

  // Initialization of top row, i.e. cells (0,j)
  for (j=0; j<=t->L; ++j)
    S_prev[j] = 0.0;
  score_MAC=-INT_MAX; i2=j2=0; btr[0][0]=STOP;

  // Dynamic programming 
  // Loop through query positions i
  for (i=1; i<=q->L; ++i) {
    if (self) {
      // If q is compared to itself, ignore cells below diagonal+SELFEXCL
      jmin = i+SELFEXCL;
      jmax = t->L;

      if (jmin>jmax)
        continue;
    }
    else {
      // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
      jmin=imax( 1, i+min_overlap-q->L);  // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq}
      jmax=imin(t->L,i-min_overlap+t->L);  // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt}
    }

    // Initialize cells
    S_curr[jmin-1] = 0.0;
    if (jmax<t->L) S_prev[jmax] = 0.0; // initialize at (i-1,jmax) if upper right triagonal is excluded due to min overlap

    // Loop through template positions j
    for (j=jmin; j<=jmax; ++j) {
      if (cell_off[i][j]) {
        S_curr[j] = -FLT_MIN;
        btr[i][j] = STOP;
        //          if (i>135 && i<140)
        //        printf("Cell off   i=%i  j=%i b:%i\n",i,j,btr[i][j]);
      }
      else {
        // Recursion

        // GAPPENALTY suppresses alignments XX-xXX IF the posterior prob for aligning (x,y) as in XXxXX
        //                                  YYy-YY                                                YYyYY
        // is larger than par.macins. For par.macins = 0.0, such insertions are always compressed and aligned on top of each other.

        // NOT the state before the first MM state)
        CALCULATE_MAX4(
               S_curr[j],
               P_MM[i][j] - par.mact,     // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
               S_prev[j-1] + P_MM[i][j] - par.mact, // P_MM[i][j] contains posterior probability
               S_prev[j] - GAPPENALTY,
               S_curr[j-1] - GAPPENALTY,
               btr[i][j]   // backtracing matrix
               );

        //if (i>36 && i<40 && j>2200 && j<2230)
        //printf("i=%i  j=%i  S[i][j]=%8.3f  MM:%7.3f  MI:%7.3f  IM:%7.3f  b:%i\n",i,j,S[i][j],S[i-1][j-1]+P_MM[i][j]-par.mact,S[i-1][j],S[i][j-1],btr[i][j]);

        // Find maximum score; global alignment: maximize only over last row and last column
        if(S_curr[j]>score_MAC && (par.loc || i==q->L)) {
          i2=i; j2=j; score_MAC=S_curr[j];
        }

//        S[i][j] = S_curr[j]; //DEBUG
      } // end else
    } //end for j
      
      // if global alignment: look for best cell in last column
      if (!par.loc && S_curr[jmax]>score_MAC) { i2=i; j2=jmax; score_MAC=S_curr[jmax]; }
      
      for (j=0; j<=t->L; ++j)
    S_prev[j] = S_curr[j];
    } // end for i
  

  // // DEBUG
  // if (!strncmp(t->name,"UP20|QED",8))
  //    {
  //      printf("\nTemplate=%-12.12s  i=%-4i j=%-4i score=%6.3f  irep=%i, Pforward=%6.3f\n",t->name,i2,j2,score_MAC,irep,Pforward);
  //     printf("\nP_MM  ");
  //     for (j=0; j<=j2; ++j) printf("%3i   ",j);
  //     printf("\n");
  //     for (i=0; i<=i2; ++i) 
  //    {
  //      printf("%2i:    ",i);
  //      for (j=0; j<=j2; ++j)
  //        printf("%5.2f ",P_MM[i][j]);
  //      printf("\n");
  //    }
  //     printf("\nScore  ");
  //     for (j=0; j<=j2; ++j) printf("%3i   ",j);
  //     printf("\n");
  //     for (i=0; i<=i2; ++i) 
  //    {
  //      printf("%2i:    ",i);
  //      for (j=0; j<=j2; ++j)
  //        printf("%5.2f ",S[i][j]);
  //      printf("\n");
  //    }
  //     printf("***************\n");
  //   }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Trace back Viterbi alignment of two profiles based on matrix btr[][]
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Backtrace(HMM* q, HMM* t) {
  // Trace back trough the matrices bXY[i][j] until first match state is found (STOP-state)

  unsigned char* csSeq = NULL;
  if (par.useCSScoring) {
    std::string id(t->name);
    csSeq = columnStateSequences[id];
  }

  int step;      // counts steps in path through 5-layered dynamic programming matrix
  int i,j;       // query and template match state indices

  InitializeBacktrace(q,t);

  // Make sure that backtracing stops when t:M1 or q:M1 is reached (Start state), e.g. sMM[i][1], or sIM[i][1] (M:MM, B:IM)
  for (i=0; i<=q->L; ++i)
    btr[i][1] = STOP;

  for (j=1; j<=t->L; ++j)
    btr[1][j] = STOP;
  
  // Back-tracing loop
  matched_cols=0; // for each MACTH (or STOP) state matched_col is incremented by 1
  step=0;         // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  //  state=MM;       // state with maximum score must be MM state  // already set at the end of Viterbi()
  i=i2; j=j2;     // last aligned pair is (i2,j2)

  // while (state!=STOP)  because STOP=0
  while (state && i>0 && j>0) {
    step++;
    states[step] = state;
    this->i[step] = i;
    this->j[step] = j;

      // Exclude cells in direct neighbourhood from all further alignments
    for (int ii=imax(i-2,1); ii<=imin(i+2,q->L); ++ii)
      cell_off[ii][j]=1;
    for (int jj=imax(j-2,1); jj<=imin(j+2,t->L); ++jj)
      cell_off[i][jj]=1;
      
    switch (state) {
      case MM: // current state is MM, previous state is btr[i][j]
        matched_cols++;
        state = (0x07 & btr[i--][j--]);
        break;
      case GD: // current state is GD
        if (0x08 & btr[i][j--]) state = MM;
        break;
      case IM:
        if (0x10 & btr[i][j--]) state = MM;
        break;
      case DG:
        if (0x20 & btr[i--][j]) state = MM;
        break;
      case MI:
        if (0x40 & btr[i--][j]) state = MM;
        break;
      default:
        fprintf(stderr,"Error in %s: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i) with template %s\n",par.argv[0],state,step,i,j,t->name);
        fprintf(stderr,"Dumping alignment and terminating:\n");
        state=0;
        v=4;
        break;
    } //end switch (state)
  } //end while (state)
 
  i1 = this->i[step];
  j1 = this->j[step];

  if (state!=0) {
    fprintf(stderr,"Error in %s: reached  (i,j)=(%i,%i) in state value %i at  at step %i  with template %s during backtracing,\n",par.argv[0],i,j,state,step,t->name);
    fprintf(stderr,"Dumping alignment and terminating:\n");
    state=0;
    v=4;
  }

  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step;
  
  // Allocate new space for alignment scores
  S    = new( float[nsteps+1]);
  S_ss = new( float[nsteps+1]);
  if (!S_ss) MemoryError("space for HMM-HMM alignments");

  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  score_ss=0.0f;
  int ssm=ssm1+ssm2;
  for (step=1; step<=nsteps; step++) {
    switch(states[step]) {
      case MM:
        i = this->i[step];
        j = this->j[step];

        S[step] = (par.useCSScoring && csSeq) ?
            columnStateScoring->substitutionScores[i][csSeq[j]] :
            Score(q->p[i], t->p[j]);

        S_ss[step] = ScoreSS(q,t,i,j,ssm);
        score_ss += S_ss[step];
        break;
      case MI: //if gap in template
      case DG:
      default: //if gap in T or Q
        S[step]=S_ss[step]=0.0f;
        break;
    }
  }

  if (ssm2>=1)
    score-=score_ss;    // subtract SS score added during alignment!!!!

  // Add contribution from correlation of neighboring columns to score
  float Scorr=0;
  if (nsteps) {
    for (step=2; step<=nsteps; step++)
      Scorr+=S[step]*S[step-1];
    for (step=3; step<=nsteps; step++)
      Scorr+=S[step]*S[step-2];
    for (step=4; step<=nsteps; step++)
      Scorr+=S[step]*S[step-3];
    for (step=5; step<=nsteps; step++)
      Scorr+=S[step]*S[step-4];
    score+=par.corr*Scorr;
  }
  
  // Set score, P-value etc.
  score_sort = score_aass = -score;
  logPval=0; Pval=1;
  if (t->mu) {
    logPvalt=logPvalue(score,t->lamda,t->mu);
    Pvalt=Pvalue(score,t->lamda,t->mu);
  }
  else {
    logPvalt=0;
    Pvalt=1;
  }
  //   printf("%-10.10s lamda=%-9f  score=%-9f  logPval=%-9g\n",name,t->lamda,score,logPvalt);
  

  //DEBUG: Print out Viterbi path and exit
  if (v>=4) {
    printf("NAME=%7.7s score=%7.3f  score_ss=%7.3f\n",name,score,score_ss);
    printf("step  Q T    i    j  state   score    bt T Q cf ss-score\n");

    for (step=nsteps; step>=1; step--) {
      switch(states[step]) {
        case MM:
          printf("%4i  %1c %1c ",step,q->seq[q->nfirst][this->i[step]],seq[nfirst][this->j[step]]);
          break;
        case GD:
        case IM:
          printf("%4i  - %1c ",step,seq[nfirst][this->j[step]]);
          break;
        case DG:
        case MI:
          printf("%4i  %1c - ",step,q->seq[q->nfirst][this->i[step]]);
          break;
      }

      printf("%4i %4i     %2i %7.2f    %2x ",this->i[step],this->j[step],(int)states[step],S[step],btr[this->i[step]][this->j[step]]);
      printf("%c %c %1i %7.2f\n",i2ss(t->ss_dssp[this->j[step]]),i2ss(q->ss_pred[this->i[step]]),q->ss_conf[this->i[step]]-1,S_ss[step]);
    }

    cerr<<"Exiting in DEBUG output loop of BacktraceViterbi\n";
    exit(1);
  }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Trace back alignment of two profiles based on matrix btr[][]
/////////////////////////////////////////////////////////////////////////////////////
void Hit::BacktraceMAC(HMM* q, HMM* t) {
  // Trace back trough the matrix b[i][j] until STOP state is found
  unsigned char* csSeq = NULL;
  if (par.useCSScoring) {
    std::string id(t->name);
    csSeq = columnStateSequences[id];
  }

  char** b=btr;  // define alias for backtracing matrix
  int step;      // counts steps in path through 5-layered dynamic programming matrix
  int i,j;       // query and template match state indices

  InitializeBacktrace(q,t);
  
  // Make sure that backtracing stops when t:M1 or q:M1 is reached (Start state), e.g. sMM[i][1], or sIM[i][1] (M:MM, B:IM)
  for (i=0; i<=q->L; ++i)
    b[i][1] = STOP;
  for (j=1; j<=t->L; ++j)
    b[1][j] = STOP;
  
  // Back-tracing loop
  // In contrast to the Viterbi-Backtracing, STOP signifies the first Match-Match state, NOT the state before the first MM state
  matched_cols=1; // for each MACTH (or STOP) state matched_col is incremented by 1
  state=MM;       // lowest state with maximum score must be match-match state
  step=0;         // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  i=i2; j=j2;     // last aligned pair is (i2,j2)

  if (b[i][j] != MM) {
    if (v>3) fprintf(stderr,"Error in %s: backtrace does not start in match-match state, but in state %i, (i,j)=(%i,%i)\n",par.argv[0],b[i][j],i,j);
    step = 1;
    this->i[step] = i;
    this->j[step] = j;
    alt_i->Push(i);
    alt_j->Push(j);
    state = STOP;
  }
  else {
    while (state!=STOP) {
      step++;
      states[step] = state = b[i][j];
      this->i[step] = i;
      this->j[step] = j;
      alt_i->Push(i);
      alt_j->Push(j);
      // Exclude cells in direct neighbourhood from all further alignments
      for (int ii=imax(i-2,1); ii<=imin(i+2,q->L); ++ii)
        cell_off[ii][j]=1;
      for (int jj=imax(j-2,1); jj<=imin(j+2,t->L); ++jj)
        cell_off[i][jj]=1;
      if (state==MM) matched_cols++;

      switch (state) {
        case MM:
          i--;
          j--;
          break;
        case IM:
          j--;
          break;
        case MI:
          i--;
          break;
        case STOP:
          break;
        default:
          fprintf(stderr,"Errorin %s: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i)\n",par.argv[0],state,step,i,j);
          state=0;
          v=4;
          break;
      } //end switch (state)
    } //end while (state)
  }

  i1 = this->i[step];
  j1 = this->j[step];
  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step;

  // Allocate new space for alignment scores
  S    = new( float[nsteps+1]);
  S_ss = new( float[nsteps+1]);
  P_posterior = new( float[nsteps+1]);

  if (!P_posterior)
    MemoryError("space for HMM-HMM alignments");

  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  score_ss = 0.0f;
  sum_of_probs = 0.0;       // number of identical residues in query and template sequence
  int ssm = ssm1 + ssm2;
  //   printf("Hit=%s\n",name); /////////////////////////////////////////////////////////////
  for (step = 1; step <= nsteps; step++) {
    switch(states[step]) {
      case MM:
        i = this->i[step];
        j = this->j[step];

        S[step] = (par.useCSScoring && csSeq) ?
            columnStateScoring->substitutionScores[i][csSeq[j]] :
            Score(q->p[i], t->p[j]);

        S_ss[step] = ScoreSS(q,t,i,j,ssm);
        score_ss += S_ss[step];
        P_posterior[step] = P_MM[this->i[step]][this->j[step]];
        // Add probability to sum of probs if no dssp states given or dssp states exist and state is resolved in 3D structure
        if (t->nss_dssp<0 || t->ss_dssp[j]>0) sum_of_probs += P_posterior[step];
        //      printf("j=%-3i  dssp=%1i  P=%4.2f  sum=%6.2f\n",j,t->ss_dssp[j],P_posterior[step],sum_of_probs); //////////////////////////
        break;
      case MI: //if gap in template
      case DG:
      default: //if gap in T or Q
        S[step] = S_ss[step] = P_posterior[step] = 0.0;
        break;
    }
  }

  //   printf("\n"); /////////////////////////////////////////////////////////////
  if (ssm2>=1)
    score-=score_ss;    // subtract SS score added during alignment!!!!

  // Add contribution from correlation of neighboring columns to score
  float Scorr=0;
  if (nsteps) {
    for (step=1; step<=nsteps-1; step++)
      Scorr+=S[step]*S[step+1];
    for (step=1; step<=nsteps-2; step++)
      Scorr+=S[step]*S[step+2];
    for (step=1; step<=nsteps-3; step++)
      Scorr+=S[step]*S[step+3];
    for (step=1; step<=nsteps-4; step++)
      Scorr+=S[step]*S[step+4];
    score+=par.corr*Scorr;
  }
  
  // Set score, P-value etc.
  score_sort = score_aass = -score;
  logPval=0;
  Pval=1;

  if (t->mu) {
    logPvalt=logPvalue(score,t->lamda,t->mu);
    Pvalt=Pvalue(score,t->lamda,t->mu);
  }
  else {
    logPvalt=0;
    Pvalt=1;
  }
  //   printf("%-10.10s lamda=%-9f  score=%-9f  logPval=%-9g\n",name,t->lamda,score,logPvalt);

  //DEBUG: Print out MAC alignment path
  if (v>=4) {
    float sum_post=0.0;
    printf("NAME=%7.7s score=%7.3f  score_ss=%7.3f\n",name,score,score_ss);
    printf("step  Q T    i    j  state   score    T Q cf ss-score   P_post Sum_post\n");
    for (step=nsteps; step>=1; step--) {
      switch(states[step]) {
        case MM:
          sum_post+=P_posterior[step];
          printf("%4i  %1c %1c ",step,q->seq[q->nfirst][this->i[step]],seq[nfirst][this->j[step]]);
          break;
        case IM:
          printf("%4i  - %1c ",step,seq[nfirst][this->j[step]]);
          break;
        case MI:
          printf("%4i  %1c - ",step,q->seq[q->nfirst][this->i[step]]);
          break;
      }
      printf("%4i %4i     %2i %7.1f    ",this->i[step],this->j[step],(int)states[step],S[step]);
      printf("%c %c  %1i  %7.1f  ",i2ss(t->ss_dssp[this->j[step]]),i2ss(q->ss_pred[this->i[step]]),q->ss_conf[this->i[step]]-1,S_ss[step]);
      printf("%7.5f  %7.2f\n",P_posterior[step],sum_post);
    }
  }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Functions that calculate probabilities
/////////////////////////////////////////////////////////////////////////////////////
void Hit::InitializeForAlignment(HMM* q, HMM* t, bool vit) {
  int i,j;

  if (vit) {
    alt_i = alt_j = NULL;
  } else {
    if (irep == 1) {
      //    if (alt_i && alt_i->Size()>0) delete alt_i;
      alt_i = new List<int>();
      //    if (alt_j && alt_j->Size()>0) delete alt_j;
      alt_j = new List<int>();
    }
  }

  // SS scoring during (ssm2>0) or after (ssm1>0) alignment? Query SS known or Template SS known?
  switch (par.ssm) {
    case 0:
      ssm1=0;
      ssm2=0;
      break;
    case 1:
      ssm2=0;  // SS scoring after alignment
      if (t->nss_dssp>=0 && q->nss_pred>=0)
        ssm1=1;
      else if (q->nss_dssp>=0 && t->nss_pred>=0)
        ssm1=2;
      else if (q->nss_pred>=0 && t->nss_pred>=0)
        ssm1=3;
      else
        ssm1=0;
      break;
    case 2:
      ssm1=0;  // SS scoring during alignment
      if (t->nss_dssp>=0 && q->nss_pred>=0)
        ssm2=1;
      else if (q->nss_dssp>=0 && t->nss_pred>=0)
        ssm2=2;
      else if (q->nss_pred>=0 && t->nss_pred>=0)
        ssm2=3;
      else
        ssm2=0;
      break;
    case 3:
      ssm2=0;  // SS scoring after alignment
      if (q->nss_pred>=0 && t->nss_pred>=0)
        ssm1=3;
      else
        ssm1=0;
      break;
    case 4:
      ssm1=0;  // SS scoring during alignment
      if (q->nss_pred>=0 && t->nss_pred>=0)
        ssm2=3;
      else
        ssm2=0;
      break;
      //     case 5:
      //       ssm2=0;  // SS scoring after alignment
      //       if (q->nss_dssp>=0 && t->nss_dssp>=0) ssm1=4; else ssm1=0;  
      //       break;
      //     case 6:
      //       ssm1=0;  // SS scoring during alignment
      //       if (q->nss_dssp>=0 && t->nss_dssp>=0) ssm2=4; else ssm2=0;
      //       break;
  }

  if (self) {
      // Cross out cells in lower diagonal for self-comparison?
      for (i=1; i<=q->L; ++i) {
        int jmax = imin(i+SELFEXCL,t->L);
        for (j=1; j<=jmax; ++j)
          cell_off[i][j]=1;   // cross out cell near diagonal
        for (j=jmax+1; j<=t->L+1; ++j)
          cell_off[i][j]=0;   // no other cells crossed out yet
      }
  }
  else {
    if (par.block_shading && !strcmp(par.block_shading_mode,"tube") && par.block_shading->Contains(t->file)) {
      // Deactivate all cells in dynamic programming matrix
      for (i=1; i<=q->L; ++i)
        for (j=1; j<=t->L; ++j)
          cell_off[i][j]=1;

      int* tmp = par.block_shading->Show(t->file);
      int counter = par.block_shading_counter->Show(t->file);

      //printf("Hit %s:\n",t->name);

      int m = 0;
      while (m < counter) {
        int i0 = tmp[m++];
        int i1 = tmp[m++];
        int j0 = tmp[m++];
        int j1 = tmp[m++];

        int d1 = imin(i0-j0,i1-j1) - par.block_shading_space;
        int d2 = imax(i0-j0,i1-j1) + par.block_shading_space;

        //printf("query: %i-%i   template: %i-%i    d1: %i  d2: %i\n",i0,i1,j0,j1,d1,d2);

        int istart = imax(1,i0-par.block_shading_space);
        int istop = imin(q->L,i1+par.block_shading_space);
        int jstart = imax(1,j0-par.block_shading_space);
        int jstop = imin(t->L,j1+par.block_shading_space);

        for (i = istart; i <= istop; ++i)
          for (j = jstart; j <= jstop; ++j)
            if ((i-j) > d1 && (i-j) < d2)
              cell_off[i][j]=0;
      }

      // int after = 0;
      // for (i=1; i<=q->L; ++i)
      //   for (j=1; j<=t->L; ++j)
      //     {
      //    if (cell_off[i][j]==1)
      //      after++;
      //     }
      // printf("cells cross out : %i\n",after);
      // printf("cells total     : %i\n",q->L*t->L);
      // printf("Ersparnis       : %4.2f %%\n",(double)after/(double)(q->L*t->L)*100.0);

      min_overlap = 0;

      // Cross out cells that are excluded by the minimum-overlap criterion
      // if (par.min_overlap==0)
      //   min_overlap = imin(60, (int)(0.333f*imin(q->L,t->L))+1); // automatic minimum overlap
      // else
      //   min_overlap = imin(par.min_overlap, (int)(0.8f*imin(q->L,t->L)));

      // for (i=0; i<min_overlap; ++i)
      //   for (j=i-min_overlap+t->L+1; j<=t->L; ++j) // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt}
      //     cell_off[i][j]=1;
      // for (i=q->L-min_overlap+1; i<=q->L; ++i)
      //   for (j=1; j<i+min_overlap-q->L; ++j)      // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq}
      //     cell_off[i][j]=1;

    }
    // Compare two different HMMs Q and T
    else {
      // Activate all cells in dynamic programming matrix
      for (i=1; i<=q->L; ++i)
        for (j=1; j<=t->L; ++j)
          cell_off[i][j]=0;   // no other cells crossed out yet

      // Cross out cells that are excluded by the minimum-overlap criterion
      if (par.min_overlap==0)
        min_overlap = imin(60, (int)(0.333f*imin(q->L,t->L))+1); // automatic minimum overlap
      else
        min_overlap = imin(par.min_overlap, (int)(0.8f*imin(q->L,t->L)));

      for (i=0; i<min_overlap; ++i)
        for (j=i-min_overlap+t->L+1; j<=t->L; ++j) // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt}
          cell_off[i][j]=1;
      for (i=q->L-min_overlap+1; i<=q->L; ++i)
        for (j=1; j<i+min_overlap-q->L; ++j)      // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq}
          cell_off[i][j]=1;

      // Cross out rows which are contained in range given by exclstr ("3-57,238-314")
      if (par.exclstr) {
        char* ptr=par.exclstr;
        int i0, i1;
        while (1) {
          i0 = abs(strint(ptr));
          i1 = abs(strint(ptr));
          if (!ptr) break;
          for (i=i0; i<=imin(i1,q->L); ++i)
            for (j=1; j<=t->L; ++j)
              cell_off[i][j]=1;
        }
      }

      // Cross out cells not contained in the range of the prefiltering in HHblits
      if (par.block_shading && par.block_shading->Contains(t->file)) {
        int* tmp = par.block_shading->Show(t->file);
        int counter = par.block_shading_counter->Show(t->file);

        int i0, i1, j0, j1;
        i0 = j0 = 1000000;
        i1 = j1 = 0;

        int m = 0;

        while (m < counter) {
          i0 = imin(tmp[m++]-par.block_shading_space,i0);
          i1 = imax(tmp[m++]+par.block_shading_space,i1);
          j0 = imin(tmp[m++]-par.block_shading_space,j0);
          j1 = imax(tmp[m++]+par.block_shading_space,j1);
        }

          //printf("Hit: %40s query: %i-%i   template: %i-%i \n",t->name,i0,i1,j0,j1);

          //printf("cell_off in query: 1-%i\n",i0);
        for (i=1; i<i0; ++i)
          for (j=1; j<=t->L; ++j)
            cell_off[i][j]=1;
          //printf("cell_off in query: %i-%i\n",i1+1,q->L);

        for (i=i1+1; i<=q->L; ++i)
          for (j=1; j<=t->L; ++j)
            cell_off[i][j]=1;
          //printf("cell_off in template: 1-%i\n",j0);

        for (j=1; j<j0; ++j)
          for (i=1; i<=q->L; ++i)
            cell_off[i][j]=1;
          //printf("cell_off in template: %i-%i\n",j1+1,t->L);
        for (j=j1+1; j<=t->L; ++j)
          for (i=1; i<=q->L; ++i)
            cell_off[i][j]=1;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Allocate memory for data of new alignment (sequence names, alignment, scores,...)
/////////////////////////////////////////////////////////////////////////////////////
void Hit::InitializeBacktrace(HMM* q, HMM* t) {
  //Copy information about template profile to hit and reset template pointers to avoid destruction
  longname = new (char[strlen(t->longname) + 1]);
  name = new (char[strlen(t->name) + 1]);
  file = new (char[strlen(t->file) + 1]);
  if (!file)
    MemoryError(
        "space for alignments with database HMMs. \nNote that all alignments have to be kept in memory");
  strcpy(longname, t->longname);
  strcpy(name, t->name);
  strcpy(fam, t->fam);
  strcpy(sfam, t->sfam);
  strcpy(fold, t->fold);
  strcpy(cl, t->cl);
  strcpy(file, t->file);

  // Allocate new space
  this->i = new (int[i2 + j2 + 2]);
  this->j = new (int[i2 + j2 + 2]);
  states = new (char[i2 + j2 + 2]);
  S = S_ss = P_posterior = NULL;

  sname = new (char*[t->n_display]);
  seq = new (char*[t->n_display]);
  if (!sname || !seq)
    MemoryError(
        "space for alignments with database HMMs.\nNote that all sequences for display have to be kept in memory");
  
  if (irep == 1) {
    // Make flat copy for first alignment of template seqs to save speed
    for (int k = 0; k < t->n_display; k++) {
      sname[k] = t->sname[k];
      seq[k] = t->seq[k];
    }
    t->dont_delete_seqs = 1;
  }
  else {
    // Make deep copy for all further alignments
    for (int k = 0; k < t->n_display; k++) {
      sname[k] = new (char[strlen(t->sname[k]) + 1]);
      seq[k] = new (char[strlen(t->seq[k]) + 1]);
      strcpy(sname[k], t->sname[k]);
      strcpy(seq[k], t->seq[k]);
    }
    if (dbfile) {
      char* ptr = dbfile;
      dbfile = new (char[strlen(ptr) + 1]);
      strcpy(dbfile, ptr);
    }
  }

  n_display = t->n_display;
  ncons = t->ncons;
  nfirst = t->nfirst;
  nss_dssp = t->nss_dssp;
  nsa_dssp = t->nsa_dssp;
  nss_pred = t->nss_pred;
  nss_conf = t->nss_conf;
  L = t->L;
  Neff_HMM = t->Neff_HMM;
  Eval = 1.0;
  logEval = 0.0;
  Pval = 1.0;
  Pvalt = 1.0;
  logPval = 0.0;
  logPvalt = 0.0;
  Probab = 1.0;
}

/////////////////////////////////////////////////////////////////////////////////////
// Some score functions 
/////////////////////////////////////////////////////////////////////////////////////

// Calculate score for a given alignment
void Hit::ScoreAlignment(HMM* q, HMM* t, int steps) {
  unsigned char* csSeq = NULL;
  if (par.useCSScoring) {
    std::string id(t->name);
    csSeq = columnStateSequences[id];
  }

  score = 0;
  for (int step = 0; step < steps; step++) {
    float substitutionScore =
        (par.useCSScoring && csSeq) ?
            columnStateScoring->substitutionScores[i[step]][csSeq[j[step]]] :
            Score(q->p[i[step]], t->p[j[step]]);

    if (v > 2) {
      cout << "Score at step " << step << "!\n";
      cout << "i: " << i[step] << "  j: " << j[step] << "   score: " << substitutionScore << "\n";
    }
    score += substitutionScore;
  }
}


//Calculate score between columns i and j of two HMMs (query and template)
inline float Hit::Score(float* qi, float* tj) {
  return fast_log2(ProbFwd(qi,tj));
}

// Calculate score between columns i and j of two HMMs (query and template)
inline float Hit::ProbFwd(float* qi, float* tj) {
  return ScalarProd20(qi,tj); //
}


// Calculate secondary structure score between columns i and j of two HMMs (query and template)
inline float Hit::ScoreSS(HMM* q, HMM* t, int i, int j, int ssm)
{
  switch (ssm) //SS scoring during alignment 
    {
    case 0: // no SS scoring during alignment 
      return 0.0;
    case 1: // t has dssp information, q has psipred information 
      return par.ssw * S73[ (int)t->ss_dssp[j]][ (int)q->ss_pred[i]][ (int)q->ss_conf[i]];
    case 2: // q has dssp information, t has psipred information 
      return par.ssw * S73[ (int)q->ss_dssp[i]][ (int)t->ss_pred[j]][ (int)t->ss_conf[j]];
    case 3: // q has dssp information, t has psipred information 
      return par.ssw * S33[ (int)q->ss_pred[i]][ (int)q->ss_conf[i]][ (int)t->ss_pred[j]][ (int)t->ss_conf[j]];
      //     case 4: // q has dssp information, t has dssp information 
      //       return par.ssw*S77[ (int)t->ss_dssp[j]][ (int)t->ss_conf[j]];
    }
  return 0.0;
}


// Calculate secondary structure score between columns i and j of two HMMs (query and template)
inline float Hit::ScoreSS(HMM* q, HMM* t, int i, int j)
{
  return ScoreSS(q,t,i,j,ssm2);
}

// /////////////////////////////////////////////////////////////////////////////////////
// //// Function for Viterbi()
// /////////////////////////////////////////////////////////////////////////////////////
// inline float max2(const float& xMM, const float& xX, char& b) 
// {
//   if (xMM>xX) { b=MM; return xMM;} else { b=SAME;  return xX;}
// }
inline float max2(const float& xMM, const float& xSAME, char& b, const unsigned char bit) 
{
  if (xMM>xSAME) {b |= bit; return xMM;} else {  /* b |= 0x00!*/ return xSAME;}
}


/////////////////////////////////////////////////////////////////////////////////////
// Functions for StochasticBacktrace()
/////////////////////////////////////////////////////////////////////////////////////

inline int pickprob2(const double& xMM, const double& xX, const int& state)
{
  if ( (xMM+xX)*frand() < xMM) return MM; else return state;
}
inline int pickprob3_GD(const double& xMM, const double& xDG, const double& xGD)
{
  double x = (xMM+xDG+xGD)*frand();
  if ( x<xMM) return MM;
  else if ( x<xMM+xDG) return DG;
  else return GD;
}
inline int pickprob3_IM(const double& xMM, const double& xMI, const double& xIM)
{
  double x = (xMM+xMI+xIM)*frand();
  if ( x<xMM) return MM;
  else if ( x<xMM+xMI) return MI;
  else return IM;
}
inline int pickprob6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI)
{
  double x = (x0+xMM+xGD+xIM+xDG+xMI)*frand();
  x-=xMM; if (x<0) return MM;
  x-=x0;  if (x<0) return STOP;
  x-=xGD; if (x<0) return GD;
  x-=xIM; if (x<0) return IM;
  if (x < xDG) return DG; else return MI;
}

inline int pickmax2(const double& xMM, const double& xX, const int& state)
{
  if (xMM > xX) return MM; else return state;
}
inline int pickmax3_GD(const double& xMM, const double& xDG, const double& xGD)
{
  char state;
  double x;
  if ( xMM>xDG) {state=MM; x=xMM;}
  else          {state=DG; x=xDG;}
  if ( xGD>x)   {state=GD; x=xGD;}
  return state;
}inline int pickmax3_IM(const double& xMM, const double& xMI, const double& xIM)
{
  char state;
  double x;
  if ( xMM>xMI) {state=MM; x=xMM;}
  else          {state=MI; x=xMI;}
  if ( xIM>x)   {state=IM; x=xIM;}
  return state;
}
inline int pickmax6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI)
{
  char state;
  double x;
  if ( x0 >xMM) {state=STOP; x=x0;}
  else          {state=MM; x=xMM;}
  if ( xGD>x)   {state=GD; x=xGD;}
  if ( xIM>x)   {state=IM; x=xIM;}
  if ( xDG>x)   {state=DG; x=xDG;}
  if ( xMI>x)   {state=MI; x=xMI;}
  return state;
}



/////////////////////////////////////////////////////////////////////////////////////
//// Functions that calculate P-values and probabilities 
/////////////////////////////////////////////////////////////////////////////////////


//// Evaluate the CUMULATIVE extreme value distribution at point x
//// p(s)ds = lamda * exp{ -exp[-lamda*(s-mu)] - lamda*(s-mu) } ds = exp( -exp(-x) - x) dx = p(x) dx
//// => P(s>S) = integral_-inf^inf {p(x) dx}  = 1 - exp{ -exp[-lamda*(S-mu)] }
inline double Pvalue(double x, double a[])
{
  //a[0]=lamda, a[1]=mu
  double h = a[0]*(x-a[1]);
  return (h>10)? exp(-h) : double(1.0)-exp( -exp(-h));
}

inline double Pvalue(float x, float lamda, float mu)
{
  double h = lamda*(x-mu);
  return (h>10)? exp(-h) : (double(1.0)-exp( -exp(-h)));
}

inline double logPvalue(float x, float lamda, float mu)
{
  double h = lamda*(x-mu);
  return (h>10)? -h : (h<-2.5)? -exp(-exp(-h)): log( ( double(1.0) - exp(-exp(-h)) ) );
}

inline double logPvalue(float x, double a[])
{
  double h = a[0]*(x-a[1]);
  return (h>10)? -h : (h<-2.5)? -exp(-exp(-h)): log( ( double(1.0) - exp(-exp(-h)) ) );
}

// Calculate probability of true positive : p_TP(score)/( p_TP(score)+p_FP(score) )
// TP: same superfamily OR MAXSUB score >=0.1
inline double Hit::CalcProbab()
{
  double s=-score_aass;
  double t;
  if (s>200) return 100.0; 
  if (par.loc) 
    {
      if (par.ssm && (ssm1 || ssm2) && par.ssw>0) 
    {
      // local with SS
      const double a=sqrt(6000.0);
      const double b=2.0*2.5;
      const double c=sqrt(0.12);
      const double d=2.0*32.0;
      t = a*exp(-s/b) + c*exp(-s/d);
    }
      else
    {
      // local no SS
      const double a=sqrt(4000.0);
      const double b=2.0*2.5;
      const double c=sqrt(0.15);
      const double d=2.0*34.0;
      t = a*exp(-s/b) + c*exp(-s/d);
    }
    }
  else
    {
      if (par.ssm>0 && par.ssw>0)
    {
      // global with SS
      const double a=sqrt(4000.0);
      const double b=2.0*3.0;
      const double c=sqrt(0.13);
      const double d=2.0*34.0;
      t = a*exp(-s/b) + c*exp(-s/d);
    }
      else
    {
      // global no SS
      const double a=sqrt(6000.0);
      const double b=2.0*2.5;
      const double c=sqrt(0.10);
      const double d=2.0*37.0;
      t = a*exp(-s/b) + c*exp(-s/d);
    }

    }

  return 100.0/(1.0+t*t); // ??? JS Jul'12
}

// Calculate Evalue, score_aass, Proba from logPval and score_ss
inline void Hit::CalcEvalScoreProbab(int N_searched, float lamda)
{
  Eval    = exp(logPval+log(N_searched));
  logEval = logPval+log(N_searched);
  // P-value = 1 - exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
  score_aass = (logPval<-10.0? logPval : log(-log(1-Pval)) )/0.45 - fmin(lamda*score_ss,fmax(0.0,0.2*(score-8.0)))/0.45 - 3.0;
  score_sort = score_aass;
  Probab = CalcProbab();
}
