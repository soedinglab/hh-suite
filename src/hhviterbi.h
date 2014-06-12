//
//  hhviterbi4.h
//
//  Created by Martin Steinegger on 19.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef HHVITERBI4_h
#define HHVITERBI4_h
#include <float.h>

#include "hhviterbimatrix.h"
#include "simd.h"
#include "hhhmmsimd.h"


class Viterbi
{
public:
    static const int VEC_SIZE=HMMSimd::VEC_SIZE;
    
    struct ViterbiResult {
        int i[VEC_SIZE];
        int j[VEC_SIZE];
        float score[VEC_SIZE];
    };
    
    
    struct BacktraceResult {
        int * i_steps;
        int * j_steps;
        char * states;
        int count;
        int matched_cols;
    };
    
    struct BacktraceScore {
        float score_ss;
        float score;
        float score_sort;
        float score_aass;
        float Pvalt;
        float logPvalt;
        float * S;
        float * S_ss;
    };
    
    Viterbi(int maxres,bool local,float penalty_gap_query,float penalty_gap_template, float correlation, int par_min_overlap, float shift);
    ~Viterbi();
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Align
    // Alignes two HMMSimd objects
    /////////////////////////////////////////////////////////////////////////////////////
    ViterbiResult Align(HMMSimd* q, HMMSimd* t,ViterbiMatrix * viterbiMatrix);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Align
    // Alignes two HMMSimd objects
    /////////////////////////////////////////////////////////////////////////////////////
    ViterbiResult AlignWithOutCellOff(HMMSimd* q, HMMSimd* t,ViterbiMatrix * viterbiMatrix);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Align with Cell Off
    // Alignes two HMMSimd objects and excludes alignments
    /////////////////////////////////////////////////////////////////////////////////////
    ViterbiResult AlignWithCellOff(HMMSimd* q, HMMSimd* t,ViterbiMatrix * viterbiMatrix);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Backtrace
    // Makes backtrace from start i, j position.
    /////////////////////////////////////////////////////////////////////////////////////
    static BacktraceResult Backtrace(ViterbiMatrix * matrix,int elem,int * start_i,int * start_j);
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    // ScoreForBacktrace
    // Computes the score from a backtrace result
    /////////////////////////////////////////////////////////////////////////////////////
    BacktraceScore  ScoreForBacktrace(HMMSimd* q_four, HMMSimd* t_four,int elem,
                                                      Viterbi::BacktraceResult *backtraceResult,
                                                      float * alignmentScore, int ssm1,int ssm2);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // ExcludeAlignment
    // Excludes the alignment from the Matrix
    // This alignment will not be considered at the next Align call
    /////////////////////////////////////////////////////////////////////////////////////
    static void ExcludeAlignment(ViterbiMatrix * matrix,HMMSimd* q_four, HMMSimd* t_four,int elem,
                                 int * i_steps, int * j_steps, int nsteps);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Set cell off for excl paramenter, ...
    /////////////////////////////////////////////////////////////////////////////////////
    static void InitializeForAlignment(HMM* q, HMM* t, ViterbiMatrix * matrix,int elem, bool self, int par_min_overlap);
    
    //Calculate score between columns i and j of two HMMs (query and template)
    static inline float Score(float* qi, float* tj)
    {
        return fast_log2(ProbFwd(qi,tj));
    }

    // Calculate score between columns i and j of two HMMs (query and template)
    static inline float ProbFwd(float* qi, float* tj)
    {
        return ScalarProd20(qi,tj); //
    }


    // Calculate secondary structure score between columns i and j of two HMMs (query and template)
    static inline float ScoreSS(const HMM* q, const HMM* t, const int i, const int j, const int ssm)
    {
        //TODO
        return 0.0;
    }

    static inline simd_float ScalarProd20Vec(simd_float* qi, simd_float* tj)
    {
        _mm_prefetch((char *) &qi[4] , _MM_HINT_T0 );
        _mm_prefetch((char *) &tj[4] , _MM_HINT_T0 );
        simd_float res0 = simdf32_mul(tj[ 0],qi[ 0]);
        simd_float res1 = simdf32_mul(tj[ 1],qi[ 1]);
        simd_float res2 = simdf32_mul(tj[ 2],qi[ 2]);
        simd_float res3 = simdf32_mul(tj[ 3],qi[ 3]);
        _mm_prefetch((char *) &qi[8] , _MM_HINT_T0 );
        _mm_prefetch((char *) &tj[8] , _MM_HINT_T0 );
        res0 = simdf32_add(simdf32_mul(tj[ 4],qi[ 4]),res0);
        res1 = simdf32_add(simdf32_mul(tj[ 5],qi[ 5]),res1);
        res2 = simdf32_add(simdf32_mul(tj[ 6],qi[ 6]),res2);
        res3 = simdf32_add(simdf32_mul(tj[ 7],qi[ 7]),res3);
        _mm_prefetch((char *) &qi[12] , _MM_HINT_T0 );
        _mm_prefetch((char *) &tj[12] , _MM_HINT_T0 );
        res0 = simdf32_add(simdf32_mul(tj[ 8],qi[ 8]),res0);
        res1 = simdf32_add(simdf32_mul(tj[ 9],qi[ 9]),res1);
        res2 = simdf32_add(simdf32_mul(tj[10],qi[10]),res2);
        res3 = simdf32_add(simdf32_mul(tj[11],qi[11]),res3);
        _mm_prefetch((char *) &qi[16] , _MM_HINT_T0 );
        _mm_prefetch((char *) &tj[16] , _MM_HINT_T0 );
        res0 = simdf32_add(simdf32_mul(tj[12],qi[12]),res0);
        res1 = simdf32_add(simdf32_mul(tj[13],qi[13]),res1);
        res2 = simdf32_add(simdf32_mul(tj[14],qi[14]),res2);
        res3 = simdf32_add(simdf32_mul(tj[15],qi[15]),res3);
        
        res0 = simdf32_add(simdf32_mul(tj[16],qi[16]),res0);
        res1 = simdf32_add(simdf32_mul(tj[17],qi[17]),res1);
        res2 = simdf32_add(simdf32_mul(tj[18],qi[18]),res2);
        res3 = simdf32_add(simdf32_mul(tj[19],qi[19]),res3);
        
        res0 = simdf32_add(res0,res1);
        res2 = simdf32_add(res2,res3);
        
        return simdf32_add(res0,res2);
        
    }

private:




    static void PrintDebug(const HMM * q,const HMM *t,Viterbi::BacktraceScore * backtraceScore,Viterbi::BacktraceResult * backtraceResult,
                           const int ssm);
    
    int maxres;
    bool local;
    float penalty_gap_query;
    float penalty_gap_template;
    float correlation;
    int par_min_overlap;
    float shift;
//    char* exclstr;
    // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
    // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete)
    // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
    // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
    // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins)
    simd_float * sMM_DG_MI_GD_IM_vec; // one vector for cache line optimization

    
};


#endif
