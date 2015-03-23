//
//  hhviterbi4.h
//
//  Created by Martin Steinegger on 19.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef HHVITERBI4_h
#define HHVITERBI4_h
#include <float.h>
#include "hhhit.h"
#include "hhviterbimatrix.h"
#include "simd.h"
#include "hhhmm.h"
#include "hhhmmsimd.h"


class Viterbi {
  public:
    static const int VEC_SIZE = HMMSimd::VEC_SIZE;

    struct ViterbiResult {
        int i[VEC_SIZE];
        int j[VEC_SIZE];
        float score[VEC_SIZE];
        ViterbiResult(){
            for(int idx = 0; idx< VEC_SIZE;idx++){
                i[idx] = -1;
                j[idx] = -1;
                score[idx] = 0.0f;
            }
        }
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

    Viterbi(int maxres, bool local, float penalty_gap_query,
        float penalty_gap_template, float correlation, int par_min_overlap,
        float shift, const int ss_mode, float ssw, const float S73[NDSSP][NSSPRED][MAXCF],
        const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF], const float S37[NSSPRED][MAXCF][NDSSP]);
    ~Viterbi();

    /////////////////////////////////////////////////////////////////////////////////////
    // Align
    // Alignes two HMMSimd objects
    /////////////////////////////////////////////////////////////////////////////////////
    ViterbiResult* Align(HMMSimd* q, HMMSimd* t, ViterbiMatrix * viterbiMatrix,
        int maxres, int ss_hmm_mode);

    /////////////////////////////////////////////////////////////////////////////////////
    // Align
    // Alignes two HMMSimd objects
    /////////////////////////////////////////////////////////////////////////////////////
    void AlignWithOutCellOff(HMMSimd* q, HMMSimd* t,
        ViterbiMatrix * viterbiMatrix, int maxres, ViterbiResult* result);

    /////////////////////////////////////////////////////////////////////////////////////
    // Align with Cell Off
    // Alignes two HMMSimd objects and excludes alignments
    /////////////////////////////////////////////////////////////////////////////////////
    void AlignWithCellOff(HMMSimd* q, HMMSimd* t, ViterbiMatrix * viterbiMatrix,
        int maxres, ViterbiResult* result);

    /////////////////////////////////////////////////////////////////////////////////////
    // Align with SS score
    // Alignes two HMMSimd objects
    /////////////////////////////////////////////////////////////////////////////////////
    void AlignWithOutCellOffAndSS(HMMSimd* q, HMMSimd* t,
            ViterbiMatrix * viterbiMatrix, int maxres, ViterbiResult* result, int ss_hmm_mode);


    /////////////////////////////////////////////////////////////////////////////////////
    // Align with Cell Off and SS score
    // Alignes two HMMSimd objects
    /////////////////////////////////////////////////////////////////////////////////////
    void AlignWithCellOffAndSS(HMMSimd* q, HMMSimd* t,
            ViterbiMatrix * viterbiMatrix, int maxres, ViterbiResult* result, int ss_hmm_mode);

    /////////////////////////////////////////////////////////////////////////////////////
    // Backtrace
    // Makes backtrace from start i, j position.
    /////////////////////////////////////////////////////////////////////////////////////
    static BacktraceResult Backtrace(ViterbiMatrix * matrix, int elem,
        int start_i[VEC_SIZE], int start_j[VEC_SIZE]);

    /////////////////////////////////////////////////////////////////////////////////////
    // ScoreForBacktrace
    // Computes the score from a backtrace result
    /////////////////////////////////////////////////////////////////////////////////////
    BacktraceScore ScoreForBacktrace(HMMSimd* q_four, HMMSimd* t_four, int elem,
        Viterbi::BacktraceResult *backtraceResult,
        float alignmentScore[VEC_SIZE], int ss_hmm_mode);

    /////////////////////////////////////////////////////////////////////////////////////
    // ExcludeAlignment
    // Excludes the alignment from the Matrix
    // This alignment will not be considered at the next Align call
    /////////////////////////////////////////////////////////////////////////////////////
    static void ExcludeAlignment(ViterbiMatrix * matrix, HMMSimd* q_four,
        HMMSimd* t_four, int elem, int * i_steps, int * j_steps, int nsteps);

    /////////////////////////////////////////////////////////////////////////////////////
    // Set cell off for excl paramenter, ...
    /////////////////////////////////////////////////////////////////////////////////////
    static void InitializeForAlignment(HMM* q, HMM* t, ViterbiMatrix * matrix,
        int elem, bool self, int par_min_overlap);


    static inline simd_float ScalarProd20Vec(simd_float* qi, simd_float* tj) {
      _mm_prefetch((char * ) &qi[4], _MM_HINT_T0);
      _mm_prefetch((char * ) &tj[4], _MM_HINT_T0);
      simd_float res0 = simdf32_mul(tj[0], qi[0]);
      simd_float res1 = simdf32_mul(tj[1], qi[1]);
      simd_float res2 = simdf32_mul(tj[2], qi[2]);
      simd_float res3 = simdf32_mul(tj[3], qi[3]);
      _mm_prefetch((char * ) &qi[8], _MM_HINT_T0);
      _mm_prefetch((char * ) &tj[8], _MM_HINT_T0);
      res0 = simdf32_add(simdf32_mul(tj[ 4],qi[ 4]), res0);
      res1 = simdf32_add(simdf32_mul(tj[ 5],qi[ 5]), res1);
      res2 = simdf32_add(simdf32_mul(tj[ 6],qi[ 6]), res2);
      res3 = simdf32_add(simdf32_mul(tj[ 7],qi[ 7]), res3);
      _mm_prefetch((char * ) &qi[12], _MM_HINT_T0);
      _mm_prefetch((char * ) &tj[12], _MM_HINT_T0);
      res0 = simdf32_add(simdf32_mul(tj[ 8],qi[ 8]), res0);
      res1 = simdf32_add(simdf32_mul(tj[ 9],qi[ 9]), res1);
      res2 = simdf32_add(simdf32_mul(tj[10],qi[10]), res2);
      res3 = simdf32_add(simdf32_mul(tj[11],qi[11]), res3);
      _mm_prefetch((char * ) &qi[16], _MM_HINT_T0);
      _mm_prefetch((char * ) &tj[16], _MM_HINT_T0);
      res0 = simdf32_add(simdf32_mul(tj[12],qi[12]), res0);
      res1 = simdf32_add(simdf32_mul(tj[13],qi[13]), res1);
      res2 = simdf32_add(simdf32_mul(tj[14],qi[14]), res2);
      res3 = simdf32_add(simdf32_mul(tj[15],qi[15]), res3);

      res0 = simdf32_add(simdf32_mul(tj[16],qi[16]), res0);
      res1 = simdf32_add(simdf32_mul(tj[17],qi[17]), res1);
      res2 = simdf32_add(simdf32_mul(tj[18],qi[18]), res2);
      res3 = simdf32_add(simdf32_mul(tj[19],qi[19]), res3);

      res0 = simdf32_add(res0, res1);
      res2 = simdf32_add(res2, res3);
      return simdf32_add(res0, res2);

    }
    
    
    void ss_score_simd(__m128i score_matrix_vec01
                        , __m128i score_matrix_vec16
                        , __m128i template_sequence
                        , __m128 ssw
                        , float * storeResult){
        const __m128i sixteen_vec  = _mm_set1_epi8(16);
        const __m128i fiveteen_vec = _mm_set1_epi8(15);
        const __m128i zero = _mm_setzero_si128();
        // create slice mask
        // Example:
        //	15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
        //                      if lt 16
        //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
        __m128i lookup_mask01=_mm_cmplt_epi8(template_sequence,sixteen_vec);
        __m128i lookup_mask16=_mm_cmpgt_epi8(template_sequence,fiveteen_vec);
        // slice index
        // Example:
        //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
        //  15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
        //                          min
        //  15	12	11	0	0	0	0	11	15	12	11	0	0	0	0	155
        __m128i lookup_index01=_mm_min_epu8(lookup_mask01,template_sequence);
        __m128i lookup_index16=_mm_min_epu8(lookup_mask16,template_sequence);
        // 2xmal array lookup
        __m128i score01 = _mm_shuffle_epi8(score_matrix_vec01, lookup_index01);
        __m128i score16 = _mm_shuffle_epi8(score_matrix_vec16, lookup_index16);
        // merge 0_15 and 16_31
        __m128i res = _mm_add_epi8(score01,score16);

        __m128i lo_16 = _mm_unpacklo_epi8(res, zero);
        __m128i hi_16 = _mm_unpackhi_epi8(res,  zero);
        __m128i in1 = _mm_unpacklo_epi16(lo_16, zero);
        __m128i in2 = _mm_unpackhi_epi16(lo_16, zero);
        __m128i in3 = _mm_unpacklo_epi16(hi_16, zero);
        __m128i in4 = _mm_unpackhi_epi16(hi_16, zero);
        __m128 flt_0_3   = _mm_cvtepi32_ps(in1);
        flt_0_3   = _mm_mul_ps(flt_0_3, ssw);

        __m128 flt_4_7   = _mm_cvtepi32_ps(in2);
        flt_4_7   = _mm_mul_ps(flt_4_7, ssw);

        __m128 flt_8_11  = _mm_cvtepi32_ps(in3);
        flt_8_11   = _mm_mul_ps(flt_8_11, ssw);

        __m128 flt_12_15 = _mm_cvtepi32_ps(in4);
        flt_12_15   = _mm_mul_ps(flt_12_15, ssw);

        _mm_storer_ps(storeResult,      flt_0_3);
        _mm_storer_ps(storeResult + 4,  flt_4_7);
        _mm_storer_ps(storeResult + 8,  flt_8_11);
        _mm_storer_ps(storeResult + 12, flt_12_15);
    }
    
    // Calculate secondary structure score between columns i and j of two HMMs (query and template)
    static inline float ScoreSS(const HMM* q, const HMM* t, const int i,
                                const int j, const float ssw, const int ssm,
                                const float S73[NDSSP][NSSPRED][MAXCF],
                                const float S37[NSSPRED][MAXCF][NDSSP],
                                const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF])
    {
        switch (ssm) //SS scoring during alignment
        {
            case HMM::NO_SS_INFORMATION: // no SS scoring during alignment
                return 0.0;
            case HMM::PRED_DSSP: // t has dssp information, q has psipred information
                return ssw * S37[ (int)q->ss_pred[i]][ (int)q->ss_conf[i]][ (int)t->ss_dssp[j]];
            case HMM::DSSP_PRED: // q has dssp information, t has psipred information
                return ssw * S73[ (int)q->ss_dssp[i]][ (int)t->ss_pred[j]][ (int)t->ss_conf[j]];
            case HMM::PRED_PRED: // q has dssp information, t has psipred information
                return ssw * S33[ (int)q->ss_pred[i]][ (int)q->ss_conf[i]][ (int)t->ss_pred[j]][ (int)t->ss_conf[j]];
        }
        return 0.0;
    }


    
    
    inline void read_scoreline_sscore(char * templateSeq,
                                      float * output){
//        const simd_i zero =  _mm_setzero_si128();
//        for(unsigned int i = 0; i < L; i +=16)
//        {
//            simd_int seq = simdi_load(templateSeq+i);
//            simd_int res = scoreLookup30(seq);
//            simd_int lo_16 = _mm_unpacklo_epi8(res, zero);
//            simd_int hi_16 = _mm_unpackhi_epi8(res, zero);
//            simd_int in1 = _mm_unpacklo_epi16(lo_16, zero);
//            simd_int in2 = _mm_unpackhi_epi16(lo_16, zero);
//            simd_int in3 = _mm_unpacklo_epi16(hi_16, zero);
//            simd_int in4 = _mm_unpackhi_epi16(hi_16, zero);
//            __m128 ou1 = _mm_cvtepi32_ps(in1);
//            __m128 ou2 = _mm_cvtepi32_ps(in2);
//            __m128 ou3 = _mm_cvtepi32_ps(in3);
//            __m128 ou4 = _mm_cvtepi32_ps(in4);
//            simdf32_store(output + i, out1);
//            simdf32_store(output + i + 4, out2);
//            simdf32_store(output + i + 8, out3);
//            simdf32_store(output + i + 12, out4);
//        }
//
    }
    
    
    inline simd_int scoreLookup30(simd_int score_matrix_vec01,
                                  simd_int score_matrix_vec16,
                                  simd_int template_sequence){
//        const simd_int sixteen_vec  = simdi_set(16);
//        const simd_int fiveteen_vec = simdi_set(15);
//
//        // create slice mask
//        // Example:
//        //	15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
//        //                      if lt 16
//        //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
//        simd_int lookup_mask01 = simdi8_gt(sixteen_vec, template_sequence); // six > tmp 0-15
//        simd_int lookup_mask16 = _mm_andnot_si128(lookup_mask01,fiveteen_vec); // tmp > five 16 - 32
//        //print_sse((char *)&lookup_mask01);
//        //print_sse((char *)&lookup_mask16);
//        // slice index
//        // Example:
//        //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
//        //  15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
//        //                          min
//        //  15	12	11	0	0	0	0	11	15	12	11	0	0	0	0	155
//        simd_int lookup_index01 = simdui8_max(lookup_mask01,template_sequence);
//        simd_int lookup_index16 = simdui8_max(lookup_mask16,template_sequence);
//        //print_sse((char *)&lookup_index01);
//        //print_sse((char *)&lookup_index16);
//
//
//        // 2xmal array lookup
//        simd_int score01 = _mm_shuffle_epi8(score_matrix_vec01, lookup_index01);
//        simd_int score16 = _mm_shuffle_epi8(score_matrix_vec16, lookup_index16);
//        //print_sse((char *)&score01);
//        //print_sse((char *)&score16);
//
//        return simdui8_adds(score01,score16);
        return simdi8_set(0);
    }

    void setSSLookup(float S73[NDSSP][NSSPRED][MAXCF], float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

private:

    void PrintDebug(const HMM * q, const HMM *t,
        Viterbi::BacktraceScore * backtraceScore,
        Viterbi::BacktraceResult * backtraceResult, const int ssm);

    bool local;
    float penalty_gap_query;
    float penalty_gap_template;
    float correlation;
    int par_min_overlap;
    int max_seq_length;
    float shift;
//    char* exclstr;
    // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match)
    // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete)
    // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
    // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
    // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins)
    simd_float * sMM_DG_MI_GD_IM_vec; // one vector for cache line optimization
    // look up of linear ss scores scaled by (S/maxXX)*255
//    simd_int *ss33_lookup;
//    simd_int *ss73_lookup;

//    float S73[NDSSP][NSSPRED][MAXCF];
    const float (*S73)[NSSPRED][MAXCF];
//    float S33[NSSPRED][MAXCF][NSSPRED][MAXCF];
    const float (*S33)[MAXCF][NSSPRED][MAXCF];
//    float S37[NSSPRED][MAXCF][NDSSP];
    const float (*S37)[MAXCF][NDSSP];

    // needed to scale SS back
//    float max33;
//    float max73;
    float *ss_score;
    // weight for ss score
    float ssw;
    // set the scoring mode for ss score (DSSP_PRED, PRED_DSSP, PRED_PRED)
    int ss_mode;
};

#endif
