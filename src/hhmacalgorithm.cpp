/*
 * hhmacalgorithm.C
 *
 *  Created on: 14/03/2014
 *      Author: Stefan Haunsberger
 *
 *  Copyright (c) 2014 -. All rights reserved.
 *
 *	Info:
 * 		Use Forward and Backward matrices to find that alignment which
 * 		maximizes the expected number of correctly aligned pairs of residues (mact=0)
 * 		or, more generally, which maximizes the expectation value of the number of
 * 		correctly aligned pairs minus (mact x number of aligned pairs)
 * 		"Correctly aligned" can be based on posterior probabilities calculated with
 * 		a local or a global version of the Forward-Backward algorithm.
 *
 */

#include "hhposteriordecoder.h"
#include "simd.h"
#include <float.h>
#include <utility>

using std::swap;

#define MAX2(vec1, vec2, vec3, res)           \
mask_gt = (simd_int)simdf32_gt(vec1,vec2);          \
index_mask = simdi_and(mask_gt,vec3);          \
res        = simdui8_max(res,index_mask);

void PosteriorDecoder::macAlgorithm(HMMSimd & q_hmm, HMMSimd & t_hmm,
		std::vector<Hit *> & hit_vec, PosteriorMatrix & p_mm,
		ViterbiMatrix & viterbi_matrix, float par_mact) {

  //  char * nam = "PF07714";
  //  int elem_idx = 0;
  //  bool eq = false;
  //  for (int elem = 0; elem < (int)hit_vec.size(); elem++) {
  //      if (!strcmp(hit_vec.at(elem)->name, nam)) {
  //          eq = true;
  //          elem_idx = elem;
  //          break;
  //      }
  //  }

      const simd_int stop_vec = simdi32_set(ViterbiMatrix::STOP); //    00000000
      const simd_int mm_vec   = simdi32_set(ViterbiMatrix::MM);   // MM 00000010
      const simd_int im_vec   = simdi32_set(ViterbiMatrix::IM);   // IM 00000100
      const simd_int mi_vec   = simdi32_set(ViterbiMatrix::MI);   // MI 00000110

      const __m128i tmp_vec   = _mm_set_epi32(0x40000000,0x00400000,0x00004000,0x00000040);//01000000010000000100000001000000

  #ifdef AVX2_SUPPORT
      const simd_int co_vec               = _mm256_inserti128_si256(_mm256_castsi128_si256(tmp_vec), tmp_vec, 1);
      const simd_int shuffle_mask_extract = _mm256_setr_epi8(0,  4,  8,  12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                                                                                                  -1, -1, -1,  -1,  0,  4,  8, 12, -1, -1, -1, -1, -1, -1, -1, -1);
      const simd_int shuffle_mask_celloff = _mm256_set_epi8(
                                                            15, 14, 13, 12,
                                                            15, 14, 13, 12,
                                                            15, 14, 13, 12,
                                                            15, 14, 13, 12,
                                                            3, 2,  1, 0,
                                                            3, 2,  1, 0,
                                                            3, 2,  1, 0,
                                                            3, 2,  1, 0);
  #else
      const simd_int co_vec = tmp_vec;
  #endif

      const simd_float mact = simdf32_set(par_mact);
      const simd_float mact_half = simdf32_mul(simdf32_set(0.5), mact);
      const simd_float zeros_float = simdf32_setzero();
      const simd_int float_min_vec = (simd_int) simdf32_set(-FLT_MAX);

      simd_int state_result = simdi32_set(0);

      simd_float score_mac;
      simd_int i2_vec = simdi32_set(0);
      simd_int j2_vec = simdi32_set(0);
      int i,j; // query and template match state indices

      // same as: s_curr[0] = simdf32_setzero();
      simdf32_store((float *)&m_s_curr[0], zeros_float);

      // Initialization of top row, i.e. cells (0,j)
      for (j = 0; j <= t_hmm.L; j++) {
          simdf32_store((float *)&m_s_prev[j], zeros_float);
      }

      score_mac = simdf32_set(-FLT_MAX);

  #ifdef AVX2_SUPPORT
      unsigned long long * sCO_MI_DG_IM_GD_MM_vec = (unsigned long long *) viterbi_matrix.getRow(0);
  #else
      unsigned int   * sCO_MI_DG_IM_GD_MM_vec   = (unsigned int *) viterbi_matrix.getRow(0);
  #endif
      // same as: hit.bMM[0][0] = STOP;
      state_result = stop_vec;
      // Set "state_result" to Viterbi matrix
      // write values back to ViterbiMatrix
  #ifdef AVX2_SUPPORT
      /* byte_result_vec        000H  000G  000F  000E   000D  000C  000B  000A */
      /* abcdefgh               0000  0000  HGFE  0000   0000  0000  0000  DCBA */
      const __m256i abcdefgh = _mm256_shuffle_epi8(state_result, shuffle_mask_extract);
      /* abcd                                            0000  0000  0000  DCBA */
      const __m128i abcd     = _mm256_castsi256_si128(abcdefgh);
      /* efgh                                            0000  0000  HGFE  0000 */
      const __m128i efgh     = _mm256_extracti128_si256(abcdefgh, 1);
      _mm_storel_epi64((__m128i*)&sCO_MI_DG_IM_GD_MM_vec[0], _mm_or_si128(abcd, efgh));
  #else
      state_result = _mm_packs_epi32(state_result, state_result);
      state_result = _mm_packus_epi16(state_result, state_result);
      int int_result = _mm_cvtsi128_si32(state_result);
      sCO_MI_DG_IM_GD_MM_vec[0] = int_result;
  #endif

      ///////////////////////////////////////////////////////////////
      // Loop through query positions i
      for (i = 1; i <= q_hmm.L; i++) {
          const simd_float * p_mm_row_ptr = p_mm.getRow(i);   // pointer to current row of the posterior probability matrix
  #ifdef AVX2_SUPPORT
          unsigned long long * sCO_MI_DG_IM_GD_MM_vec = (unsigned long long *) viterbi_matrix.getRow(i);
  #else
          unsigned int   * sCO_MI_DG_IM_GD_MM_vec   = (unsigned int *) viterbi_matrix.getRow(i);
  #endif
          const simd_int curr_pos_i = simdi32_set(i);

          /////////////////////////////////////////////////////////////
          // Loop through template positions j
          for (j = 1; j <= t_hmm.L; j++) {

              const simd_int curr_pos_j = simdi32_set(j);
              // Find maximum score; global alignment: maxize only over last row and last column
              const bool findMaxInnerLoop = (m_local || i == q_hmm.L);

  //      CALCULATE_MAX4(
  //              S_curr[j],
  //              fpow2(hit.P_MM[i][j]) - par.mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
  //              S_prev[j-1] + fpow2(hit.P_MM[i][j]) - par.mact, // hit.P_MM[i][j] contains log2-posterior probability
  //              S_prev[j] - 0.5 * par.mact,  // gap penalty prevents alignments such as this: XX--xxXX
  //              S_curr[j-1] - 0.5 * par.mact,  //                                               YYyy--YY
  //              val //  hit.bMM[i][j]   backtracing matrix
  //      );
  //CALCULATE_MAX4(max, var1, var2, var3, var4, varb)
  //  if (var1>var2) { max=var1; varb=STOP;}
  //  else           { max=var2; varb=MM;};
  //  if (var3>max)  { max=var3; varb=MI;};
  //  if (var4>max)  { max=var4; varb=IM;};

              // check p_mm values -> OK
  //          if (eq && viterbi_matrix.getCellOff(i, j, elem_idx) == 0 && j <= 259) {
  //              fprintf(stdout, "%2.20f\n", ((float*)&p_mm_row_ptr[j])[elem_idx]);
  //          }

              const simd_float s_prev_j_min_1 = simdf32_load((float *)&m_s_prev[j-1]);
              const simd_float s_prev_j               = simdf32_load((float *)&m_s_prev[j]);
              const simd_float s_curr_j_min_1 = simdf32_load((float *)&m_s_curr[j-1]);
              const simd_float p_mm_i_j           = simdf32_load((float *)&p_mm_row_ptr[j]);

              simd_int mask_gt;
              simd_float s_mm_i_j = simdf32_set(0);
              simd_int index_mask = simdi32_set(0);
              // set to zero == STOP
              state_result = stop_vec;
              // same as: powf(2, p_mm.getSingleValue(i, j, elem)) - par.mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
              const simd_float term_1 = simdf32_sub(simdf32_fpow2(p_mm_i_j), mact);
              // same as: S_prev[j-1] + powf(2, p_mm.getSingleValue(i, j, elem)) - par.mact, // hit.P_MM[i][j] contains log2-posterior probability
              const simd_float term_2 = simdf32_add(s_prev_j_min_1, term_1);

              // if term1 > term2 { 1 = STOP } -> max = term1 (set as default)
              s_mm_i_j = term_1;
              // if term2 > term1 { 2 = MM } -> max = term2, otherwise keep zero
              mask_gt = (simd_int)simdf32_gt(term_2, term_1);
              state_result = simdi_and(mask_gt, mm_vec);
              s_mm_i_j = simdf32_max(s_mm_i_j, term_2);

              // same as: S_prev[j] - 0.5 * par.mact,  // gap penalty prevents alignments such as this: XX--xxXX
              const simd_float term_3 = simdf32_sub(s_prev_j, mact_half);
              // if term3 > s_curr[j] { 3 = MI } -> max = term3
              MAX2(term_3, s_mm_i_j, mi_vec, state_result);
              s_mm_i_j = simdf32_max(s_mm_i_j, term_3);

              // same as: S_curr[j-1] - 0.5 * par.mact,  //
              const simd_float term_4 = simdf32_sub(s_curr_j_min_1, mact_half);
              // if term4 > s_curr[j] { 4 = IM } -> max = term4
              MAX2(term_4, s_mm_i_j, im_vec, state_result);
              s_mm_i_j = simdf32_max(s_mm_i_j, term_4);

  //          if (eq && viterbi_matrix.getCellOff(i, j, elem_idx) == 0 && j <= 259) {
                  // check terms ->
  //              fprintf(stdout, "%2.20f,%2.20f,%2.20f,%2.20f\n", ((float *)&term_1)[elem_idx], ((float *)&term_2)[elem_idx], ((float *)&term_3)[elem_idx], ((float *)&term_4)[elem_idx]);
  //              fprintf(stdout, "%2.20f,%2.20f,%2.20f", s_)
                  // check score ->
  //              fprintf(stdout, "%2.20f\n", ((float *)&s_mm_i_j)[elem_idx]);
  //          }

              // Cell of logic
              // if (cell_off[i][j])
              //shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
              //because 10000000000000000000000000000000 = -2147483648 kills cmplt
  #ifdef AVX2_SUPPORT
              simd_int matrix_vec    = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
              matrix_vec             = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
  #else
              simd_int matrix_vec    = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
  #endif
              simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
              simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
              simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec, float_min_vec); // inverse
              s_mm_i_j = simdf32_add(s_mm_i_j, cell_off_float_min_vec);    // add the cell off vec to sMM_i_j. Set -FLT_MAX to cell off

              // Store current score value
              simdf32_store((float *)&m_s_curr[j], s_mm_i_j);

  //          if (eq) {
  //              fprintf(stdout, "%2.20f,%i\n", ((float *)&s_curr[j])[elem_idx], ((int *)&state_result)[elem_idx]);
  //          }

  //          if (eq && j <= 259) {
  //              fprintf(stdout, "%2.20f\n", ((float *)&s_mm_i_j)[elem_idx]);
  //              fprintf(stdout, "%2.20f\n", ((float *)&score_mac)[elem_idx]);
  //              fprintf(stdout, "state:%i\n", ((int*)&state_result)[elem_idx]);
  //          }

              // Set "state_result" to Viterbi matrix
              // write values back to ViterbiMatrix
  #ifdef AVX2_SUPPORT
              /* byte_result_vec        000H  000G  000F  000E   000D  000C  000B  000A */
              /* abcdefgh               0000  0000  HGFE  0000   0000  0000  0000  DCBA */
              const __m256i abcdefgh = _mm256_shuffle_epi8(state_result, shuffle_mask_extract);
              /* abcd                                            0000  0000  0000  DCBA */
              const __m128i abcd     = _mm256_castsi256_si128(abcdefgh);
              /* efgh                                            0000  0000  HGFE  0000 */
              const __m128i efgh     = _mm256_extracti128_si256(abcdefgh, 1);
              _mm_storel_epi64((__m128i*)&sCO_MI_DG_IM_GD_MM_vec[j], _mm_or_si128(abcd, efgh));
  #else
              state_result = _mm_packs_epi32(state_result, state_result);
              state_result = _mm_packus_epi16(state_result, state_result);
              int int_result = _mm_cvtsi128_si32(state_result);
              sCO_MI_DG_IM_GD_MM_vec[j] = int_result;
  #endif

              if (findMaxInnerLoop){
                  // new score is higer
                  // output
                  //  0   0   0   MAX
                  simd_int lookup_mask_hi = (simd_int) simdf32_gt(s_mm_i_j, score_mac);
                  // old score is higher
                  // output
                  //  MAX MAX MAX 0
                  simd_int lookup_mask_lo = (simd_int) simdf32_lt(s_mm_i_j, score_mac);

                  simd_int new_j_pos_hi = simdi_and(lookup_mask_hi, curr_pos_j);
                  simd_int old_j_pos_lo = simdi_and(lookup_mask_lo, j2_vec);
                  j2_vec = simdi32_add(new_j_pos_hi, old_j_pos_lo);
                  simd_int new_i_pos_hi = simdi_and(lookup_mask_hi, curr_pos_i);
                  simd_int old_i_pos_lo = simdi_and(lookup_mask_lo, i2_vec);
                  i2_vec = simdi32_add(new_i_pos_hi, old_i_pos_lo);

                  score_mac = simdf32_max(s_mm_i_j, score_mac);
              }


              // if global alignment: look for best cell in last column
            if (!m_local){

                // new score is higer
                  // output
                  //  0   0   0   MAX
                  simd_int lookup_mask_hi = (simd_int) simdf32_gt(s_mm_i_j, score_mac);
                  // old score is higher
                  // output
                  //  MAX MAX MAX 0
                  simd_int lookup_mask_lo = (simd_int) simdf32_lt(s_mm_i_j, score_mac);

                simd_int new_j_pos_hi = simdi_and(lookup_mask_hi, curr_pos_j);
                simd_int old_j_pos_lo = simdi_and(lookup_mask_lo, j2_vec);
                j2_vec = simdi32_add(new_j_pos_hi,old_j_pos_lo);
                simd_int new_i_pos_hi = simdi_and(lookup_mask_hi, curr_pos_i);
                simd_int old_i_pos_lo = simdi_and(lookup_mask_lo, i2_vec);
                i2_vec = simdi32_add(new_i_pos_hi, old_i_pos_lo);

                score_mac = simdf32_max(s_mm_i_j, score_mac);
            }

            // check score_mac
  //        if (eq && j <= 259) {
  //              fprintf(stdout, "%2.20f\n", ((float *)&s_mm_i_j)[elem_idx]);
  //              fprintf(stdout, "%2.20f\n", ((float *)&score_mac)[elem_idx]);
  //          fprintf(stdout, "state:%i\n", ((int*)&state_result)[elem_idx]);
  //        }

          }   // end for j

          // Swap pointers so that previous becomes current and vice versa
          swap(m_s_curr, m_s_prev);

      }   // end for i

  //  if (eq) {
  ////                fprintf(stdout, "%2.20f\n", ((float *)&s_mm_i_j)[elem_idx]);
  //      fprintf(stdout, "[i,j] values: [%i,%i]\n", ((int *)&i2_vec)[elem_idx], ((int *)&j2_vec)[elem_idx]);
  //  }

      // Set found i2 and j2 value to hit hit objects respectively
      for (int elem = 0; elem < (int)hit_vec.size(); elem++) {
          hit_vec.at(elem)->i2 = ((int *)&i2_vec)[elem];
          hit_vec.at(elem)->j2 = ((int *)&j2_vec)[elem];
      }

  }
