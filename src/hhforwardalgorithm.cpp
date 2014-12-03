/*
 * hhforwardalgorithm.C
 *
 *  Created on: 16/04/2014
 *      Author: Stefan Haunsberger
 *
 *	Copyright (c) 2014 -. All rights reserved.
 *
 *	Info:
 *		Compare two HMMs with Forward Algorithm in log-space with SIMD instructions
 *
 */


#include "hhposteriordecoder.h"
#include "simd.h"
#include <float.h>
#include <utility>

using std::swap;

void PosteriorDecoder::forwardAlgorithm(HMMSimd & q_hmm, HMMSimd & t_hmm,
		std::vector<Hit *> & hit_vec, PosteriorMatrix & p_mm,
		ViterbiMatrix & viterbi_matrix, float shift) {
    
    //	char * nam = "PF07714";
    //	int elem_idx = 0;
    //	bool eq = false;
    //	for (int elem = 0; elem < (int)hit_vec.size(); elem++) {
    //		if (!strcmp(hit_vec.at(elem)->name, nam)) {
    //			eq = true;
    //			elem_idx = elem;
    //			break;
    //		}
    //	}
    
	int i,j;      //query and template match state indices
	const unsigned int n_elem = hit_vec.size();
    
	// Initialize SIMD variable for previous probability
	simd_float mm_prev_j_1 = simdf32_set(-FLT_MAX);
	simd_float im_prev_j_1 = simdf32_set(-FLT_MAX);
	simd_float mi_prev_j_1 = simdf32_set(-FLT_MAX);
	simd_float gd_prev_j_1 = simdf32_set(-FLT_MAX);
	simd_float dg_prev_j_1 = simdf32_set(-FLT_MAX);
    
	simd_float mm_prev_j = simdf32_set(-FLT_MAX);
	simd_float im_prev_j = simdf32_set(-FLT_MAX);
	simd_float mi_prev_j = simdf32_set(-FLT_MAX);
	simd_float gd_prev_j = simdf32_set(-FLT_MAX);
	simd_float dg_prev_j = simdf32_set(-FLT_MAX);
    
	// Initialize SIMD variable for current probability
	simd_float mm_curr_j_1 = simdf32_set(-FLT_MAX);
	simd_float im_curr_j_1 = simdf32_set(-FLT_MAX);
	simd_float gd_curr_j_1 = simdf32_set(-FLT_MAX);
    
	simd_float mm_curr_j = simdf32_set(-FLT_MAX);
	simd_float im_curr_j = simdf32_set(-FLT_MAX);
	simd_float mi_curr_j = simdf32_set(-FLT_MAX);
	simd_float gd_curr_j = simdf32_set(-FLT_MAX);
	simd_float dg_curr_j = simdf32_set(-FLT_MAX);
    simd_float m_p_min = (m_local ? simdf32_set(0.0f) : simdf32_set(-FLT_MAX));
	simd_float p_match;						// vector for storing emission probability
	simd_float transition_prob;			// vector for transition probabilities
    
	simd_float * p_mm_row_ptr = NULL;
    
        simd_int * t_length = t_hmm.lengths; 
	simdi_store(m_t_lengths_le, simdi32_add(*t_length, simdi32_set(1)));	// add one because of le comparison
	simdi_store(m_t_lengths_ge, simdi32_sub(*t_length, simdi32_set(1)));	// add one because of ge comparison
    
	const simd_float score_offset_shift = simdf32_set(shift);
    
	const simd_int float_min_int_vec = (simd_int)simdf32_set(-FLT_MAX);
	const simd_float float_min_vec = simdf32_set(-FLT_MAX);
    
    
    
	// for global alignment
	simd_float * p_last_col = malloc_simd_float(q_hmm.L * sizeof(simd_float));	// store forward probability in last column
    
	simd_int j_vec = simdi32_set(0);
    
	const __m128i tmp_vec = _mm_set_epi32(0x40000000,0x00400000,0x00004000,0x00000040);//01000000010000000100000001000000
	simd_int matrix_vec;
    
#ifdef AVX2
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
    
        simdf32_store((float *)m_p_forward, float_min_vec);

	// Store probabilities to array-vector variables respectively
	simdf32_store((float *)&m_mm_prev[0], float_min_vec);
	simdf32_store((float *)&m_im_prev[0], float_min_vec);
	simdf32_store((float *)&m_gd_prev[0], float_min_vec);
	simdf32_store((float *)&m_mi_prev[0], float_min_vec);
	simdf32_store((float *)&m_dg_prev[0], float_min_vec);
    
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize m_xx_prev (representing i = 1) and posterior matrix at p_mm[1][j]
	i = 1;
	p_mm_row_ptr = p_mm.getRow(i);
    
	const unsigned int start_pos_tr_i = i * 7;
	// load query transitions M2I and I2I (just needed for this initialization)
	const simd_float q_tr_1_m2i = simdf32_load((float *) (q_hmm.tr+start_pos_tr_i+1));
	const simd_float q_tr_1_i2i = simdf32_load((float *) (q_hmm.tr+start_pos_tr_i));
    
#ifdef AVX2
	unsigned long long * sCO_MI_DG_IM_GD_MM_vec = (unsigned long long *) viterbi_matrix.getRow(i);
#else
	unsigned int   * sCO_MI_DG_IM_GD_MM_vec   = (unsigned int *) viterbi_matrix.getRow(i);
#endif
    
	for (j = 1; j <= t_hmm.L; ++j) {
		j_vec = simdi32_set(j);
		const unsigned int start_pos_tr_j_1 = (j-1) * 7;
        
#ifdef AVX2
		matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
		matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
		matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
#endif
		simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
		simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
		// Check if all vector elements are cell off
		if (simd_hmax((const float *)&res_eq_co_vec, n_elem) == 0) {
            
			// load template transitions
			const simd_float t_m2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 2));	// M2M
			const simd_float t_m2d = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 3));	// M2D
			const simd_float t_d2d = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 5));	// D2D
            
			// Load probabilities into vector variables respectively
			mm_prev_j_1 = simdf32_load((float *)&m_mm_prev[j-1]);
			im_prev_j_1 = simdf32_load((float *)&m_im_prev[j-1]);
			gd_prev_j_1 = simdf32_load((float *)&m_gd_prev[j-1]);
            
			// same as:
			//	F_MM_prev[j] = flog2(hit.ProbFwd(q.p[1],t.p[j])) + par.shift;
			p_match = simdf32_flog2(Viterbi::ScalarProd20Vec((simd_float *) q_hmm.p[i],(simd_float *) t_hmm.p[j]));
			p_match = simdf32_add(p_match, score_offset_shift);
			mm_prev_j = p_match;
            
			// same as: F_MI_prev[j] = F_DG_prev[j] = -FLT_MAX;
			mi_prev_j = dg_prev_j = float_min_vec;
            
			// same as:
			//	F_IM_prev[j] = flog2_sum_fpow2(
			//								 F_MM_prev[j-1] + q.tr[1][M2I] + t.tr[j-1][M2M],	// MM -> IM
			//								 F_IM_prev[j-1] + q.tr[1][I2I] + t.tr[j-1][M2M]	// IM -> IM
			//							 );
			im_prev_j = simd_flog2_sum_fpow2(
                                             simdf32_add(mm_prev_j_1, simdf32_add(q_tr_1_m2i, t_m2m)),	// MM -> IM
                                             simdf32_add(im_prev_j_1, simdf32_add(q_tr_1_i2i, t_m2m))	// IM -> IM
                                             );
			// same as:
			//	F_GD_prev[j] = flog2_sum_fpow2(
			//										F_MM_prev[j-1] + t.tr[j-1][M2D],	// DG -> MM
			//										F_GD_prev[j-1] + t.tr[j-1][D2D]		// GD -> GD
			//									);
			gd_prev_j = simd_flog2_sum_fpow2(
                                             simdf32_add(mm_prev_j_1, t_m2d),	// DG -> MM
                                             simdf32_add(gd_prev_j_1, t_d2d)		// GD -> GD
                                             );
            
			// Cell of logic
			// if (cell_off[i][j])
			//shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
			//because 10000000000000000000000000000000 = -2147483648 kills cmplt
            
			// Cell off just for selected elements
			simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec, float_min_int_vec); // inverse
			mm_prev_j = simdf32_add(mm_prev_j, cell_off_float_min_vec);    // add the cell off vec to sMM_i_j. Set -FLT_MAX to cell off
			im_prev_j = simdf32_add(im_prev_j, cell_off_float_min_vec);
			gd_prev_j = simdf32_add(gd_prev_j, cell_off_float_min_vec);
            
		} else {
			mm_prev_j = im_prev_j = gd_prev_j = mi_prev_j = dg_prev_j = float_min_vec;
		}	// end if celloff
        
		// Store probabilities to array-vectors variables respectively
		simdf32_store((float *)&m_mm_prev[j], mm_prev_j);
		simdf32_store((float *)&m_im_prev[j], im_prev_j);
		simdf32_store((float *)&m_gd_prev[j], gd_prev_j);
		simdf32_store((float *)&m_mi_prev[j], mi_prev_j);
		simdf32_store((float *)&m_dg_prev[j], dg_prev_j);
        
		// Set values to m_mm_prev as long as j < t.L or rather replace with -FLT_MAX, respectively
		simdf32_store((float *)&p_mm_row_ptr[j], mm_prev_j);
        
		// Add probability to p_forward
		simdf32_store((float *)m_p_forward, simd_flog2_sum_fpow2(*m_p_forward, mm_prev_j));
        
        //		if (eq) {
        //			float * mm = (float *)&mm_prev_j;
        //			float * mi = (float *)&mi_prev_j;
        //			float * dg = (float *)&dg_prev_j;
        //			float * im = (float *)&im_prev_j;
        //			float * gd = (float *)&gd_prev_j;
        //			fprintf(stdout, "1,co:%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", viterbi_matrix.getCellOff(i, j, elem_idx), mm[elem_idx], gd[elem_idx], im[elem_idx], dg[elem_idx], mi[elem_idx]);
        //		}
        
	}
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Loop through query positions i
    
	const simd_float t_m2i_jmin = simdf32_load((float *) (t_hmm.tr + m_jmin * 7 + 1));
	const simd_float t_i2i_jmin = simdf32_load((float *) (t_hmm.tr + m_jmin * 7));
	j_vec = simdi32_set(m_jmin);
    
	for (i = 2; i <= q_hmm.L; i++) {
		simd_float * p_mm_row_ptr = p_mm.getRow(i);	// pointer to current row of the posterior probability matrix
		const unsigned int start_pos_tr_i = i * 7;
		const unsigned int start_pos_tr_i_1 = (i - 1) * 7;
		simdf32_store((float *)&p_mm_row_ptr[0], float_min_vec);
		j = m_jmin;
		// Read vector for cell off logic
#ifdef AVX2
		unsigned long long * sCO_MI_DG_IM_GD_MM_vec = (unsigned long long *) viterbi_matrix.getRow(i);
#else
		unsigned int   * sCO_MI_DG_IM_GD_MM_vec   = (unsigned int *) viterbi_matrix.getRow(i);
#endif
        
		// Cell of logic
		// if (cell_off[i][j])
		//	Replace values by -FLT_MAX if cell_off
		//
		//shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
		//because 10000000000000000000000000000000 = -2147483648 kills cmplt
#ifdef AVX2
		matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
		matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
		matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
#endif
		simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
		simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
        
        
        // load query transitions
        const simd_float q_i2i = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i			));	// I2I
        const simd_float q_m2i = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i 	 + 1));	// M2I
        const simd_float q_m2m = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i_1 + 2));	// M2M
        const simd_float q_m2d = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i_1 + 3));	// M2D
        const simd_float q_d2m = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i_1 + 4));	// D2M
        const simd_float q_d2d = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i_1 + 5));	// D2D
        const simd_float q_i2m = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i_1 + 6));	// I2M
        
        // Load probabilities into vector variables respectively
        mm_prev_j = simdf32_load((float *)&m_mm_prev[j]);
        mi_prev_j = simdf32_load((float *)&m_mi_prev[j]);
        dg_prev_j = simdf32_load((float *)&m_dg_prev[j]);
        
		// Check if all vector elements are cell off
		if (simd_hmax((const float *)&res_eq_co_vec, n_elem) == 0) {
            
			// same as: F_MM_curr[j_min] = flog2(hit.ProbFwd(q.p[i],t.p[j_min])) + par.shift;
			mm_curr_j = simdf32_add(
                                    simdf32_flog2(Viterbi::ScalarProd20Vec((simd_float *) q_hmm.p[i],(simd_float *) t_hmm.p[j])),
                                    score_offset_shift
                                    );
			// Same as:
			//	F_MI_curr[j_min] = flog2_sum_fpow2(
			//												F_MM_prev[j_min] + q.tr[i-1][M2M] + t.tr[j_min][M2I], // MI -> MM
			//												F_MI_prev[j_min] + q.tr[i-1][M2M] + t.tr[j_min][I2I]  // MI -> MI
			//											);
			mi_curr_j = simd_flog2_sum_fpow2(
                                             simdf32_add(mm_prev_j, simdf32_add(q_m2m, t_m2i_jmin)),	// MI -> MM
                                             simdf32_add(mi_prev_j, simdf32_add(q_m2m, t_i2i_jmin))	// MI -> MI
                                             );
			// Same as:
			//	F_DG_curr[j_min] = flog2_sum_fpow2(
			//										 F_MM_prev[j_min] + q.tr[i-1][M2D], // DG -> MM
			//										 F_DG_prev[j_min] + q.tr[i-1][D2D]	// DG -> DG
			//									 );
			dg_curr_j = simd_flog2_sum_fpow2(
                                             simdf32_add(mm_prev_j, q_m2d),	// DG -> MM
                                             simdf32_add(dg_prev_j, q_d2d)		// DG -> DG
                                             );
            
			// same as: F_IM_curr[j_min] = F_GD_curr[j_min] = -FLT_MAX;
			im_curr_j = gd_curr_j = simdf32_set(-FLT_MAX);
            
			// Cell off just for selected elements
			simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec,float_min_int_vec); // inverse
			mm_curr_j = simdf32_add(mm_curr_j,cell_off_float_min_vec);    // add the cell off vec to sMM_i_j. Set -FLT_MAX to cell off
			mi_curr_j = simdf32_add(mi_curr_j,cell_off_float_min_vec);
			dg_curr_j = simdf32_add(dg_curr_j,cell_off_float_min_vec);
		} else {
            mm_curr_j = mi_curr_j = im_curr_j = dg_curr_j = gd_curr_j = float_min_vec;
		}	// end if celloff
        
		// Store probabilities to array-vectors variables respectively
		simdf32_store((float *)&m_mm_curr[j], mm_curr_j);
		simdf32_store((float *)&m_mi_curr[j], mi_curr_j);
		simdf32_store((float *)&m_im_curr[j], im_curr_j);
		simdf32_store((float *)&m_dg_curr[j], dg_curr_j);
		simdf32_store((float *)&m_gd_curr[j], gd_curr_j);
        
		// Set values to matrix where j <= t.L or rather replace with -FLT_MAX, respectively
		simdf32_store((float *)&p_mm_row_ptr[j], mm_curr_j);
        
        //		if (eq && viterbi_matrix.getCellOff(i, j, 2) == 0) {
        //		if (eq) {
        //			float * mm = (float *)&mm_curr_j;
        //			float * mi = (float *)&mi_curr_j;
        //			float * dg = (float *)&dg_curr_j;
        //			float * im = (float *)&im_curr_j;
        //			float * gd = (float *)&gd_curr_j;
        //			fprintf(stdout, "2,co:%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", viterbi_matrix.getCellOff(i, j, elem_idx), mm[elem_idx], gd[elem_idx], im[elem_idx], dg[elem_idx], mi[elem_idx]);
        //		}
        
        //		if (par.loc) {
        //			calcLocalPForward(j_vec, mm_curr_j);
        //		}
        
		// Add probability to p_forward
		simdf32_store((float *)m_p_forward, simd_flog2_sum_fpow2(*m_p_forward, mm_curr_j));
        
		///////////////////////////////////////////////////////////////////////////
		// Loop through template positions j
		for (j = m_jmin + 1; j <= t_hmm.L; j++) {
            
			j_vec = simdi32_set(j);
            
			const unsigned int start_pos_tr_j = j * 7;
			const unsigned int start_pos_tr_j_1 = (j-1) * 7;
            
#ifdef AVX2
			matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
			matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
			matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
#endif
			simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
			simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
            
			// Check if all vector elements are cell off
			if (simd_hmax((const float *)&res_eq_co_vec, n_elem) == 0) {
                
                // Load template transitions
                const simd_float t_i2i = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j			));	// I2I
                const simd_float t_m2i = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j 	 + 1));	// M2I
                const simd_float t_m2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 2));	// M2M
                const simd_float t_m2d = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 3));	// M2D
                const simd_float t_d2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 4));	// D2M
                const simd_float t_d2d = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 5));	// D2D
                const simd_float t_i2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j_1 + 6));	// I2M
                
                // Load probabilities into vector variables respectively
                mm_prev_j_1 = simdf32_load((float *)&m_mm_prev[j-1]);
                gd_prev_j_1 = simdf32_load((float *)&m_gd_prev[j-1]);
                im_prev_j_1 = simdf32_load((float *)&m_im_prev[j-1]);
                dg_prev_j_1 = simdf32_load((float *)&m_dg_prev[j-1]);
                mi_prev_j_1 = simdf32_load((float *)&m_mi_prev[j-1]);
                
                mm_prev_j = simdf32_load((float *)&m_mm_prev[j]);
                dg_prev_j = simdf32_load((float *)&m_dg_prev[j]);
                mi_prev_j = simdf32_load((float *)&m_mi_prev[j]);
                
                mm_curr_j_1 = simdf32_load((float *)&m_mm_curr[j-1]);
                gd_curr_j_1 = simdf32_load((float *)&m_gd_curr[j-1]);
                im_curr_j_1 = simdf32_load((float *)&m_im_curr[j-1]);
                
                
                // same as:
                // 	F_MM_curr[j] = flog2(hit.ProbFwd(q.p[i],t.p[j])) + par.shift +
                // 									flog2_sum_fpow2(
                //											m_p_min_scalar,
                //											F_MM_prev[j-1] + q.tr[i-1][M2M] + t.tr[j-1][M2M],	// BB -> MM (BB = Begin/Begin, for local alignment)
                //											F_GD_prev[j-1] + q.tr[i-1][M2M] + t.tr[j-1][D2M],	// GD -> MM
                //											F_IM_prev[j-1] + q.tr[i-1][I2M] + t.tr[j-1][M2M],	// IM -> MM
                //											F_DG_prev[j-1] + q.tr[i-1][D2M] + t.tr[j-1][M2M],	// DG -> MM
                //											F_MI_prev[j-1] + q.tr[i-1][M2M] + t.tr[j-1][I2M]	// MI -> MM
                //										);
                p_match = simdf32_add(
                                      simdf32_flog2(Viterbi::ScalarProd20Vec((simd_float *) q_hmm.p[i],(simd_float *) t_hmm.p[j])),
                                      score_offset_shift
                                      );
                transition_prob = simd_flog2_sum_fpow2(
                                                       m_p_min,
                                                       simdf32_add(mm_prev_j_1, simdf32_add(q_m2m, t_m2m)),	// BB -> MM (BB = Begin/Begin, for local alignment)
                                                       simdf32_add(gd_prev_j_1, simdf32_add(q_m2m, t_d2m)),	// GD -> MM
                                                       simdf32_add(im_prev_j_1, simdf32_add(q_i2m, t_m2m)),	// IM -> MM
                                                       simdf32_add(dg_prev_j_1, simdf32_add(q_d2m, t_m2m)),	// DG -> MM
                                                       simdf32_add(mi_prev_j_1, simdf32_add(q_m2m, t_i2m))		// MI -> MM
                                                       );
                mm_curr_j = simdf32_add(p_match, transition_prob);
                // same as:
                //	F_GD_curr[j] = flog2_sum_fpow2(
                //										F_MM_curr[j-1] + t.tr[j-1][M2D],		// GD -> MM
                //										F_GD_curr[j-1] + t.tr[j-1][D2D]			// GD -> GD
                //								  );
                gd_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(mm_curr_j_1, t_m2d),	// GD -> MM
                                                 simdf32_add(gd_curr_j_1, t_d2d)		// GD -> GD
                                                 );
                // same as:
                // F_IM_curr[j] = flog2_sum_fpow2(
                //										F_MM_curr[j-1] + q.tr[i][M2I] + t.tr[j-1][M2M],	// MM -> IM
                //										F_IM_curr[j-1] + q.tr[i][I2I] + t.tr[j-1][M2M]	// IM -> IM
                //								 );
                im_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(mm_curr_j_1, simdf32_add(q_m2i, t_m2m)),	// MM -> IM
                                                 simdf32_add(im_curr_j_1, simdf32_add(q_i2i, t_m2m))		// IM -> IM
                                                 );
                //	F_DG_curr[j] = flog2_sum_fpow2(
                //										F_MM_prev[j] + q.tr[i-1][M2D],	// DG -> MM
                //										F_DG_prev[j] + q.tr[i-1][D2D]		// DG -> DG
                //									);
                dg_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(mm_prev_j, q_m2d),	// DG -> MM
                                                 simdf32_add(dg_prev_j, q_d2d)		// DG -> DG
                                                 );
                //	F_MI_curr[j] = flog2_sum_fpow2(
                //										F_MM_prev[j] + q.tr[i-1][M2M] + t.tr[j][M2I],		// MI -> MM
                //										F_MI_prev[j] + q.tr[i-1][M2M] + t.tr[j][I2I]		// MI -> MI
                //	 							 );
                mi_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(mm_prev_j, simdf32_add(q_m2m, t_m2i)),	// MI -> MM
                                                 simdf32_add(mi_prev_j, simdf32_add(q_m2m, t_i2i))		// MI -> MI
                                                 );
                
                // Cell of logic
                // if (cell_off[i][j])
                //	Replace values by -FLT_MAX if cell_off
                //
                //shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
                //because 10000000000000000000000000000000 = -2147483648 kills cmplt
                
                // Cell off just for selected elements
                simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec,float_min_int_vec); // inverse
                mm_curr_j = simdf32_add(mm_curr_j,cell_off_float_min_vec);    // add the cell off vec to sMM_i_j. Set -FLT_MAX to cell off
                mi_curr_j = simdf32_add(mi_curr_j,cell_off_float_min_vec);
                dg_curr_j = simdf32_add(dg_curr_j,cell_off_float_min_vec);
                im_curr_j = simdf32_add(im_curr_j,cell_off_float_min_vec);
                gd_curr_j = simdf32_add(gd_curr_j,cell_off_float_min_vec);
                
			} else {
				mm_curr_j = mi_curr_j = im_curr_j = dg_curr_j = gd_curr_j = float_min_vec;
			}	// end if celloff
            
			// Store probabilities to array-vectors variables respectively
			simdf32_store((float *)&m_mm_curr[j], mm_curr_j);
			simdf32_store((float *)&m_mi_curr[j], mi_curr_j);
			simdf32_store((float *)&m_dg_curr[j], dg_curr_j);
			simdf32_store((float *)&m_im_curr[j], im_curr_j);
			simdf32_store((float *)&m_gd_curr[j], gd_curr_j);
            
			// Set values to matrix where j <= t.L or rather replace with -FLT_MAX, respectively
			simdf32_store((float *)&p_mm_row_ptr[j], mm_curr_j);
            
			// Add probability to p_forward
			simdf32_store((float *)m_p_forward, simd_flog2_sum_fpow2(*m_p_forward, mm_curr_j));
            
			if (!m_local) {
				setGlobalColumnPForward(p_last_col, j_vec, i, mm_curr_j);
			}
            
            //			if (eq && viterbi_matrix.getCellOff(i, j, 2) == 0) {
            //			if (eq) {
            //				float * mm = (float *)&mm_curr_j;
            //				float * mi = (float *)&mi_curr_j;
            //				float * dg = (float *)&dg_curr_j;
            //				float * im = (float *)&im_curr_j;
            //				float * gd = (float *)&gd_curr_j;
            //				fprintf(stdout, "3,co:%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", viterbi_matrix.getCellOff(i, j, 2), mm[2], gd[2], im[2], dg[2], mi[2]);
            //			}
		}	// end for j
        
		// Swap pointers so that previous becomes current and vice versa
		swap(m_mm_prev, m_mm_curr);
		swap(m_mi_prev, m_mi_curr);
		swap(m_im_prev, m_im_curr);
		swap(m_dg_prev, m_dg_curr);
		swap(m_gd_prev, m_gd_curr);
        
	}	// end for i
    
    
	if (!m_local) {		// global alignment
		// sum over last column
		simdf32_store((float *)m_p_forward, simdf32_set(-FLT_MAX));
        
		for (i = 1; i < q_hmm.L; i++) {
			simdf32_store((float *)m_p_forward,  simd_flog2_sum_fpow2(
                                                *m_p_forward,
                                                p_last_col[i]
                                                ));
		}
        
		// sum over last row
		simd_float values;
		const simd_float * p_mm_row_ptr = p_mm.getRow(q_hmm.L);
		for (j = 1; j <= t_hmm.L; j++) {
			j_vec = simdi32_set(j);
			const simd_float mask_lt = (simd_float)simdi32_lt(j_vec, *m_t_lengths_le);
			values = simdf32_or(simdf32_andnot(mask_lt, *m_p_forward),
                                simdf32_and(mask_lt, p_mm_row_ptr[j]));
			simdf32_store((float *)m_p_forward, simd_flog2_sum_fpow2(*m_p_forward, values));
		}
        
	}

	for (i = 1; i <= q_hmm.L; i++) {
		m_forward_profile[i] = simdf32_set(0);
		for (j = 1; j <= t_hmm.L; j++) {
			m_forward_profile[i] = simdf32_add(m_forward_profile[i], simdf32_fpow2(simdf32_sub(p_mm.getValue(i, j), *m_p_forward)));
		}
	}
    
	for (int elem = 0; elem < (int)hit_vec.size(); elem++) {
		hit_vec.at(elem)->Pforward = ((float *)m_p_forward)[elem];
	}

	free(p_last_col);
  }

void PosteriorDecoder::setGlobalColumnPForward(simd_float * column,	const simd_int & j_vec, const int i_count, const simd_float & values) {

	const simd_float mask_lt = (simd_float)simdi32_lt(j_vec, *m_t_lengths_le);
	// andnot: 	to keep values that are already in column
	// and:			to add the new values
	simd_float col_i_vec = simdf32_load((float *)&column[i_count]);
	col_i_vec = simdf32_or(simdf32_andnot(mask_lt, col_i_vec),
														 simdf32_and(mask_lt, values));
	simdf32_store((float *)&column[i_count], col_i_vec);

}

