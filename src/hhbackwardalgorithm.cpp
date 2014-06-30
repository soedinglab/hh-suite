/*
 * hhbackwardalgorithm.C
 *
 *  Created on: Apr 16, 2014
 *      Author: Stefan Haunsberger
 *
 *  Copyright (c) 2014 -. All rights reserved.
 *
 *	Info:
 *    Compare two HMMs with Backward Algorithm in log-space with SIMD instructions
 *
 */

#include "hhposteriordecoder.h"
#include "simd.h"
#include <float.h>
#include <utility>

using std::swap;

void PosteriorDecoder::backwardAlgorithm(HMMSimd & q_hmm, HMMSimd & t_hmm, std::vector<Hit *> & hit_vec, 
                                         PosteriorMatrix & p_mm,
                                         ViterbiMatrix & viterbi_matrix, float shift) {
    
	int i,j;      //query and template match state indices
	simdi_store(m_t_lengths_le, simdi32_add(*t_hmm.lengths, simdi32_set(1)));	// add one because of lt comparison
    
	const unsigned int n_elem = hit_vec.size();
    
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
    
	const simd_float score_offset_shift = simdf32_set(shift);
    
	simd_float p_match;						// vector for storing emission probability
	const simd_int float_min_int_vec = (simd_int)simdf32_set(-FLT_MAX);
	const simd_float float_min_vec = simdf32_set(-FLT_MAX);
    
	// Initialize SIMD variable for previous probability
	simd_float mm_prev_j_1 = simdf32_set(-FLT_MAX);
    
	simd_float mm_prev_j = simdf32_set(-FLT_MAX);
	simd_float mi_prev_j = simdf32_set(-FLT_MAX);
	simd_float dg_prev_j = simdf32_set(-FLT_MAX);
    
	// Initialize SIMD variable for current probability
	simd_float mm_curr_j_1 = simdf32_set(-FLT_MAX);
	simd_float im_curr_j_1 = simdf32_set(-FLT_MAX);
	simd_float mi_curr_j_1 = simdf32_set(-FLT_MAX);
	simd_float gd_curr_j_1 = simdf32_set(-FLT_MAX);
	simd_float dg_curr_j_1 = simdf32_set(-FLT_MAX);
    
	simd_float mm_curr_j = simdf32_set(-FLT_MAX);
	simd_float im_curr_j = simdf32_set(-FLT_MAX);
	simd_float mi_curr_j = simdf32_set(-FLT_MAX);
	simd_float gd_curr_j = simdf32_set(-FLT_MAX);
	simd_float dg_curr_j = simdf32_set(-FLT_MAX);
    
	simd_float p_mm_i_j = simdf32_set(-FLT_MAX);
	simd_float p_mm_i_j_1 = simdf32_set(-FLT_MAX);
    simd_float m_p_min = (m_local ? simdf32_set(0.0f) : simdf32_set(-FLT_MAX));

	for(int i = 1; i <= q_hmm.L; i++) {
		m_backward_profile[i] = simdf32_set(0);
	}
    
	simd_float * p_mm_row_ptr = NULL;
    
	// Variables for cell off logic
	const __m128i tmp_vec = _mm_set_epi32(0x40000000,0x00400000,0x00004000,0x00000040);//01000000010000000100000001000000
    
	simd_int matrix_vec;
    
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
    
    //	m_mm_prev[0] = m_mi_prev[0] = m_dg_prev[0] = simdf32_set(-FLT_MAX);
    
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize m_xx_prev (representing i = q.length) and posterior matrix at p_mm[q.L][j] = bottom row
	i = q_hmm.L;
#ifdef AVX2_SUPPORT
	unsigned long long * sCO_MI_DG_IM_GD_MM_vec = (unsigned long long *) viterbi_matrix.getRow(i);
#else
	unsigned int   * sCO_MI_DG_IM_GD_MM_vec   = (unsigned int *) viterbi_matrix.getRow(i);
#endif
	p_mm_row_ptr = p_mm.getRow(i);
	for (j = t_hmm.L; j >= 1; j--) {
		// Load probability in i,j
		p_mm_i_j = simdf32_load((float*)&p_mm_row_ptr[j]);
		// Only set 0.0f where j <= t.l respectively (otherwise -FLT_MAX)
		mm_prev_j = simdf32_setzero();
		// Set values to m_mm_prev as long as j < t.L or rather replace with -FLT_MAX, respectively
		p_mm_i_j = simdf32_sub(p_mm_i_j, *m_p_forward);
        
		// Cell of logic
		// if (cell_off[i][j])
		//shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
		//because 10000000000000000000000000000000 = -2147483648 kills cmplt
#ifdef AVX2_SUPPORT
		matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
		matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
		matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
#endif
		simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
		simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
		simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec, float_min_int_vec); // inverse
		mm_prev_j = simdf32_add(mm_prev_j, cell_off_float_min_vec);    // add the cell off vec. Set -FLT_MAX to cell off
		p_mm_i_j = simdf32_add(p_mm_i_j, cell_off_float_min_vec);    // add the cell off vec. Set -FLT_MAX to cell off
        
		// Store values to arrays respectively
		simdf32_store((float *)&m_mm_prev[j], mm_prev_j);
		simdf32_store((float *)&m_mi_prev[j], float_min_vec);
		simdf32_store((float *)&m_dg_prev[j], float_min_vec);
        
		simdf32_store((float *)&p_mm_row_ptr[j], p_mm_i_j);
        
        //		if (eq) {
        //			float * mm = (float *)&m_mm_prev[j];
        //			float * mi = (float *)&m_mi_prev[j];
        //			float * dg = (float *)&m_dg_prev[j];
        //			fprintf(stdout, "1,co:%i,%2.20f,%2.20f,%2.20f\n", viterbi_matrix.getCellOff(i, j, elem_idx), mm[elem_idx], dg[elem_idx], mi[elem_idx]);
        //		}
        
	}
    
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//	Backward probability calculation
	// 	Backward algorithm
	// 	Loop through query positions i

	for (i = q_hmm.L - 1; i >= 1; i--) {
		p_mm_row_ptr = p_mm.getRow(i);	// pointer to a row of the posterior probability matrix
		const unsigned int start_pos_tr_i = i * 7;	// start position for indexing query transitions
        
#ifdef AVX2_SUPPORT
        unsigned long long * sCO_MI_DG_IM_GD_MM_vec = (unsigned long long *) viterbi_matrix.getRow(i);
#else
        unsigned int   * sCO_MI_DG_IM_GD_MM_vec   = (unsigned int *) viterbi_matrix.getRow(i);
#endif
        
		j = t_hmm.L;
        
        // load query transitions
		const simd_float q_i2i = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i		));	// I2I
		const simd_float q_m2i = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i + 1));	// M2I
		const simd_float q_m2m = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i + 2));	// M2M
		const simd_float q_m2d = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i + 3));	// M2D
		const simd_float q_d2m = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i + 4));	// D2M
		const simd_float q_d2d = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i + 5));	// D2D
		const simd_float q_i2m = simdf32_load((float *) (q_hmm.tr + start_pos_tr_i + 6));	// I2M
        
		mm_curr_j = simdf32_setzero();
        
		p_mm_i_j = simdf32_load((float *)&p_mm_row_ptr[j]);
		p_mm_i_j = simdf32_sub(p_mm_i_j, *m_p_forward);
		// Cell of logic
		// if (cell_off[i][j])
		//shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
		//because 10000000000000000000000000000000 = -2147483648 kills cmplt
#ifdef AVX2_SUPPORT
		matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
		matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
		matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
#endif
		simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
		simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
		simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec, float_min_int_vec); // inverse
		p_mm_i_j = simdf32_add(p_mm_i_j,cell_off_float_min_vec);    // add the cell off vec. Set -FLT_MAX to cell off
		mm_curr_j = simdf32_add(mm_curr_j,cell_off_float_min_vec);    // add the cell off vec. Set -FLT_MAX to cell off
        
		// Store values to arrays respectively
		simdf32_store((float *)&m_mm_curr[j], mm_curr_j);
		simdf32_store((float *)&m_im_curr[j], float_min_vec);
		simdf32_store((float *)&m_mi_curr[j], float_min_vec);
		simdf32_store((float *)&m_dg_curr[j], float_min_vec);
		simdf32_store((float *)&m_gd_curr[j], float_min_vec);
        
		simdf32_store((float *)&p_mm_row_ptr[j], p_mm_i_j);
        
		///////////////////////////////////////////////////////////////////////////
		// Loop through template positions j
		for (j = t_hmm.L - 1; j >= m_jmin; j--) {
			// set counter vector
			const simd_int j_vec_prev = simdi32_set(j+1);
            
			const unsigned int start_pos_tr_j = j * 7;
            
			// Load values into current(j+1) vectors
			mm_curr_j_1 = simdf32_load((float *)&m_mm_curr[j+1]);
			mi_curr_j_1 = simdf32_load((float *)&m_mi_curr[j+1]);
			im_curr_j_1 = simdf32_load((float *)&m_im_curr[j+1]);
			dg_curr_j_1 = simdf32_load((float *)&m_dg_curr[j+1]);
			gd_curr_j_1 = simdf32_load((float *)&m_gd_curr[j+1]);
            
			// If j >= t.L replace value by 0.0f otherwise keep the value, respectively
			simd_float mask_gt = (simd_float)simdi32_gt(j_vec_prev, *m_t_lengths_ge);
			mm_curr_j_1 = simdf32_andnot(mask_gt, mm_curr_j_1);
            
			p_mm_i_j_1 = simdf32_load((float *)&p_mm_row_ptr[j+1]);
			// Cell of logic
			// if (cell_off[i][j])
			//shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
			//because 10000000000000000000000000000000 = -2147483648 kills cmplt
#ifdef AVX2_SUPPORT
			matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j+1]>>1);
			matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
			matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j+1]>>1);
#endif
			simd_int cell_off_vec  = simdi_and(matrix_vec,co_vec);
			simd_int res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
			simd_float  cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec, float_min_int_vec); // inverse
			mm_curr_j_1 = simdf32_add(mm_curr_j_1,cell_off_float_min_vec);  // add the cell off vec. Set -FLT_MAX to cell off
			im_curr_j_1 = simdf32_add(im_curr_j_1,cell_off_float_min_vec);  // add the cell off vec. Set -FLT_MAX to cell off
			mi_curr_j_1 = simdf32_add(mi_curr_j_1,cell_off_float_min_vec);  // add the cell off vec. Set -FLT_MAX to cell off
			dg_curr_j_1 = simdf32_add(dg_curr_j_1,cell_off_float_min_vec);  // add the cell off vec. Set -FLT_MAX to cell off
			gd_curr_j_1 = simdf32_add(gd_curr_j_1,cell_off_float_min_vec);  // add the cell off vec. Set -FLT_MAX to cell off
			p_mm_i_j_1  = simdf32_add(p_mm_i_j_1, cell_off_float_min_vec);  // add the cell off vec. Set -FLT_MAX to cell off
            
            
			// Store values to arrays respectively
			simdf32_store((float *)&m_mm_curr[j+1], mm_curr_j_1);
			simdf32_store((float *)&m_im_curr[j+1], im_curr_j_1);
			simdf32_store((float *)&m_mi_curr[j+1], mi_curr_j_1);
			simdf32_store((float *)&m_dg_curr[j+1], dg_curr_j_1);
			simdf32_store((float *)&m_gd_curr[j+1], gd_curr_j_1);
            
			// Store probability value to posterior matrix
			simdf32_store((float *)&p_mm_row_ptr[j+1], p_mm_i_j_1);
            
            //			if (eq && j <= ((int *)&t_hmm.lengths)[elem_idx]) {
            //			if (eq && j < 259) {
            ////				float * mm = (float *)&m_mm_curr[j+1];
            ////				float * mi = (float *)&m_mi_curr[j+1];
            ////				float * dg = (float *)&m_dg_curr[j+1];
            ////				float * im = (float *)&m_im_curr[j+1];
            ////				float * gd = (float *)&m_gd_curr[j+1];
            ////				fprintf(stdout, "2,co:%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", viterbi_matrix.getCellOff(i, j+1, elem_idx), mm[elem_idx], gd[elem_idx], im[elem_idx], dg[elem_idx], mi[elem_idx]);
            //				fprintf(stdout, "p_mm: %2.20f\n", ((float *)&p_mm_row_ptr[j+1])[elem_idx]);
            //			}
            
			// Cell of logic
			// if (cell_off[i][j])
			//shift   10000000100000001000000010000000 -> 01000000010000000100000001000000
			//because 10000000000000000000000000000000 = -2147483648 kills cmplt
#ifdef AVX2_SUPPORT
			matrix_vec = _mm256_set1_epi64x(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
			matrix_vec = _mm256_shuffle_epi8(matrix_vec,shuffle_mask_celloff);
#else
			matrix_vec = simdi32_set(sCO_MI_DG_IM_GD_MM_vec[j]>>1);
#endif
			cell_off_vec  = simdi_and(matrix_vec,co_vec);
			res_eq_co_vec = simdi32_lt(cell_off_vec,co_vec); // shift is because signed can't be checked here
            
			// Check if all vector elements are cell off
			if (simd_hmax((const float *)&res_eq_co_vec, n_elem) == 0) {
                
				// load template transitions
				const simd_float t_i2i = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j		));	// I2I
				const simd_float t_m2i = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j + 1));	// M2I
				const simd_float t_m2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j + 2));	// M2M
				const simd_float t_m2d = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j + 3));	// M2D
				const simd_float t_d2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j + 4));	// D2M
				const simd_float t_d2d = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j + 5));	// D2D
				const simd_float t_i2m = simdf32_load((float *) (t_hmm.tr + start_pos_tr_j + 6));	// I2M
                
				mm_prev_j_1 = simdf32_load((float *)&m_mm_prev[j+1]);
				dg_prev_j = simdf32_load((float *)&m_dg_prev[j]);
				mi_prev_j = simdf32_load((float *)&m_mi_prev[j]);
                
				p_match = simdf32_add(
                                      mm_prev_j_1, simdf32_add(
                                                               simdf32_flog2(Viterbi::ScalarProd20Vec((simd_float *) q_hmm.p[i+1],(simd_float *) t_hmm.p[j+1])),
                                                               score_offset_shift)
                                      );
				mm_curr_j = simd_flog2_sum_fpow2(
                                                 m_p_min,																							// MM -> EE (End/End, for local alignment)
                                                 simdf32_add(p_match,     simdf32_add(q_m2m, t_m2m)),	// MM -> MM
                                                 simdf32_add(gd_curr_j_1, 									   t_m2d),	// MM -> GD (q->tr[i][M2M] is already contained in GD->MM)
                                                 simdf32_add(im_curr_j_1, simdf32_add(q_m2i, t_m2m)),	// MM -> IM
                                                 simdf32_add(dg_prev_j, 						   q_m2d				),	// MM -> DG (t->tr[j][M2M] is already contained in DG->MM)
                                                 simdf32_add(mi_prev_j,   simdf32_add(q_m2m, t_m2i))	  // MM -> MI
                                                 );
				gd_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(p_match, 		 simdf32_add(q_m2m, t_d2m)),	// GD -> MM
                                                 simdf32_add(gd_curr_j_1, 									  t_d2d)		// DG -> DG
                                                 );
				im_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(p_match, 		 simdf32_add(q_i2m, t_m2m)),	// IM -> MM
                                                 simdf32_add(im_curr_j_1, simdf32_add(q_i2i, t_m2m))		// IM -> IM
                                                 );
				dg_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(p_match, 	 simdf32_add(q_d2m, t_m2m)),	// DG -> MM
                                                 simdf32_add(dg_prev_j, 						 q_d2d				)		// DG -> DG
                                                 );
				mi_curr_j = simd_flog2_sum_fpow2(
                                                 simdf32_add(p_match, 	 simdf32_add(q_m2m, t_i2m)),	// MI -> MM
                                                 simdf32_add(mi_prev_j, simdf32_add(q_m2m, t_i2i))		// MI -> MI
                                                 );
                
				cell_off_float_min_vec = (simd_float) simdi_andnot(res_eq_co_vec, float_min_int_vec); // inverse
				mm_curr_j = simdf32_add(mm_curr_j, cell_off_float_min_vec);    // add the cell off vec. Set -FLT_MAX to cell off
				mi_curr_j = simdf32_add(mi_curr_j, cell_off_float_min_vec);
				dg_curr_j = simdf32_add(dg_curr_j, cell_off_float_min_vec);
				im_curr_j = simdf32_add(im_curr_j, cell_off_float_min_vec);
				gd_curr_j = simdf32_add(gd_curr_j, cell_off_float_min_vec);
                
				p_mm_i_j = simdf32_load((float *)&p_mm_row_ptr[j]);
				p_mm_i_j = simdf32_sub(simdf32_add(p_mm_i_j, mm_curr_j), *m_p_forward);
				p_mm_i_j = simdf32_add(p_mm_i_j, cell_off_float_min_vec);
                
			} else {
				mm_curr_j = mi_curr_j = im_curr_j = dg_curr_j = gd_curr_j = p_mm_i_j = float_min_vec;
			} // end if celloff
            
			// Store values to arrays respectively
			simdf32_store((float *)&m_mm_curr[j], mm_curr_j);
			simdf32_store((float *)&m_im_curr[j], im_curr_j);
			simdf32_store((float *)&m_mi_curr[j], mi_curr_j);
			simdf32_store((float *)&m_dg_curr[j], dg_curr_j);
			simdf32_store((float *)&m_gd_curr[j], gd_curr_j);
            
			// Store probability value to posterior matrix
			simdf32_store((float *)&p_mm_row_ptr[j], p_mm_i_j);
            
            //			if (eq && j <= 259) {
            ////				float * mm = (float *)&m_mm_curr[j+1];
            ////				float * mi = (float *)&m_mi_curr[j+1];
            ////				float * dg = (float *)&m_dg_curr[j+1];
            ////				float * im = (float *)&m_im_curr[j+1];
            ////				float * gd = (float *)&m_gd_curr[j+1];
            ////				fprintf(stdout, "2,co:%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", viterbi_matrix.getCellOff(i, j+1, elem_idx), mm[elem_idx], gd[elem_idx], im[elem_idx], dg[elem_idx], mi[elem_idx]);
            //				fprintf(stdout, "p_mm: %2.20f\n", ((float *)&p_mm_row_ptr[j])[elem_idx]);
            //			}
            
			simd_float actual_backward = simdf32_mul(Viterbi::ScalarProd20Vec((simd_float *) q_hmm.p[i],(simd_float *) t_hmm.p[j]),
					simdf32_fpow2(simdf32_sub(simdf32_add(score_offset_shift, mm_curr_j), *m_p_forward)));
			m_backward_profile[i] = simdf32_add(m_backward_profile[i], actual_backward);
		}	// end for j
        
		// Swap pointers so that previous becomes current and vice versa
		swap(m_mm_prev, m_mm_curr);
		swap(m_mi_prev, m_mi_curr);
		swap(m_im_prev, m_im_curr);
		swap(m_dg_prev, m_dg_curr);
		swap(m_gd_prev, m_gd_curr);
        
	}	// end for i
}

