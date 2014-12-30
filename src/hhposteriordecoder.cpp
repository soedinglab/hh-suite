/*
 * hhposteriordecoder.C
 *
 *  Created on: 19/02/2014
 *      Author: Stefan Haunsberger
 *
 *  Copyright (c) 2014 -. All rights reserved.
 *
 *  Info:
 *		This class contains the needed algorithms to perform the
 *		MAC algorithm:
 *			+ Forward (scalar linear),
 *			+ Backward (scalar linear),
 *			+ MAC (scalar linear) and
 *			+ MAC backtrace (scalar).
 *
 *		The realign method is called by the posterior consumer thread.
 *		It prepares all needed matrices and parameters. This includes the
 *			- exclusion of previously found alignments
 *				(found by Viterbi or MAC) and the current Viterbi alignment
 *			- memorize selected hit values (Eval, Pval,...) that are restored
 *				after the posterior decoding, MAC and backtrace.
 *
 *		After the preparation step the algorithms in the above represented
 *		order are computed.
 *
 */

#include "hhposteriordecoder.h"
#include "util.h"
#include <float.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
// This method checks if the counter (j_vec) is lower than t.L+1 respectively (j <= t.L)
//	If this is the case the value is kept, if not it is replaced by -FLT_MAX and the new vector is returned.
//	Because there is no lower-equal t.L has been incremented by one (--> m_t_lengths_le = t.L + 1)
///////////////////////////////////////////////////////////////////////////////////////////////////
#define SetMinFloat_GT(j_vec, values, result) \
	min_float_vec = simdf32_set(-FLT_MAX);	\
	mask_lt = (simd_float)simdi32_lt(j_vec, m_t_lengths_le);		\
	result = simdf_or(simdf_and(mask_lt, values),																\
									simdf_andnot(mask_lt, min_float_vec));	// set to -FLT_MAX


PosteriorDecoder::PosteriorDecoder(int maxres, bool local, int q_length) :
	m_max_res(maxres),
	m_local(local),
	m_q_length(q_length) {

	m_jmin = 1;
	m_mm_prev = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_gd_prev = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_dg_prev = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_im_prev = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_mi_prev = (float*)  malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));

	m_mm_curr = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_gd_curr = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_dg_curr = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_im_curr = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_mi_curr = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));

	m_s_curr = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));
	m_s_prev = (float*) malloc_simd_float((m_max_res + 2 ) * sizeof(simd_float));

	p_last_col = (float*) malloc_simd_float(q_length * sizeof(simd_float));

	//m_p_min = (m_local ? simdf32_set(0.0f) : simdf32_set(-FLT_MAX));
	m_p_min_scalar = (m_local ? 1.0f : 0.0);

	m_p_forward = malloc_simd_float( sizeof(float));

	m_backward_profile = (float*) malloc_simd_float((q_length + 1) * sizeof(simd_float));
	m_forward_profile =  (float*) malloc_simd_float((q_length + 1) * sizeof(simd_float));
	scale = new double[m_max_res + 2];
	m_temp_hit = new Hit;
}

PosteriorDecoder::~PosteriorDecoder() {
	free(m_mm_prev);
	free(m_gd_prev);
	free(m_dg_prev);
	free(m_im_prev);
	free(m_mi_prev);

	free(m_mm_curr);
	free(m_gd_curr);
	free(m_dg_curr);
	free(m_im_curr);
	free(m_mi_curr);
	free(m_s_curr);
	free(m_s_prev);
	free(p_last_col);
	free(m_p_forward);

	free(m_backward_profile);
	free(m_forward_profile);
	delete [] scale;
	delete m_temp_hit;
}


/////////////////////////////////////////////////////////////////////////////////////
// Realign hits: compute F/B/MAC and MAC-backtrace algorithms
/////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::realign(HMM &q, HMM &t, Hit &hit,
		PosteriorMatrix &p_mm, ViterbiMatrix &viterbi_matrix,
		int par_min_overlap, float shift, float mact, float corr) {

	HMM & curr_q_hmm = q;
	HMM & curr_t_hmm = t;
	memorizeHitValues(hit);
	initializeForAlignment(curr_q_hmm, curr_t_hmm, hit, viterbi_matrix, 0, t.L, par_min_overlap);
	
	forwardAlgorithm(curr_q_hmm, curr_t_hmm, hit, p_mm, viterbi_matrix, shift, 0);
	//std::cout << hit->score << hit[elem]->Pforward << std::endl;

	backwardAlgorithm(curr_q_hmm, curr_t_hmm, hit, p_mm, viterbi_matrix, shift, 0);
	macAlgorithm(curr_q_hmm, curr_t_hmm, hit, p_mm, viterbi_matrix, mact, 0);
	backtraceMAC(curr_q_hmm, curr_t_hmm, p_mm, viterbi_matrix, 0, hit, corr);
	restoreHitValues(hit);
	writeProfilesToHits(curr_q_hmm, curr_t_hmm, p_mm, hit);

}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Prepare template and hit for the forthcoming alignment computation
//	- Eventually initialize some selected transition probabilities
//	- Exclude previously found alignments (including MAC and Viterbi)
//			--> Initialization of cell off matrix
/////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::initializeForAlignment(HMM &q, HMM &t, Hit &hit, ViterbiMatrix &celloff_matrix,
											const int elem, const int t_max_L, int par_min_overlap) {

	// First alignment of this pair of HMMs?
	t.tr[0][M2M] = 1.0f;
	t.tr[0][M2D] = t.tr[0][M2I] = 0.0f;
	t.tr[0][I2M] = t.tr[0][I2I] = 0.0f;
	t.tr[0][D2M] = t.tr[0][D2D] = 0.0f;
	t.tr[t.L][M2M] = 1.0f;
	t.tr[t.L][M2D] = t.tr[t.L][M2I] = 0.0f;
	t.tr[t.L][I2M] = t.tr[t.L][I2I] = 0.0f;
	t.tr[t.L][D2M] = 1.0f;
	t.tr[t.L][D2D] = 0.0f;
	//    if (alt_i && alt_i->Size()>0) delete alt_i;
	if(hit.alt_i) {
	  delete hit.alt_i;
	}
	hit.alt_i = new std::vector<int>();
	//    if (alt_j && alt_j->Size()>0) delete alt_j;
	if(hit.alt_j) {
		delete hit.alt_j;
	}
	hit.alt_j = new std::vector<int>();
	hit.realign_around_viterbi = true;

	// Call Viterbi - InitializeForAlignment
	Viterbi::InitializeForAlignment(&q, &t, &celloff_matrix, elem, hit.self, par_min_overlap);

	// Mask out the Viterbi alignment of the current hit
	if (hit.realign_around_viterbi) {
		maskViterbiAlignment(q.L, t.L, celloff_matrix, elem, hit);
	}

	// Mask out the outstanding matrix elements (t_hmm_vec_L - t_L)
	for (int i = 0; i <= q.L; i++) {
		for (int j = t.L + 1; j <= t_max_L; j++) {
			celloff_matrix.setCellOff(i, j, elem, true);
		}
	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Exclude Viterbi alignment of a hit
///////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::maskViterbiAlignment(const int q_length, const int t_length,
		ViterbiMatrix &celloff_matrix, const int elem, Hit const &hit) const {

	int i, j;

	// Switch off all cells (CellOff=true) except the upper left rectangle above (i1,j1)
	// and the cells in the lower right rectangle below (i2,j2), which are set to CellOff=false:
	//         j1              j2
	// 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1
	// 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1
	// 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 <-i1
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 <-i2
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
	// 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
	for (i = 1; i <= q_length; ++i)
		for (j = 1; j <= t_length; ++j)
			celloff_matrix.setCellOff(i, j, elem,   !((i < hit.i1 && j < hit.j1) || (i > hit.i2 && j > hit.j2))    );
	// Now switch on all cells (CellOff=false) in vicinity of the Viterbi path
	for (int step = hit.nsteps; step >= 1; step--) {
		for (i = imax(1, hit.i[step] - FWD_BKW_PATHWITDH); i <= imin(q_length, hit.i[step] + FWD_BKW_PATHWITDH); ++i)
			celloff_matrix.setCellOff(i, hit.j[step], elem, false);
	}
	for (int step = hit.nsteps; step >= 1; step--) {
		for (j = imax(1, hit.j[step] - FWD_BKW_PATHWITDH); j <= imin(t_length, hit.j[step] + FWD_BKW_PATHWITDH); ++j)
			celloff_matrix.setCellOff(hit.i[step], j, elem, false);
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Mask previous found alternative MAC alignments
///////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::excludeMACAlignment(const int q_length, const int t_length, ViterbiMatrix & celloff_matrix, const int elem,
											PosteriorDecoder::MACBacktraceResult & alignment) {

	int i,j;
	if (alignment.alt_i && alignment.alt_j) {
		alignment.alt_i->clear();
		alignment.alt_j->clear();
		for(size_t q = 0; q < alignment.alt_i->size(); q++) { //TODO: does not make sense
			i = alignment.alt_i->at(q);
			j = alignment.alt_j->at(q);

			for (int ii = imax(i - 2, 1); ii <= imin(i + 2, q_length); ++ii)
				celloff_matrix.setCellOff(ii, j, elem, true);
			for (int jj = imax(j - 2, 1); jj <= imin(j + 2, t_length); ++jj)
				celloff_matrix.setCellOff(i, jj, elem, true);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Memorize values that are going to be restored after computation
///////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::memorizeHitValues(Hit &curr_hit) {
	m_temp_hit->score      = curr_hit.score;
	m_temp_hit->score_ss   = curr_hit.score_ss;
	m_temp_hit->score_aass = curr_hit.score_aass;
	m_temp_hit->score_sort = curr_hit.score_sort;
	m_temp_hit->Pval       = curr_hit.Pval;
	m_temp_hit->Pvalt      = curr_hit.Pvalt;
	m_temp_hit->logPval    = curr_hit.logPval;
	m_temp_hit->logPvalt   = curr_hit.logPvalt;
	m_temp_hit->Eval       = curr_hit.Eval;
	m_temp_hit->logEval    = curr_hit.logEval;
	m_temp_hit->Probab     = curr_hit.Probab;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Restore the current hit with Viterbi scores, probabilities etc. of hit_cur
///////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::restoreHitValues(Hit &curr_hit) {
	curr_hit.score = m_temp_hit->score;
	curr_hit.score_ss = m_temp_hit->score_ss;
	curr_hit.score_aass = m_temp_hit->score_aass;
	curr_hit.score_sort = m_temp_hit->score_sort;
	curr_hit.Pval = m_temp_hit->Pval;
	curr_hit.Pvalt = m_temp_hit->Pvalt;
	curr_hit.logPval = m_temp_hit->logPval;
	curr_hit.logPvalt = m_temp_hit->logPvalt;
	curr_hit.Eval = m_temp_hit->Eval;
	curr_hit.logEval = m_temp_hit->logEval;
	curr_hit.Probab = m_temp_hit->Probab;
}

void PosteriorDecoder::printVector(float * vec) {

	for (int i = 0; i < VEC_SIZE; i++) {
//		printf("[%i],\t[%i],\t[%i],\t[%i]\n", (int*)&vec[i]);
		printf("[%10.5f],\t", vec[i]);
	}
	printf("\n");

}

void PosteriorDecoder::printVector(simd_float * vec) {

	for (int i = 0; i < VEC_SIZE; i++) {
//		printf("[%i],\t[%i],\t[%i],\t[%i]\n", (int*)&vec[i]);
		printf("[%10.5f],\t", ((float*)vec)[i]);
	}
	printf("\n");

}
void PosteriorDecoder::printVector(simd_int * vec) {

	for (int i = 0; i < VEC_SIZE; i++) {
//		printf("[%i],\t[%i],\t[%i],\t[%i]\n", (int*)&vec[i]);
		printf("[%i],\t", ((int*)vec)[i]);
	}
	printf("\n");

}

