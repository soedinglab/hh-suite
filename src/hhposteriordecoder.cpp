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
 *			+ Forward (SIMD),
 *			+ Backward (SIMD),
 *			+ MAC (SIMD) and
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

	m_mm_prev = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_gd_prev = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_dg_prev = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_im_prev = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_mi_prev = malloc_simd_float(m_max_res * sizeof(simd_float));

	m_mm_curr = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_gd_curr = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_dg_curr = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_im_curr = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_mi_curr = malloc_simd_float(m_max_res * sizeof(simd_float));

	m_s_curr = malloc_simd_float(m_max_res * sizeof(simd_float));
	m_s_prev = malloc_simd_float(m_max_res * sizeof(simd_float));

	p_last_col = malloc_simd_float(q_length * sizeof(simd_float));

	m_p_min = (m_local ? simdf32_set(0.0f) : simdf32_set(-FLT_MAX));
	m_p_min_scalar = (m_local ? 0.0f : -FLT_MAX);

	m_p_forward = simdf32_set(-FLT_MAX);
	m_t_lengths_le = simdi32_set(0);
	m_t_lengths_ge = simdi32_set(0);

	for (int elem = 0; elem < VEC_SIZE; elem++) {
		m_temp_hit_vec.push_back(new Hit);
	}

}

PosteriorDecoder::~PosteriorDecoder() {
	// TODO Auto-generated destructor stub
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

	for (int elem = 0; elem < VEC_SIZE; elem++) {
//	  m_temp_hit_vec.at(elem)->Delete();
	  delete(m_temp_hit_vec.at(elem));
	}
	m_temp_hit_vec.clear();

}


/////////////////////////////////////////////////////////////////////////////////////
// Realign hits: compute F/B/MAC and MAC-backtrace algorithms
/////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::realign(HMMSimd & q_hmm, HMMSimd & t_hmm,
		std::vector<Hit*> & hit_vec, PosteriorMatrix & p_mm,
		ViterbiMatrix & viterbi_matrix,
		std::vector<std::vector<PosteriorDecoder::MACBacktraceResult *> * > & alignment_exclusion_vec, int par_min_overlap, float shift, float mact, float corr) {

	simd_int min_overlap;

	// For DEBUGGING
//	char * nam = "PF07714";
//	int elem_idx = 0;
//	bool eq = false;

	HMM & curr_q_hmm = *(q_hmm.GetHMM(0));
	int num_t = (int)hit_vec.size();
	// Iterate over vector elements
	for (int elem = 0; elem < num_t; elem++) {
		// DEBUGGING OPTION
//		if (!strcmp(hit_vec.at(elem)->name, nam)) {
//			eq = true;
//			elem_idx = elem;
//		}

		HMM * curr_t_hmm = t_hmm.GetHMM(elem);

		// Store some hit values temporary
		memorizeHitValues(hit_vec.at(elem), elem);

		// Clear current Viterbi path and also exclude previous found MAC alignments (from earlier alignments)
		initializeForAlignment(curr_q_hmm, *curr_t_hmm, hit_vec.at(elem), viterbi_matrix, elem,
				*alignment_exclusion_vec.at(elem), t_hmm.L, par_min_overlap);
	}

	// Compute SIMD Forward algorithm
	forwardAlgorithm(q_hmm, t_hmm, hit_vec, p_mm, viterbi_matrix, shift);

//	for (int elem = 0; elem < num_t; elem++) {
//		printf("%s: %20.20f\n", hit_vec.at(elem)->file, hit_vec.at(elem)->Pforward);
//	}

	// Compute SIMD Backward algorithm
	backwardAlgorithm(q_hmm, t_hmm, p_mm, viterbi_matrix, shift);


	// Compute SIMD MAC algorithm
	macAlgorithm(q_hmm, t_hmm, hit_vec, p_mm, viterbi_matrix, min_overlap, mact);

	for (int elem = 0; elem < num_t; elem++) {
//		HMM * curr_t_hmm = t_hmm.GetHMM(elem);
		// Perform scalar MAC backtrace
		backtraceMAC(curr_q_hmm, *t_hmm.GetHMM(elem), p_mm, viterbi_matrix, elem, *hit_vec.at(elem), corr);
		// Restore selected values
		restoreHitValues(*hit_vec.at(elem), elem);
	}

	// Print posterior matrix values of selected alignment
//	if (eq) {
//		for (int i = 0; i < q_hmm.L; i++) {
//			for (int j = 0; j < t_hmm.GetHMM(elem_idx)->L; j++) {
//				printf("%i,%i,%20.20f\n", i, j, p_mm.getSingleValue(i, j, elem_idx));
//			}
//		}
//	}

}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Prepare template and hit for the forthcoming alignment computation
//	- Eventually initialize some selected transition probabilities
//	- Exclude previously found alignments (including MAC and Viterbi)
//			--> Initialization of cell off matrix
/////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::initializeForAlignment(HMM & q, HMM & t, Hit* hit, ViterbiMatrix & celloff_matrix, const int elem,
		std::vector<PosteriorDecoder::MACBacktraceResult *> & alignment_to_exclude, const int t_max_L, int par_min_overlap) {

	// First alignment of this pair of HMMs?
	t.tr[0][M2M] = 0.0;
	t.tr[0][M2D] = t.tr[0][M2I] = -FLT_MAX;
	t.tr[0][I2M] = t.tr[0][I2I] = -FLT_MAX;
	t.tr[0][D2M] = t.tr[0][D2D] = -FLT_MAX;
	t.tr[t.L][M2M] = 0.0;
	t.tr[t.L][M2D] = t.tr[t.L][M2I] = -FLT_MAX;
	t.tr[t.L][I2M] = t.tr[t.L][I2I] = -FLT_MAX;
	t.tr[t.L][D2M] = 0.0;
	t.tr[t.L][D2D] = -FLT_MAX;
	//    if (alt_i && alt_i->Size()>0) delete alt_i;
	hit->alt_i = new std::vector<int>();
	//    if (alt_j && alt_j->Size()>0) delete alt_j;
	hit->alt_j = new std::vector<int>();

	hit->realign_around_viterbi = true;

	// Call Viterbi - InitializeForAlignment
	Viterbi::InitializeForAlignment(&q, &t, &celloff_matrix, elem, hit->self, par_min_overlap);

	for (int idx = 0; idx < (int) alignment_to_exclude.size(); idx++) {
		// Mask out previous found MAC alignments
		excludeMACAlignment(q.L, t.L, celloff_matrix, elem, *alignment_to_exclude.at(idx));
	}

	// Mask out the Viterbi alignment of the current hit
	if (hit->realign_around_viterbi) {
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
					ViterbiMatrix & celloff_matrix, const int elem, const Hit * hit) const {

	int i, j;

	// Cross out regions
	for (i = 1; i <= q_length; ++i)
		for (j = 1; j <= t_length; ++j)
			if (!((i < hit->i1 && j < hit->j1) || (i > hit->i2 && j > hit->j2)))
				celloff_matrix.setCellOff(i, j, elem, true);

//	maskViterbiAlignment(q_length, t_length, celloff_matrix, elem, alignment);
	// Clear Viterbi path
	for (int step = hit->nsteps; step >= 1; step--) {
		int path_width = 40;
		for (i = imax(1, hit->i[step] - path_width); i <= imin(q_length, hit->i[step] + path_width); ++i)
			celloff_matrix.setCellOff(i, hit->j[step], elem, false);

		for (j = imax(1, hit->j[step] - path_width); j <= imin(t_length, hit->j[step] + path_width); ++j)
			celloff_matrix.setCellOff(hit->i[step], j, elem, false);
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
		for(int q = 0; q < alignment.alt_i->size(); q++) {
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
void PosteriorDecoder::memorizeHitValues(Hit * curr_hit, const int i) {

	m_temp_hit_vec.at(i)->score      = curr_hit->score;
	m_temp_hit_vec.at(i)->score_ss   = curr_hit->score_ss;
	m_temp_hit_vec.at(i)->score_aass = curr_hit->score_aass;
	m_temp_hit_vec.at(i)->score_sort = curr_hit->score_sort;
	m_temp_hit_vec.at(i)->Pval       = curr_hit->Pval;
	m_temp_hit_vec.at(i)->Pvalt      = curr_hit->Pvalt;
	m_temp_hit_vec.at(i)->logPval    = curr_hit->logPval;
	m_temp_hit_vec.at(i)->logPvalt   = curr_hit->logPvalt;
	m_temp_hit_vec.at(i)->Eval       = curr_hit->Eval;
	m_temp_hit_vec.at(i)->logEval    = curr_hit->logEval;
	m_temp_hit_vec.at(i)->Probab     = curr_hit->Probab;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Restore the current hit with Viterbi scores, probabilities etc. of hit_cur
///////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoder::restoreHitValues(Hit& curr_hit, const int i) {
	curr_hit.score = m_temp_hit_vec.at(i)->score;
	curr_hit.score_ss = m_temp_hit_vec.at(i)->score_ss;
	curr_hit.score_aass = m_temp_hit_vec.at(i)->score_aass;
	curr_hit.score_sort = m_temp_hit_vec.at(i)->score_sort;
	curr_hit.Pval = m_temp_hit_vec.at(i)->Pval;
	curr_hit.Pvalt = m_temp_hit_vec.at(i)->Pvalt;
	curr_hit.logPval = m_temp_hit_vec.at(i)->logPval;
	curr_hit.logPvalt = m_temp_hit_vec.at(i)->logPvalt;
	curr_hit.Eval = m_temp_hit_vec.at(i)->Eval;
	curr_hit.logEval = m_temp_hit_vec.at(i)->logEval;
	curr_hit.Probab = m_temp_hit_vec.at(i)->Probab;
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

