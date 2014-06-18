/*
 * hhposteriordecoder.h
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

#ifndef HHPOSTERIORDECODER_H_
#define HHPOSTERIORDECODER_H_
#include <map>
#include <vector>
#include <stddef.h>

#include "hhhmmsimd.h"
#include "hhviterbimatrix.h"
#include "hhposteriormatrix.h"
#include "hhviterbi.h"
using std::map;

class PosteriorDecoder {
public:

	static const int VEC_SIZE = HMMSimd::VEC_SIZE;

	struct MACBacktraceResult {
		std::vector<int> * alt_i;
		std::vector<int> * alt_j;
	};


	PosteriorDecoder(int maxres, bool local, int q_length);

	virtual ~PosteriorDecoder();

	/////////////////////////////////////////////////////////////////////////////////////
	// Realign hits: compute F/B/MAC and MAC-backtrace algorithms SIMD
	/////////////////////////////////////////////////////////////////////////////////////
	void realign(HMMSimd& q_hmm, HMMSimd& t_hmm, std::vector<Hit*>& hit_vec,
			PosteriorMatrix& p_mm, ViterbiMatrix& viterbi_matrix,
			std::vector<std::vector<MACBacktraceResult *> * > & alignment_exclusion_vec, int par_min_overlap, float shift, float mact, float corr);

private:

	simd_float * m_mm_prev;
	simd_float * m_gd_prev;
	simd_float * m_dg_prev;
	simd_float * m_im_prev;
	simd_float * m_mi_prev;

	simd_float * m_mm_curr;
	simd_float * m_gd_curr;
	simd_float * m_dg_curr;
	simd_float * m_im_curr;
	simd_float * m_mi_curr;

	simd_float * m_s_curr;		// MAC scores - current
	simd_float * m_s_prev;		// MAC scores - previous
	simd_float * p_last_col;

	simd_float * m_backward_profile;
	simd_float * m_forward_profile;


//	PosteriorSharedVariables m_column_vars;

	simd_float m_p_min;    // used to distinguish between SW and NW algorithms in maximization
	float m_p_min_scalar;    // used to distinguish between SW and NW algorithms in maximization

//	std::vector<Hit *> m_temp_hit_vec;	// temporary used hit objects for computation

	const int m_max_res;
	const bool m_local;				// local alignment
	const int m_q_length;			// query length

	int m_jmin;

	simd_float m_p_forward;
	simd_int m_t_lengths_le;	// for comparison le
	simd_int m_t_lengths_ge;	// for comparison ge

	std::vector<Hit *> m_temp_hit_vec;

	void forwardAlgorithm(HMMSimd & q_hmm, HMMSimd & t_hmm, std::vector<Hit *> & hit_vec, PosteriorMatrix & p_mm, ViterbiMatrix & viterbi_matrix, float shift);
	void backwardAlgorithm(HMMSimd & q_hmm, HMMSimd & t_hmm, std::vector<Hit *> & hit_vec, PosteriorMatrix & p_mm, ViterbiMatrix & viterbi_matrix, float shift);
	void macAlgorithm(HMMSimd & q_hmm, HMMSimd & t_hmm, std::vector<Hit *> & hit_vec, PosteriorMatrix & p_mm,
			ViterbiMatrix & viterbi_matrix, float par_mact);
	void backtraceMAC(HMM & q, HMM & t, PosteriorMatrix & p_mm, ViterbiMatrix & backtrace_matrix, const int elem, Hit & hit, float corr);
	void writeProfilesToHits(HMM & q, HMM & t, PosteriorMatrix & p_mm, const int elem, Hit & hit);
	void initializeBacktrace(HMM & t, Hit & hit);

	void initializeForAlignment(HMM & q, HMM & t, Hit * hit, ViterbiMatrix & viterbi_matrix, const int elem,
			std::vector<MACBacktraceResult *> & alignment_to_exclude,
			const int t_max_L, int par_min_overlap);
	void maskViterbiAlignment(const int q_length, const int t_length, ViterbiMatrix & celloff_matrix,
																									const int elem, const Hit * hit) const;
	void memorizeHitValues(Hit * curr_hit, const int i);
	void restoreHitValues(Hit & curr_hit, const int i);

	void setGlobalColumnPForward(simd_float * column, const simd_int & j_vec, const int i_count, const simd_float & values);

	void excludeMACAlignment(const int q_length, const int t_length, ViterbiMatrix &celloff_matrix, const int elem,
			MACBacktraceResult & alignment);

	void printVector(simd_float * vec);
	void printVector(simd_int * vec);
	void printVector(float * vec);

};

#endif /* HHPOSTERIORDECODER_H_ */
