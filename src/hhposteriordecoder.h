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
#include <algorithm>

#include "hhhmmsimd.h"
#include "hhviterbimatrix.h"
#include "hhposteriormatrix.h"
#include "hhviterbi.h"
using std::map;

struct MACTriple {
  int i;
  int j;
  float value;
};

bool compareIndices(const MACTriple &a, const MACTriple &b);

class PosteriorDecoder {
public:

	struct MACBacktraceResult {
		std::vector<int> * alt_i;
		std::vector<int> * alt_j;
		MACBacktraceResult(std::vector<int> * alt_i, std::vector<int> * alt_j):
				alt_i(alt_i), alt_j(alt_j) {};
	};

	PosteriorDecoder(int maxres, bool local, int q_length, const float ssw,
			         const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF],
					 const float S37[NSSPRED][MAXCF][NDSSP]);

	virtual ~PosteriorDecoder();

	/////////////////////////////////////////////////////////////////////////////////////
	// Realign hits: compute F/B/MAC and MAC-backtrace algorithms SIMD
	/////////////////////////////////////////////////////////////////////////////////////
	void realign(HMM &q, HMM &t, Hit &hit, PosteriorMatrix &p_mm, ViterbiMatrix &viterbi_matrix,
				 std::vector<PosteriorDecoder::MACBacktraceResult> alignment_to_exclude, char * exclstr,
				 char* template_exclstr, int par_min_overlap, float shift, float mact, float corr);
	void excludeMACAlignment(const int q_length, const int t_length, ViterbiMatrix &celloff_matrix, const int elem,
			MACBacktraceResult & alignment);

private:

	struct PosteriorMatrixCol {
		double mm;
		double gd;
		double im;
		double dg;
		double mi;
	};

	PosteriorMatrixCol * m_prev;
	PosteriorMatrixCol * m_curr;

	//	sec. structure data
	float ssw;
	//  SCORE_ALIGNMENT SCORE_BACKTRACE
	int ss_mode;
	//    float S73[NDSSP][NSSPRED][MAXCF];
	const float (*S73)[NSSPRED][MAXCF];
	//    float S33[NSSPRED][MAXCF][NSSPRED][MAXCF];
	const float (*S33)[MAXCF][NSSPRED][MAXCF];
	//    float S37[NSSPRED][MAXCF][NDSSP];
	const float (*S37)[MAXCF][NDSSP];


	double * m_s_curr;		// MAC scores - current
	double * m_s_prev;		// MAC scores - previous
	double * p_last_col;

	float m_back_forward_matrix_threshold;
	std::vector<MACTriple> m_backward_entries;
	std::vector<MACTriple> m_forward_entries;

	double * scale;

//	PosteriorSharedVariables m_column_vars;

//	simd_float m_p_min;    // used to distinguish between SW and NW algorithms in maximization

//	std::vector<Hit *> m_temp_hit;	// temporary used hit objects for computation

	const int m_max_res;
	const bool m_local;				// local alignment
	const int m_q_length;			// query length

	simd_float * m_p_forward;

	Hit * m_temp_hit;

	void forwardAlgorithm(HMM & q_hmm, HMM & t_hmm, Hit & hit_vec, PosteriorMatrix & p_mm,
			ViterbiMatrix & viterbi_matrix, float shift, const int elem);
	void backwardAlgorithm(HMM & q_hmm,HMM & t_hmm, Hit & hit_vec, PosteriorMatrix & p_mm,
			ViterbiMatrix & viterbi_matrix, float shift, const int elem);
	void macAlgorithm(HMM & q_hmm, HMM & t_hmm, Hit & hit_vec, PosteriorMatrix & p_mm,
			ViterbiMatrix & viterbi_matrix, float par_mact, const int elem);
	void backtraceMAC(HMM & q, HMM & t, PosteriorMatrix & p_mm, ViterbiMatrix & backtrace_matrix, const int elem, Hit & hit, float corr);
	void writeProfilesToHits(HMM &q, HMM &t, PosteriorMatrix &p_mm, Hit &hit);
	void initializeBacktrace(HMM & t, Hit & hit);

	void initializeForAlignment(HMM &q, HMM &t, Hit &hit, ViterbiMatrix &viterbi_matrix, const int elem, const int t_max_L, int par_min_overlap);
    void maskViterbiAlignment(const int q_length, const int t_length, ViterbiMatrix &celloff_matrix,
			const int elem, Hit const &hit) const;
	void memorizeHitValues(Hit & curr_hit);
	void restoreHitValues(Hit &curr_hit);

	void printVector(simd_float * vec);
	void printVector(simd_int * vec);
	void printVector(float * vec);

	void exclude_regions(char *exclstr, HMM &q_hmm, HMM &t_hmm, ViterbiMatrix &viterbiMatrix);
        void exclude_template_regions(char* exclstr, HMM & q_hmm, HMM & t_hmm, ViterbiMatrix& viterbiMatrix);
};

#endif /* HHPOSTERIORDECODER_H_ */
