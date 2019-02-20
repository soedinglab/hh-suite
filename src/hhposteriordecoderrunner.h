/*
 * hhposteriordecoderunner.h
 *
 *  Created on: 22/03/2014
 *      Author: Stefan Haunsberger
 */
#ifndef HHPOSTERIORDECODERRUNNER_H_
#define HHPOSTERIORDECODERRUNNER_H_

#include "hhposteriordecoder.h"
#include "hhposteriormatrix.h"
#include "hhviterbimatrix.h"
#include "hhfunc.h"

class PosteriorDecoderRunner {
public:
	PosteriorDecoderRunner(PosteriorMatrix **posterior_matrices, ViterbiMatrix **backtrace_matrix,
						   const int n_threadsconst, float ssw, const float S73[NDSSP][NSSPRED][MAXCF],
						   const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF],
						   const float S37[NSSPRED][MAXCF][NDSSP]);
	virtual ~PosteriorDecoderRunner();

	void executeComputation(HMM &q, std::vector<Hit *> hits, Parameters &par,
			const float qsc, float *pb, const float S[20][20], const float Sim[20][20], const float R[20][20]);

private:

	//    float S73[NDSSP][NSSPRED][MAXCF];
	const float (*S73)[NSSPRED][MAXCF];
	//    float S33[NSSPRED][MAXCF][NSSPRED][MAXCF];
	const float (*S33)[MAXCF][NSSPRED][MAXCF];
	//    float S37[NSSPRED][MAXCF][NDSSP];
	const float (*S37)[MAXCF][NDSSP];


	// map-key: irep of hits (1,2,...); vector: contains hits respectively to their irep value
	PosteriorMatrix** m_posterior_matrices;
	ViterbiMatrix ** m_backtrace_matrix;		// ViterbiMatrix used as backtrace and celloff matrix
	const int m_n_threads;			// Number of threads used to process m_worker_queue

	std::vector<PosteriorDecoder*> * initializeConsumerThreads(char loc, size_t max_target_size, size_t query_size,
			           										   const float ssw, const float S73[NDSSP][NSSPRED][MAXCF],
															   const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF],
															   const float S37[NSSPRED][MAXCF][NDSSP]);
	void initializeQueryHMMTransitions(HMM & q);
	void cleanupThread(std::vector<PosteriorDecoder*> * threads);
};

#endif /* HHPOSTERIORDECODERRUNNER_H_ */
