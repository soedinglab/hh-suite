/*
 * hhposteriordecoderunner.h
 *
 *  Created on: 22/03/2014
 *      Author: Stefan Haunsberger
 */

#include "hhposteriordecoderinputdata.h"
#include "hhposteriordecoder.h"
#include "hhposteriormatrix.h"
#include "hhviterbimatrix.h"
#include "hhfunc.h"

#ifndef HHPOSTERIORDECODERRUNNER_H_
#define HHPOSTERIORDECODERRUNNER_H_

class PosteriorDecoderRunner {
public:
	PosteriorDecoderRunner(PosteriorDecoderRunnerInputData &data_container, HMMSimd &q_vec,
			PosteriorMatrix** posterior_matrices, ViterbiMatrix** backtrace_matrix, const int n_threads);
	virtual ~PosteriorDecoderRunner();

	void executeComputation(Parameters& par, const float qsc, float* pb, const float S[20][20],
	    const float Sim[20][20], const float R[20][20]);

protected:
	void * run();

private:

	std::vector<HHblitsDatabase* > & m_databases;
	std::vector<HHDatabaseEntry* > & m_dbfiles;
	// map-key: irep of hits (1,2,...); vector: contains hits respectively to their irep value
	std::map<short int, std::vector<Hit *> > & m_alignments;
	HMMSimd & m_q_simd;				// Four equal query HMMs
	const int m_n_t_hmms;			// Number of hits that are realigned by F/B/MAC (all alignments)
	const int m_t_maxres;				// maximum template resolution (t.L + 2)
	PosteriorMatrix** m_posterior_matrices;
	ViterbiMatrix ** m_backtrace_matrix;		// ViterbiMatrix used as backtrace and celloff matrix
	const int m_n_threads;			// Number of threads used to process m_worker_queue

	std::vector<PosteriorDecoder*> * initializeConsumerThreads(char loc);
	void reinitializeSelectedMembers();
	void calcNElementsToSubmit();
	void initializeQueryHMMTransitions(HMM & q);
	void cleanupThread(std::vector<PosteriorDecoder*> * threads);
	void mergeThreadResults(int irep_counter,
			std::map<std::string, std::vector<PosteriorDecoder::MACBacktraceResult *> *> & alignments_to_exclude);

};

#endif /* HHPOSTERIORDECODERRUNNER_H_ */
