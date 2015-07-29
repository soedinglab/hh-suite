/*
 * hhposteriordecoderunner.cpp
 *
 *  Created on: 22/03/2014
 *      Author: Stefan Haunsberger
 *
 *	Copyright (c) 2014 -. All rights reserved.
 *
 *	Info:
 *		This class inherits from the pthread class. It is used to read the hits from the input map and
 *		creates a worker queue item that contains: a hit-vector, template-HMM-vector and an alignments-to-exclude vector.
 *		This worker queue item is released into the workerQueue and signals waiting consumer threads
 *		to consume items. On the other side the consumer thread passes items through the second queue that are then reused.
 *
 */

#include "hhposteriordecoderrunner.h"

#include "util.h"

int compareIrep(const void * a, const void * b) {
    Hit pa = *(Hit*)a;
    Hit pb = *(Hit*)b;
    return (pa.irep < pb.irep);
}

PosteriorDecoderRunner::PosteriorDecoderRunner( PosteriorMatrix **posterior_matrices,
        ViterbiMatrix **backtrace_matrix, const int n_threads, const float ssw,
        const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF],
        const float S37[NSSPRED][MAXCF][NDSSP])
        : m_posterior_matrices(posterior_matrices),
          m_backtrace_matrix(backtrace_matrix),
          m_n_threads(n_threads), S73(S73), S33(S33), S37(S37) {}

PosteriorDecoderRunner::~PosteriorDecoderRunner() {
}

void PosteriorDecoderRunner::executeComputation(HMM &q, std::vector<Hit *>  hits, Parameters &par,
        const float qsc, float *pb, const float S[20][20], const float Sim[20][20], const float R[20][20]) {

    HMM * q_hmm = &q;    // Initialize single query HMM for PrepareTemplateHMM
    // algorithm performs in linear space
    q_hmm->Log2LinTransitionProbs(1.0);
    // Initialize selected query transitions
    initializeQueryHMMTransitions(*q_hmm);
    // prepeare data structure
    size_t target_max_length = 0;
    std::map<std::string, std::vector<Hit *> > alignments_map;
    for(size_t i = 0; i < hits.size(); i++){
        alignments_map[hits[i]->entry->getName()].push_back(hits[i]);
        target_max_length = std::max(target_max_length, (size_t) hits[i]->L);
    }
    // sort each std::vector<Hit *> by irep
    std::vector<std::vector<Hit *> > alignment;
    for (std::map<std::string, std::vector<Hit *> >::iterator alignment_vec =
            alignments_map.begin(); alignment_vec != alignments_map.end(); alignment_vec++) {
        std::sort(alignment_vec->second.begin(), alignment_vec->second.end(), compareIrep);
        alignment.push_back(alignment_vec->second);
    }

    // Routine to start consumer threads
    std::vector<PosteriorDecoder *> *threads = initializeConsumerThreads(par.loc, target_max_length, q.L, par.ssw, S73, S33, S37);
    // create one hmm for each threads
    HMM **t_hmm = new HMM *[m_n_threads];
    for (int i = 0; i < m_n_threads; i++) {
        t_hmm[i] = new HMM(MAXSEQDIS, par.maxres);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Iterate over all alignment vectors.
    // Each vector contains all alternative alignments for one Target
#pragma omp parallel for schedule(static)
    for (size_t idx = 0; idx < alignment.size(); idx++) {
        std::vector<Hit *> hits = alignment[idx];
        // find next free worker thread
        int current_thread_id = 0;
#ifdef OPENMP
     current_thread_id = omp_get_thread_num();
#endif
        PosteriorDecoder * decoder = threads->at(current_thread_id);
        std::vector<PosteriorDecoder::MACBacktraceResult> alignment_to_exclude;
        for(size_t idb = 0; idb < hits.size(); idb++){
            Hit *hit_cur = hits.at(idb);
            int format_tmp = 0;
            //char wg = 0;
            if(idb == 0){ // just read in the first time (less IO/CPU usage)
                hit_cur->entry->getTemplateHMM(par, par.wg, qsc, format_tmp, pb, S, Sim, t_hmm[current_thread_id]);
                PrepareTemplateHMM(par, q_hmm, t_hmm[current_thread_id], format_tmp, true, pb, R);
            }

            //TODO: par.ssw_realign not used???
            // start realignment process
            decoder->realign(*q_hmm, *t_hmm[current_thread_id],
                    *hit_cur, *m_posterior_matrices[current_thread_id],
                    *m_backtrace_matrix[current_thread_id], alignment_to_exclude, par.exclstr, par.min_overlap, par.shift, par.mact, par.corr);
            // add result to exclution paths (needed to align 2nd, 3rd, ... best alignment)
            alignment_to_exclude.push_back(PosteriorDecoder::MACBacktraceResult(hit_cur->alt_i, hit_cur->alt_j));
        } // end idb
        // remove all backtrace paths
        alignment_to_exclude.clear();
    }    // end - alignment vector

    for (int i = 0; i < m_n_threads; i++) {
        delete t_hmm[i];
    }
    delete[] t_hmm;

    cleanupThread(threads);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize consumer threads that process submitted items from shared worker queue
//	- Memory of the PosteriorMatrix is freed by the consumer threads respectively
///////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<PosteriorDecoder *> *PosteriorDecoderRunner::initializeConsumerThreads(char loc,
                                                size_t max_target_size, size_t query_size,
                                                const float ssw, const float S73[NDSSP][NSSPRED][MAXCF],
                                                const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF],
                                                const float S37[NSSPRED][MAXCF][NDSSP]) {
    std::vector<PosteriorDecoder *> *threads = new std::vector<PosteriorDecoder *>;
    for (int thread_id = 0; thread_id < m_n_threads; thread_id++) {
        PosteriorDecoder *thread = new PosteriorDecoder(max_target_size, loc, query_size, ssw, S73, S33, S37);
        threads->push_back(thread);
    }
    return threads;
}


void PosteriorDecoderRunner::cleanupThread(
        std::vector<PosteriorDecoder *> *threads) {
    for (size_t i = 0; i < threads->size(); i++) {
        delete threads->at(i);
    }
    threads->clear();
    delete threads;
}

void PosteriorDecoderRunner::initializeQueryHMMTransitions(HMM &q) {
    q.tr[0][M2D] = q.tr[0][M2I] = 0.0f;
    q.tr[0][I2M] = q.tr[0][I2I] = 0.0f;
    q.tr[0][D2M] = q.tr[0][D2D] = 0.0f;
    q.tr[q.L][M2M] = 1.0f;
    q.tr[q.L][M2D] = q.tr[q.L][M2I] = 0.0f;
    q.tr[q.L][I2M] = q.tr[q.L][I2I] = 0.0f;
    q.tr[q.L][D2M] = 1.0f;
    q.tr[q.L][D2D] = 0.0f;
}

