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
#include <float.h>

PosteriorDecoderRunner::PosteriorDecoderRunner(
    PosteriorDecoderRunnerInputData & data_container, HMMSimd & q_simd,
    PosteriorMatrix** posterior_matrices,
    ViterbiMatrix** backtrace_matrix, const int n_threads) :
        m_dbfiles(data_container.dbentries), m_alignments(
        data_container.alignments), m_q_simd(q_simd), m_n_t_hmms(
        data_container.n_t_hmms_to_align), m_t_maxres(data_container.t_maxres), m_posterior_matrices(
        posterior_matrices), m_backtrace_matrix(backtrace_matrix), m_n_threads(
        n_threads) {
}

PosteriorDecoderRunner::~PosteriorDecoderRunner() {
//	std::cout << "Destruct PosteriorDecoderConsumerThread" << std::endl;

}

void PosteriorDecoderRunner::executeComputation(Parameters& par, const float qsc, float* pb, const float S[20][20],
    const float Sim[20][20], const float R[20][20]) {

  HMM * q_hmm = m_q_simd.GetHMM(0);	// Initialize single query HMM for PrepareTemplateHMM

  int irep_counter = 0;

  // Alignments that are excluded for the next irep-run
  std::map<std::string, std::vector<PosteriorDecoder::MACBacktraceResult *>* > alignments_to_exclude;

  // Initialize selected query transitions
  initializeQueryHMMTransitions(*m_q_simd.GetHMM(0));

  // Routine to start consumer threads
  std::vector<PosteriorDecoder*> * threads = initializeConsumerThreads(par.loc);

  HMMSimd* t_hmm_simd[m_n_threads];
  for (int i = 0; i < m_n_threads; i++) {
    t_hmm_simd[i] = new HMMSimd(par.maxres);
  }

  HMM t_hmm[(HMMSimd::VEC_SIZE * m_n_threads)];


  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read all alignments
  //	Vector contains a map with the file-name as key and a hit object
  //		The index counter of the vector therefore is the irep number
  for (std::map<short int, std::vector<Hit *> >::iterator outer_map =
      m_alignments.begin(); outer_map != m_alignments.end(); outer_map++) {
    irep_counter++;
    // Vector contains hits
    std::vector<Hit *> & hits = outer_map->second;
    #pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int idb = 0; idb < hits.size(); idb +=
        HMMSimd::VEC_SIZE) {

      // find next free worker thread
      const int current_thread_id = omp_get_thread_num();
      const int current_t_index = (current_thread_id * HMMSimd::VEC_SIZE);

      std::vector<HMM *> templates_to_align;
      std::vector<Hit *> hit_items;
      std::vector<std::vector<PosteriorDecoder::MACBacktraceResult *> *> alignment_exclusions;

      // read in alignment
      int maxResElem = imin((hits.size()) - (idb), HMMSimd::VEC_SIZE);
      for (int i = 0; i < maxResElem; i++) {
        Hit* hit_cur = hits.at(idb+i);
        hit_items.push_back(hit_cur);

        int format_tmp = 0;
        char wg = 0;
        hit_cur->entry->getTemplateHMM(par, wg, qsc, format_tmp, pb, S, Sim, &t_hmm[current_t_index + i]);

        if (alignments_to_exclude.find(std::string(hit_cur->file)) != alignments_to_exclude.end()) {
          alignment_exclusions.push_back(alignments_to_exclude.at(std::string(hit_cur->file)));
        }
        else {
          std::vector< PosteriorDecoder::MACBacktraceResult *>* vec = new std::vector<PosteriorDecoder::MACBacktraceResult *>();
          alignment_exclusions.push_back(vec);
          alignments_to_exclude[std::string(hit_cur->file)] = vec;
        }


        PrepareTemplateHMM(par, q_hmm, &t_hmm[current_t_index + i], format_tmp, pb, R);
        templates_to_align.push_back(&t_hmm[current_t_index + i]);
      }
      t_hmm_simd[current_thread_id]->MapHMMVector(templates_to_align);
      // start next job
      threads->at(current_thread_id)->realign(m_q_simd, *t_hmm_simd[current_thread_id], hit_items,
          *m_posterior_matrices[current_thread_id], *m_backtrace_matrix[current_thread_id], alignment_exclusions, par.min_overlap, par.shift, par.mact, par.corr);
    } // idb loop

    mergeThreadResults(irep_counter, alignments_to_exclude);
  }	// end - outer map

  std::map<std::string, std::vector<PosteriorDecoder::MACBacktraceResult *>* >::iterator it;
  for(it = alignments_to_exclude.begin(); it != alignments_to_exclude.end(); it++) {
    for(size_t i = 0; i < (*it).second->size(); i++) {
      delete (*it).second->at(i);
    }
    (*it).second->clear();
    delete (*it).second;
  }
  alignments_to_exclude.clear();

  for (int i = 0; i < m_n_threads; i++) {
    delete t_hmm_simd[i];
  }

  cleanupThread(threads);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//	Add all alignments from the current 'irep' to the map where all alignments are stored that
//		have to be excluded in the next run (previous found MAC alignments).
///////////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorDecoderRunner::mergeThreadResults(int irep_counter,
    std::map<std::string, std::vector<PosteriorDecoder::MACBacktraceResult *> *> & alignments_to_exclude) {

  std::vector<Hit *> hits = m_alignments.at(irep_counter);
  for (int hit_elem = 0; hit_elem < (int) hits.size(); hit_elem++) {
    Hit & hit = *hits.at(hit_elem);
    // Initialize MACBacktraceResult object
    PosteriorDecoder::MACBacktraceResult * mac_btr =
        new PosteriorDecoder::MACBacktraceResult;
    mac_btr->alt_i = hit.alt_i;
    mac_btr->alt_j = hit.alt_j;

    // Add mac and viterbi backtrace result to alignments_to_exclude
    if(alignments_to_exclude.find(std::string(hit.file)) == alignments_to_exclude.end()) {
      alignments_to_exclude[std::string(hit.file)] = new std::vector<PosteriorDecoder::MACBacktraceResult *>();
    }
    alignments_to_exclude[std::string(hit.file)]->push_back(mac_btr);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize consumer threads that process submitted items from shared worker queue
//	- Memory of the PosteriorMatrix is freed by the consumer threads respectively
///////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<PosteriorDecoder*> * PosteriorDecoderRunner::initializeConsumerThreads(char loc) {

  std::vector<PosteriorDecoder *> * threads = new std::vector<PosteriorDecoder *>;
  for (int thread_id = 0; thread_id < m_n_threads; thread_id++) {
    PosteriorDecoder* thread = new PosteriorDecoder(m_t_maxres, loc, m_q_simd.L);
    threads->push_back(thread);
  }

  return threads;
}



void PosteriorDecoderRunner::cleanupThread(
    std::vector<PosteriorDecoder*> * threads) {

  for(size_t i = 0; i < threads->size(); i++) {
    delete threads->at(i);
  }
  threads->clear();

  delete threads;

}


void PosteriorDecoderRunner::initializeQueryHMMTransitions(HMM & q) {
  q.tr[0][M2D] = q.tr[0][M2I] = -FLT_MAX;
  q.tr[0][I2M] = q.tr[0][I2I] = -FLT_MAX;
  q.tr[0][D2M] = q.tr[0][D2D] = -FLT_MAX;
  q.tr[q.L][M2M] = 0.0;
  q.tr[q.L][M2D] = q.tr[q.L][M2I] = -FLT_MAX;
  q.tr[q.L][I2M] = q.tr[q.L][I2I] = -FLT_MAX;
  q.tr[q.L][D2M] = 0.0;
  q.tr[q.L][D2D] = -FLT_MAX;
}

