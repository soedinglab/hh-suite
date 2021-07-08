#include "hhviterbirunner.h"

#ifdef OPENMP
#include <omp.h>
#endif

void ViterbiConsumerThread::clear() {
    hits.clear();
    excludeAlignments.clear();
}

void ViterbiConsumerThread::align(int maxres, int nseqdis, const float smin, const char ssm) {

    int consensus_ss_hmm_mode = 0xFF;
    for(size_t i = 0; i < maxres; i++){
        consensus_ss_hmm_mode &=  HMM::computeScoreSSMode(q_simd->GetHMM(0), t_hmm_simd->GetHMM(i));
    }
    // The following code solves the problem if more than 1 bit is set in "consensus_ss_hmm_mode".
    // It will pick the best possible mode
    int ss_hmm_mode = (consensus_ss_hmm_mode & HMM::PRED_DSSP);
    ss_hmm_mode = (ss_hmm_mode == 0) ? consensus_ss_hmm_mode & HMM::DSSP_PRED : 0;
    ss_hmm_mode = (ss_hmm_mode == 0) ? consensus_ss_hmm_mode & HMM::PRED_PRED : 0;

    Viterbi::ViterbiResult* viterbiResult = viterbiAlgo->Align(q_simd, t_hmm_simd, viterbiMatrix, maxres, ss_hmm_mode);
    for (int elem = 0; elem < maxres; elem++) {
        HMM * curr_t_hmm = t_hmm_simd->GetHMM(elem);
        HMM * curr_q_hmm = q_simd->GetHMM(elem);
        Viterbi::BacktraceResult backtraceResult = Viterbi::Backtrace(viterbiMatrix,
                                                                      elem, viterbiResult->i, viterbiResult->j);

        Viterbi::BacktraceScore backtraceScore = viterbiAlgo->ScoreForBacktrace(
                                                                                q_simd, t_hmm_simd, elem, &backtraceResult, viterbiResult->score, ss_hmm_mode);

        // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
        Hit hit_cur;
        hit_cur.lastrep = (backtraceScore.score <= smin) ? 1 : 0;

        hit_cur.initHitFromHMM(curr_q_hmm, curr_t_hmm, nseqdis, ssm);

        hit_cur.realign_around_viterbi = false;
        hit_cur.score = backtraceScore.score;
        hit_cur.score_ss = backtraceScore.score_ss;
        hit_cur.score_aass = backtraceScore.score_aass;
        hit_cur.S = backtraceScore.S;
        hit_cur.S_ss = backtraceScore.S_ss;

        hit_cur.i1 = backtraceResult.i_steps[backtraceResult.count];
        hit_cur.j1 = backtraceResult.j_steps[backtraceResult.count];

        hit_cur.i = backtraceResult.i_steps;

        hit_cur.j = backtraceResult.j_steps;
        hit_cur.nsteps = backtraceResult.count;

        hit_cur.states = backtraceResult.states;

        hit_cur.matched_cols = backtraceResult.matched_cols;

        hit_cur.i2 = viterbiResult->i[elem];
        hit_cur.j2 = viterbiResult->j[elem];

        hit_cur.entry = curr_t_hmm->entry;
        
//                        std::cout << "Thread: " << thread_id << std::endl;
//       HH_LOG(INFO) << string_format ("%d %-12.12s  %-12.12s   irep=%-2i  score=%6.2f ss_scor=%6.2f i=%d j=%d nstep=%d ssm_mode=%d t_ss_pred=%d t_ss_dssp=%d",elem, hit_cur.name,hit_cur.fam,hit_cur.irep,hit_cur.score, hit_cur.score_ss,viterbiResult->i[elem], viterbiResult->j[elem], hit_cur.nsteps, ss_hmm_mode, t_hmm_simd->GetHMM(elem)->nss_pred, t_hmm_simd->GetHMM(elem)->nss_dssp) << std::endl;
//                        printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit_cur.file,hit_cur.fam,hit_cur.irep,backtraceScore.score);
        hits.push_back(hit_cur); // insert hit at beginning of list (last repeats first!) Deep Copy of hit_cur
    }

    delete viterbiResult;
}



std::vector<Hit> ViterbiRunner::alignment(Parameters& par, HMMSimd * q_simd,
    std::vector<HHEntry*> dbfiles, const float qsc, float* pb,
    const float S[20][20], const float Sim[20][20], const float R[20][20], const int ssm_mode,
    const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF],
    const float S37[NSSPRED][MAXCF][NDSSP]) {

    HMM * q = q_simd->GetHMM(0);
    // Initialize memory
    std::vector<HMM*> t_hmm;
    for(size_t i = 0; i < VECSIZE_FLOAT * thread_count; i++) {
      HMM* t = new HMM(MAXSEQDIS, par.maxres);
      t_hmm.push_back(t);
    }

    HMMSimd** t_hmm_simd = new HMMSimd*[thread_count];
    std::vector<ViterbiConsumerThread *> threads;
    for (int thread_id = 0; thread_id < thread_count; thread_id++) {
        t_hmm_simd[thread_id] = new HMMSimd(par.maxres);
        ViterbiConsumerThread * thread = new ViterbiConsumerThread(thread_id, par, q_simd, t_hmm_simd[thread_id],viterbiMatrix[thread_id], ssm_mode, S73, S33, S37);
        threads.push_back(thread);
    }

    std::vector<Hit> ret_hits;
    std::vector<HHEntry*> dbfiles_to_align;
    std::map<std::string, std::vector<Viterbi::BacktraceResult> > excludeAlignments;
    // For all the databases comming through prefilter
    std::copy(dbfiles.begin(), dbfiles.end(), std::back_inserter(dbfiles_to_align));

    // loop to detect second/thrid/... best alignemtns
    for (int alignment = 0; alignment < par.altali; alignment++) {
        HH_LOG(INFO) << "Alternative alignment: " << alignment << std::endl;
        unsigned int allElementToAlignCount = dbfiles_to_align.size();
        unsigned int seqBlockSize = allElementToAlignCount;

        if(alignment == 0 && par.early_stopping_filter){
            seqBlockSize = 2000;
        }

        for(unsigned int seqJunkStart = 0; seqJunkStart <  allElementToAlignCount; seqJunkStart += seqBlockSize ){
            //sort by length to improve performance.
            //desc sort (for better utilisation ofthreads)
            unsigned int seqJunkSize = imin(allElementToAlignCount - (seqJunkStart), seqBlockSize);
            sort(dbfiles_to_align.begin() + seqJunkStart,
                 dbfiles_to_align.begin() + (seqJunkStart + seqJunkSize),
                 HHDatabaseEntryCompare());

            // read in data for thread
#pragma omp parallel for schedule(dynamic, 1) num_threads(thread_count)
            for (unsigned int idb = seqJunkStart; idb < (seqJunkStart + seqJunkSize); idb +=VECSIZE_FLOAT) {
                int current_thread_id = 0;
                #ifdef OPENMP
                    current_thread_id = omp_get_thread_num();
                #endif
                const int current_t_index = (current_thread_id *VECSIZE_FLOAT);

                std::vector<HMM *> templates_to_align;

                // read in alignment
                int maxResElem = imin((seqJunkStart + seqJunkSize) - (idb),
                                     VECSIZE_FLOAT);
                for (int i = 0; i < maxResElem; i++) {
                    HHEntry* entry = dbfiles_to_align.at(idb + i);
//                    if(strcmp(entry->getName(), "A0A075AHE7") == 0){
//                        i -= 1;
//                        std::cout << "##### ALIGN=" << entry->getName() << " i=" << i << " maxRes" << maxResElem << std::endl;
//
//                    }
                    int format_tmp = 0;
                    char wg = 1; // performance reason
                    entry->getTemplateHMM(par, wg, qsc, format_tmp, pb, S, Sim, t_hmm[current_t_index + i]);
                    t_hmm[current_t_index + i]->entry = entry;

                    PrepareTemplateHMM(par, q, t_hmm[current_t_index + i], format_tmp, false, pb, R);
                    templates_to_align.push_back(t_hmm[current_t_index + i]);

                }
                t_hmm_simd[current_thread_id]->MapHMMVector(templates_to_align);
                exclude_alignments(maxResElem, q_simd, t_hmm_simd[current_thread_id],
                                   excludeAlignments, viterbiMatrix[current_thread_id]);


                if(par.exclstr) {
                  // Mask excluded regions
                  exclude_regions(par.exclstr, maxResElem, q_simd, t_hmm_simd[current_thread_id], viterbiMatrix[current_thread_id]);
                }

                if(par.template_exclstr) {
                  // Mask excluded regions
                  exclude_template_regions(par.template_exclstr, maxResElem, q_simd, t_hmm_simd[current_thread_id], viterbiMatrix[current_thread_id]);
                }

                // start next job
                threads[current_thread_id]->align(maxResElem, par.nseqdis, par.smin, par.ssm);
            } // idb loop
            // merge thread results
            // search hits for next alignment
            HH_LOG(INFO) << (seqJunkStart + seqJunkSize) <<  " alignments done" << std::endl;

            merge_thread_results(ret_hits, dbfiles_to_align, excludeAlignments, threads, alignment, par.smin);
            for (unsigned int thread = 0; thread < threads.size(); thread++) {
                threads[thread]->clear();
            }

            if ( alignment == 0  && par.early_stopping_filter )
            {
                float early_stopping_sum = calculateEarlyStop(par, q, ret_hits, seqJunkStart);
                float filter_cutoff = seqJunkSize * par.filter_thresh;

                if( early_stopping_sum < filter_cutoff){
                    HH_LOG(INFO) << "Stop after DB-HHM: " << (seqJunkStart + seqJunkSize) << " because early stop  "
                    << early_stopping_sum << " < filter cutoff " << filter_cutoff << "\n";
                    break; // stop junk loop and just find alternative alignments
                }
            }
        } // junk loop
        // earse first elements. These are the elements from alignment run before,
        // new elements are after  + elementToAlignCount
        dbfiles_to_align.erase(dbfiles_to_align.begin(), dbfiles_to_align.begin() + allElementToAlignCount);

    }  // Alignment loop

    // clean memory
    for (int thread_id = 0; thread_id < thread_count; thread_id++) {
        delete t_hmm_simd[thread_id];
        delete threads[thread_id];
    }
    threads.clear();
    delete[] t_hmm_simd;

    for(size_t i = 0; i <VECSIZE_FLOAT * thread_count; i++) {
      delete t_hmm[i];
    }
    t_hmm.clear();

    return ret_hits;
}


float ViterbiRunner::calculateEarlyStop(Parameters& par, HMM * q, std::vector<Hit> &all_hits,
                                        unsigned int startPos){
    float early_stop_result = 0.0;
    for (unsigned int hit = startPos; hit < all_hits.size(); hit++) {
        Hit current_hit = all_hits[hit];
        float q_len = log(q->L) / LOG1000;
        float hit_len = log(current_hit.L) / LOG1000;
        float q_neff = q->Neff_HMM / 10.0;
        float hit_neff = current_hit.Neff_HMM / 10.0;
        float lamda = lamda_NN( q_len, hit_len, q_neff, hit_neff );
        float mu    =    mu_NN( q_len, hit_len, q_neff, hit_neff );
        current_hit.logPval = logPvalue(current_hit.score, lamda, mu);
        float alpha = 0;
        float log_Pcut = log(par.prefilter_evalue_thresh / par.dbsize);
        float log_dbsize = log(par.dbsize);

        if (par.prefilter)
            alpha = par.alphaa + par.alphab * (hit_neff - 1) * (1 - par.alphac * (q_neff - 1));

        current_hit.Eval = exp(current_hit.logPval + log_dbsize + (alpha * log_Pcut));
        current_hit.logEval = current_hit.logPval + log_dbsize + (alpha * log_Pcut);

        // Rolling average: replace oldest data point at par.filter_counter by newest one
        float eval_normalized = 1.0/(1.0+current_hit.Eval);
        early_stop_result += eval_normalized;

//        printf("%s Score %4.2f E-val: %4.2f eval_normalized: %4.2f  1/(1+Eval): %4.2f\n", current_hit.name, current_hit.score, current_hit.Eval, eval_normalized, early_stop_result);
//        printf("QLen: %4.2f q_neff: %4.2f  lambda: %4.2f mu: %4.2f\n",q_len, q_neff, lamda, mu);
//        printf("log_Pcut: %4.2f log_dbsize: %4.2f alpha: %4.2f   \n",log_Pcut, log_dbsize, alpha);



    }
    return early_stop_result;
}

void ViterbiRunner::merge_thread_results(std::vector<Hit> &all_hits,
                                         std::vector<HHEntry*> &dbfiles_to_align,
                                         std::map<std::string, std::vector<Viterbi::BacktraceResult> > &excludeAlignments,
                                         std::vector<ViterbiConsumerThread *> &threads, int alignment, const float smin) {
    for (unsigned int thread = 0; thread < threads.size(); thread++) {
        ViterbiConsumerThread * current_thread = threads[thread];
        for (unsigned int hit = 0; hit < current_thread->hits.size(); hit++) {
            Hit current_hit = current_thread->hits[hit];
            current_hit.irep = (alignment + 1);
            //        printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",current_hit.name,current_hit.fam,                                                                       current_hit.irep,current_hit.score);
            all_hits.push_back(current_hit);
            if (current_hit.score > smin) { // add to next alignmentif score for previous hit is better than SMIN
                dbfiles_to_align.push_back(current_hit.entry);
                Viterbi::BacktraceResult backtraceResult;
                backtraceResult.i_steps = current_hit.i;
                backtraceResult.j_steps = current_hit.j;
                backtraceResult.count = current_hit.nsteps;
                excludeAlignments[std::string(current_hit.entry->getName())].push_back(
                                                                                         backtraceResult);
            }
        }
    }
}

void ViterbiRunner::exclude_alignments(int maxResElem, HMMSimd* q_simd,
                                       HMMSimd* t_hmm_simd,
                                       std::map<std::string, std::vector<Viterbi::BacktraceResult> > &excludeAlignments,
                                       ViterbiMatrix* viterbiMatrix) {
    for (int elem = 0; elem < maxResElem; elem++) {
        HMM * curr_t_hmm = t_hmm_simd->GetHMM(elem);
        if (excludeAlignments.find(std::string(curr_t_hmm->entry->getName())) != excludeAlignments.end() && excludeAlignments[std::string(curr_t_hmm->entry->getName())].size() > 0) {
            std::vector<Viterbi::BacktraceResult> to_exclude = excludeAlignments[std::string(curr_t_hmm->entry->getName())];
            for (unsigned int i = 0; i < to_exclude.size(); i++) {
                Viterbi::BacktraceResult backtraceResult = to_exclude[i];
                Viterbi::ExcludeAlignment(viterbiMatrix, q_simd, t_hmm_simd, elem,
                                          backtraceResult.i_steps, backtraceResult.j_steps,
                                          backtraceResult.count);
            }
        }
    }
}

void ViterbiRunner::exclude_regions(char* exclstr, int maxResElem, HMMSimd* q_hmm_simd, HMMSimd* t_hmm_simd, ViterbiMatrix* viterbiMatrix) {
  char* ptr = exclstr;
  while (true) {
    int i0 = abs(strint(ptr));
    int i1 = abs(strint(ptr));

    if (!ptr) break;

    for (int elem = 0; elem < maxResElem; elem++) {
      HMM * curr_t_hmm = t_hmm_simd->GetHMM(elem);
      HMM * curr_q_hmm = q_hmm_simd->GetHMM(elem);
      for (int i = i0; i <= std::min(i1, curr_q_hmm->L); ++i) {
        for (int j=1; j <= curr_t_hmm->L; ++j) {
          viterbiMatrix->setCellOff(i, j, elem, true);
        }
      }
    }
  }
}

void ViterbiRunner::exclude_template_regions(char* exclstr, int maxResElem, HMMSimd* q_hmm_simd, HMMSimd* t_hmm_simd, ViterbiMatrix* viterbiMatrix) {
  char* ptr = exclstr;
  while (true) {
    int j0 = abs(strint(ptr));
    int j1 = abs(strint(ptr));

    if (!ptr) break;

    for (int elem = 0; elem < maxResElem; elem++) {
      HMM * curr_t_hmm = t_hmm_simd->GetHMM(elem);
      HMM * curr_q_hmm = q_hmm_simd->GetHMM(elem);

      for (int j=j0; j <= std::min(j1, curr_t_hmm->L); ++j) {
        for (int i = 1; i <= curr_q_hmm->L; ++i) {
          viterbiMatrix->setCellOff(i, j, elem, true);
        }
      }
    }
  }
}


