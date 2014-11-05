/*
 * hhprefilter.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: meiermark
 */

#include "hhprefilter.h"

namespace hh {

#define SWAP(tmp, arg1, arg2) tmp = arg1; arg1 = arg2; arg2 = tmp;

  Prefilter::Prefilter(const char* cs_library,
      FFindexDatabase* cs219_database) {
    // Prepare column state lib (context size =1 )
    FILE* fin = fopen(cs_library, "r");
    if (!fin)
      OpenFileError(cs_library, __FILE__, __LINE__, __func__);

    cs_lib = new cs::ContextLibrary<cs::AA>(fin);
    fclose(fin);

    cs::TransformToLin(*cs_lib);

    init_prefilter(cs219_database);
  }

  Prefilter::~Prefilter() {
    free(length);
    free(first);

    for (size_t n = 0; n < num_dbs; n++)
      delete[] dbnames[n];
    free(dbnames);

    delete cs_lib;
  }

  int Prefilter::swStripedByte(unsigned char *querySeq, int queryLength,
      unsigned char *dbSeq, int dbLength, unsigned short gapOpen,
      unsigned short gapExtend, simd_int *pvHLoad, simd_int *pvHStore,
      simd_int *pvE, unsigned short bias) {
    int i, j;

    unsigned int cmp;
    int element_count = (VECSIZE_INT * 4);

    int iter = (queryLength + (element_count - 1)) / element_count; // width of bands in query and score matrix = hochgerundetes LQ/16

    simd_int *pv;

    simd_int vE, vF, vH;

    simd_int vMaxScore;
    simd_int vBias;
    simd_int vGapOpen;
    simd_int vGapExtend;

    simd_int vTemp;
    simd_int vZero;

    simd_int *pvScore;

    simd_int *pvQueryProf = (simd_int*) querySeq;

    /* Load the bias to all elements of a constant */
    vBias = simdi8_set(bias);

    /* Load gap opening penalty to all elements of a constant */
    vGapOpen = simdi8_set(gapOpen);

    /* Load gap extension penalty to all elements of a constant */
    vGapExtend = simdi8_set(gapExtend);

    vMaxScore = simdi_setzero();
    vZero = simdi_setzero();

    /* Zero out the storage vector */
    for (i = 0; i < iter; ++i) {
      simdi_store(pvE + i, vMaxScore);
      simdi_store(pvHStore + i, vMaxScore);
    }

    for (i = 0; i < dbLength; ++i) {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;

      /* zero out F. */
      vF = simdi_setzero();

      /* load the next h value */
      vH = simdi_load(pvHStore + iter - 1);
      vH = simdi8_shiftl(vH, 1);

      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;

      for (j = 0; j < iter; ++j) {
        /* load values of vF and vH from previous row (one unit up) */
        vE = simdi_load(pvE + j);

        /* add score to vH */
        vH = simdui8_adds(vH, *(pvScore++));
        vH = simdui8_subs(vH, vBias);

        /* Update highest score encountered this far */
        vMaxScore = simdui8_max(vMaxScore, vH);

        /* get max from vH, vE and vF */
        vH = simdui8_max(vH, vE);
        vH = simdui8_max(vH, vF);

        /* save vH values */
        simdi_store(pvHStore + j, vH);

        /* update vE value */
        vH = simdui8_subs(vH, vGapOpen);
        vE = simdui8_subs(vE, vGapExtend);
        vE = simdui8_max(vE, vH);

        /* update vF value */
        vF = simdui8_subs(vF, vGapExtend);
        vF = simdui8_max(vF, vH);

        /* save vE values */
        simdi_store(pvE + j, vE);

        /* load the next h value */
        vH = simdi_load(pvHLoad + j);
      }

      /* reset pointers to the start of the saved data */
      j = 0;
      vH = simdi_load(pvHStore);

      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */
      vF = simdi8_shiftl(vF, 1);
      vTemp = simdui8_subs(vH, vGapOpen);
      vTemp = simdui8_subs(vF, vTemp);
      vTemp = simdi8_eq(vTemp, vZero);
      cmp = simdi8_movemask(vTemp);
#ifdef AVX2
      while (cmp != 0xffffffff)
#else
      while (cmp != 0xffff)
#endif
      {
        vE = simdi_load(pvE + j);

        vH = simdui8_max(vH, vF);

        /* save vH values */
        simdi_store(pvHStore + j, vH);

        /*  update vE incase the new vH value would change it */
        vH = simdui8_subs(vH, vGapOpen);
        vE = simdui8_max(vE, vH);
        simdi_store(pvE + j, vE);

        /* update vF value */
        vF = simdui8_subs(vF, vGapExtend);

        ++j;
        if (j >= iter) {
          j = 0;
          vF = simdi8_shiftl(vF, 1);
        }

        vH = simdi_load(pvHStore + j);

        vTemp = simdui8_subs(vH, vGapOpen);
        vTemp = simdui8_subs(vF, vTemp);
        vTemp = simdi8_eq(vTemp, vZero);
        cmp = simdi8_movemask(vTemp);
      }
    }

    /* find largest score in the vMaxScore vector */
    int score = simd_hmax((unsigned char *) &vMaxScore, element_count);

    /* return largest score */
    return score;
  }

  int Prefilter::ungapped_sse_score(const unsigned char* query_profile,
      const int query_length, const unsigned char* db_sequence,
      const int dbseq_length, const unsigned char score_offset,
      simd_int* workspace) {
    int i; // position in query bands (0,..,W-1)
    int j; // position in db sequence (0,..,dbseq_length-1)
    int element_count = (VECSIZE_INT * 4);
    const int W = (query_length + (element_count - 1)) / element_count; // width of bands in query and score matrix = hochgerundetes LQ/16

    simd_int *p;
    simd_int S;              // 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15)
    simd_int Smax = simdi_setzero();
    simd_int Soffset; // all scores in query profile are shifted up by Soffset to obtain pos values
    simd_int *s_prev, *s_curr; // pointers to Score(i-1,j-1) and Score(i,j), resp.
    simd_int *qji;             // query profile score in row j (for residue x_j)
    simd_int *s_prev_it, *s_curr_it;
    simd_int *query_profile_it = (simd_int *) query_profile;
    simd_int Zero = simdi_setzero();

    // Load the score offset to all 16 unsigned byte elements of Soffset
    Soffset = simdi8_set(score_offset);

    // Initialize  workspace to zero
    for (i = 0, p = workspace; i < 2 * W; ++i)
      simdi_store(p++, Zero);

    s_curr = workspace;
    s_prev = workspace + W;

    for (j = 0; j < dbseq_length; ++j) // loop over db sequence positions
        {

      // Get address of query scores for row j
      qji = query_profile_it + db_sequence[j] * W;

      // Load the next S value
      S = simdi_load(s_curr + W - 1);
      S = simdi8_shiftl(S, 1);

      // Swap s_prev and s_curr, smax_prev and smax_curr
      SWAP(p, s_prev, s_curr);

      s_curr_it = s_curr;
      s_prev_it = s_prev;

      for (i = 0; i < W; ++i) // loop over query band positions
          {
        // Saturated addition and subtraction to score S(i,j)
        S = simdui8_adds(S, *(qji++)); // S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset)
        S = simdui8_subs(S, Soffset);       // S(i,j) = max(0, S(i,j) - Soffset)
        simdi_store(s_curr_it++, S);       // store S to s_curr[i]
        Smax = simdui8_max(Smax, S);       // Smax(i,j) = max(Smax(i,j), S(i,j))

        // Load the next S and Smax values
        S = simdi_load(s_prev_it++);
      }
    }
    int score = simd_hmax((unsigned char *) &Smax, element_count);

    /* return largest score */
    return score;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
// Pull out all names from prefilter db file and copy into dbfiles_new for full HMM-HMM comparison
///////////////////////////////////////////////////////////////////////////////////////////////////
  void Prefilter::init_no_prefiltering(FFindexDatabase* cs219_database,
      std::vector<std::pair<int, std::string>>& prefiltered_entries) {
    ffindex_index_t* db_index = cs219_database->db_index;

    for (size_t n = 0; n < db_index->n_entries; n++) {
      ffindex_entry_t* entry = ffindex_get_entry_by_index(db_index, n);

      prefiltered_entries.push_back(
          std::make_pair<int, std::string>(entry->length,
              std::string(entry->name)));
    }

    HH_LOG(LogLevel::INFO) << "Searching " << prefiltered_entries.size()
        << " database HHMs without prefiltering" << std::endl;
  }

  void Prefilter::init_selected(FFindexDatabase* cs219_database,
      std::vector<std::string> templates,
      std::vector<std::pair<int, std::string>>& prefiltered_entries) {

    ffindex_index_t* db_index = cs219_database->db_index;

    for (size_t n = 0; n < templates.size(); n++) {
      ffindex_entry_t* entry = ffindex_get_entry_by_name(db_index, const_cast<char *>(templates[n].c_str()));

      prefiltered_entries.push_back(
          std::make_pair<int, std::string>(entry->length,
              std::string(entry->name)));
    }
  }

//////////////////////////////////////////////////////////////
// Reading in column state sequences for prefiltering
//////////////////////////////////////////////////////////////
  void Prefilter::init_prefilter(FFindexDatabase* cs219_database) {
    // Set up variables for prefiltering
    num_dbs = cs219_database->db_index->n_entries;
    first = (unsigned char**) mem_align(16, num_dbs * sizeof(unsigned char*));
    length = (int*) mem_align(16, num_dbs * sizeof(int));
    dbnames = (char**) mem_align(16, num_dbs * sizeof(char*));
    for (size_t n = 0; n < num_dbs; n++) {
      ffindex_entry_t* entry = ffindex_get_entry_by_index(
          cs219_database->db_index, n);
      first[n] = (unsigned char*) ffindex_get_data_by_entry(
          cs219_database->db_data, entry);
      length[n] = entry->length - 1;
      dbnames[n] = new char[strlen(entry->name) + 1];
      strcpy(dbnames[n], entry->name);
    }

    //check if cs219 format is new binary format
    checkCSFormat(5);

    HH_LOG(LogLevel::INFO) << "Searching " << num_dbs
        << " column state sequences." << std::endl;
  }

  void Prefilter::checkCSFormat(size_t nr_checks) {
    for (size_t n = 0; n < std::min(nr_checks, num_dbs); n++) {
      if (first[n][0] == '>') {
        nr_checks--;
      }
    }

    if (nr_checks == 0) {
      std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      std::cerr << "\tYour cs database is in an old format!" << std::endl;
      std::cerr << "\tThis format is no longer supportet!" << std::endl;
      std::cerr << "\tCorrespond to the user manual!" << std::endl;
      exit(1);
    }
  }

////////////////////////////////////////////////////////////////////////
// Prepare query profile for prefitering
////////////////////////////////////////////////////////////////////////
  void Prefilter::stripe_query_profile(HMM* q_tmp,
      const int prefilter_score_offset, const int prefilter_bit_factor,
      const int W, unsigned char* qc) {
    int LQ = q_tmp->L;
    float** query_profile = NULL;
    int a, h, i, j, k;

    // Build query profile with 219 column states
    query_profile = new float*[LQ + 1];
    for (i = 0; i < LQ + 1; ++i)
      query_profile[i] = (float*) malloc_simd_int(NUMCOLSTATES * sizeof(float));

    const cs::ContextLibrary<cs::AA>& lib = *cs_lib;

    // log (S(i,k)) = log ( SUM_a p(i,a) * p(k,a) / f(a) )   k: column state, i: pos in ali, a: amino acid
    for (i = 0; i < LQ; ++i)
      for (k = 0; k < NUMCOLSTATES; ++k) {
        float sum = 0;
        for (a = 0; a < 20; ++a)
          sum += ((q_tmp->p[i][a] * lib[k].probs[0][a]) / q_tmp->pav[a]);
        query_profile[i + 1][k] = sum;

      }

    /////////////////////////////////////////
    // Stripe query profile with chars
    int element_count = (VECSIZE_INT * 4);

    for (a = 0; a < NUMCOLSTATES; ++a) {
      h = a * W * element_count;
      for (i = 0; i < W; ++i) {
        j = i;
        for (k = 0; k < element_count; ++k) {
          if (j >= LQ)
            qc[h] = (unsigned char) prefilter_score_offset;
          else {
            float dummy = flog2(query_profile[j + 1][a])
                * prefilter_bit_factor + prefilter_score_offset + 0.5;
            if (dummy > 255.0)
              qc[h] = 255;
            else if (dummy < 0)
              qc[h] = 0;
            else
              qc[h] = (unsigned char) dummy; // 1/3 bits & make scores >=0 everywhere
          }
          ++h;
          j += W;
        }
      }
    }

    // Add extra ANY-state (220'th state)
    h = NUMCOLSTATES * W * element_count;
    for (i = 0; i < W; ++i) {
      j = i;
      for (k = 0; k < element_count; ++k) {
        if (j >= LQ)
          qc[h] = (unsigned char) prefilter_score_offset;
        else
          qc[h] = (unsigned char) (prefilter_score_offset - 1);
        h++;
        j += W;
      }
    }

    for (i = 0; i < LQ + 1; ++i)
      free(query_profile[i]);
    delete[] query_profile;
  }


////////////////////////////////////////////////////////////////////////
// Main prefilter function
////////////////////////////////////////////////////////////////////////
  void Prefilter::prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits,
      const int threads, const int prefilter_gap_open,
      const int prefilter_gap_extend, const int prefilter_score_offset,
      const int prefilter_bit_factor, const double prefilter_evalue_thresh,
      const double prefilter_evalue_coarse_thresh,
      const int preprefilter_smax_thresh, const int min_prefilter_hits, const int maxnumdb,
      const float R[20][20],
      std::vector<std::pair<int, std::string>>& new_prefilter_hits,
      std::vector<std::pair<int, std::string>>& old_prefilter_hits) {

    Hash<char>* doubled = new Hash<char>;
    doubled->New(16381, 0);

    int element_count = (VECSIZE_INT * 4);
    //W = (LQ+15) / 16;   // band width = hochgerundetes LQ/16
    int W = (q_tmp->L + (element_count - 1)) / element_count;
    // query profile (states + 1 because of ANY char)
    unsigned char* qc = (unsigned char*)malloc_simd_int((NUMCOLSTATES+1)*(q_tmp->L+element_count)*sizeof(unsigned char));
    stripe_query_profile(q_tmp, prefilter_score_offset, prefilter_bit_factor, W, qc);

    simd_int ** workspace = new simd_int *[threads];

    std::vector<std::pair<int, int> > first_prefilter;
    std::vector<std::pair<double, int> > hits;

    int count_dbs = 0;
    int gap_init = prefilter_gap_open + prefilter_gap_extend;
    int gap_extend = prefilter_gap_extend;
    int LQ = q_tmp->L;
    const float log_qlen = flog2(LQ);
    const double factor = (double) num_dbs * LQ;

    for (int i = 0; i < threads; i++)
      workspace[i] = (simd_int*) malloc_simd_int(
          3 * (LQ + element_count) * sizeof(char));

#pragma omp parallel for schedule(static)
    // Loop over all database sequences
    for (size_t n = 0; n < num_dbs; n++) {
      int thread_id = omp_get_thread_num();

      // Perform search step
      int score = ungapped_sse_score(qc, LQ, first[n], length[n],
          prefilter_score_offset, workspace[thread_id]);

      score = score
          - (int) (prefilter_bit_factor * (log_qlen + flog2(length[n])));

#pragma omp critical
      first_prefilter.push_back(std::pair<int, int>(score, n));
    }
    //filter after calculation of ungapped sse score to include at least min_prefilter_hits
    std::vector<std::pair<int, int> >::iterator it;

    sort(first_prefilter.begin(), first_prefilter.end());
    std::reverse(first_prefilter.begin(), first_prefilter.end());

    std::vector<std::pair<int, int> >::iterator first_prefilter_begin_erase =
        first_prefilter.end();
    std::vector<std::pair<int, int> >::iterator first_prefilter_end_erase =
        first_prefilter.end();
    count_dbs = 0;
    for (it = first_prefilter.begin(); it < first_prefilter.end(); it++) {
      if (count_dbs >= min_prefilter_hits
          && (*it).first <= preprefilter_smax_thresh) {
        first_prefilter_begin_erase = it;
        break;
      }
      else {
        count_dbs++;
      }
    }

    first_prefilter.erase(first_prefilter_begin_erase,
        first_prefilter_end_erase);

    HH_LOG(LogLevel::INFO)
        << "HMMs passed 1st prefilter (gapless profile-profile alignment)  : "
        << count_dbs << std::endl;

#pragma omp parallel for schedule(static)
    // Loop over all database sequences
//  for (int n = 0; n < count_dbs; n++) {
    for (size_t i = 0; i < first_prefilter.size(); i++) {
      int thread_id = omp_get_thread_num();

      int n = first_prefilter[i].second;

      // Perform search step
      int score = swStripedByte(qc, LQ, first[n], length[n], gap_init,
          gap_extend, workspace[thread_id], workspace[thread_id] + W,
          workspace[thread_id] + 2 * W, prefilter_score_offset);

      double evalue = factor * length[n] * fpow2(-score / prefilter_bit_factor);

      if (evalue < prefilter_evalue_coarse_thresh) {
#pragma omp critical
        hits.push_back(std::pair<double, int>(evalue, n));
      }
    }

    //filter after calculation of evalues to include at least min_prefilter_hits
    sort(hits.begin(), hits.end());

    std::vector<std::pair<double, int> >::iterator second_prefilter_begin_erase =
        hits.end();
    std::vector<std::pair<double, int> >::iterator second_prefilter_end_erase =
        hits.end();
    std::vector<std::pair<double, int> >::iterator it2;

    count_dbs = 0;
    for (it2 = hits.begin(); it2 < hits.end(); it2++) {
      if (count_dbs >= min_prefilter_hits
          && (*it2).first > prefilter_evalue_thresh) {
        second_prefilter_begin_erase = it2;
        break;
      }
      else {
        count_dbs++;
      }
    }

    hits.erase(second_prefilter_begin_erase, second_prefilter_end_erase);

    count_dbs = 0;

    for (it2 = hits.begin(); it2 < hits.end(); it2++) {
      // Add hit to dbfiles
      count_dbs++;
      char db_name[NAMELEN];
      strcpy(db_name, dbnames[(*it2).second]);

      char name[NAMELEN];
      RemoveExtension(name, db_name);

      if (!doubled->Contains(db_name)) {
        doubled->Add(db_name);

        std::pair<int, std::string> result;
        result.first = length[(*it2).second];
        result.second = std::string(db_name);

        // check, if DB was searched in previous rounds
        strcat(name, "__1");  // irep=1

        if (previous_hits->Contains(name)) {
          old_prefilter_hits.push_back(result);
        }
        else {
          new_prefilter_hits.push_back(result);
        }
      }
      if (count_dbs >= maxnumdb)
      {
        HH_LOG(LogLevel::WARNING)
        << "Number of hits passing 2nd prefilter (reduced from " << hits.size() << " to allowed maximum of " << maxnumdb << ").\n"
        <<"You can increase the allowed maximum using the -maxfilt <max> option.\n";
        break;
      }
    }


      
      
    // Free memory
    free(qc);
    for (int i = 0; i < threads; i++)
      free(workspace[i]);
    delete[] workspace;
    if (doubled)
      delete doubled;
  }

} /* namespace hh */
