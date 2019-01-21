/*
 * hhprefilter.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: meiermark
 */

#include "hhprefilter.h"
#include "ssw.h"

namespace hh {

namespace {

int ungapped_sse_score(simd_int* query_profile, const int query_length,
                       const char* db_sequence, const int dbseq_length,
                       const unsigned char score_offset, simd_int* workspace) {
  int element_count = (VECSIZE_INT * 4);
  const int W = (query_length + (element_count - 1)) /
                element_count;  // width of bands in query and score matrix =
                                // hochgerundetes LQ/16

  simd_int Smax = simdi_setzero();
  simd_int Zero = simdi_setzero();

  // All scores in query profile are shifted up by Soffset to obtain pos values.
  // Load the score offset to all 16 unsigned byte elements of Soffset.
  simd_int Soffset = simdi8_set(score_offset);

  // Initialize  workspace to zero
  simd_int* p = workspace;
  for (int i = 0; i < 2 * W; ++i) {
    simdi_store(p++, Zero);
  }

  // Pointers to Score(i-1,j-1) and Score(i,j), resp.
  simd_int* s_prev = workspace + W;
  simd_int* s_curr = workspace;

  for (int j = 0; j < dbseq_length; ++j)  // loop over db sequence positions
  {
    // Get address of query scores for row j (for residue x_j).
    simd_int* qji = query_profile + db_sequence[j] * W;

    // Load the next S value - 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15).
    simd_int S = simdi_load(s_curr + W - 1);
    S = simdi8_shiftl(S, 1);

    // Swap s_prev and s_curr, smax_prev and smax_curr
    std::swap(s_prev, s_curr);

    simd_int* s_curr_it = s_curr;
    simd_int* s_prev_it = s_prev;

    for (int i = 0; i < W; ++i) {  // loop over query band positions
      // Saturated addition and subtraction to score S(i,j)
      S = simdui8_adds(S,
                       *(qji++));  // S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset)
      S = simdui8_subs(S, Soffset);  // S(i,j) = max(0, S(i,j) - Soffset)
      simdi_store(s_curr_it++, S);   // store S to s_curr[i]
      Smax = simdui8_max(Smax, S);   // Smax(i,j) = max(Smax(i,j), S(i,j))

      // Load the next S and Smax values
      S = simdi_load(s_prev_it++);
    }
  }

  /* return largest score */
  return simd_hmax((unsigned char*)&Smax, element_count);
}

}  // namespace

///////////////////////////////////////////////////////////////////////////////////////////////////
// Pull out all names from prefilter db file and copy into dbfiles_new for full
// HMM-HMM comparison
///////////////////////////////////////////////////////////////////////////////////////////////////
void Prefilter::init_no_prefiltering(
    FFindexDatabase* query_database,
    std::vector<std::pair<int, std::string>>& prefiltered_entries) {
  ffindex_index_t* db_index = query_database->db_index;

  for (int n = 0; n < db_index->n_entries; n++) {
    ffindex_entry_t* entry = ffindex_get_entry_by_index(db_index, n);

    prefiltered_entries.emplace_back(std::make_pair<int, std::string>(
        entry->length, std::string(entry->name)));
  }

  HH_LOG(INFO) << "Searching " << prefiltered_entries.size()
               << " database HHMs without prefiltering" << std::endl;
}

void Prefilter::init_selected(
    FFindexDatabase* cs219_database, std::vector<std::string> templates,
    std::vector<std::pair<int, std::string> >& prefiltered_entries) {
  ffindex_index_t* db_index = cs219_database->db_index;

  for (const auto& template_it : templates) {
    ffindex_entry_t* entry = ffindex_get_entry_by_name(
        db_index, const_cast<char*>(template_it.c_str()));

    prefiltered_entries.emplace_back(std::make_pair<int, std::string>(
        entry->length, std::string(entry->name)));
  }
}

//////////////////////////////////////////////////////////////
// Reading in column state sequences for prefiltering
//////////////////////////////////////////////////////////////
void Prefilter::init_prefilter(FFindexDatabase* cs219_database) {
  // Set up variables for prefiltering
  for (int n = 0; n < cs219_database->db_index->n_entries; n++) {
    ffindex_entry_t* entry =
        ffindex_get_entry_by_index(cs219_database->db_index, n);
    DBEntry db_entry;
    db_entry.first = ffindex_get_data_by_entry(cs219_database->db_data, entry);
    db_entry.length = entry->length - 1;
    db_entry.name = entry->name;
    dbs_.push_back(db_entry);
  }

  // check if cs219 format is new binary format
  checkCSFormat(5);

  HH_LOG(INFO) << "Searching " << dbs_.size() << " column state sequences."
               << std::endl;
}

void Prefilter::checkCSFormat(size_t nr_checks) {
  for (int n = 0; n < std::min(nr_checks, dbs_.size()); n++) {
    if (dbs_[n].first[0] == '>') {
      nr_checks--;
    }
  }

  if (nr_checks == 0) {
    HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__
                  << ":" << std::endl;
    HH_LOG(ERROR) << "\tYour cs database is in an old format!" << std::endl;
    HH_LOG(ERROR) << "\tThis format is no longer supportet!" << std::endl;
    HH_LOG(ERROR) << "\tCorrespond to the user manual!" << std::endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////
// Prepare query profile for prefitering
////////////////////////////////////////////////////////////////////////
void Prefilter::stripe_query_profile(HMM* q_tmp,
                                     const int prefilter_score_offset,
                                     const int prefilter_bit_factor,
                                     const int W, unsigned char* qc) {
  // Build query profile with 219 column states
  int LQ = q_tmp->L;
  std::vector<float*> query_profile(LQ + 1);
  for (int i = 0; i < LQ + 1; ++i)
    query_profile[i] = (float*)malloc_simd_int(NUMCOLSTATES * sizeof(float));

  const cs::ContextLibrary<cs::AA>& lib = *cs_lib_;

  // log (S(i,k)) = log ( SUM_a p(i,a) * p(k,a) / f(a) )   k: column state, i:
  // pos in ali, a: amino acid
  for (int i = 0; i < LQ; ++i)
    for (int k = 0; k < NUMCOLSTATES; ++k) {
      float sum = 0;
      for (int a = 0; a < 20; ++a)
        sum += ((q_tmp->p[i][a] * lib[k].probs[0][a]) / q_tmp->pav[a]);
      query_profile[i + 1][k] = sum;
    }

  /////////////////////////////////////////
  // Stripe query profile with chars
  int element_count = (VECSIZE_INT * 4);

  for (int a = 0; a < NUMCOLSTATES; ++a) {
    int h = a * W * element_count;
    for (int i = 0; i < W; ++i) {
      int j = i;
      for (int k = 0; k < element_count; ++k) {
        if (j >= LQ)
          qc[h] = static_cast<unsigned char>(prefilter_score_offset);
        else {
          float dummy = flog2(query_profile[j + 1][a]) * prefilter_bit_factor +
                        prefilter_score_offset + 0.5;
          if (dummy > 255.0) {
            qc[h] = 255;
          } else if (dummy < 0) {
            qc[h] = 0;
          } else {
            // 1/3 bits & make scores >=0 everywhere
            qc[h] = static_cast<unsigned char>(dummy);
          }
        }
        ++h;
        j += W;
      }
    }
  }

  // Add extra ANY-state (220'th state)
  int h = NUMCOLSTATES * W * element_count;
  for (int i = 0; i < W; ++i) {
    int j = i;
    for (int k = 0; k < element_count; ++k) {
      if (j >= LQ)
        qc[h] = static_cast<unsigned char>(prefilter_score_offset);
      else
        qc[h] = static_cast<unsigned char>(prefilter_score_offset - 1);
      h++;
      j += W;
    }
  }

  for (int i = 0; i < LQ + 1; ++i) {
    free(query_profile[i]);
  }
}

Prefilter::Prefilter(const char* cs_library, FFindexDatabase* cs219_database) {
  // Prepare column state lib (context size =1 )
  FILE* fin = fopen(cs_library, "r");
  if (!fin) OpenFileError(cs_library, __FILE__, __LINE__, __func__);

  cs_lib_ = std::make_unique<cs::ContextLibrary<cs::AA>>(fin);
  fclose(fin);

  cs::TransformToLin(*cs_lib_);

  init_prefilter(cs219_database);
}

////////////////////////////////////////////////////////////////////////
// Main prefilter function
////////////////////////////////////////////////////////////////////////
void Prefilter::prefilter_db(
    HMM* q_tmp, Hash<Hit>* previous_hits, const int threads,
    const int prefilter_gap_open, const int prefilter_gap_extend,
    const int prefilter_score_offset, const int prefilter_bit_factor,
    const double prefilter_evalue_thresh,
    const double prefilter_evalue_coarse_thresh,
    const int preprefilter_smax_thresh, const int min_prefilter_hits,
    const int maxnumdb, const float R[20][20],
    std::vector<std::pair<int, std::string>>& new_prefilter_hits,
    std::vector<std::pair<int, std::string>>& old_prefilter_hits) {
  std::unique_ptr<Hash<char>> doubled = std::make_unique<Hash<char>>();
  doubled->New(16381, 0);

  int element_count = (VECSIZE_INT * 4);
  int W = (q_tmp->L + (element_count - 1)) / element_count;
  // query profile (states + 1 because of ANY char)
  int profile_element_count = (NUMCOLSTATES + 1) * (q_tmp->L + element_count);
  simd_int* qc = malloc_simd_int(profile_element_count * sizeof(unsigned char));
  stripe_query_profile(q_tmp, prefilter_score_offset, prefilter_bit_factor, W,
                       reinterpret_cast<unsigned char*>(qc));

  std::vector<simd_int*> workspaces(threads);

  std::vector<std::pair<int, int>> first_prefilter;
  std::vector<std::pair<double, int>> hits;

  int count_dbs = 0;
  int gap_init = prefilter_gap_open + prefilter_gap_extend;
  int gap_extend = prefilter_gap_extend;
  int LQ = q_tmp->L;
  const float log_qlen = flog2(LQ);
  const double factor = static_cast<double>(dbs_.size() * LQ);

  s_profile* ssw_profile =
      ssw_hhblits_init(reinterpret_cast<__m128i*>(qc), LQ, NUMCOLSTATES + 1, prefilter_score_offset);

  for (auto& workspace : workspaces) {
    workspace =
        (simd_int*)malloc_simd_int(3 * (LQ + element_count) * sizeof(char));
  }

#pragma omp parallel for schedule(static)
  // Loop over all database sequences
  for (size_t n = 0; n < dbs_.size(); n++) {
    const auto& db_entry = dbs_[n];
    int thread_id = 0;
#ifdef OPENMP
    thread_id = omp_get_thread_num();
#endif
    // Perform search step
    int score =
        ungapped_sse_score(qc, LQ, db_entry.first, db_entry.length,
                           prefilter_score_offset, workspaces[thread_id]);

    score = score -
            (int)(prefilter_bit_factor * (log_qlen + flog2(db_entry.length)));

#pragma omp critical
    first_prefilter.push_back(std::pair<int, int>(score, n));
  }
  // filter after calculation of ungapped sse score to include at least
  // min_prefilter_hits
  std::vector<std::pair<int, int>>::iterator it;

  sort(first_prefilter.begin(), first_prefilter.end());
  std::reverse(first_prefilter.begin(), first_prefilter.end());

  std::vector<std::pair<int, int>>::iterator first_prefilter_begin_erase =
      first_prefilter.end();
  std::vector<std::pair<int, int>>::iterator first_prefilter_end_erase =
      first_prefilter.end();
  count_dbs = 0;
  for (it = first_prefilter.begin(); it < first_prefilter.end(); it++) {
    if (count_dbs >= min_prefilter_hits &&
        (*it).first <= preprefilter_smax_thresh) {
      first_prefilter_begin_erase = it;
      break;
    } else {
      count_dbs++;
    }
  }

  first_prefilter.erase(first_prefilter_begin_erase, first_prefilter_end_erase);

  HH_LOG(INFO)
      << "HMMs passed 1st prefilter (gapless profile-profile alignment)  : "
      << count_dbs << std::endl;

#pragma omp parallel for schedule(static)
  // Loop over all database sequences
  for (size_t i = 0; i < first_prefilter.size(); i++) {
    // int thread_id = 0;
#ifdef OPENMP
    // thread_id = omp_get_thread_num();
#endif

    int n = first_prefilter[i].second;
    const auto& db_entry = dbs_[n];

    s_align* align_result =
        ssw_align(ssw_profile, reinterpret_cast<const uint8_t*>(db_entry.first),
                  db_entry.length, gap_init, gap_extend, 0, 0, 0, 15);

    if (align_result == nullptr) {
      HH_LOG(ERROR) << "SSW algorithm failed.";
      continue;
    }
    int score = align_result->score1;
    align_destroy(align_result);

    double evalue =
        factor * db_entry.length * fpow2(-score / prefilter_bit_factor);

    if (evalue < prefilter_evalue_coarse_thresh) {
#pragma omp critical
      hits.push_back(std::pair<double, int>(evalue, n));
    }
  }

  // filter after calculation of evalues to include at least
  // min_prefilter_hits
  sort(hits.begin(), hits.end());

  std::vector<std::pair<double, int>>::iterator second_prefilter_begin_erase =
      hits.end();
  std::vector<std::pair<double, int>>::iterator second_prefilter_end_erase =
      hits.end();
  std::vector<std::pair<double, int>>::iterator it2;

  count_dbs = 0;
  for (it2 = hits.begin(); it2 < hits.end(); it2++) {
    if (count_dbs >= min_prefilter_hits &&
        (*it2).first > prefilter_evalue_thresh) {
      second_prefilter_begin_erase = it2;
      break;
    } else {
      count_dbs++;
    }
  }

  hits.erase(second_prefilter_begin_erase, second_prefilter_end_erase);

  count_dbs = 0;

  for (it2 = hits.begin(); it2 < hits.end(); it2++) {
    // Add hit to dbfiles
    count_dbs++;
    const auto& db_name = dbs_[(*it2).second].name;

    char name[NAMELEN];
    RemoveExtension(name, (char*)db_name.c_str());

    if (!doubled->Contains(db_name.c_str())) {
      doubled->Add((char*)db_name.c_str());

      std::pair<int, std::string> result;
      result.first = dbs_[(*it2).second].length;
      result.second = db_name;

      // check, if DB was searched in previous rounds
      std::stringstream ss_tmp;
      ss_tmp << name << "__" << 1;

      if (previous_hits->Contains((char*)ss_tmp.str().c_str())) {
        old_prefilter_hits.push_back(result);
      } else {
        new_prefilter_hits.push_back(result);
      }
    }
    if (count_dbs >= maxnumdb) {
      HH_LOG(WARNING) << "Number of hits passing 2nd prefilter (reduced from "
                      << hits.size() << " to allowed maximum of " << maxnumdb
                      << ").\n"
                      << "You can increase the allowed maximum using the "
                         "-maxfilt <max> option.\n";
      break;
    }
  }

  // Free memory
  free(qc);
  for (auto& workspace : workspaces) {
    free(workspace);
  }
  init_destroy(ssw_profile);
}
} /* namespace hh */
