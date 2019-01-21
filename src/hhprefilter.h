/*
 * hhprefilter.h
 *
 *  Created on: Apr 17, 2014
 *      Author: meiermark
 */

#ifndef HHPREFILTER_H_
#define HHPREFILTER_H_

namespace hh {
class Prefilter;
}


#ifdef OPENMP
#include <omp.h>
#endif

#include <memory>
#include <string>
#include <vector>

#include "hash.h"
#include "hhhit.h"
#include "hhhmm.h"

namespace hh {

struct DBEntry {
  // Sequence names in prefilter db file.
  std::string name;

  // Pointer to first letter of next sequence in db_data.
  char* first;

  // Length of next sequence.
  int length;
};

// HHSuite Prefilter, applying the Smith-Waterman algorithm to a
// protein sequence database. 

class Prefilter {
 public:
  // Initialize the Prefilter. This initializer will open a file handle and read
  // the data from disk, and perform some pre-initialization work. Initializing
  // an instance of this class should be considered an I/O-blocking non-trivial
  // operation.
  //
  // cs_library: The file location of the ContextLibrary that should be used.
  // cs219_database: The FFindexDatabase to filter.
  Prefilter(const char* cs_library, FFindexDatabase* cs219_database);

  // Reads the FFindexDatabase entries, populating their length and names into
  // `prefiltered_entries`. `prefiltered_entries` will be populated with pairs
  // of `<entry length>`, `<entry name>`.
  static void init_no_prefiltering(
      FFindexDatabase* cs219_database,
      std::vector<std::pair<int, std::string>>& prefiltered_entries);

  // Reads the FFindexDatabase entries from `cs219_database` corresponding to
  // the entry names in `templates`, populating their length and names into
  // `prefiltered_entries`. `prefiltered_entries` will be populated with pairs
  // of `<entry length>`, `<entry name>`.
  static void init_selected(
      FFindexDatabase* cs219_database, std::vector<std::string> templates,
      std::vector<std::pair<int, std::string>>& prefiltered_entries);

  // Pre-filter the database.
  void prefilter_db(
      HMM* q_tmp, Hash<Hit>* previous_hits, const int threads,
      const int prefilter_gap_open, const int prefilter_gap_extend,
      const int prefilter_score_offset, const int prefilter_bit_factor,
      const double prefilter_evalue_thresh,
      const double prefilter_evalue_coarse_thresh,
      const int preprefilter_smax_thresh, const int min_prefilter_hits,
      const int maxnumdb, const float R[20][20],
      std::vector<std::pair<int, std::string>>& new_prefilter_hits,
      std::vector<std::pair<int, std::string>>& old_prefilter_hits);

 private:
  // array containing all sequence names in prefilter db file.
  std::vector<DBEntry> dbs_;

  void stripe_query_profile(HMM* q_tmp, const int prefilter_score_offset,
                            const int prefilter_bit_factor, const int W,
                            unsigned char* qc);
  std::unique_ptr<cs::ContextLibrary<cs::AA>> cs_lib_;

  void init_prefilter(FFindexDatabase* cs219_database);
  void checkCSFormat(size_t nr_checks);

  const int NUMCOLSTATES = cs::AS219::kSize;
};

} /* namespace hh */

#endif /* HHPREFILTER_H_ */
