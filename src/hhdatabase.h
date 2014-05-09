/*
 * HHDatabase.h
 *
 *  Created on: Apr 7, 2014
 *      Author: meiermark
 */


#ifndef HHDATABASE_H_
#define HHDATABASE_H_

class HHDatabaseEntry;
class FFindexDatabase;


extern "C" {
#include <ffindex.h>
}

#include "hhutil.h"
#include "hhprefilter.h"
#include "hash.h"
#include "hhhit.h"



class FFindexDatabase {
  public:
    FFindexDatabase(char* data_filename, char* index_filename, int superId, bool isCompressed);
    virtual ~FFindexDatabase();

    ffindex_index_t* db_index = NULL;
    char* db_data;
    char* data_filename;

    bool isCompressed;

    int id;
    int superId;

  private:
    size_t data_size;
    FILE* db_data_fh;
};

class HHDatabase {
  public:
    HHDatabase();
    virtual ~HHDatabase();

  protected:
    void buildDatabaseName(const char* base, const char* extension, const char* suffix, char* databaseName);
};

class HHsearchDatabase : HHDatabase {
  public:
    HHsearchDatabase(char* base);
    ~HHsearchDatabase();
    FFindexDatabase* database;
    char* basename;

  private:
};

class HHblitsDatabase : HHDatabase {
  public:
    HHblitsDatabase(const char* base);
    ~HHblitsDatabase();

    void initPrefilter(const char* cs_library);
    void initNoPrefilter(std::vector<HHDatabaseEntry*>& new_prefilter_hits);
    void prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits,
			const int threads, const int prefilter_gap_open, const int prefilter_gap_extend,
			const int prefilter_score_offset, const int prefilter_bit_factor,
			const double prefilter_evalue_thresh, const double prefilter_evalue_coarse_thresh,
			const int preprefilter_smax_thresh, const int min_prefilter_hits, const float R[20][20],
			std::vector<HHDatabaseEntry*>& new_entries, std::vector<HHDatabaseEntry*>& old_entries);

    char* basename;

    int id;

    FFindexDatabase* cs219_database = NULL;

    FFindexDatabase* a3m_database = NULL;
    FFindexDatabase* hhm_database = NULL;

    bool use_compressed = false;
    FFindexDatabase* ca3m_database = NULL;
    FFindexDatabase* sequence_database = NULL;
    FFindexDatabase* header_database = NULL;

  private:
    void getEntriesFromNames(std::vector<std::string>& names, std::vector<HHDatabaseEntry*>& entries);
    bool checkAndBuildCompressedDatabase(const char* base);

    hh::Prefilter* prefilter;
};

struct HHDatabaseEntry {
	FFindexDatabase* ffdatabase;
	ffindex_entry_t* entry;
};

#endif /* HHDATABASE_H_ */
