/*
 * HHDatabase.h
 *
 *  Created on: Apr 7, 2014
 *      Author: meiermark
 */

#ifndef HHDATABASE_H_
#define HHDATABASE_H_

class HHEntry;
class HHDatabaseEntry;
class Alignment;
class Prefilter;

#include <cstdlib>

#include "ffindexdatabase.h"
#include "hhutil.h"
#include "hash.h"
#include "hhhit.h"
#include "log.h"
#include "hhalignment.h"

class HHDatabase {
  public:
    HHDatabase();
    virtual ~HHDatabase();

    static bool checkDatabaseConflicts(const char* base);

protected:
    static void buildDatabaseName(const char* base, const char* extension,
                                  const char* suffix, char* databaseName);
};

class HHblitsDatabase: HHDatabase {
  public:
    HHblitsDatabase(const char* base, bool initCs219 = true);
    ~HHblitsDatabase();

    void initPrefilter(const std::string& cs_library);
    void initNoPrefilter(std::vector<HHEntry*>& new_prefilter_hits);
    void initSelected(std::vector<std::string>& selected_templates,
        std::vector<HHEntry*>& new_entries);

    void prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits, const int threads,
        const int prefilter_gap_open, const int prefilter_gap_extend,
        const int prefilter_score_offset, const int prefilter_bit_factor,
        const double prefilter_evalue_thresh,
        const double prefilter_evalue_coarse_thresh,
        const int preprefilter_smax_thresh, const int min_prefilter_hits, const int maxnumdb,
        const float R[20][20], std::vector<HHEntry*>& new_entries,
        std::vector<HHEntry*>& old_entries);

    char* basename;

    FFindexDatabase* cs219_database;

    FFindexDatabase* a3m_database;
    FFindexDatabase* hhm_database;

    FFindexDatabase* query_database;

    bool use_compressed;
    FFindexDatabase* ca3m_database;
    FFindexDatabase* sequence_database;
    FFindexDatabase* header_database;

  private:
    void getEntriesFromNames(std::vector<std::pair<int, std::string> >& names,
        std::vector<HHEntry*>& entries);
    bool checkAndBuildCompressedDatabase(const char* base);

    Prefilter* prefilter;
};

class HHEntry {
  public:
    int sequence_length;

    HHEntry(int sequence_length);
    virtual ~HHEntry();

    virtual void getTemplateA3M(Parameters& par, float* pb, const float S[20][20],
        const float Sim[20][20], Alignment& tali) {};
    virtual void getTemplateHMM(Parameters& par, char use_global_weights, const float qsc, int& format,
        float* pb, const float S[20][20], const float Sim[20][20], HMM* t) {};

    virtual char* getName() {return NULL;};

  protected:
    void getTemplateHMM(FILE* inf, char* name, Parameters& par, char use_global_weights,
        const float qsc, int& format, float* pb, const float S[20][20],
        const float Sim[20][20], HMM* t);
};

class HHDatabaseEntry : public HHEntry {
  public:
    HHDatabaseEntry(int sequence_length, HHblitsDatabase* hhdatabase, FFindexDatabase* ffdatabase, ffindex_entry_t* entry);
    ~HHDatabaseEntry();

    void getTemplateA3M(Parameters& par, float* pb, const float S[20][20],
        const float Sim[20][20], Alignment& tali);
    void getTemplateHMM(Parameters& par, char use_global_weights, const float qsc, int& format,
        float* pb, const float S[20][20], const float Sim[20][20], HMM* t);

    char* getName();

  private:
    HHblitsDatabase* hhdatabase;
    FFindexDatabase* ffdatabase;
    ffindex_entry_t* entry;
};

class HHFileEntry : public HHEntry {
  public:
    HHFileEntry(const char* file, int sequence_length);
    ~HHFileEntry();

    void getTemplateA3M(Parameters& par, float* pb, const float S[20][20],
        const float Sim[20][20], Alignment& tali);
    void getTemplateHMM(Parameters& par, char use_global_weights, const float qsc, int& format,
        float* pb, const float S[20][20], const float Sim[20][20], HMM* t);

    char* getName();

  private:
    char* file;
};

struct HHDatabaseEntryCompare {
  bool operator()(const HHEntry* l, const HHEntry* r) {
    return (*l).sequence_length > (*r).sequence_length;
  }
};

int getMaxTemplateLength(std::vector<HHEntry*>& entries);

#endif /* HHDATABASE_H_ */
