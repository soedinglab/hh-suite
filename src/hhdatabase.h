/*
 * HHDatabase.h
 *
 *  Created on: Apr 7, 2014
 *      Author: meiermark
 */


#ifndef HHDATABASE_H_
#define HHDATABASE_H_

class HHDatabaseEntry;

extern "C" {
#include <ffindex.h>
}

#include "hhutil.h"



class FFindexDatabase {
  public:
    FFindexDatabase(char* data_filename, char* index_filename);
    virtual ~FFindexDatabase();

    ffindex_index_t* db_index = NULL;
    char* db_data;
    char* data_filename;

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
    HHblitsDatabase(char* base);
    ~HHblitsDatabase();

    char* basename;

    FFindexDatabase* cs219_database = NULL;

    FFindexDatabase* a3m_database = NULL;
    FFindexDatabase* hhm_database = NULL;

    bool use_compressed = false;
    FFindexDatabase* ca3m_database = NULL;
    FFindexDatabase* sequence_database = NULL;
    FFindexDatabase* header_database = NULL;

  private:
    bool checkAndBuildCompressedDatabase(char* base);
    char* base;
};

class HHDatabaseEntry {
  public:
    HHDatabaseEntry(ffindex_entry_t* entry, HHsearchDatabase* db, FFindexDatabase* ffdb);
    virtual ~HHDatabaseEntry();

    HHsearchDatabase* database;
    FFindexDatabase* ffdatabase;
    ffindex_entry_t* entry;
  private:
};

#endif /* HHDATABASE_H_ */
