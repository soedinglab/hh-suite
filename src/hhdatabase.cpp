/*
 * HHDatabase.cpp
 *
 *  Created on: Apr 7, 2014
 *      Author: meiermark
 */

#include "hhdatabase.h"

#include "hhdecl.h"

FFindexDatabase::FFindexDatabase(char* data_filename, char* index_filename) {
  this->data_filename = new char [strlen(data_filename) + 1];
  strcpy(this->data_filename, data_filename);

  db_data_fh = fopen(data_filename, "r");
  FILE* db_index_fh = fopen(index_filename, "r");

  if (db_data_fh == NULL) {
    OpenFileError(data_filename, __FILE__, __LINE__, __func__);
  }

  if (db_index_fh == NULL) {
    OpenFileError(data_filename, __FILE__, __LINE__, __func__);
  }

  size_t ca3m_data_size = CountLinesInFile(index_filename);

  db_index = ffindex_index_parse(db_index_fh, ca3m_data_size);

  if (db_index == NULL) {
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    std::cerr << "\tcould not read index file" << index_filename << ". Is the file empty or corrupted?\n";
    exit(1);
  }

  db_data = ffindex_mmap_data(db_data_fh, &data_size);
}

FFindexDatabase::~FFindexDatabase() {
  delete[] db_index;
  fclose(db_data_fh);
}

HHDatabase::HHDatabase() {

}

HHDatabase::~HHDatabase() {
}

void HHDatabase::buildDatabaseName(const char* base, const char* extension,
    const char* suffix, char* databaseName) {
  strcpy(databaseName, base);
  if(strlen(extension) != 0)
    strcat(databaseName, "_");
  strcat(databaseName, extension);
  strcat(databaseName, suffix);
}

HHsearchDatabase::HHsearchDatabase(char* base) {
  char index_filename[NAMELEN];
  char data_filename[NAMELEN];

  basename = new char[strlen(base)+1];
  strcpy(basename, base);

  buildDatabaseName(base, "", ".ffdata", data_filename);
  buildDatabaseName(base, "", ".ffindex", index_filename);

  database = new FFindexDatabase(data_filename, index_filename);
}

HHsearchDatabase::~HHsearchDatabase() {
  delete[] basename;
  delete database;
}

HHblitsDatabase::HHblitsDatabase(char* base) {
  char cs219_index_filename[NAMELEN];
  char cs219_data_filename[NAMELEN];

  buildDatabaseName(base, "cs219", ".ffdata", cs219_data_filename);
  buildDatabaseName(base, "cs219", ".ffindex", cs219_index_filename);

  cs219_database = new FFindexDatabase(cs219_data_filename,
      cs219_index_filename);

  if (!checkAndBuildCompressedDatabase(base)) {
    char a3m_index_filename[NAMELEN];
    char a3m_data_filename[NAMELEN];

    char hhm_index_filename[NAMELEN];
    char hhm_data_filename[NAMELEN];

    buildDatabaseName(base, "a3m", ".ffdata", a3m_data_filename);
    buildDatabaseName(base, "a3m", ".ffindex", a3m_data_filename);

    buildDatabaseName(base, "hhm", ".ffdata", hhm_data_filename);
    buildDatabaseName(base, "hhm", ".ffindex", hhm_data_filename);

    a3m_database = new FFindexDatabase(a3m_data_filename, a3m_index_filename);
    hhm_database = new FFindexDatabase(hhm_data_filename, hhm_index_filename);
  }
}

HHblitsDatabase::~HHblitsDatabase() {
}


bool HHblitsDatabase::checkAndBuildCompressedDatabase(char* base) {
  char ca3m_index_filename[NAMELEN];
  char ca3m_data_filename[NAMELEN];

  char sequence_index_filename[NAMELEN];
  char sequence_data_filename[NAMELEN];

  char header_index_filename[NAMELEN];
  char header_data_filename[NAMELEN];

  buildDatabaseName(base, "ca3m", ".ffdata", ca3m_data_filename);
  buildDatabaseName(base, "ca3m", ".ffindex", ca3m_data_filename);

  buildDatabaseName(base, "sequence", ".ffdata", sequence_data_filename);
  buildDatabaseName(base, "sequence", ".ffindex", sequence_index_filename);

  buildDatabaseName(base, "header", ".ffdata", header_data_filename);
  buildDatabaseName(base, "header", ".ffindex", header_index_filename);

  if (file_exists(ca3m_index_filename) && file_exists(ca3m_data_filename)
      && file_exists(header_index_filename) && file_exists(header_data_filename)
      && file_exists(sequence_index_filename)
      && file_exists(sequence_data_filename)) {
    use_compressed = true;

    ca3m_database = new FFindexDatabase(ca3m_data_filename,
        ca3m_index_filename);
    sequence_database = new FFindexDatabase(sequence_data_filename,
        sequence_index_filename);
    header_database = new FFindexDatabase(header_data_filename,
        header_index_filename);
  }
  else {
    use_compressed = false;
  }

  return use_compressed;
}

HHDatabaseEntry::HHDatabaseEntry(ffindex_entry* entry, HHsearchDatabase* db, FFindexDatabase* ffdb) {
  this->entry = entry;
  this->database = db;
  this->ffdatabase = ffdb;
}


HHDatabaseEntry::~HHDatabaseEntry() {
}


