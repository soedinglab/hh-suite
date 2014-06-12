/*
 * HHDatabase.cpp
 *
 *  Created on: Apr 7, 2014
 *      Author: meiermark
 */

#include "hhdatabase.h"

#include <stddef.h>
#include <sys/mman.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "hhalignment.h"
#include "hhdecl.h"
#include "hhhmm.h"
#include "util-inl.h"

FFindexDatabase::FFindexDatabase(char* data_filename, char* index_filename,
    int superId, bool isCompressed) {
  static int runaway;
  this->id = runaway++;

  this->superId = superId;

  this->isCompressed = isCompressed;

  this->data_filename = new char[strlen(data_filename) + 1];
  strcpy(this->data_filename, data_filename);

  db_data_fh = fopen(data_filename, "r");
  if (db_data_fh == NULL) {
    OpenFileError(data_filename, __FILE__, __LINE__, __func__);
  }

  FILE* db_index_fh = fopen(index_filename, "r");
  if (db_index_fh == NULL) {
    OpenFileError(index_filename, __FILE__, __LINE__, __func__);
  }

  size_t ca3m_data_size = CountLinesInFile(index_filename);
  db_index = ffindex_index_parse(db_index_fh, ca3m_data_size);

  fclose(db_index_fh);

  if (db_index == NULL) {
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << std::endl;
    std::cerr << "\tcould not read index file" << index_filename
        << ". Is the file empty or corrupted?\n";
    exit(1);
  }

  db_data = ffindex_mmap_data(db_data_fh, &data_size);
}

FFindexDatabase::~FFindexDatabase() {
  delete[] data_filename;
  munmap(db_data, data_size);
  free(db_index);
  fclose(db_data_fh);
}

HHDatabase::HHDatabase() {

}

HHDatabase::~HHDatabase() {
}

void HHDatabase::buildDatabaseName(const char* base, const char* extension,
    const char* suffix, char* databaseName) {
  strcpy(databaseName, base);
  if (strlen(extension) != 0)
    strcat(databaseName, "_");
  strcat(databaseName, extension);
  strcat(databaseName, suffix);
}

HHsearchDatabase::HHsearchDatabase(char* base) {
  char index_filename[NAMELEN];
  char data_filename[NAMELEN];

  basename = new char[strlen(base) + 1];
  strcpy(basename, base);

  buildDatabaseName(base, "", ".ffdata", data_filename);
  buildDatabaseName(base, "", ".ffindex", index_filename);

  //TODO:
  database = new FFindexDatabase(data_filename, index_filename, 0, false);
}

HHsearchDatabase::~HHsearchDatabase() {
  delete[] basename;
  delete database;
}

HHblitsDatabase::HHblitsDatabase(const char* base) {
  static int runaway;
  id = runaway++;

  basename = new char[strlen(base) + 1];
  strcpy(basename, base);

  char cs219_index_filename[NAMELEN];
  char cs219_data_filename[NAMELEN];

  buildDatabaseName(base, "cs219", ".ffdata", cs219_data_filename);
  buildDatabaseName(base, "cs219", ".ffindex", cs219_index_filename);

  cs219_database = new FFindexDatabase(cs219_data_filename,
      cs219_index_filename, id, use_compressed);

  if (!checkAndBuildCompressedDatabase(base)) {
    char a3m_index_filename[NAMELEN];
    char a3m_data_filename[NAMELEN];

    char hhm_index_filename[NAMELEN];
    char hhm_data_filename[NAMELEN];

    buildDatabaseName(base, "a3m", ".ffdata", a3m_data_filename);
    buildDatabaseName(base, "a3m", ".ffindex", a3m_index_filename);

    buildDatabaseName(base, "hhm", ".ffdata", hhm_data_filename);
    buildDatabaseName(base, "hhm", ".ffindex", hhm_index_filename);

    a3m_database = new FFindexDatabase(a3m_data_filename, a3m_index_filename,
        id, use_compressed);
    hhm_database = new FFindexDatabase(hhm_data_filename, hhm_index_filename,
        id, use_compressed);
  }

  prefilter = NULL;
}

HHblitsDatabase::~HHblitsDatabase() {
  delete[] basename;
  delete cs219_database;

  if (use_compressed) {
    delete ca3m_database;
    delete sequence_database;
    delete header_database;
  }
  else {
    delete a3m_database;
    delete hhm_database;
  }

  if (prefilter) {
    delete prefilter;
  }
}

void HHblitsDatabase::initPrefilter(const char* cs_library) {
  prefilter = new hh::Prefilter(cs_library, cs219_database);
}

void HHblitsDatabase::initNoPrefilter(
    std::vector<HHDatabaseEntry*>& new_entries) {
  std::vector<std::pair<int, std::string>> new_entry_names;
  hh::Prefilter::init_no_prefiltering(cs219_database, new_entry_names);

  getEntriesFromNames(new_entry_names, new_entries);
}

void HHblitsDatabase::prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits,
    const int threads, const int prefilter_gap_open,
    const int prefilter_gap_extend, const int prefilter_score_offset,
    const int prefilter_bit_factor, const double prefilter_evalue_thresh,
    const double prefilter_evalue_coarse_thresh,
    const int preprefilter_smax_thresh, const int min_prefilter_hits,
    const float R[20][20], std::vector<HHDatabaseEntry*>& new_entries,
    std::vector<HHDatabaseEntry*>& old_entries) {

  std::vector<std::pair<int, std::string>> prefiltered_new_entry_names;
  std::vector<std::pair<int, std::string>> prefiltered_old_entry_names;

  prefilter->prefilter_db(q_tmp, previous_hits, threads, prefilter_gap_open,
      prefilter_gap_extend, prefilter_score_offset, prefilter_bit_factor,
      prefilter_evalue_thresh, prefilter_evalue_coarse_thresh,
      preprefilter_smax_thresh, min_prefilter_hits, R,
      prefiltered_new_entry_names, prefiltered_old_entry_names);

  getEntriesFromNames(prefiltered_new_entry_names, new_entries);
  getEntriesFromNames(prefiltered_old_entry_names, old_entries);
}

void HHblitsDatabase::getEntriesFromNames(std::vector<std::pair<int, std::string>>& hits,
    std::vector<HHDatabaseEntry*>& entries) {
  for (size_t i = 0; i < hits.size(); i++) {
    ffindex_entry_t* entry;

    if (use_compressed) {
      entry = ffindex_get_entry_by_name(ca3m_database->db_index, const_cast<char*>(hits[i].second.c_str()));
      if (entry == NULL) {
        //TODO: error
        HH_LOG(LogLevel::WARNING) << "warning: could not fetch entry from compressed a3m!"
            << std::endl;
        HH_LOG(LogLevel::WARNING) << "\tentry: " << hits[i].second << std::endl;
        HH_LOG(LogLevel::WARNING) << "\tdb: " << ca3m_database->data_filename << std::endl;
        continue;
      }

      HHDatabaseEntry* hhentry = new HHDatabaseEntry();
      hhentry->entry = entry;
      hhentry->ffdatabase = ca3m_database;
      hhentry->sequence_length = hits[i].first;

      entries.push_back(hhentry);
    }
    else {
      entry = ffindex_get_entry_by_name(hhm_database->db_index,
          const_cast<char*>(hits[i].second.c_str()));

      if (entry != NULL) {
        HHDatabaseEntry* hhentry = new HHDatabaseEntry();
        hhentry->entry = entry;
        hhentry->ffdatabase = hhm_database;
        hhentry->sequence_length = hits[i].first;
        entries.push_back(hhentry);
        continue;
      }

      entry = ffindex_get_entry_by_name(a3m_database->db_index,
          const_cast<char*>(hits[i].second.c_str()));

      if (entry != NULL) {
        HHDatabaseEntry* hhentry = new HHDatabaseEntry();
        hhentry->entry = entry;
        hhentry->ffdatabase = a3m_database;
        hhentry->sequence_length = hits[i].first;
        entries.push_back(hhentry);
      }
      else {
        //TODO: error
        HH_LOG(LogLevel::WARNING) << "warning: could not fetch entry from a3m or hhm!" << std::endl;
        HH_LOG(LogLevel::WARNING) << "\tentry: " << hits[i].second << std::endl;
        HH_LOG(LogLevel::WARNING) << "\ta3m_db: " << a3m_database->data_filename << std::endl;
        HH_LOG(LogLevel::WARNING) << "\thhm_db: " << hhm_database->data_filename << std::endl;
        continue;
      }
    }
  }
}

bool HHblitsDatabase::checkAndBuildCompressedDatabase(const char* base) {
  char ca3m_index_filename[NAMELEN];
  char ca3m_data_filename[NAMELEN];

  char sequence_index_filename[NAMELEN];
  char sequence_data_filename[NAMELEN];

  char header_index_filename[NAMELEN];
  char header_data_filename[NAMELEN];

  buildDatabaseName(base, "ca3m", ".ffdata", ca3m_data_filename);
  buildDatabaseName(base, "ca3m", ".ffindex", ca3m_index_filename);

  buildDatabaseName(base, "sequence", ".ffdata", sequence_data_filename);
  buildDatabaseName(base, "sequence", ".ffindex", sequence_index_filename);

  buildDatabaseName(base, "header", ".ffdata", header_data_filename);
  buildDatabaseName(base, "header", ".ffindex", header_index_filename);

  if (file_exists(ca3m_index_filename) && file_exists(ca3m_data_filename)
      && file_exists(header_index_filename) && file_exists(header_data_filename)
      && file_exists(sequence_index_filename)
      && file_exists(sequence_data_filename)) {
    use_compressed = true;

    ca3m_database = new FFindexDatabase(ca3m_data_filename, ca3m_index_filename,
        id, use_compressed);
    sequence_database = new FFindexDatabase(sequence_data_filename,
        sequence_index_filename, id, use_compressed);
    header_database = new FFindexDatabase(header_data_filename,
        header_index_filename, id, use_compressed);
  }
  else {
    use_compressed = false;
  }

  return use_compressed;
}


void getTemplateHMM(Parameters& par, HHDatabaseEntry& entry, std::vector<HHblitsDatabase*>& dbs,
    char use_global_weights, int& format, float* pb, const float S[20][20], const float Sim[20][20], HMM* t) {
    if (entry.ffdatabase->isCompressed) {
        Alignment tali;

        char* data = ffindex_get_data_by_entry(entry.ffdatabase->db_data,
                entry.entry);

        if (data == NULL) {
            std::cerr << "Could not fetch data for a3m " << entry.entry->name
                    << "!" << std::endl;
            exit(4);
        }

        HHblitsDatabase* db = getHHblitsDatabase(entry, dbs);
        if (db == NULL) {
            //TODO throw error
            std::cerr << "this should not happen!!!!" << std::endl;
            exit(0);
        }

        tali.ReadCompressed(entry.entry, data, db->sequence_database->db_index,
                db->sequence_database->db_data, db->header_database->db_index,
                db->header_database->db_data, par.mark, par.maxcol);

        tali.Compress(entry.entry->name, par.cons, par.maxres, par.maxcol,
                par.M, par.Mgaps);

        tali.N_filtered = tali.Filter(par.max_seqid_db, S, par.coverage_db,
                par.qid_db, par.qsc_db, par.Ndiff_db);
        t->name[0] = t->longname[0] = t->fam[0] = '\0';
        tali.FrequenciesAndTransitions(t, use_global_weights, par.mark,
                par.cons, par.showcons, par.maxres, pb, Sim);

        format = 0;
    } else {
        FILE* dbf = ffindex_fopen_by_entry(entry.ffdatabase->db_data,
                entry.entry);

        if (dbf != NULL) {
            char line[LINELEN];
            if (!fgetline(line, LINELEN, dbf)) {
                std::cerr << "this should not happen!" << std::endl;
                //TODO: throw error
            }

            while (strscn(line) == NULL)
                fgetline(line, LINELEN, dbf); // skip lines that contain only white space

            if (!strncmp(line, "HMMER3", 6))      // read HMMER3 format
                    {
                format = 1;
                t->ReadHMMer3(dbf, par.showcons, pb, entry.entry->name);
                par.hmmer_used = true;
            } else if (!strncmp(line, "HMMER", 5))      // read HMMER format
                    {
                format = 1;
                t->ReadHMMer(dbf, par.showcons, pb, entry.entry->name);
                par.hmmer_used = true;
            } else if (!strncmp(line, "HH", 2))    // read HHM format
                    {
                char path[NAMELEN];
                Pathname(path, entry.entry->name);

                format = 0;
                t->Read(dbf, par.maxcol, par.nseqdis, pb, path);

            }
            //TODO: old hhm format discarded
            // read a3m alignment
            else if (line[0] == '#' || line[0] == '>') {
                Alignment tali;
                tali.Read(dbf, entry.entry->name, par.mark, par.maxcol,
                        par.nseqdis, line);
                tali.Compress(entry.entry->name, par.cons, par.maxres,
                        par.maxcol, par.M, par.Mgaps);
                //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
                tali.N_filtered = tali.Filter(par.max_seqid_db, S,
                        par.coverage_db, par.qid_db, par.qsc_db, par.Ndiff_db);
                t->name[0] = t->longname[0] = t->fam[0] = '\0';
                tali.FrequenciesAndTransitions(t, use_global_weights, par.mark,
                        par.cons, par.showcons, par.maxres, pb, Sim);
                format = 0;
            } else {
                std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
                        << __func__ << ":" << std::endl;
                std::cerr << "\tunrecognized HMM file format in \'"
                        << entry.entry->name << "\'. \n";
                std::cerr << "Context:\n'" << line << "\n";
                fgetline(line, LINELEN, dbf);
                std::cerr << line << "\n";
                fgetline(line, LINELEN, dbf);
                std::cerr << line << "'\n";
                exit(1);
            }

            fclose(dbf);
            return;
        }
    }
}

HHblitsDatabase* getHHblitsDatabase(HHDatabaseEntry& entry, std::vector<HHblitsDatabase*>& dbs) {
    for (std::vector<HHblitsDatabase*>::size_type i = 0; i < dbs.size(); i++) {
        if (dbs[i]->id == entry.ffdatabase->superId) {
            return dbs[i];
        }
    }

    return NULL;
}


