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
                                 bool isCompressed) {

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

HHblitsDatabase::HHblitsDatabase(const char* base) {
  use_compressed = false;
  basename = new char[strlen(base) + 1];
  strcpy(basename, base);

  char cs219_index_filename[NAMELEN];
  char cs219_data_filename[NAMELEN];

  buildDatabaseName(base, "cs219", ".ffdata", cs219_data_filename);
  buildDatabaseName(base, "cs219", ".ffindex", cs219_index_filename);

  cs219_database = new FFindexDatabase(cs219_data_filename,
                                       cs219_index_filename, use_compressed);

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
                                       use_compressed);
    hhm_database = new FFindexDatabase(hhm_data_filename, hhm_index_filename,
                                       use_compressed);
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
  } else {
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

void HHblitsDatabase::initNoPrefilter(std::vector<HHEntry*>& new_entries) {
  std::vector<std::pair<int, std::string> > new_entry_names;
  hh::Prefilter::init_no_prefiltering(cs219_database, new_entry_names);

  getEntriesFromNames(new_entry_names, new_entries);
}

void HHblitsDatabase::initSelected(std::vector<std::string>& selected_templates,
                                   std::vector<HHEntry*>& new_entries) {

  std::vector<std::pair<int, std::string> > new_entry_names;
  hh::Prefilter::init_selected(cs219_database, selected_templates,
                               new_entry_names);

  getEntriesFromNames(new_entry_names, new_entries);
}

void HHblitsDatabase::prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits,
                                   const int threads,
                                   const int prefilter_gap_open,
                                   const int prefilter_gap_extend,
                                   const int prefilter_score_offset,
                                   const int prefilter_bit_factor,
                                   const double prefilter_evalue_thresh,
                                   const double prefilter_evalue_coarse_thresh,
                                   const int preprefilter_smax_thresh,
                                   const int min_prefilter_hits,
                                   const int maxnumbdb, const float R[20][20],
                                   std::vector<HHEntry*>& new_entries,
                                   std::vector<HHEntry*>& old_entries) {

  std::vector<std::pair<int, std::string> > prefiltered_new_entry_names;
  std::vector<std::pair<int, std::string> > prefiltered_old_entry_names;

  prefilter->prefilter_db(q_tmp, previous_hits, threads, prefilter_gap_open,
                          prefilter_gap_extend, prefilter_score_offset,
                          prefilter_bit_factor, prefilter_evalue_thresh,
                          prefilter_evalue_coarse_thresh,
                          preprefilter_smax_thresh, min_prefilter_hits,
                          maxnumbdb, R, prefiltered_new_entry_names,
                          prefiltered_old_entry_names);

  getEntriesFromNames(prefiltered_new_entry_names, new_entries);
  getEntriesFromNames(prefiltered_old_entry_names, old_entries);
}

void HHblitsDatabase::getEntriesFromNames(
    std::vector<std::pair<int, std::string> >& hits,
    std::vector<HHEntry*>& entries) {
  for (size_t i = 0; i < hits.size(); i++) {
    ffindex_entry_t* entry;

    if (use_compressed) {
      entry = ffindex_get_entry_by_name(
          ca3m_database->db_index, const_cast<char*>(hits[i].second.c_str()));
      if (entry == NULL) {
        //TODO: error
        HH_LOG(WARNING)
            << "warning: could not fetch entry from compressed a3m!"
            << std::endl;
        HH_LOG(WARNING) << "\tentry: " << hits[i].second << std::endl;
        HH_LOG(WARNING) << "\tdb: " << ca3m_database->data_filename
                                  << std::endl;
        continue;
      }

      HHEntry* hhentry = new HHDatabaseEntry(hits[i].first, this, ca3m_database,
                                             entry);
      entries.push_back(hhentry);
    } else {
      entry = ffindex_get_entry_by_name(
          hhm_database->db_index, const_cast<char*>(hits[i].second.c_str()));

      if (entry != NULL) {
        HHEntry* hhentry = new HHDatabaseEntry(hits[i].first, this,
                                               hhm_database, entry);
        entries.push_back(hhentry);
        continue;
      }

      entry = ffindex_get_entry_by_name(
          a3m_database->db_index, const_cast<char*>(hits[i].second.c_str()));

      if (entry != NULL) {
        HHEntry* hhentry = new HHDatabaseEntry(hits[i].first, this,
                                               a3m_database, entry);
        entries.push_back(hhentry);
      } else {
        //TODO: error
        HH_LOG(WARNING)
            << "warning: could not fetch entry from a3m or hhm!" << std::endl;
        HH_LOG(WARNING) << "\tentry: " << hits[i].second << std::endl;
        HH_LOG(WARNING) << "\ta3m_db: " << a3m_database->data_filename
                                  << std::endl;
        HH_LOG(WARNING) << "\thhm_db: " << hhm_database->data_filename
                                  << std::endl;
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
                                        use_compressed);
    sequence_database = new FFindexDatabase(sequence_data_filename,
                                            sequence_index_filename,
                                            use_compressed);
    header_database = new FFindexDatabase(header_data_filename,
                                          header_index_filename,
                                          use_compressed);
  } else {
    use_compressed = false;
  }

  return use_compressed;
}

HHEntry::HHEntry(int sequence_length) {
  this->sequence_length = sequence_length;
}

HHEntry::~HHEntry() {
}

HHDatabaseEntry::HHDatabaseEntry(int sequence_length,
                                 HHblitsDatabase* hhdatabase,
                                 FFindexDatabase* ffdatabase,
                                 ffindex_entry_t* entry)
    : HHEntry(sequence_length) {
  this->hhdatabase = hhdatabase;
  this->ffdatabase = ffdatabase;
  this->entry = entry;
}

HHDatabaseEntry::~HHDatabaseEntry() {
}

void HHDatabaseEntry::getTemplateHMM(Parameters& par, char use_global_weights,
                                     const float qsc, int& format, float* pb,
                                     const float S[20][20],
                                     const float Sim[20][20], HMM* t) {
  if (ffdatabase->isCompressed) {
    Alignment tali;

    char* data = ffindex_get_data_by_entry(ffdatabase->db_data, entry);

    if (data == NULL) {
      std::cerr << "Could not fetch data for a3m " << entry->name << "!"
                << std::endl;
      exit(4);
    }

    tali.ReadCompressed(entry, data, hhdatabase->sequence_database->db_index,
                        hhdatabase->sequence_database->db_data,
                        hhdatabase->header_database->db_index,
                        hhdatabase->header_database->db_data, par.mark,
                        par.maxcol);

    tali.Compress(entry->name, par.cons, par.maxres, par.maxcol, par.M,
                  par.Mgaps);

    tali.N_filtered = tali.Filter(par.max_seqid_db, S, par.coverage_db,
                                  par.qid_db, qsc, par.Ndiff_db);
    t->name[0] = t->longname[0] = t->fam[0] = '\0';
    tali.FrequenciesAndTransitions(t, use_global_weights, par.mark, par.cons,
                                   par.showcons, par.maxres, pb, Sim);

    format = 0;
  } else {
    FILE* dbf = ffindex_fopen_by_entry(ffdatabase->db_data, entry);
    char* name = new char[strlen(entry->name) + 1];
    strcpy(name, entry->name);
    HHEntry::getTemplateHMM(dbf, name, par, use_global_weights, qsc, format, pb,
                            S, Sim, t);
    fclose(dbf);
    delete[] name;
  }
}

void HHDatabaseEntry::getTemplateA3M(Parameters& par, float* pb,
                                     const float S[20][20],
                                     const float Sim[20][20], Alignment& tali) {
  if (hhdatabase->use_compressed) {
    ffindex_entry_t* entry = ffindex_get_entry_by_name(
        hhdatabase->ca3m_database->db_index, this->entry->name);

    if (entry == NULL) {
      HH_LOG(ERROR) << "Could not fetch entry " << this->entry->name
                              << " from compressed hhblits database!"
                              << std::endl;
      exit(1);
    }

    char* data = ffindex_get_data_by_entry(hhdatabase->ca3m_database->db_data,
                                           entry);

    if (data == NULL) {
      std::cerr << "Could not fetch data for a3m " << entry->name << "!"
                << std::endl;
      exit(4);
    }

    tali.ReadCompressed(entry, data, hhdatabase->sequence_database->db_index,
                        hhdatabase->sequence_database->db_data,
                        hhdatabase->header_database->db_index,
                        hhdatabase->header_database->db_data, par.mark,
                        par.maxcol);
  } else {
    FILE* dbf = ffindex_fopen_by_name(hhdatabase->a3m_database->db_data,
                                      hhdatabase->a3m_database->db_index,
                                      entry->name);

    if (dbf == NULL) {
      std::cerr << std::endl << "Error: opening A3M " << entry->name
                << std::endl;
      exit(4);
    }

    char line[LINELEN];
    if (!fgetline(line, LINELEN, dbf)) {
      std::cerr << "this should not happen!" << std::endl;
      //TODO: throw error
    }

    while (strscn(line) == NULL)
      fgetline(line, LINELEN, dbf);  // skip lines that contain only white space

    tali.Read(dbf, entry->name, par.mark, par.maxcol, par.nseqdis, line);
    fclose(dbf);
  }

  tali.Compress(entry->name, par.cons, par.maxres, par.maxcol, par.M,
                par.Mgaps);
}

void HHEntry::getTemplateHMM(FILE* dbf, char* name, Parameters& par,
                             char use_global_weights, const float qsc,
                             int& format, float* pb, const float S[20][20],
                             const float Sim[20][20], HMM* t) {
  if (dbf != NULL) {
    char line[LINELEN];
    if (!fgetline(line, LINELEN, dbf)) {
      std::cerr << "this should not happen!" << std::endl;
      //TODO: throw error
    }

    while (strscn(line) == NULL)
      fgetline(line, LINELEN, dbf);  // skip lines that contain only white space

    // read HMMER3 format
    if (!strncmp(line, "HMMER3", 6)) {
      format = 1;
      t->ReadHMMer3(dbf, par.showcons, pb, name);
      par.hmmer_used = true;
    }
    // read HMMER format
    else if (!strncmp(line, "HMMER", 5)) {
      format = 1;
      t->ReadHMMer(dbf, par.showcons, pb, name);
      par.hmmer_used = true;
    }
    // read HHM format
    else if (!strncmp(line, "HH", 2)) {
      char path[NAMELEN];
      Pathname(path, name);

      format = 0;
      t->Read(dbf, par.maxcol, par.nseqdis, pb, path);
    }
    //TODO: old hhm format discarded
    // read a3m alignment
    else if (line[0] == '#' || line[0] == '>') {
      Alignment tali;
      tali.Read(dbf, name, par.mark, par.maxcol, par.nseqdis, line);
      tali.Compress(name, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
      //              qali.FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);
      tali.N_filtered = tali.Filter(par.max_seqid_db, S, par.coverage_db,
                                    par.qid_db, qsc, par.Ndiff_db);
      t->name[0] = t->longname[0] = t->fam[0] = '\0';
      tali.FrequenciesAndTransitions(t, use_global_weights, par.mark, par.cons,
                                     par.showcons, par.maxres, pb, Sim);
      format = 0;
    } else {
      std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
                << __func__ << ":" << std::endl;
      std::cerr << "\tunrecognized HMM file format in \'" << name << "\'. \n";
      std::cerr << "Context:\n'" << line << "\n";
      fgetline(line, LINELEN, dbf);
      std::cerr << line << "\n";
      fgetline(line, LINELEN, dbf);
      std::cerr << line << "'\n";
      exit(1);
    }
  }
}

char* HHDatabaseEntry::getName() {
  return entry->name;
}

HHFileEntry::HHFileEntry(char* file)
    : HHEntry(MAXRES) {
  this->file = new char[strlen(file) + 1];
  strcpy(this->file, file);
}

HHFileEntry::~HHFileEntry() {
  delete[] file;
}

void HHFileEntry::getTemplateHMM(Parameters& par, char use_global_weights,
                                 const float qsc, int& format, float* pb,
                                 const float S[20][20], const float Sim[20][20],
                                 HMM* t) {

  FILE * dbf = fopen(file, "r");
  if(dbf == NULL) {
    //TODO: throw error
    HH_LOG(ERROR) << "Template File does not exist: " << file << std::endl;
    exit(1);
  }

  HHEntry::getTemplateHMM(dbf, file, par, use_global_weights, qsc, format, pb,
                          S, Sim, t);
  fclose(dbf);
}

void HHFileEntry::getTemplateA3M(Parameters& par, float* pb,
                                 const float S[20][20], const float Sim[20][20],
                                 Alignment& tali) {

  char line[LINELEN];
  HMM* t = new HMM(MAXSEQDIS, par.maxres);

  FILE* inf = fopen(file, "r");

  if (!fgetline(line, LINELEN, inf)) {
    HH_LOG(ERROR) << "Error in " << __FILE__ << ":" << __LINE__
                            << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\t" << file << " is empty!\n";
    exit(4);
  }

  while (strscn(line) == NULL)
    fgetline(line, LINELEN, inf);  // skip lines that contain only white space

  // Is infile a HMMER file?
  if (!strncmp(line, "HMMER", 5)) {
    // Uncomment this line to allow HMMER2/HMMER3 models as queries:
    HH_LOG(ERROR)
        << "Error: Use of HMMER format as input will result in severe loss of sensitivity!\n";
  }
  // ... or is it an hhm file?
  else if (!strncmp(line, "NAME", 4) || !strncmp(line, "HH", 2)) {
    char path[NAMELEN];
    Pathname(path, file);

    HH_LOG(INFO) << "Query file is in HHM format\n";

    // Rewind to beginning of line and read query hhm file
    rewind(inf);
    t->Read(inf, par.maxcol, par.nseqdis, pb, path);

    Alignment ali_tmp;
    ali_tmp.GetSeqsFromHMM(t);
    ali_tmp.Compress(file, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
    tali = ali_tmp;
  }
  // ... or is it an alignment file
  else if (line[0] == '#' || line[0] == '>') {
    if (par.calibrate) {
      HH_LOG(ERROR) << "Error in " << __FILE__ << ":" << __LINE__
                              << ": " << __func__ << ":" << std::endl;
      HH_LOG(ERROR) << "\tonly HHM files can be calibrated.\n";
      HH_LOG(ERROR)
          << "\tBuild an HHM file from your alignment with hhmake and rerun with the hhm file"
          << std::endl;
      exit(1);
    }

    Alignment ali_tmp;

    // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
    ali_tmp.Read(inf, file, par.mark, par.maxcol, par.nseqdis, line);

    // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
    // and store marked sequences in name[k] and seq[k]
    ali_tmp.Compress(file, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);

    tali = ali_tmp;
  } else {
    HH_LOG(ERROR) << "Error in " << __FILE__ << ":" << __LINE__
                            << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\tunrecognized input file format in \'" << file
                            << "\'\n";
    HH_LOG(ERROR) << "\tline = " << line << "\n";
    exit(1);
  }

  fclose(inf);

  delete t;
}

char* HHFileEntry::getName() {
  return file;
}

int getMaxTemplateLength(std::vector<HHEntry*>& entries) {
  int max_template_length = 0;

  for (size_t i = 0; i < entries.size(); i++) {
    max_template_length = std::max(max_template_length,
                                   entries[i]->sequence_length);
  }

  return max_template_length;
}

