/*
 * hhblits_mpi.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include <stdio.h>
#include <sys/mman.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "hhdecl.h"
#include "hhblits.h"

extern "C" {
#include <ffindex.h>
}

struct OutputFFIndex {
    char base[NAMELEN];
    FILE* data_fh;
    FILE* index_fh;
    size_t offset;
    size_t number_entries;
    void (*print)(HHblits&, std::stringstream&);

    void close() {
      fclose(data_fh);
      fclose(index_fh);
    }

    void saveOutput(HHblits& hhblits, char* name) {
      std::stringstream out;
      print(hhblits, out);

      std::string tmp = out.str();
      ffindex_insert_memory(data_fh, index_fh, &offset,
          const_cast<char*>(tmp.c_str()), tmp.size(), name);

      fflush(data_fh);
      fflush(index_fh);
      number_entries++;
    }

    void sort() {
      /* Sort the index entries and write back */
      char index_filename[NAMELEN];
      snprintf(index_filename, FILENAME_MAX, "%s.ffindex", base);

      rewind(index_fh);
      ffindex_index_t* index = ffindex_index_parse(index_fh, number_entries);
      fclose(index_fh);

      if (index == NULL) {
        HH_LOG(ERROR) << "Could not read index from " << index_filename << " for sorting!" << std::endl;
        return;
      }

      ffindex_sort_index_file(index);

      index_fh = fopen(index_filename, "w");

      if (index_fh == NULL) {
        HH_LOG(ERROR) << "Could not open " << index_filename << " for sorting!" << std::endl;
        return;
      }

      ffindex_write(index, index_fh);
      free(index);
    }
};


void makeOutputFFIndex(char* par, void (*print)(HHblits&, std::stringstream&),
    std::vector<OutputFFIndex>& outDatabases) {
  if (*par) {
    OutputFFIndex db;

    strcpy(db.base, par);
    db.offset = 0;
    db.print = print;
    db.number_entries = 0;

    char data_filename_out_rank[NAMELEN];
    char index_filename_out_rank[NAMELEN];

    snprintf(data_filename_out_rank, FILENAME_MAX, "%s.ffdata", par);
    snprintf(index_filename_out_rank, FILENAME_MAX, "%s.ffindex", par);

    db.data_fh = fopen(data_filename_out_rank, "w+");
    db.index_fh = fopen(index_filename_out_rank, "w+");

    if (db.data_fh == NULL) {
      HH_LOG(ERROR) << "Could not open datafile " << data_filename_out_rank << "!" << std::endl;
      return;
    }

    if (db.index_fh == NULL) {
      HH_LOG(ERROR) << "Could not open indexfile " << index_filename_out_rank << "!" << std::endl;
      return;
    }

    outDatabases.push_back(db);
  }
}


void checkOutput(Parameters& par) {
  if (!*par.outfile) {
    RemoveExtension(par.outfile, par.infile);
    strcat(par.outfile, "_hhr");
    HH_LOG(INFO) << "Search results will be written to " << par.outfile << "!\n";
  }
}


int main(int argc, char **argv) {
  Parameters par;
  HHblits::ProcessAllArguments(argc, argv, par);

  char data_filename[NAMELEN];
  char index_filename[NAMELEN];

  strcpy(data_filename, par.infile);
  strcat(data_filename, ".ffdata");

  strcpy(index_filename, par.infile);
  strcat(index_filename, ".ffindex");

  FILE *data_file = fopen(data_filename, "r");
  FILE *index_file = fopen(index_filename, "r");

  if (data_file == NULL) {
    HH_LOG(ERROR) << "Input data file " << data_filename << " does not exist!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (index_file == NULL) {
    HH_LOG(ERROR) << "Input index file " << index_filename << " does not exist!" << std::endl;
    exit(EXIT_FAILURE);
  }

  //init input ffindex
  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  size_t number_input_index_lines = CountLinesInFile(index_filename);
  ffindex_index_t* index = ffindex_index_parse(index_file, number_input_index_lines);
  if (index == NULL) {
    HH_LOG(ERROR) << "Could not parse index from " << index_filename << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<OutputFFIndex> outputDatabases;
  makeOutputFFIndex(par.outfile, &HHblits::writeHHRFile, outputDatabases);
  makeOutputFFIndex(par.scorefile, &HHblits::writeScoresFile, outputDatabases);
  makeOutputFFIndex(par.pairwisealisfile, &HHblits::writePairwiseAlisFile, outputDatabases);
  makeOutputFFIndex(par.alitabfile, &HHblits::writeAlitabFile, outputDatabases);
  makeOutputFFIndex(par.psifile, &HHblits::writePsiFile, outputDatabases);
  makeOutputFFIndex(par.hhmfile, &HHblits::writeHMMFile, outputDatabases);
  makeOutputFFIndex(par.alnfile, &HHblits::writeA3MFile, outputDatabases);
  makeOutputFFIndex(par.matrices_output_file, &HHblits::writeMatricesFile, outputDatabases);

  std::vector<HHblitsDatabase*> databases;
  HHblits::prepareDatabases(par, databases);

  //need to save for cleanup in the end
  int threads = par.threads;

#ifdef OPENMP
  omp_set_num_threads(threads);
#endif

  //no openmp parallelization in hhblits methods

  HHblits* hhblits_instances[MAXBINS];
  par.threads = 1;

  for(int i = 0; i < threads; i++) {
    hhblits_instances[i] = new HHblits(par, databases);
  }

  size_t range_start = 0;
  size_t range_end = index->n_entries;
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t entry_index = range_start; entry_index < range_end; entry_index++) {
    ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);

    if (entry == NULL) {
      HH_LOG(WARNING) << "Could not open entry " << entry_index << " from input ffindex!" << std::endl;
      continue;
    }

    int bin = 0;
#ifdef OPENMP
    bin = omp_get_thread_num();
    omp_set_num_threads(1);
#endif

    FILE* inf = ffindex_fopen_by_entry(data, entry);
    if(inf == NULL) {
      HH_LOG(WARNING) << "Could not open input entry (" << entry->name << ")!" << std::endl;
      continue;
    }

    HH_LOG(INFO) << "Thread " << bin << "\t" << entry->name << std::endl;
    hhblits_instances[bin]->run(inf, entry->name);

    #pragma omp critical
    {
      for (size_t i = 0; i < outputDatabases.size(); i++) {
        outputDatabases[i].saveOutput(*hhblits_instances[bin], entry->name);
      }
    }

    hhblits_instances[bin]->Reset();
  }

  munmap(data, data_size);
  free(index);

  fclose(data_file);
  fclose(index_file);

  for(int i = 0; i < threads; i++) {
    delete hhblits_instances[i];
  }

  for (size_t i = 0; i < databases.size(); i++) {
    delete databases[i];
  }
  databases.clear();

  for (size_t i = 0; i < outputDatabases.size(); i++) {
    outputDatabases[i].sort();
    outputDatabases[i].close();
  }
}

