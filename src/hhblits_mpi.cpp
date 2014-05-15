/*
 * hhblits_mpi.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include <mpi.h>
#include <stdio.h>
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
    }

    void merge(const int mpi_num_procs) {
      char data_filename[FILENAME_MAX];
      char index_filename[FILENAME_MAX];

      snprintf(data_filename, FILENAME_MAX, "%s.ffdata", base);
      snprintf(index_filename, FILENAME_MAX, "%s.ffindex", base);

      FILE *data_file = NULL, *index_file = NULL;
      size_t offset = 0;

      ffindex_index_open(data_filename, index_filename, "w", &data_file, &index_file, &offset);
      //TODO: check error

      /* Append other ffindexes */
      for (int i = 0; i < mpi_num_procs; i++) {
        char data_file_name_to_add[FILENAME_MAX];
        char index_file_name_to_add[FILENAME_MAX];

        snprintf(data_file_name_to_add, FILENAME_MAX, "%s.%d.ffdata", base, i);
        snprintf(index_file_name_to_add, FILENAME_MAX, "%s.%d.ffindex", base,
            i);

        FILE* data_file_to_add = fopen(data_file_name_to_add, "r");
        if (data_file_to_add == NULL) {
          //TODO: throw error
        }

        FILE* index_file_to_add = fopen(index_file_name_to_add, "r");
        if (index_file_to_add == NULL) {
          //TODO: throw error
        }

        size_t data_size;
        char *data_to_add = ffindex_mmap_data(data_file_to_add, &data_size);
        ffindex_index_t* index_to_add = ffindex_index_parse(index_file_to_add,
            0);

        ffindex_insert_ffindex(data_file, index_file, &offset, data_to_add,
            index_to_add);
        //TODO: free memory???
      }

      fclose(data_file);

      /* Sort the index entries and write back */
      rewind(index_file);
      ffindex_index_t* index = ffindex_index_parse(index_file, 0);
      if (index == NULL) {
        //TODO: throw error
      }
      fclose(index_file);
      ffindex_sort_index_file(index);
      index_file = fopen(index_filename, "w");
      if (index_file == NULL) {
        //TODO: throw error
      }
      ffindex_write(index, index_file);

      fclose(index_file);
      //TODO: free memory???
    }
};

void makeOutputFFIndex(char* par, const int mpi_rank,
    void (*print)(HHblits&, std::stringstream&),
    std::vector<OutputFFIndex>& outDatabases) {
  if (*par) {
    OutputFFIndex db;

    strcpy(db.base, par);
    db.offset = 0;
    db.print = print;

    char data_filename_out_rank[NAMELEN];
    char index_filename_out_rank[NAMELEN];

    snprintf(data_filename_out_rank, FILENAME_MAX, "%s.%d.ffdata", par,
        mpi_rank);
    snprintf(index_filename_out_rank, FILENAME_MAX, "%s.%d.ffindex", par,
        mpi_rank);

    db.data_fh = fopen(data_filename_out_rank, "w+");
    db.index_fh = fopen(index_filename_out_rank, "w+");

    if (db.data_fh == NULL) {
      std::cerr << "could not open datafile " << data_filename_out_rank
          << std::endl;
      return;
    }

    if (db.index_fh == NULL) {
      std::cerr << "could not open indexfile " << index_filename_out_rank
          << std::endl;
      return;
    }

    outDatabases.push_back(db);
  }
}

int main(int argc, char **argv) {
  int mpi_error, mpi_rank, mpi_num_procs;

  mpi_error = MPI_Init(&argc, &argv);
  mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  mpi_error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);

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
    std::cerr << "input data file " << data_filename << " does not exist!"
        << std::endl;
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  if (index_file == NULL) {
    std::cerr << "input index file " << index_filename << " does not exist!"
        << std::endl;
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  //init input ffindex
  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  ffindex_index_t* index = ffindex_index_parse(index_file, 0);
  if (index == NULL) {
    std::cerr << "Could not parse index from " << index_filename << std::endl;
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  std::vector<OutputFFIndex> outputDatabases;
  makeOutputFFIndex(par.outfile, mpi_rank, &HHblits::writeHHRFile,
      outputDatabases);
  makeOutputFFIndex(par.scorefile, mpi_rank, &HHblits::writeScoresFile,
      outputDatabases);
  makeOutputFFIndex(par.pairwisealisfile, mpi_rank,
      &HHblits::writePairwiseAlisFile, outputDatabases);
  makeOutputFFIndex(par.alitabfile, mpi_rank, &HHblits::writeAlitabFile,
      outputDatabases);
  makeOutputFFIndex(par.reduced_outfile, mpi_rank,
      &HHblits::writeReducedHHRFile, outputDatabases);
  makeOutputFFIndex(par.psifile, mpi_rank, &HHblits::writePsiFile,
      outputDatabases);
  makeOutputFFIndex(par.hhmfile, mpi_rank, &HHblits::writeHMMFile,
      outputDatabases);
  makeOutputFFIndex(par.alnfile, mpi_rank, &HHblits::writeA3MFile,
      outputDatabases);

//	initOutputFFDatabase(par.alisbasename, mpi_rank, print_alis, alis_data_fh, alis_index_fh);

  size_t batch_size, range_start, range_end;
  if (index->n_entries >= mpi_num_procs)
    batch_size = index->n_entries / mpi_num_procs;
  else
    batch_size = 0;

  range_start = mpi_rank * batch_size;
  range_end = range_start + batch_size;

  HHblits hhblits(par);

  // Foreach entry
  if (batch_size > 0) {
    for (size_t entry_index = range_start; entry_index < range_end;
        entry_index++) {
      ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);
      if (entry == NULL) {
        continue;
      }

      hhblits.Reset();

      FILE* inf = ffindex_fopen_by_entry(data, entry);
      hhblits.run(inf, entry->name);

      for (size_t i = 0; i < outputDatabases.size(); i++) {
        outputDatabases[i].saveOutput(hhblits, entry->name);
      }

//			if (print_alis)
//				saveAlisOutput(hhblits.getAlis(), alis_data_fh, alis_data_fh,
//						entry->name, alis_index_offset);
    }
  }

  ssize_t left_over = index->n_entries - (batch_size * mpi_num_procs);
  if (mpi_rank < left_over) {
    size_t left_over_entry_index = (batch_size * mpi_num_procs) + mpi_rank;
    ffindex_entry_t* entry = ffindex_get_entry_by_index(index,
        left_over_entry_index);

    hhblits.Reset();

    FILE* inf = ffindex_fopen_by_entry(data, entry);
    hhblits.run(inf, entry->name);

    for (size_t i = 0; i < outputDatabases.size(); i++) {
      outputDatabases[i].saveOutput(hhblits, entry->name);
    }

//		if (print_alis)
//			saveAlisOutput(hhblits.getAlis(), alis_data_fh, alis_data_fh,
//					entry->name, alis_index_offset);
  }

  fclose(data_file);
  fclose(index_file);

  for (size_t i = 0; i < outputDatabases.size(); i++) {
    outputDatabases[i].close();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    for (size_t i = 0; i < outputDatabases.size(); i++) {
      outputDatabases[i].merge(mpi_num_procs);
    }
  }

  MPI_Finalize();
}

