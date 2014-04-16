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
#include <ffindex.h>     // fast index-based database reading
}

void initOutputFFDatabase(char* par, int mpi_rank, bool& print,
    char* data_filename_out_rank, char* index_filename_out_rank,
    FILE*& data_file_out, FILE*& index_file_out) {
  if (*par) {
    data_filename_out_rank = new char[FILENAME_MAX];
    index_filename_out_rank = new char[FILENAME_MAX];

    snprintf(data_filename_out_rank, FILENAME_MAX, "%s.%d.ffdata", par, mpi_rank);
    snprintf(index_filename_out_rank, FILENAME_MAX, "%s.%d.ffindex", par, mpi_rank);

    data_file_out = fopen(data_filename_out_rank, "w+");
    index_file_out = fopen(index_filename_out_rank, "w+");

    if (data_file_out == NULL) {
      std::cerr << "could not open datafile " << data_filename_out_rank << std::endl;
      //TODO: throw error
      print = false;
      return;
    }

    if (index_file_out == NULL) {
      std::cerr << "could not open indexfile " << index_filename_out_rank << std::endl;
      //TODO: throw error
      print = false;
      return;
    }

    print = true;
  }
}

void saveOutput(std::stringstream& out, FILE* data_file_out,
    FILE* index_file_out, char* name, size_t& offset) {
  std::string tmp = out.str();
  ffindex_insert_memory(data_file_out, index_file_out, &offset,
      const_cast<char*>(tmp.c_str()), tmp.size(), name);
}

void saveAlisOutput(std::map<int, Alignment>& alis, FILE* data_file_out,
    FILE* index_file_out, char* name, size_t& offset) {
  std::map<int, Alignment>::iterator it;
  for (it = alis.begin(); it != alis.end(); it++) {
    std::stringstream ss_tmp;
    ss_tmp << name << "_" << (*it).first << ".a3m";
    std::string id = ss_tmp.str();

    std::stringstream content;
    (*it).second.WriteToFile(content, "a3m");

    std::string content_str = content.str();

    ffindex_insert_memory(data_file_out, index_file_out, &offset,
        const_cast<char*>(content_str.c_str()), content_str.length(),
        const_cast<char*>(id.c_str()));
  }
}

void closeFH(bool print, FILE* data_file_out, FILE* index_file_out) {
  if (print) {
    fclose(data_file_out);
    fclose(index_file_out);
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

  bool print_hhr = false, print_alis = false, print_scores = false,
      print_pairwise_alis = false, print_alitab = false, print_reduced_hhr =
          false, print_psi = false, print_hmm = false, print_a3m = false;
  char* hhr_data_file, *alis_data_file, *scores_data_file,
      *pairwise_alis_data_file, *alitab_data_file, *reduced_hhr_data_file,
      *psi_data_file, *hmm_data_file, *a3m_data_file;
  char* hhr_index_file, *alis_index_file, *scores_index_file,
      *pairwise_alis_index_file, *alitab_index_file, *reduced_hhr_index_file,
      *psi_index_file, *hmm_index_file, *a3m_index_file;
  FILE* hhr_data_fh = NULL, *alis_data_fh = NULL, *scores_data_fh = NULL,
      *pairwise_alis_data_fh = NULL, *alitab_data_fh = NULL,
      *reduced_hhr_data_fh = NULL, *psi_data_fh = NULL, *hmm_data_fh = NULL,
      *a3m_data_fh = NULL;
  FILE* hhr_index_fh = NULL, *alis_index_fh = NULL, *scores_index_fh = NULL,
      *pairwise_alis_index_fh = NULL, *alitab_index_fh = NULL,
      *reduced_hhr_index_fh = NULL, *psi_index_fh = NULL, *hmm_index_fh = NULL,
      *a3m_index_fh = NULL;
  size_t hhr_index_offset = 0, alis_index_offset = 0, scores_index_offset = 0,
      pairwise_alis_index_offset = 0, alitab_index_offset = 0,
      reduced_hhr_index_offset = 0, psi_index_offset = 0, hmm_index_offset = 0,
      a3m_index_offset = 0;

  initOutputFFDatabase(par.outfile, mpi_rank, print_hhr, hhr_data_file,
      hhr_index_file, hhr_data_fh, hhr_index_fh);
  initOutputFFDatabase(par.alisbasename, mpi_rank, print_alis, alis_data_file,
      alis_index_file, alis_data_fh, alis_index_fh);
  initOutputFFDatabase(par.scorefile, mpi_rank, print_scores, scores_data_file,
      scores_index_file, scores_data_fh, scores_index_fh);
  initOutputFFDatabase(par.pairwisealisfile, mpi_rank, print_pairwise_alis,
      pairwise_alis_data_file, pairwise_alis_index_file, pairwise_alis_data_fh,
      pairwise_alis_index_fh);
  initOutputFFDatabase(par.alitabfile, mpi_rank, print_alitab, alitab_data_file,
      alitab_index_file, alitab_data_fh, alitab_index_fh);
  initOutputFFDatabase(par.reduced_outfile, mpi_rank, print_reduced_hhr,
      reduced_hhr_data_file, reduced_hhr_index_file, reduced_hhr_data_fh,
      reduced_hhr_index_fh);
  initOutputFFDatabase(par.psifile, mpi_rank, print_psi, psi_data_file,
      psi_index_file, psi_data_fh, psi_index_fh);
  initOutputFFDatabase(par.hhmfile, mpi_rank, print_hmm, hmm_data_file,
      hmm_index_file, hmm_data_fh, hmm_index_fh);
  initOutputFFDatabase(par.alnfile, mpi_rank, print_a3m, a3m_data_file,
      a3m_index_file, a3m_data_fh, a3m_index_fh);

  HHblits hhblits(par);

  size_t batch_size, range_start, range_end;
  if (index->n_entries >= mpi_num_procs)
    batch_size = index->n_entries / mpi_num_procs;
  else
    batch_size = 0;

  range_start = mpi_rank * batch_size;
  range_end = range_start + batch_size;

  // Foreach entry
  if (batch_size > 0) {
    for (size_t entry_index = range_start; entry_index < range_end;
        entry_index++) {
      ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);
      if (entry == NULL) {
        continue;
      }

      FILE* inf = ffindex_fopen_by_entry(data, entry);
      hhblits.Reset();
      hhblits.run(inf, entry->name);

      if (print_hhr) {
        std::stringstream out;
        hhblits.writeHHRFile(out);
        saveOutput(out, hhr_data_fh, hhr_index_fh,
            entry->name, hhr_index_offset);
      }
      if (print_scores) {
        std::stringstream out;
        hhblits.writeScoresFile(out);
        saveOutput(out, scores_data_fh, scores_index_fh,
            entry->name, scores_index_offset);
      }
      if (print_pairwise_alis) {
        std::stringstream out;
        hhblits.writePairwiseAlisFile(par.outformat, out);
        saveOutput(out,
            pairwise_alis_data_fh, pairwise_alis_index_fh, entry->name,
            pairwise_alis_index_offset);
      }
      if (print_alitab) {
        std::stringstream out;
        hhblits.writeAlitabFile(out);
        saveOutput(out, alitab_data_fh, alitab_index_fh,
            entry->name, alitab_index_offset);
      }
      if (print_reduced_hhr) {
        std::stringstream out;
        hhblits.writeReducedHHRFile(out);
        saveOutput(out, reduced_hhr_data_fh,
            reduced_hhr_index_fh, entry->name, reduced_hhr_index_offset);
      }
      if (print_psi) {
        std::stringstream out;
        hhblits.writePsiFile(out);
        saveOutput(out, psi_data_fh, psi_index_fh,
            entry->name, psi_index_offset);
      }
      if (print_hmm) {
        std::stringstream out;
        hhblits.writeHMMFile(out);
        saveOutput(out, hmm_data_fh, hmm_index_fh,
            entry->name, hmm_index_offset);
      }
      if (print_a3m) {
        std::stringstream out;
        hhblits.writeA3MFile(out);
        saveOutput(out, a3m_data_fh, a3m_index_fh,
            entry->name, a3m_index_offset);
      }
      if (print_alis)
        saveAlisOutput(hhblits.getAlis(), alis_data_fh, alis_data_fh,
            entry->name, alis_index_offset);
    }
  }

  ssize_t left_over = index->n_entries - (batch_size * mpi_num_procs);
  if (mpi_rank < left_over) {
    size_t left_over_entry_index = (batch_size * mpi_num_procs) + mpi_rank;
    ffindex_entry_t* entry = ffindex_get_entry_by_index(index,
        left_over_entry_index);

    FILE* inf = ffindex_fopen_by_entry(data, entry);
    hhblits.Reset();
    hhblits.run(inf, entry->name);

    if (print_hhr) {
      std::stringstream out;
      hhblits.writeHHRFile(out);
      saveOutput(out, hhr_data_fh, hhr_index_fh,
          entry->name, hhr_index_offset);
    }
    if (print_scores) {
      std::stringstream out;
      hhblits.writeScoresFile(out);
      saveOutput(out, scores_data_fh, scores_index_fh,
          entry->name, scores_index_offset);
    }
    if (print_pairwise_alis) {
      std::stringstream out;
      hhblits.writePairwiseAlisFile(par.outformat, out);
      saveOutput(out,
          pairwise_alis_data_fh, pairwise_alis_index_fh, entry->name,
          pairwise_alis_index_offset);
    }
    if (print_alitab) {
      std::stringstream out;
      hhblits.writeAlitabFile(out);
      saveOutput(out, alitab_data_fh, alitab_index_fh,
          entry->name, alitab_index_offset);
    }
    if (print_reduced_hhr) {
      std::stringstream out;
      hhblits.writeReducedHHRFile(out);
      saveOutput(out, reduced_hhr_data_fh,
          reduced_hhr_index_fh, entry->name, reduced_hhr_index_offset);
    }
    if (print_psi) {
      std::stringstream out;
      hhblits.writePsiFile(out);
      saveOutput(out, psi_data_fh, psi_index_fh,
          entry->name, psi_index_offset);
    }
    if (print_hmm) {
      std::stringstream out;
      hhblits.writeHMMFile(out);
      saveOutput(out, hmm_data_fh, hmm_index_fh,
          entry->name, hmm_index_offset);
    }
    if (print_a3m) {
      std::stringstream out;
      hhblits.writeA3MFile(out);
      saveOutput(out, a3m_data_fh, a3m_index_fh,
          entry->name, a3m_index_offset);
    }
    if (print_alis)
      saveAlisOutput(hhblits.getAlis(), alis_data_fh, alis_data_fh,
          entry->name, alis_index_offset);
  }

  fclose(data_file);
  fclose(index_file);

  closeFH(print_hhr, hhr_data_fh, hhr_index_fh);
  closeFH(print_scores, scores_data_fh, scores_index_fh);
  closeFH(print_pairwise_alis, pairwise_alis_data_fh, pairwise_alis_index_fh);
  closeFH(print_alitab, alitab_data_fh, alitab_index_fh);
  closeFH(print_reduced_hhr, reduced_hhr_data_fh, reduced_hhr_index_fh);
  closeFH(print_psi, psi_data_fh, psi_index_fh);
  closeFH(print_hmm, hmm_data_fh, hmm_index_fh);
  closeFH(print_a3m, a3m_data_fh, a3m_index_fh);
  closeFH(print_alis, alis_data_fh, alis_index_fh);

  MPI_Barrier(MPI_COMM_WORLD );

////TODO merge all specified ffindex databases
//// merge FFindexes in master
//if (data_filename_out != NULL && mpi_rank == 0) {
//	char* merge_command = malloc(FILENAME_MAX * 5);
//	for (int i = 0; i < mpi_num_procs; i++) {
//		snprintf(merge_command, FILENAME_MAX,
//				"ffindex_build -as %s %s -d %s.%d -i %s.%d", data_filename_out,
//				index_filename_out, data_filename_out, i, index_filename_out,
//				i);
//		//puts(merge_command);
//		if (system(merge_command) == 0) {
//			snprintf(merge_command, FILENAME_MAX, "%s.%d", data_filename_out,
//					i);
//			unlink(merge_command);
//			snprintf(merge_command, FILENAME_MAX, "%s.%d", index_filename_out,
//					i);
//			unlink(merge_command);
//		}
//	}
//}

  MPI_Finalize();

}

