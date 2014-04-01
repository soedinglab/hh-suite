/*
 * hhblits_mpi.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */



#include <mpi.h>
#include <stdio.h>


extern "C" {
  #include <ffindex.h>     // fast index-based database reading
}

int main(int argc, char **argv) {
  int mpi_error,
      mpi_rank,
      mpi_num_procs;

  mpi_error = MPI_Init(&argc, &argv);
  mpi_error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  mpi_error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);

  //TODO read parameters
  char *data_filename  = NULL;
  char *index_filename = NULL;
  char *data_filename_out  = NULL;
  char *index_filename_out = NULL;

  FILE *data_file  = fopen(data_filename,  "r");
  FILE *index_file = fopen(index_filename, "r");

  if( data_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], data_filename);  exit(EXIT_FAILURE); }
  if(index_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], index_filename);  exit(EXIT_FAILURE); }

  //TODO: must be done for all outputs
  FILE *data_file_out = NULL, *index_file_out = NULL;
  if(data_filename_out != NULL && index_filename_out != NULL)
  {
    char* data_filename_out_rank  = malloc(FILENAME_MAX);
    char* index_filename_out_rank = malloc(FILENAME_MAX);
    snprintf( data_filename_out_rank, FILENAME_MAX, "%s.%d", data_filename_out,  mpi_rank);
    snprintf(index_filename_out_rank, FILENAME_MAX, "%s.%d", index_filename_out, mpi_rank);
    data_file_out  = fopen(data_filename_out_rank,  "w+");
    index_file_out = fopen(index_filename_out_rank, "w+");

    if( data_file_out == NULL) { fferror_print(__FILE__, __LINE__, argv[0], data_filename_out);  exit(EXIT_FAILURE); }
    if(index_file_out == NULL) { fferror_print(__FILE__, __LINE__, argv[0], index_filename_out);  exit(EXIT_FAILURE); }
  }

  //init input dataframe
  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  ffindex_index_t* index = ffindex_index_parse(index_file, 0);
  if(index == NULL)
  {
    fferror_print(__FILE__, __LINE__, "ffindex_index_parse", index_filename);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  size_t batch_size, range_start, range_end;
  if(index->n_entries >= mpi_num_procs)
    batch_size = index->n_entries / mpi_num_procs;
  else
    batch_size = 0;
  range_start = mpi_rank * batch_size;
  range_end = range_start + batch_size;

  size_t offset = 0;
  // Foreach entry
  if(batch_size > 0) {
    for(size_t entry_index = range_start; entry_index < range_end; entry_index++)
    {
      ffindex_entry_t* entry = ffindex_get_entry_by_index(index, entry_index);
      if(entry == NULL) {
        continue;
      }

      //TODO run hhblits class...
      //TODO for each specified result type
        //TODO get results of hhblits
        //TODO print results to output ffindeces
    }
  }

  ssize_t left_over = index->n_entries - (batch_size * mpi_num_procs);
  if(mpi_rank < left_over) {
    size_t left_over_entry_index = (batch_size * mpi_num_procs) + mpi_rank;
    ffindex_entry_t* entry = ffindex_get_entry_by_index(index, left_over_entry_index);
    if(entry == NULL) {
      perror(entry->name);
    }

    //TODO run hhblits class...
    //TODO for each specified result type
      //TODO get results of hhblits
      //TODO print results to output ffindeces
  }

  fclose(data_file_out);
  fclose(index_file_out);

  MPI_Barrier(MPI_COMM_WORLD);

  //TODO merge all specified ffindex databases
  // merge FFindexes in master
  if(data_filename_out != NULL && mpi_rank == 0) {
    char* merge_command  = malloc(FILENAME_MAX * 5);
    for(int i = 0; i < mpi_num_procs; i++) {
      snprintf( merge_command, FILENAME_MAX, "ffindex_build -as %s %s -d %s.%d -i %s.%d",
                data_filename_out, index_filename_out, data_filename_out, i, index_filename_out, i);
      //puts(merge_command);
      if(system(merge_command) == 0)
      {
        snprintf(merge_command, FILENAME_MAX, "%s.%d", data_filename_out, i);
        unlink(merge_command);
        snprintf(merge_command, FILENAME_MAX, "%s.%d", index_filename_out, i);
        unlink(merge_command);
      }
    }
  }

  MPI_Finalize();

}

