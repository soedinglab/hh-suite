/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "cstranslate_mpi_app.h"

namespace cs {
  template<class Abc>
  int CSTranslateMpiApp<Abc>::Run() {
    ignore_signal(SIGPIPE);

    const bool isCa3m = this->opts_.informat == "ca3m";

    std::string input_data_file = this->opts_.infile + ".ffdata";
    std::string input_index_file = this->opts_.infile + ".ffindex";

    // required ffindex files for ca3m case
    size_t input_header_data_size = 0;
    FILE *input_header_data_fh = NULL, *input_header_index_fh = NULL;
    std::string input_header_data_file, input_header_index_file;
    char *input_header_data = NULL;
    ffindex_index_t *input_header_index = NULL;

    size_t input_sequence_data_size = 0;
    FILE *input_sequence_data_fh = NULL, *input_sequence_index_fh = NULL;
    std::string input_sequence_data_file, input_sequence_index_file;
    char *input_sequence_data = NULL;
    ffindex_index_t *input_sequence_index = NULL;

    if (isCa3m) {
      // infile has to be the ffindex basepath with no suffices
      input_data_file = this->opts_.infile + "_ca3m.ffdata";
      input_index_file = this->opts_.infile + "_ca3m.ffindex";

      input_header_data_file = this->opts_.infile + "_header.ffdata";
      input_header_index_file = this->opts_.infile + "_header.ffindex";

      input_header_data_fh = fopen(input_header_data_file.c_str(), "r");
      input_header_index_fh = fopen(input_header_index_file.c_str(), "r");

      if (input_header_data_fh == NULL) {
        fprintf(this->out_, "Could not open ffindex input data file! (%s)!\n", input_header_data_file.c_str());
        exit(1);
      }
      if (input_header_index_fh == NULL) {
        fprintf(this->out_, "Could not open ffindex input index file! (%s)!\n", input_header_index_file.c_str());
        exit(1);
      }

      input_header_data = ffindex_mmap_data(input_header_data_fh, &input_header_data_size);
      size_t header_entries = ffcount_lines(input_header_index_file.c_str());
      input_header_index = ffindex_index_parse(input_header_index_fh, header_entries);

      if (input_header_index == NULL) {
        fprintf(this->out_, "Input index could not be mapped!\n");
        exit(1);
      }

      input_sequence_data_file = this->opts_.infile + "_sequence.ffdata";
      input_sequence_index_file = this->opts_.infile + "_sequence.ffindex";

      input_sequence_data_fh = fopen(input_sequence_data_file.c_str(), "r");
      input_sequence_index_fh = fopen(input_sequence_index_file.c_str(), "r");

      if (input_sequence_data_fh == NULL) {
        fprintf(this->out_, "Could not open ffindex input data file! (%s)!\n", input_sequence_data_file.c_str());
        exit(1);
      }
      if (input_sequence_index_fh == NULL) {
        fprintf(this->out_, "Could not open ffindex input index file! (%s)!\n", input_sequence_index_file.c_str());
        exit(1);
      }

      input_sequence_data = ffindex_mmap_data(input_sequence_data_fh, &input_sequence_data_size);
      size_t sequence_entries = ffcount_lines(input_sequence_index_file.c_str());
      input_sequence_index = ffindex_index_parse(input_sequence_index_fh, sequence_entries);

      if (input_sequence_index == NULL) {
        fprintf(this->out_, "Input index could not be mapped!\n");
        exit(1);
      }
    }

    FILE *input_data_fh = fopen(input_data_file.c_str(), "r");
    FILE *input_index_fh = fopen(input_index_file.c_str(), "r");

    if (input_data_fh == NULL) {
      fprintf(this->out_, "Could not open ffindex input data file! (%s)!\n", input_data_file.c_str());
      exit(1);
    }
    if (input_index_fh == NULL) {
      fprintf(this->out_, "Could not open ffindex input index file! (%s)!\n", input_index_file.c_str());
      exit(1);
    }

    size_t input_data_size;
    char *input_data = ffindex_mmap_data(input_data_fh, &input_data_size);
    size_t entries = ffcount_lines(input_index_file.c_str());
    ffindex_index_t *input_index = ffindex_index_parse(input_index_fh, entries);

    if (input_index == NULL) {
      fprintf(this->out_, "Input index could not be mapped!\n");
      exit(1);
    }

    //prepare output ffindex cs219 database
    std::string data_filename_out = this->opts_.outfile + ".ffdata";
    std::string index_filename_out = this->opts_.outfile + ".ffindex";
    std::string log_filename_out = this->opts_.outfile + ".log";

    int mpq_status = MPQ_Init(this->argc_, this->argv_, input_index->n_entries);
    if (mpq_status == MPQ_SUCCESS) {

      if (MPQ_rank == MPQ_MASTER) {
        MPQ_Master(1);
      } else {
        FILE *data_file_out = NULL, *index_file_out = NULL, *log_file_out = NULL;
        if (MPQ_rank != MPQ_MASTER) {
          char data_filename_out_rank[FILENAME_MAX];
          snprintf(data_filename_out_rank, FILENAME_MAX, "%s.%d", data_filename_out.c_str(), MPQ_rank);

          data_file_out = fopen(data_filename_out_rank, "w+");
          if (data_file_out == NULL) {
            fprintf(this->out_, "Could not open ffindex output file! (%s)!\n", data_filename_out_rank);
            exit(1);
          }

          char index_filename_out_rank[FILENAME_MAX];
          snprintf(index_filename_out_rank, FILENAME_MAX, "%s.%d", index_filename_out.c_str(), MPQ_rank);

          index_file_out = fopen(index_filename_out_rank, "w+");
          if (index_file_out == NULL) {
            fprintf(this->out_, "Could not open ffindex output file! (%s)!\n", index_filename_out_rank);
            exit(1);
          }

          char log_filename_out_rank[FILENAME_MAX];
          snprintf(log_filename_out_rank, FILENAME_MAX, "%s.%d",
                   log_filename_out.c_str(), MPQ_rank);

          log_file_out = fopen(log_filename_out_rank, "w+");
          if (log_file_out == NULL) {
            fprintf(this->out_, "Could not open ffindex output file! (%s)!\n", log_filename_out_rank);
            exit(1);
          }
        }

        this->SetupPseudocountEngine();
        this->SetupAbstractStateEngine();

        // Create emission functor needed for translation into abstract states
        Emission<Abc> emission(1, this->opts_.weight_as, 1.0);

        CSTranslateMpiEnvironment<Abc> *env = new CSTranslateMpiEnvironment<Abc>
            (input_index, input_data, input_sequence_index,
             input_sequence_data, input_header_index,
             input_header_data, data_file_out, index_file_out, emission, 0,
             isCa3m, false, log_file_out, *this);

        MPQ_Worker(CSTranslateMpiPayload<Abc>, env);
        if (env->log_file) {
          fclose(env->log_file);
        }
        if (env->index_file_out) {
          fclose(env->index_file_out);
        }
        if (env->data_file_out) {
          fclose(env->data_file_out);
        }
        free(env);
      }
      MPQ_Finalize();

      if (MPQ_rank == MPQ_MASTER) {
        ffmerge_splits(data_filename_out.c_str(), index_filename_out.c_str(), MPQ_size, 1);
      }
    } else {
      if (mpq_status == MPQ_ERROR_NO_WORKERS) {
        fprintf(stderr, "MPQ_Init: Needs at least one worker process.\n");
      } else if (mpq_status == MPQ_ERROR_TOO_MANY_WORKERS) {
        fprintf(stderr, "MPQ_Init: Too many worker processes.\n");
      }
      exit(EXIT_FAILURE);
    }

    munmap(input_index->index_data, input_index->index_data_size);
    free(input_index);

    munmap(input_data, input_data_size);

    fclose(input_index_fh);
    fclose(input_data_fh);

    if (isCa3m) {
      munmap(input_sequence_index->index_data, input_sequence_index->index_data_size);
      free(input_sequence_index);

      munmap(input_sequence_data, input_sequence_data_size);

      fclose(input_sequence_data_fh);
      fclose(input_sequence_index_fh);


      munmap(input_header_index->index_data, input_header_index->index_data_size);
      free(input_header_index);

      munmap(input_header_data, input_header_data_size);

      fclose(input_header_data_fh);
      fclose(input_header_index_fh);
    }

    return EXIT_SUCCESS;
  }
}


int main(int argc, char *argv[]) {
  return cs::CSTranslateMpiApp<cs::AA>().main(argc, argv, stdout, "cstranslate_mpi");
}
