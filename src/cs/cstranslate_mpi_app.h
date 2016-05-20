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

#include "cstranslate_app.h"

#include <unistd.h>
#include <mpi.h>

extern "C" {
#include <mpq/mpq.h>
}

namespace cs {
  // MPQ passes a temporary environment object into the worker, we can avoid this by just using class itself
  class MPQWrapper {
  public:
    void Worker () {
      int message[3];
      while (1) {
        MPI_Recv(message, 3, MPI_INT, MPQ_MASTER, MPQ_TAG_JOB, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (message[0] == MPQ_MSG_RELEASE) {
          break;
        }

        Payload(message[1], message[2]);

        message[0] = MPQ_MSG_DONE;
        message[1] = 0;
        message[2] = 0;

        MPI_Send(message, 3, MPI_INT, MPQ_MASTER, MPQ_TAG_DONE, MPI_COMM_WORLD);
      }
    }

  private:
    virtual void Payload(const size_t start, const size_t end) = 0;

    // copied from MPQ, keep in sync
    enum {
      MPQ_TAG_JOB,
      MPQ_TAG_DONE
    };

    enum {
      MPQ_MSG_RELEASE,
      MPQ_MSG_JOB,
      MPQ_MSG_DONE
    };
  };

  template<class Abc>
  class CSTranslateMpiApp : public CSTranslateApp<Abc>, private MPQWrapper {
  public:
    virtual int Run() {
      std::string input_data_file = this->opts_.infile + ".ffdata";
      std::string input_index_file = this->opts_.infile + ".ffindex";

      const bool isCa3m = this->opts_.informat == "ca3m";
      if (isCa3m) {
        // infile has to be the ffindex basepath with no suffices
        input_data_file = this->opts_.infile + "_ca3m.ffdata";
        input_index_file = this->opts_.infile + "_ca3m.ffindex";
      }

      FFindexDatabase input(const_cast<char *>(input_data_file.c_str()),
                            const_cast<char *>(input_index_file.c_str()), isCa3m);

      //prepare output ffindex cs219 database
      std::string data_filename_out[2];
      std::string index_filename_out[2];
      if (this->opts_.both) {
        data_filename_out[0] = this->opts_.outfile + "_binary.ffdata";
        index_filename_out[0] = this->opts_.outfile + "_binary.ffindex";
        data_filename_out[1] = this->opts_.outfile + "_plain.ffdata";
        index_filename_out[1] = this->opts_.outfile + "_plain.ffindex";
      } else {
        data_filename_out[0] = this->opts_.outfile + ".ffdata";
        index_filename_out[0] = this->opts_.outfile + ".ffindex";
      }

      int mpq_status = MPQ_Init(this->argc_, this->argv_, input.db_index->n_entries);
      if (mpq_status == MPQ_SUCCESS) {
        if (MPQ_rank == MPQ_MASTER) {
          MPQ_Master(1);
        } else {
          this->SetupEmissions();
          this->SetupPseudocountEngine();
          this->SetupAbstractStateEngine();

          this->input_index = input.db_index;
          this->input_data = input.db_data;

          FFindexDatabase *header_db = NULL;
          FFindexDatabase *sequence_db = NULL;

          if (isCa3m) {
            std::string input_header_data_file = this->opts_.infile + "_header.ffdata";
            std::string input_header_index_file = this->opts_.infile + "_header.ffindex";
            header_db = new FFindexDatabase(const_cast<char *>(input_header_data_file.c_str()),
                                            const_cast<char *>(input_header_index_file.c_str()), false);


            std::string input_sequence_data_file = this->opts_.infile + "_sequence.ffdata";
            std::string input_sequence_index_file = this->opts_.infile + "_sequence.ffindex";
            sequence_db = new FFindexDatabase(const_cast<char *>(input_sequence_data_file.c_str()),
                                              const_cast<char *>(input_sequence_index_file.c_str()), false);

          }

          this->input_header_index = header_db ? header_db->db_index : NULL;
          this->input_header_data = header_db ? header_db->db_data : NULL;
          this->input_sequence_index = sequence_db ? sequence_db->db_index : NULL;
          this->input_sequence_data = sequence_db ? sequence_db->db_data : NULL;

          this->data_file_out[0] = openWrite(data_filename_out[0].c_str());
          this->index_file_out[0] = openWrite(index_filename_out[0].c_str());
          this->offset[0] = 0;

          if (this->opts_.both) {
            this->data_file_out[1] = openWrite(data_filename_out[1].c_str());
            this->index_file_out[1] = openWrite(index_filename_out[1].c_str());
            this->offset[1] = 0;
          }

          std::string log_filename_out = this->opts_.outfile + ".log";
          this->log_file = openWrite(log_filename_out.c_str());

          Worker();

          if (this->log_file) {
            int fd = fileno(this->log_file);
            fflush(this->log_file);
            fsync(fd);
            fclose(this->log_file);
          }

          for (size_t i = 0; i < (this->opts_.both ? 2 : 1); i++) {
            if (this->index_file_out[i]) {
              int fd = fileno(this->index_file_out[i]);
              fflush(this->index_file_out[i]);
              fsync(fd);
              fclose(this->index_file_out[i]);
            }

            if (this->data_file_out[i]) {
              int fd = fileno(this->data_file_out[i]);
              fflush(this->data_file_out[i]);
              fsync(fd);
              fclose(this->data_file_out[i]);
            }
          }

          if (isCa3m) {
            delete sequence_db;
            delete header_db;
          }
        }
        MPQ_Finalize();

        if (MPQ_rank == MPQ_MASTER) {
          ffmerge_splits(data_filename_out[0].c_str(), index_filename_out[0].c_str(), MPQ_size, 1);
          if (this->opts_.both) {
            ffmerge_splits(data_filename_out[1].c_str(), index_filename_out[1].c_str(), MPQ_size, 1);
          }
        }
      } else {
        if (mpq_status == MPQ_ERROR_NO_WORKERS) {
          fprintf(stderr, "MPQ_Init: Needs at least one worker process.\n");
        } else if (mpq_status == MPQ_ERROR_TOO_MANY_WORKERS) {
          fprintf(stderr, "MPQ_Init: Too many worker processes.\n");
        }
        exit(EXIT_FAILURE);
      }

      return EXIT_SUCCESS;
    };

    void Payload(const size_t start, const size_t end) {
      for (size_t entry_index = start; entry_index < end; entry_index++) {
        ffindex_entry_t *entry = ffindex_get_entry_by_index(this->input_index, entry_index);

        if (entry == NULL) {
          fprintf(this->log_file, "Could not open entry %zu from input ffindex!\n", entry_index);
          continue;
        }

        if (this->opts_.verbose) {
          fprintf(this->log_file, "Processing entry: %s\n", entry->name);
        }

        std::ostringstream a3m_buffer;
        std::string a3m_string;
        FILE *inf;
        if (this->opts_.informat == "ca3m") {
          char *entry_data = ffindex_get_data_by_entry(this->input_data, entry);

          compressed_a3m::extract_a3m(entry_data, entry->length,
                                      this->input_sequence_index, this->input_sequence_data,
                                      this->input_header_index, this->input_header_data,
                                      &a3m_buffer);

          a3m_string = a3m_buffer.str();

          inf = fmemopen(static_cast<void *>(const_cast<char *>(a3m_string.c_str())), a3m_string.length(), "r");
        } else {
          inf = ffindex_fopen_by_entry(this->input_data, entry);
        }

        if (inf == NULL) {
          fprintf(this->log_file, "Could not open input entry (%s)!\n", entry->name);
          continue;
        }

        string header;
        CountProfile<Abc> profile;  // input profile we want to translate
        try {
          this->ReadProfile(inf, header, profile);
        } catch (const Exception &e) {
          fprintf(this->log_file, "Could not read entry: %s, Message: %s\n", entry->name, e.what());
          continue;
        }

        size_t profile_counts_length = profile.counts.length();

        CountProfile<AS219> as_profile(profile_counts_length);  // output profile
        this->Translate(profile, as_profile);

        // Prepare abstract sequence in AS219 format
        Sequence<AS219> as_seq(profile_counts_length);
        as_seq.set_header(header);
        this->BuildSequence(as_profile, profile_counts_length, as_seq);

        std::stringstream out_buffer[2];
        if (this->opts_.outformat == "seq") {
          if(this->opts_.both) {
            this->WriteStateSequence(as_seq, out_buffer[0], true);
            this->WriteStateSequence(as_seq, out_buffer[1], false);
          } else {
            this->WriteStateSequence(as_seq, out_buffer[0], this->opts_.binary);
          }
        } else {
          this->WriteStateProfile(as_profile, out_buffer[0]);
        }
        std::string out_string = out_buffer[0].str();
        ffindex_insert_memory(this->data_file_out[0], this->index_file_out[0],
                              &(this->offset[0]), const_cast<char *>(out_string.c_str()),
                              out_string.size(), entry->name);

        if(this->opts_.both) {
          out_string = out_buffer[1].str();
          ffindex_insert_memory(this->data_file_out[1], this->index_file_out[1],
                                &(this->offset[1]), const_cast<char *>(out_string.c_str()),
                                out_string.size(), entry->name);
        }
        // FIXME: we are leaking inf, but if we fclose we get weird crashes
        //fclose(inf);

        if (entry_index % 1000 == 0) {
          fflush(this->data_file_out[0]);
          fflush(this->index_file_out[0]);

          if(this->opts_.both) {
            fflush(this->data_file_out[1]);
            fflush(this->index_file_out[1]);
          }
        }

      }
    };

    // Parses command line options.
    virtual void ParseOptions(GetOpt_pp &ops) {
      ops >> Option('i', "infile", this->opts_.infile, this->opts_.infile);
      ops >> Option('o', "outfile", this->opts_.outfile, this->opts_.outfile);
      ops >> Option('I', "informat", this->opts_.informat, this->opts_.informat);
      ops >> Option('O', "outformat", this->opts_.outformat, this->opts_.outformat);
      ops >> Option('M', "match-assign", this->opts_.match_assign, this->opts_.match_assign);
      ops >> Option('x', "pc-admix", this->opts_.pc_admix, this->opts_.pc_admix);
      ops >> Option('c', "pc-ali", this->opts_.pc_ali, this->opts_.pc_ali);
      ops >> Option('A', "alphabet", this->opts_.alphabetfile, this->opts_.alphabetfile);
      ops >> Option('D', "context-data", this->opts_.modelfile, this->opts_.modelfile);
      ops >> Option('p', "pc-engine", this->opts_.pc_engine, this->opts_.pc_engine);
      ops >> Option('w', "weight", this->opts_.weight_as, this->opts_.weight_as);
      ops >> OptionPresent('b', "binary", this->opts_.binary);
      ops >> OptionPresent('f', "ffindex", this->opts_.ffindex);
      ops >> OptionPresent('2', "both", this->opts_.both);
      ops >> Option('v', "verbose", this->opts_.verbose, this->opts_.verbose);

      this->opts_.Validate();

      if (this->opts_.informat == "auto")
        this->opts_.informat = GetFileExt(this->opts_.infile);
      if (this->opts_.pc_engine == "auto" && !this->opts_.modelfile.empty())
        this->opts_.pc_engine = GetFileExt(this->opts_.modelfile);
    };

    // Prints options summary to stream.
    virtual void PrintOptions() const {
      fprintf(this->out_, "  %-30s %s\n", "-i, --infile <ffindex>",
              "Input ffindex with alignments or sequences");
      fprintf(this->out_, "  %-30s %s\n", "-o, --outfile <ffindex>",
              "Output ffindex for generated abstract state sequence");
      fprintf(this->out_, "  %-30s %s (def=%s)\n", "-I, --informat prf|seq|fas|...",
              "Input format: prf, seq, fas, a2m, a3m or ca3m", this->opts_.informat.c_str());
      fprintf(this->out_, "  %-30s %s (def=%s)\n", "-O, --outformat seq|prf", "Outformat: abstract state sequence or profile",
              this->opts_.outformat.c_str());
      fprintf(this->out_, "  %-30s %s\n", "-M, --match-assign [0:100]",
              "Make all FASTA columns with less than X% gaps match columns");
      fprintf(this->out_, "  %-30s %s\n", "", "(def: make columns with residue in first sequence match columns)");
      fprintf(this->out_, "  %-30s %s (def=off)\n", "-A, --alphabet <file>",
              "Abstract state alphabet consisting of exactly 219 states");
      fprintf(this->out_, "  %-30s %s (def=off)\n", "-D, --context-data <file>",
              "Add context-specific pseudocounts using given context-data");
      // fprintf(this->out_, "  %-30s %s (def=%s)\n", "-p, --pc-engine lib|crf", "Specify engine for pseudocount generation", this->opts_.pc_engine.c_str());
      fprintf(this->out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]",
              "Pseudocount admix for context-specific pseudocounts", this->opts_.pc_admix);
      fprintf(this->out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[",
              "Constant in pseudocount calculation for alignments", this->opts_.pc_ali);
      fprintf(this->out_, "  %-30s %s (def=%-.2f)\n", "-w, --weight [0,inf[",
              "Weight of abstract state column in emission calculation", this->opts_.weight_as);
      fprintf(this->out_, "  %-30s %s (def=off)\n", "-b, --binary", "Write binary instead of character sequence");
      fprintf(this->out_, "  %-30s %s (def=off)\n", "-2, --both", "Write both binary and plain character sequence");
    };

    // Prints usage banner to stream.
    virtual void PrintUsage() const {
      fputs("Usage: cstranslate_mpi -i <ffindex_in> -o <ffindex_out> -A <alphabetlib> [options]\n", this->out_);
    };

  private:
    FILE* openWrite(const char* path) {
      char out_rank[FILENAME_MAX];
      snprintf(out_rank, FILENAME_MAX, "%s.%d", path, MPQ_rank);

      FILE* out = fopen(out_rank, "w+");
      if (out == NULL) {
        fprintf(this->out_, "Could not open ffindex output file! (%s)!\n", out_rank);
        exit(1);
      }
      return out;
    };

    ffindex_index_t *input_index;
    char *input_data;

    ffindex_index_t *input_sequence_index;
    char *input_sequence_data;

    ffindex_index_t *input_header_index;
    char *input_header_data;

    FILE *data_file_out[2];
    FILE *index_file_out[2];
    size_t offset[2];

    FILE *log_file;
  };
}
