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

#include <sys/mman.h>
#include <csignal>
#include "cstranslate_app.h"

extern "C" {
#include <mpq/mpq.h>
}

namespace cs {
  void ignore_signal(int signal) {
    struct sigaction handler;
    handler.sa_handler = SIG_IGN;
    sigemptyset(&handler.sa_mask);
    handler.sa_flags = 0;
    sigaction(signal, &handler, NULL);
  }

  template<class Abc>
  class CSTranslateMpiApp : public CSTranslateApp<Abc> {
  public:
    virtual int Run();

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
    };

    // Prints usage banner to stream.
    virtual void PrintUsage() const {
      fputs("Usage: cstranslate_mpi -i <ffindex_in> -o <ffindex_out> -A <alphabetlib> [options]\n", this->out_);
    };
  };

  template<class Abc>
  struct CSTranslateMpiEnvironment {
    CSTranslateMpiEnvironment(ffindex_index_t *input_index, char *input_data, ffindex_index_t *input_sequence_index,
                              char *input_sequence_data, ffindex_index_t *input_header_index, char *input_header_data,
                              FILE *data_file_out, FILE *index_file_out, const Emission<Abc> &emission, size_t offset,
                              bool verbose, bool isCa3m, FILE *log_file, CSTranslateMpiApp<Abc> &that)
        : input_index(input_index),
          input_data(input_data),
          input_sequence_index(input_sequence_index),
          input_sequence_data(input_sequence_data),
          input_header_index(input_header_index),
          input_header_data(input_header_data),
          data_file_out(data_file_out),
          index_file_out(index_file_out),
          emission(emission),
          offset(offset),
          verbose(verbose),
          isCa3m(isCa3m),
          log_file(log_file),
          that(that) { }

    ffindex_index_t *input_index;
    char *input_data;

    ffindex_index_t *input_sequence_index;
    char *input_sequence_data;

    ffindex_index_t *input_header_index;
    char *input_header_data;

    FILE *data_file_out;
    FILE *index_file_out;

    const Emission<Abc> &emission;
    size_t offset;
    const bool verbose;
    const bool isCa3m;
    FILE *log_file;

    CSTranslateMpiApp<Abc> &that;
  };

  template<class Abc>
  int CSTranslateMpiPayload(void *environment, const size_t start, const size_t end) {
    CSTranslateMpiEnvironment<Abc> *data = static_cast<CSTranslateMpiEnvironment<Abc> *>(environment);

    int status = EXIT_SUCCESS;
    for (size_t entry_index = start; entry_index < end; entry_index++) {
      ffindex_entry_t *entry = ffindex_get_entry_by_index(data->input_index, entry_index);

      if (entry == NULL) {
        fprintf(data->log_file, "Could not open entry %zu from input ffindex!\n", entry_index);
        continue;
      }

      if (data->verbose) {
        fprintf(data->log_file, "Processing entry: %s\n", entry->name);
      }

      std::ostringstream output;
      std::string tmpOut;
      FILE *inf;
      if (data->isCa3m) {
        char *entry_data = ffindex_get_data_by_entry(data->input_data, entry);

        compressed_a3m::extract_a3m(entry_data, entry->length, data->input_sequence_index, data->input_sequence_data,
                                    data->input_header_index,
                                    data->input_header_data, &output);

        tmpOut = output.str();

        inf = fmemopen(static_cast<void *>(const_cast<char *>(tmpOut.c_str())), tmpOut.length(), "r");
      } else {
        inf = ffindex_fopen_by_entry(data->input_data, entry);
      }

      if (inf == NULL) {
        fprintf(data->log_file, "Could not open input entry (%s)!\n", entry->name);
        continue;
      }

      string header;
      CountProfile<Abc> profile;  // input profile we want to translate
      try {
        data->that.ReadProfile(inf, header, profile);
      } catch (const Exception &e) {
        fprintf(data->log_file, "Could not read entry: %s, Message: %s\n", entry->name, e.what());
        continue;
      }

      size_t profile_counts_length = profile.counts.length();

      CountProfile<AS219> as_profile(profile_counts_length);  // output profile
      data->that.Translate(profile, data->emission, as_profile);

      // Prepare abstract sequence in AS219 format
      Sequence<AS219> as_seq(profile_counts_length);
      as_seq.set_header(header);
      data->that.BuildSequence(as_profile, profile_counts_length, as_seq);

      std::stringstream out_buffer;


      data->that.WriteStateProfile(as_profile, out_buffer);

      std::string out_string = out_buffer.str();

      ffindex_insert_memory(data->data_file_out, data->index_file_out,
                            &(data->offset), const_cast<char *>(out_string.c_str()),
                            out_string.size(), entry->name);

      // FIXME: we are leaking inf, but if we fclose we get weird crashes
      //fclose(inf);
    }
    return status;
  }
}
