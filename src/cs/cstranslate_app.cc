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

#include "cs.h"
#include "alignment-inl.h"
#include "application.h"
#include "blosum_matrix.h"
#include "context_library.h"
#include "crf_pseudocounts-inl.h"
#include "crf-inl.h"
#include "count_profile-inl.h"
#include "getopt_pp.h"
#include "library_pseudocounts-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "pssm.h"
#include "sequence-inl.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <sstream>

extern "C" {
  #include <ffindex.h>     // fast index-based database reading
  #include <ffutil.h>
}


using namespace GetOpt;
using std::string;
using std::vector;

namespace cs {

struct CSTranslateAppOptions {
    static const int kAssignMatchColsByQuery = -1;

    CSTranslateAppOptions() { Init(); }
    virtual ~CSTranslateAppOptions() {}

    // Set csbuild default parameters
    void Init() {
        informat         = "auto";
        outformat        = "seq";
        pc_admix         = 0.90;
        pc_ali           = 12.0;
        pc_engine        = "auto";
        match_assign     = kAssignMatchColsByQuery;
        weight_center    = 1.6;
        weight_decay     = 0.85;
        weight_as        = 1000.0;
        binary           = false;
        ffindex          = false;
        verbose          = true;
    }

    // Validates the parameter settings and throws exception if needed.
    void Validate() {
        if (infile.empty()) throw Exception("No input file provided!");
        if (alphabetfile.empty()) throw Exception("No abstract states provided!");
    }

    // The input alignment file with training data.
    string infile;
    // The output file.
    string outfile;
    // The file to which the output should be appended.
    string appendfile;
    // Input file with context profile library or HMM for generating pseudocounts
    string modelfile;
    // Input file with profile library to be used as abstract state alphabet
    string alphabetfile;
    // Input file format
    string informat;
    // Output file format (abstract state sequence or abstract state profile)
    string outformat;
    // Overall pseudocount admixture
    double pc_admix;
    // Constant in pseudocount calculation for alignments
    double pc_ali;
    // Pseudocount engine
    string pc_engine;
    // Match column assignment for FASTA alignments
    int match_assign;
    // Weight of central column in multinomial emission
    double weight_center;
    // Exponential decay of window weights
    double weight_decay;
    // Weight in emission calculation of abstract states
    double weight_as;
    // Output binary sequence
    bool binary;
    // verbose output
    bool verbose;
    // ffindex
    bool ffindex;
};  // CSTranslateAppOptions


template<class Abc>
class CSTranslateApp : public Application {
  private:
    // Runs the csbuild application.
    virtual int Run();
    // Parses command line options.
    virtual void ParseOptions(GetOpt_pp& ops);
    // Prints options summary to stream.
    virtual void PrintOptions() const;
    // Prints short application description.
    virtual void PrintBanner() const;
    // Prints usage banner to stream.
    virtual void PrintUsage() const;
    // Writes abstract state sequence to outfile
    void WriteStateSequence(const Sequence<AS219>& seq, string outfile, bool append = false) const;
    // Writes abstract state profile to outfile
    void WriteStateProfile(const CountProfile<AS219>& prof, string outfile, bool append = false) const;

    void WriteStateSequence(const Sequence<AS219>& seq, std::stringstream& ss);
    void WriteStateProfile(const CountProfile<AS219>& prof, std::stringstream& ss);

    void ReadProfile(FILE* fin, string& header, CountProfile<Abc>& profile);
    void Translate(CountProfile<Abc>& profile, const Emission<Abc>& emission, CountProfile<cs::AS219>& as_profile);
    void BuildSequence(CountProfile<AS219>& as_profile, size_t profile_counts_length, Sequence<cs::AS219>& as_seq);

    // Parameter wrapper
    CSTranslateAppOptions opts_;
    // Profile library with abstract states
    scoped_ptr< ContextLibrary<Abc> > as_lib_;
    // Profile library for context pseudocounts
    scoped_ptr< ContextLibrary<Abc> > pc_lib_;
    // CRF for CRF context pseudocounts
    scoped_ptr< Crf<Abc> > pc_crf_;
    // Pseudocount engine
    scoped_ptr< Pseudocounts<Abc> > pc_;
};  // class CSTranslateApp



template<class Abc>
void CSTranslateApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('i', "infile", opts_.infile, opts_.infile);
    ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
    ops >> Option('a', "appendfile", opts_.appendfile, opts_.appendfile);
    ops >> Option('I', "informat", opts_.informat, opts_.informat);
    ops >> Option('O', "outformat", opts_.outformat, opts_.outformat);
    ops >> Option('M', "match-assign", opts_.match_assign, opts_.match_assign);
    ops >> Option('x', "pc-admix", opts_.pc_admix, opts_.pc_admix);
    ops >> Option('c', "pc-ali", opts_.pc_ali, opts_.pc_ali);
    ops >> Option('A', "alphabet", opts_.alphabetfile, opts_.alphabetfile);
    ops >> Option('D', "context-data", opts_.modelfile, opts_.modelfile);
    ops >> Option('p', "pc-engine", opts_.pc_engine, opts_.pc_engine);
    ops >> Option('w', "weight", opts_.weight_as, opts_.weight_as);
    ops >> OptionPresent('b', "binary", opts_.binary);
    ops >> OptionPresent('f', "ffindex", opts_.ffindex);
    ops >> Option('v', "verbose", opts_.verbose, opts_.verbose);

    opts_.Validate();

    if (opts_.outfile.compare("stdout") == 0)
        opts_.verbose = false;

    if (opts_.outfile.empty() && opts_.appendfile.empty())
        opts_.outfile = GetBasename(opts_.infile, false) + ".as";
    if (opts_.informat == "auto")
        opts_.informat = GetFileExt(opts_.infile);
    if (opts_.pc_engine == "auto" && !opts_.modelfile.empty())
        opts_.pc_engine = GetFileExt(opts_.modelfile);
}

template<class Abc>
void CSTranslateApp<Abc>::PrintBanner() const {
    fputs("Translate a sequence/alignment into an abstract state alphabet.\n", out_);
}

template<class Abc>
void CSTranslateApp<Abc>::PrintUsage() const {
    fputs("Usage: cstranslate -i <infile> -A <alphabetlib> [options]\n", out_);
}

template<class Abc>
void CSTranslateApp<Abc>::PrintOptions() const {
    fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
            "Input file with alignment or sequence");
    fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "Output file for generated abstract state sequence (def: <infile>.as)");
    fprintf(out_, "  %-30s %s\n", "-a, --append <file>", "Append generated abstract state sequence to this file");
    fprintf(out_, "  %-30s %s (def=%s)\n", "-I, --informat prf|seq|fas|...", "Input format: prf, seq, fas, a2m, or a3m", opts_.informat.c_str());
    fprintf(out_, "  %-30s %s (def=%s)\n", "-O, --outformat seq|prf", "Outformat: abstract state sequence or profile", opts_.outformat.c_str());
    fprintf(out_, "  %-30s %s\n", "-M, --match-assign [0:100]", "Make all FASTA columns with less than X% gaps match columns");
    fprintf(out_, "  %-30s %s\n", "", "(def: make columns with residue in first sequence match columns)");
    fprintf(out_, "  %-30s %s (def=off)\n", "-A, --alphabet <file>", "Abstract state alphabet consisting of exactly 219 states");
    fprintf(out_, "  %-30s %s (def=off)\n", "-D, --context-data <file>", "Add context-specific pseudocounts using given context-data");
    // fprintf(out_, "  %-30s %s (def=%s)\n", "-p, --pc-engine lib|crf", "Specify engine for pseudocount generation", opts_.pc_engine.c_str());
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-x, --pc-admix [0,1]", "Pseudocount admix for context-specific pseudocounts", opts_.pc_admix);
    fprintf(out_, "  %-30s %s (def=%-.1f)\n", "-c, --pc-ali [0,inf[", "Constant in pseudocount calculation for alignments", opts_.pc_ali);
    fprintf(out_, "  %-30s %s (def=%-.2f)\n", "-w, --weight [0,inf[", "Weight of abstract state column in emission calculation", opts_.weight_as);
    fprintf(out_, "  %-30s %s (def=off)\n", "-b, --binary", "Write binary instead of character sequence");
    fprintf(out_, "  %-30s %s (def=off)\n", "-f, --ffindex", "Read from -i <ffindex>, write to -o <ffindex>; enables openmp if possible");
}

template<class Abc>
void CSTranslateApp<Abc>::WriteStateSequence(const Sequence<AS219>& seq, string outfile, bool append) const {
    FILE* fout;
    if (outfile.compare("stdout") == 0)
        fout = stdout;
    else
        fout = fopen(outfile.c_str(), append ? "a" : "w");
    if (!fout) throw Exception("Can't %s to file '%s'!", append ? "append" : "write", outfile.c_str());
    if (opts_.binary) {
        for (size_t i = 0; i < seq.length(); ++i) {
            fputc((char)seq[i], fout);
        }
    } else {
        seq.Write(fout);
    }
    fclose(fout);
    if (opts_.verbose)
        fprintf(out_, "%s abstract state sequence to %s\n", append ? "Appended" : "Wrote", outfile.c_str());
}

template<class Abc>
void CSTranslateApp<Abc>::WriteStateSequence(const Sequence<AS219>& seq, std::stringstream& ss) {
  if (opts_.binary) {
      for (size_t i = 0; i < seq.length(); ++i) {
        ss.put((char)seq[i]);
      }
  } else {
      seq.Write(ss);
  }
}


template<class Abc>
void CSTranslateApp<Abc>::WriteStateProfile(const CountProfile<AS219>& prof, string outfile, bool append) const {
    FILE* fout;
    if (outfile.compare("stdout") == 0)
        fout = stdout;
    else
        fout = fopen(outfile.c_str(), append ? "a" : "w");
    if (!fout) throw Exception("Can't %s to output file '%s'!", append ? "append" : "write", outfile.c_str());
    prof.Write(fout);
    fclose(fout);
    if (opts_.verbose)
        fprintf(out_, "%s abstract state count profile to %s\n", append ? "Appended" : "Wrote", outfile.c_str());
}


template<class Abc>
void CSTranslateApp<Abc>::WriteStateProfile(const CountProfile<AS219>& prof, std::stringstream& ss) {
    prof.Write(ss);
}



char GetMatchSymbol(double pp) {
    char rv = '=';
    if (pp > 0.8) rv = '|';
    else if (pp > 0.6) rv = '+';
    else if (pp > 0.4) rv = '.';
    else if (pp > 0.2) rv = '-';
    return rv;
}

inline int GetConfidence(double pp) {
    return static_cast<int>(floor((pp - DBL_EPSILON) * 10));
}

template<class Abc>
void CSTranslateApp<Abc>::ReadProfile(FILE* fin, string& header, CountProfile<Abc>& profile) {
    if (opts_.informat == "prf") {  // read count profile from infile
        profile = CountProfile<Abc>(fin);;
        if (profile.name.empty()) header = GetBasename(opts_.infile, false);
        else header = profile.name;

        if (pc_) {
            if (opts_.verbose)
                fprintf(out_, "Adding cs-pseudocounts (admix=%.2f) ...\n", opts_.pc_admix);
            CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
            profile.counts = pc_->AddTo(profile, admix);
            Normalize(profile.counts, profile.neff);
        }

    } else if (opts_.informat == "seq") {  // build profile from sequence
        Sequence<Abc> seq(fin);
        header = seq.header();
        profile = CountProfile<Abc>(seq);

        if (pc_) {
            if (opts_.verbose)
                fprintf(out_, "Adding cs-pseudocounts (admix=%.2f) ...\n", opts_.pc_admix);
            ConstantAdmix admix(opts_.pc_admix);
            profile.counts = pc_->AddTo(seq, admix);
        }

    } else {  // build profile from alignment
        AlignmentFormat f = AlignmentFormatFromString(opts_.informat);
        Alignment<Abc> ali(fin, f);
        header = ali.name();

        if (f == FASTA_ALIGNMENT) {
            if (opts_.match_assign == CSTranslateAppOptions::kAssignMatchColsByQuery)
                ali.AssignMatchColumnsBySequence(0);
            else
                ali.AssignMatchColumnsByGapRule(opts_.match_assign);
        }
        profile = CountProfile<Abc>(ali);

        if (pc_) {
            if (opts_.verbose)
                fprintf(out_, "Adding cs-pseudocounts (admix=%.2f) ...\n", opts_.pc_admix);
            CSBlastAdmix admix(opts_.pc_admix, opts_.pc_ali);
            profile.counts = pc_->AddTo(profile, admix);
            Normalize(profile.counts, profile.neff);
        }
    }
    fclose(fin);  // close input file
}

template<class Abc>
void CSTranslateApp<Abc>::Translate(CountProfile<Abc>& profile, const Emission<Abc>& emission, CountProfile<cs::AS219>& as_profile) {
    for (size_t i = 0; i < as_profile.length(); ++i)
        CalculatePosteriorProbs(*as_lib_, emission, profile, i, as_profile.counts[i]);
    as_profile.name = GetBasename(opts_.infile, false);
    as_profile.name = as_profile.name.substr(0, as_profile.name.length() - 1);
}

template<class Abc>
void CSTranslateApp<Abc>::BuildSequence(CountProfile<AS219>& as_profile, size_t profile_counts_length, Sequence<cs::AS219>& as_seq) {
    // We also build an abstract state sequence, either just for pretty printing or
    // even because this is the actual output that the user wants
    for (size_t i = 0; i < profile_counts_length; ++i) {
        // Find state with maximal posterior prob and assign it to as_seq[i]
        size_t k_max = 0;
        double p_max = as_profile.counts[i][0];
        for (size_t k = 1; k < AS219::kSize; ++k) {
            if (as_profile.counts[i][k] > p_max) {
                k_max = k;
                p_max = as_profile.counts[i][k];
            }
        }
        as_seq[i] = k_max;
    }
}



template<class Abc>
int CSTranslateApp<Abc>::Run() {
    // Create emission functor needed for translation into abstract states
    Emission<Abc> emission(1, opts_.weight_as, 1.0);

    // Setup pseudocount engine
    if (!opts_.modelfile.empty() && opts_.pc_engine == "lib") {
        if (opts_.verbose)
            fprintf(out_, "Reading context library for pseudocounts from %s ...\n",
                GetBasename(opts_.modelfile).c_str());
        FILE* fin = fopen(opts_.modelfile.c_str(), "r");
        if (!fin)
            throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
        pc_lib_.reset(new ContextLibrary<Abc>(fin));
        TransformToLog(*pc_lib_);
        fclose(fin);
        pc_.reset(new LibraryPseudocounts<Abc>(*pc_lib_, opts_.weight_center,
                                               opts_.weight_decay));

    } else if (!opts_.modelfile.empty() && opts_.pc_engine == "crf") {
        fprintf(out_, "Reading CRF for pseudocounts from %s ...\n",
                GetBasename(opts_.modelfile).c_str());
        FILE* fin = fopen(opts_.modelfile.c_str(), "r");
        if (!fin)
            throw Exception("Unable to read file '%s'!", opts_.modelfile.c_str());
        pc_crf_.reset(new Crf<Abc>(fin));
        fclose(fin);
        pc_.reset(new CrfPseudocounts<Abc>(*pc_crf_));
    }

    // Setup abstract state engine
    if (opts_.verbose)
        fprintf(out_, "Reading abstract state alphabet from %s ...\n",
            GetBasename(opts_.alphabetfile).c_str());
    FILE* fin = fopen(opts_.alphabetfile.c_str(), "r");
    if (!fin) throw Exception("Unable to read file '%s'!",
                              opts_.alphabetfile.c_str());
    as_lib_.reset(new ContextLibrary<Abc>(fin));
    if (as_lib_->size() != AS219::kSize)
        throw Exception("Abstract state alphabet should have %zu states but actually "
                        "has %zu states!", AS219::kSize, as_lib_->size());
    if (static_cast<int>(as_lib_->wlen()) != 1)
        throw Exception("Abstract state alphabet should have a window length of %zu "
                        "but actually has %zu!", 1, as_lib_->wlen());
    TransformToLog(*as_lib_);
    fclose(fin);


    if(!opts_.ffindex) {
      if (strcmp(opts_.infile.c_str(), "stdin") == 0)
        fin = stdin;
      else
        fin = fopen(opts_.infile.c_str(), "r");
      if (!fin)
        throw Exception("Unable to read input file '%s'!", opts_.infile.c_str());

      string header;
      CountProfile<Abc> profile;  // input profile we want to translate
      ReadProfile(fin, header, profile);

      size_t profile_counts_lenght = profile.counts.length();

      // Prepare abstract sequence in AS219 format
      Sequence<AS219> as_seq(profile_counts_lenght);
      as_seq.set_header(header);

      // Translate count profile into abstract state count profile (Neff is one)
      if (opts_.verbose)
        fputs("Translating count profile to abstract state alphabet AS219 ...\n", out_);

      CountProfile<AS219> as_profile(profile_counts_lenght);  // output profile
      Translate(profile, emission, as_profile);

      BuildSequence(as_profile, profile_counts_lenght, as_seq);

      // Build pseudo-alignment for output
      const size_t nseqs = 5;
      const size_t header_width = 4;
      const size_t width = 100;
      vector<string> ali(nseqs, "");
      vector<string> labels(nseqs, "");
      labels[0] = "Pos";
      labels[1] = "Cons";
      labels[3] = "AS219";
      labels[4] = "Conf";

      BlosumMatrix sm;
      ali[1] = ConservationSequence(profile, sm);

      for (size_t i = 0; i < profile.counts.length(); ++i) {
        ali[0].append(strprintf("%d", (i + 1) % 10));
        ali[2].push_back(GetMatchSymbol(as_profile.counts[i][as_seq[i]]));
        ali[3].push_back(as_seq.chr(i));
        ali[4].append(strprintf("%d", GetConfidence(as_profile.counts[i][as_seq[i]])));
      }

      // Print pseudo alignment in blocks if verbose
      if (opts_.verbose){
        fputc('\n', out_);  // blank line before alignment
        while (!ali.front().empty()) {
          for (size_t k = 0; k < nseqs; ++k) {
            string label = labels[k];
            label += string(header_width - label.length() + 1, ' ');
            fputs(label.c_str(), out_);
            fputc(' ', out_);  // separator between header and sequence

            size_t len = MIN(width, ali[k].length());
            fputs(ali[k].substr(0, len).c_str(), out_);
            fputc('\n', out_);
            ali[k].erase(0, len);
          }
          fputc('\n', out_);  // blank line after each block
        }
      }

	    // Write abstract-state sequence or profile to outfile
	    if (opts_.outformat == "seq") {
	        if (!opts_.outfile.empty()) {
	            WriteStateSequence(as_seq, opts_.outfile, false);
	        }
	        if (!opts_.appendfile.empty()) {
	            WriteStateSequence(as_seq, opts_.appendfile, true);
	        }
	    } else {
	        if (!opts_.outfile.empty()) {
	            WriteStateProfile(as_profile, opts_.outfile, false);
	        }
	        if (!opts_.appendfile.empty()) {
	            WriteStateProfile(as_profile, opts_.appendfile, true);
	        }
	    }
    }
    else {
      std::string input_data_file = opts_.infile + ".ffdata";
      std::string input_index_file = opts_.infile + ".ffindex";

      FILE *input_data_fh  = fopen(input_data_file.c_str(), "r");
      FILE *input_index_fh = fopen(input_index_file.c_str(), "r");

      if (input_data_fh == NULL) {
        LOG(ERROR) << "Could not open ffindex input data file! (" << input_data_file << ")!" << std::endl;
        exit(1);
      }

      if(input_index_fh == NULL) {
        LOG(ERROR) << "Could not open ffindex input index file! (" << input_index_file << ")!" << std::endl;
        exit(1);
      }

      size_t input_offset;
      char* input_data = ffindex_mmap_data(input_data_fh, &input_offset);
      size_t entries = ffcount_lines(input_index_file.c_str());
      ffindex_index_t* input_index = ffindex_index_parse(input_index_fh, entries);

      if(input_index == NULL) {
        LOG(ERROR) << "Input index could not be loaded!" << std::endl;
        exit(1);
      }

      //prepare output ffindex cs219 database
      std::string output_data_file = opts_.outfile+".ffdata";
      std::string output_index_file = opts_.outfile+".ffindex";

      FILE *output_data_fh  = fopen(output_data_file.c_str(), "w");
      FILE *output_index_fh = fopen(output_index_file.c_str(), "w");

      if (output_data_fh == NULL) {
        LOG(ERROR) << "Could not open ffindex output data file! (" << output_data_file << ")!" << std::endl;
        exit(1);
      }

      if(output_index_fh == NULL) {
        LOG(ERROR) << "Could not open ffindex output index file! (" << output_index_file << ")!" << std::endl;
        exit(1);
      }

      size_t output_offset = 0;


      size_t input_range_start = 0;
      size_t input_range_end = input_index->n_entries;

      // Foreach entry
      #pragma omp parallel for shared(input_index, input_data, output_data_fh, output_index_fh, output_offset)
      for(size_t entry_index = input_range_start; entry_index < input_range_end; entry_index++) {
        ffindex_entry_t* entry = ffindex_get_entry_by_index(input_index, entry_index);

        if (entry == NULL) {
          LOG(WARNING) << "Could not open entry " << entry_index << " from input ffindex!" << std::endl;
          continue;
        }

        FILE* inf = ffindex_fopen_by_entry(input_data, entry);
        if(inf == NULL) {
          LOG(WARNING) << "Could not open input entry (" << entry->name << ")!" << std::endl;
          continue;
        }

        string header;
        CountProfile<Abc> profile;  // input profile we want to translate
        ReadProfile(inf, header, profile);

        size_t profile_counts_length = profile.counts.length();

        CountProfile<AS219> as_profile(profile_counts_length);  // output profile
        Translate(profile, emission, as_profile);

        // Prepare abstract sequence in AS219 format
        Sequence<AS219> as_seq(profile_counts_length);
        as_seq.set_header(header);
        BuildSequence(as_profile, profile_counts_length, as_seq);

        std::stringstream out_buffer;

        if (opts_.outformat == "seq") {
          WriteStateSequence(as_seq, out_buffer);
        } else {
          WriteStateProfile(as_profile, out_buffer);
        }

        std::string out_string = out_buffer.str();
        #pragma omp critical
        {
          ffindex_insert_memory(output_data_fh, output_index_fh, &output_offset, const_cast<char*>(out_string.c_str()), out_string.size(), entry->name);
        }

        fclose(inf);
      }
    }

    return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
    return cs::CSTranslateApp<cs::AA>().main(argc, argv, stdout, "cstranslate");
}
