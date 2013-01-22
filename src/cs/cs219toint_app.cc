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
#include "application.h"
#include "sequence-inl.h"
#include "as.h"

using namespace GetOpt;
using std::string;

namespace cs {

struct CS219toIntAppOptions {
    CS219toIntAppOptions() {}
    virtual ~CS219toIntAppOptions() {}

    // The input sequence file.
    string infile;
    // The output file.
    string outfile;
}; 


template<class Abc>
class CS219toIntApp : public Application {
  private:
    // Runs the application.
    virtual int Run();
    // Parses command line options.
    virtual void ParseOptions(GetOpt_pp& ops);
    // Prints options summary to stream.
    virtual void PrintOptions() const;
    // Prints short application description.
    virtual void PrintBanner() const;
    // Prints usage banner to stream.
    virtual void PrintUsage() const;

    // Parameter wrapper
    CS219toIntAppOptions opts_;
};  // class CS219toIntApp

template<class Abc>
void CS219toIntApp<Abc>::ParseOptions(GetOpt_pp& ops) {
    ops >> Option('i', "infile", opts_.infile, opts_.infile);
    ops >> Option('o', "outfile", opts_.outfile, opts_.outfile);
}

template<class Abc>
void CS219toIntApp<Abc>::PrintBanner() const {
    fputs("Translate a cs219 sequence into an integer sequence.\n", out_);
}

template<class Abc>
void CS219toIntApp<Abc>::PrintUsage() const {
    fputs("Usage: cs219 -i <infile> [options]\n", out_);
}

template<class Abc>
void CS219toIntApp<Abc>::PrintOptions() const {
    fprintf(out_, "  %-30s %s\n", "-i, --infile <file>",
            "Input file with sequence");
    fprintf(out_, "  %-30s %s\n", "-o, --outfile <file>", "Output file for integer sequence (def: stdout)");
}

template<class Abc>
int CS219toIntApp<Abc>::Run() {
    FILE* fin;
    if (opts_.infile.empty() || opts_.infile.compare("stdin") == 0 || opts_.infile.compare("-") == 0) {
        fin = stdin;
    } else {
        fin = fopen(opts_.infile.c_str(), "r");
        if (!fin)
            throw Exception("Unable to read input file '%s'!", opts_.infile.c_str());
    }
    Sequence<Abc> seq(fin);
    fclose(fin);

    FILE* fout;
    if (opts_.outfile.empty() || opts_.outfile.compare("stdout") == 0) {
        fout = stdout;
    } else {
        fout = fopen(opts_.outfile.c_str(), "w");
        if (!fout)
            throw Exception("Unable to write to file '%s'!", opts_.outfile.c_str());
    }
    // fprintf(fout, ">%s\n", seq.header().c_str());
    for (size_t i = 0; i < seq.length(); ++i) {
        fprintf(fout, "%c", (char)seq[i]);
    }
    fclose(fout);
    return 0;
}

}  // namespace cs

int main(int argc, char* argv[]) {
    return cs::CS219toIntApp<cs::AS219>().main(argc, argv, stdout, "cs219toint");
}
