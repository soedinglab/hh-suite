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
#include "getopt_pp.h"

using std::string;

namespace cs {

Application* Application::instance_;

const char* Application::kVersionNumber = "2.2.26";

const char* Application::kCopyright =
  "Copyright (c) 2010-2012 Andreas Biegert, Christof Angermueller, Johannes Soeding, and LMU Munich";

Application::Application()
    : log_level_(Log::to_string(Log::from_int(LOG_MAX_LEVEL))),
      log_fp_(stderr) {
  // Register the application instance
  if (instance_)
    throw Exception("Second instance of Application is prohibited");
  instance_ = this;
}

Application::~Application() {
  if (log_fp_ && log_fp_ != stderr)
    fclose(log_fp_);
}

int Application::main(int argc, char* argv[], FILE* fout, const string& name) {
  // save original arg parameters for e.g. mpi
  argc_ = argc;
  argv_ = argv;

  int status = 0;
  out_      = fout;
  app_name_ = name;
  log_file_ = "stderr";

  // Prepare command line parsing
  GetOpt_pp options(argc, argv, Include_Environment);
  options.exceptions_all();

  try{
    try {
      // Print usage?
      if (argc < 2 || argv[1][0] == '?' ||
          options >> OptionPresent(' ', std::string("help"))) {
        PrintHelp();
        return 0;
      }

#ifdef LOGGING
      // Process logging options
      options >> Option(' ', "loglevel", log_level_, log_level_);
      Log::reporting_level() = Log::from_string(log_level_);
      options >> Option(' ', "logfile", log_file_, log_file_);
      if (log_file_.empty() || log_file_ == "stderr") log_fp_ = stderr;
      else log_fp_ = fopen(log_file_.c_str(), "w");
      Log::stream() = log_fp_;
#endif

      // Let subclasses parse the command line options
      ParseOptions(options);

      // Run application
      status = Run();

    } catch(const OptionNotFoundEx& e) {
      throw Exception("Missing command line option!");
    } catch(const ArgumentNotFoundEx& e) {
      throw Exception("Missing argument for command line option!");
    } catch(const TooManyOptionsEx& e) {
      throw Exception("Invalid command line option!");
    } catch(const TooManyArgumentsEx& e) {
      throw Exception("Too many arguments for command line option!");
    } catch(const GetOptEx& e) {
      throw Exception("Error parsing command line options!");
    }

  } catch(const std::exception& e) {
    LOG(ERROR) << e.what();
    fprintf(stderr, "\nERROR: %s\n", e.what());
    return 1;
  }

  return status;
}

void Application::PrintHelp() const {
  fprintf(out_, "%s version %s\n", app_name_.c_str(), kVersionNumber);
  PrintBanner();
  fprintf(out_, "%s\n\n", kCopyright);
  PrintUsage();
  fputs("\nOptions:\n", out_);
  PrintOptions();

#ifdef LOGGING
  fprintf(out_, "  %-30s %s (def=%s)\n", "    --loglevel <level>",
          "Maximal reporting level for logging", log_level_.c_str());
  fprintf(out_, "  %-30s %s (def=%s)\n", "    --logfile <file>",
          "Output file for logging", log_file_.c_str());
#endif
}

}  // namespace cs
