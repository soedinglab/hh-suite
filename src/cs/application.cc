// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "application.h"
#include "getopt_pp.h"

using std::string;

namespace cs {

Application* Application::instance_;

const char* Application::kVersionNumber = "2.1.2";

const char* Application::kCopyright =
  "Copyright (c) 2010 Andreas Biegert, Johannes Soding, and LMU Munich";

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
  int status = 0;
  out_      = fout;
  app_name_ = name;
  log_file_ = "stderr";

  // Prepare command line parsing
  GetOpt_pp options(argc, argv, Include_Environment);
  options.exceptions_all();

  try {
    // Print usage?
    if (argc < 2 || argv[1][0] == '?' ||
        options >> OptionPresent(' ', std::string("help"))) {
      PrintHelp();
      return 1;
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

  } catch(const std::exception& e) {
    LOG(ERROR) << e.what();
    fprintf(fout, "\n%s\n", e.what());
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
