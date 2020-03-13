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

#ifndef CS_APPLICATION_H_
#define CS_APPLICATION_H_

#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "globals.h"
#include "getopt_pp.h"

using namespace GetOpt;

namespace cs {

// Basic (abstract) application class.
// Defines the high level behavior of an application. A new application is
// written by deriving a class from Application and writing an implementation
// of the Run() and maybe some other (like ParseOptions() etc.) methods.
class Application {
 public:
  // Register the application instance.
  Application();

  // Clean up the application settings, flush the diagnostic stream.
  virtual ~Application();

  // Main function (entry point) for the application.
  int main(int argc,                 // argc in a regular main
           char* argv[],             // argv in a regular main
           FILE* fout,               // output stream
           const std::string& name   // application name
           );

 protected:
  // Copyright string for usage output.
  static const char* kCopyright;

  // Runs the application and return exit code. To be implemented by derived
  // classes.
  virtual int Run() = 0;
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp& /* options */) {};
  // Prints options summary to stream.
  virtual void PrintOptions() const {};
  // Prints usage banner to stream.
  virtual void PrintUsage() const {};
  // Prints short application description.
  virtual void PrintBanner() const {};
  // Prints copyright notification.
  void PrintHelp() const;

  static Application* instance_;       // current app. instance
  std::string         app_name_;       // application name
  std::string         log_level_;      // log reporting level
  std::string         log_file_;       // name of logfile
  FILE*               log_fp_;         // file pointer to logfile
  FILE*               out_;            // file pointer to output stream

  int argc_;
  char** argv_;
};

}  // namespace cs

#endif  // CS_APPLICATION_H_
