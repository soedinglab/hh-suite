/*
  Copyright 2009 Andreas Biegert

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

#ifndef CS_LOG_H_
#define CS_LOG_H_

#include <stdio.h>
#include <sys/time.h>

#include <iostream>
#include <sstream>
#include <string>

namespace cs {

enum LogLevel { ERROR, WARNING, INFO, DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4 };

class Log {
 public:
  Log() {}
  virtual ~Log();
  std::ostringstream& get(LogLevel level = INFO);

  static LogLevel& reporting_level();
  static std::string to_string(LogLevel log_level);
  static LogLevel from_string(const std::string& log_level);
  static LogLevel from_int(int log_level);
  static FILE*& stream() {
    static FILE* p_stream = stderr;
    return p_stream;
  }

 protected:
  std::ostringstream os;
  LogLevel level;

 private:
  Log(const Log&);
  Log& operator =(const Log&);
};

#ifndef LOG_MAX_LEVEL
#define LOG_MAX_LEVEL DEBUG4
#endif

#ifdef LOGGING
#define LOG(level)                                              \
  if (level > LOG_MAX_LEVEL) ;                                  \
  else if (level > Log::reporting_level() || !Log::stream()) ; \
  else Log().get(level)
#else
#define LOG(level)                                               \
  if (true) ;                                                    \
  else if (level > Log::reporting_level() || !Log::stream()) ;   \
  else Log().get(level)
#endif

}  // namespace cs

#endif  // CS_LOG_H_
