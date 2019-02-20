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

#include "log.h"

#include <stdio.h>
#include <sys/time.h>

#include <string>
#include <sstream>

namespace cs {

inline std::string now_time() {
  char buffer[11];
  time_t now = time(0);
  strftime(buffer, sizeof(buffer), "%X", localtime(&now));
  struct timeval tv;
  gettimeofday(&tv, 0);
  char result[100] = {0};
  sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
  return result;
}

std::ostringstream& Log::get(LogLevel log_level) {
  level = log_level;
  std::string level_str(to_string(level));
  os << level_str << std::string(9 - level_str.length(), ' ')
     << now_time() << "   ";
  os << std::string(level > DEBUG ? (level - DEBUG) : 0, '\t');
  return os;
}

Log::~Log() {
  if (os.str().find('\n') == std::string::npos) {
    os << std::endl;
    fprintf(stream(), "%s", os.str().c_str());
  } else {
    const std::string margin(
        "\t\t\t" + std::string(level > DEBUG ? (level - DEBUG) : 0, '\t'));
    std::string s(os.str());
    if (*s.rbegin() == '\n') s.erase(s.begin() + s.length() - 1);

    std::string::size_type i = 0;
    std::string::size_type j = s.find('\n');
    while (j != std::string::npos) {
      fprintf(stream(), "%s\n%s", s.substr(i, j-i).c_str(), margin.c_str());
      i = ++j;
      j = s.find('\n', j);

      if (j == std::string::npos)
        fprintf(stream(), "%s\n", s.substr(i, s.length()).c_str());
    }
  }
  fflush(stream());
}

LogLevel& Log::reporting_level() {
  static LogLevel log_level = DEBUG4;
  return log_level;
}

std::string Log::to_string(LogLevel log_level) {
  static const char* const buffer[] =
    {"ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
  return buffer[log_level];
}

LogLevel Log::from_string(const std::string& log_level) {
  if (log_level == "DEBUG4")
    return DEBUG4;
  if (log_level == "DEBUG3")
    return DEBUG3;
  if (log_level == "DEBUG2")
    return DEBUG2;
  if (log_level == "DEBUG1")
    return DEBUG1;
  if (log_level == "DEBUG")
    return DEBUG;
  if (log_level == "INFO")
    return INFO;
  if (log_level == "WARNING")
    return WARNING;
  if (log_level == "ERROR")
    return ERROR;
  Log().get(WARNING) << "Unknown logging level '" << log_level
                     << "'. Using INFO level as default.";
  return INFO;
}

LogLevel Log::from_int(int log_level) {
  if (log_level == 7)
    return DEBUG4;
  if (log_level == 6)
    return DEBUG3;
  if (log_level == 5)
    return DEBUG2;
  if (log_level == 4)
    return DEBUG1;
  if (log_level == 3)
    return DEBUG;
  if (log_level == 2)
    return INFO;
  if (log_level == 1)
    return WARNING;
  if (log_level == 0)
    return ERROR;
  Log().get(WARNING) << "Unknown logging level '" << log_level
                     << "'. Using INFO level as default.";
  return INFO;
}

}  // namespace cs


