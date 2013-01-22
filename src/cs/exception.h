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

#ifndef CS_EXCEPTION_H_
#define CS_EXCEPTION_H_

#include <exception>
#include <string>
#include <stdarg.h>

// This is a simple but useful class for construction of standard exceptions
// The constructor is implemented to take variadic arguments
// in the style of printf for simple message construction
class Exception : public std::exception {
 public:
  Exception(const std::string &m) : msg(m), c(0) {}

  Exception(const char *str, ...) : c(0) {
    char *buffer = new char[1024];
    va_list ap;
    va_start(ap, str);
    vsprintf(buffer, str, ap);
    va_end(ap);
    msg = buffer;
    delete [] buffer;
  }

  Exception(const int type, const char *str, ...) : c(type) {
    char *buffer = new char[1024];
    va_list ap;
    va_start(ap, str);
    vsprintf(buffer, str, ap);
    va_end(ap);
    msg = buffer;
    delete [] buffer;
  }

  virtual ~Exception() throw() {};

  virtual const char* what() const throw() { return msg.c_str(); }

  virtual int code() const { return c; }

 private:
  std::string msg;
  int c;
};

#endif  // CS_EXCEPTION_H_
