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

#ifndef CS_IO_H_
#define CS_IO_H_

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cs {

// Removes the newline and other control characters at the end of a string (if
//  present) and returns the new length of the string (-1 if str is NULL)
inline int chomp(char* str) {
  if (!str) return -1;
  int l = 0;
  for (l = strlen(str) - 1; l >= 0 && str[l] < 32; --l)
    /* do nothing */;
  str[++l] = '\0';
  return l;
}

// Emulates the ifstream::getline method; similar to fgets(str,maxlen,FILE*),
// but removes the newline at the end and returns NULL if at end of file or read
// error
inline char* fgetline(char* str, int maxlen, FILE* file) {
  if (!fgets(str, maxlen, file)) return NULL;
  if (chomp(str) + 1 >= maxlen)  // if line is cut after maxlen characters...
    while (fgetc(file) != '\n')  // ... read in rest of line
      /* do nothing */;
  return(str);
}

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets ptr to NULL
inline int strtoi(const char*& ptr) {
  int i;
  const char* ptr0 = ptr;
  if (!ptr) return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9')) ++ptr;
  if (*ptr == '\0') {
    ptr = NULL;
    return INT_MIN;
  }
  if (ptr > ptr0 && *(ptr-1) == '-') i = -atoi(ptr); else i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9') ++ptr;
  return i;
}

// Same as strtoi, but interpretes '*' as default
inline int strastoi(const char*& ptr, int deflt = INT_MAX) {
  int i;
  const char* ptr0 = ptr;
  if (!ptr) return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9') && *ptr != '*') ++ptr;
  if (*ptr == '\0') {
    ptr = NULL;
    return INT_MIN;
  }
  if (*ptr == '*') {
    ++ptr;
    return deflt;
  }
  if (ptr > ptr0 &&  *(ptr-1) == '-') i = -atoi(ptr); else i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9') ++ptr;
  return i;
}

// Returns leftmost double in ptr and sets the pointer to first non-whitespace char
// after the double. If no integer is found, returns zero and sets ptr to NULL
inline double strtof(const char*& ptr) {
  double rv = 0.0;
  if (!ptr) return 0.0;
  while (*ptr != '\0' && isspace(*ptr)) ++ptr;
  if (*ptr == '\0') {
    ptr = NULL;
    return 0.0;
  }
  rv = atof(ptr);
  while (!isspace(*ptr)) ++ptr;
  return rv;
}

// Returns pointer to first non-white-space character in str OR to NULL if none
// found
inline const char* strscn(const char* str) {
  if (!str) return NULL;
  const char* ptr = str;
  while (*ptr != '\0' && isspace(*ptr)) ++ptr;
  return (*ptr == '\0') ? NULL : ptr;
}

// Parse serialization record and return integer value following label 'str' in
// line read from file pointer 'fp'.
inline int ReadInt(const char* line,
                   const char* label,
                   const char* errmsg = NULL) {
  int rv = INT_MIN;
  if (strstr(line, label)) {
    const char* ptr = line + strlen(label);
    rv = strtoi(ptr);
  } else if (errmsg) {
    throw Exception(errmsg);
  }
  return rv;
}

// Parse serialization record and return double value following label 'label' in
// line read from file pointer 'fp'.
inline double ReadDouble(const char* line,
                         const char* label,
                         const char* errmsg = NULL) {
  double rv = DBL_MIN;
  if (strstr(line, label)) {
    rv = atof(line + strlen(label));
  } else if (errmsg) {
    throw Exception(errmsg);
  }
  return rv;
}

// Parse serialization record and return string following label 'label' in
// line read from file pointer 'fp'.
inline std::string ReadString(const char* line,
                              const char* label,
                              const char* errmsg = NULL) {
  std::string rv;
  if (strstr(line, label)) {
    const char* ptr = strscn(line + strlen(label));
    rv = ptr;
  } else if (errmsg) {
    throw Exception(errmsg);
  }
  return rv;
}

// Parse serialization record and return bool value following label 'str' in
// line read from file pointer 'fp'.
inline bool ReadBool(const char* line,
                     const char* label,
                     const char* errmsg = NULL) {
  bool rv = false;
  if (strstr(line, label)) {
    const char* ptr = line + strlen(label);
    if (strchr(ptr, 'T') != NULL || strchr(ptr, '1') != NULL)
      rv = true;
    else if (strchr(ptr, 'F') != NULL || strchr(ptr, '0') != NULL)
      rv = false;
    else if (errmsg)
      throw Exception(errmsg);
  } else if (errmsg) {
    throw Exception(errmsg);
  }
  return rv;
}

// Returns true iff next non-blank line in file 'fp' contains string 'id'.
inline bool StreamStartsWith(FILE* fp, const char* id) {
  char buffer[KB];
  while (fgetline(buffer, KB, fp)) if (strscn(buffer)) break;
  return strstr(buffer, id) == buffer;
}

// Reads all serialized records of class type 'T' from stream 'fin' into vector.
template<class T>
void ReadAll(FILE* fin, std::vector<T>& vec, int max = -1) {
  int n = 0;
  while (!feof(fin) && (max == -1 || n < max)) {
    // Parse next record
    vec.push_back(T(fin));
    // Check for EOF
    int c = getc(fin);
    if (c == EOF) break;
    ungetc(c, fin);
    ++n;
  }
}

// Writes all serialized records of class type 'T' from vector to stream 'fout'
template<class T>
void WriteAll(const std::vector<T>& vec, FILE* fout) {
  typedef typename std::vector<T>::const_iterator ConstIter;
  for (ConstIter iter = vec.begin(); iter != vec.end(); ++iter)
    iter->Write(fout);
}

}  // namespace cs

#endif  // CS_IO_H_
