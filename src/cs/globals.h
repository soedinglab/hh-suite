// Copyright 2009, Andreas Biegert

#ifndef CS_GLOBALS_H_
#define CS_GLOBALS_H_

#include <limits>

namespace cs {

const int KB = 1024;
const int MB = KB * KB;
const int GB = KB * KB * KB;
const int kMaxInt = 0x7FFFFFFF;
const int kMinInt = -kMaxInt - 1;
const int kScale  = 1000; // scaling factor for serialization

const int kCharSize     = sizeof(char);
const int kShortSize    = sizeof(short);
const int kIntSize      = sizeof(int);
const int kFloatSize    = sizeof(float);
const int kDoubleSize   = sizeof(double);
const int kPointerSize  = sizeof(void*);

#ifdef _WIN32
const char kDirSep = '\\';
#else
const char kDirSep = '/';
#endif

// A macro to disallow the evil copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName)      \
  TypeName(const TypeName&);                    \
  void operator=(const TypeName&)

// A macro to disallow all the implicit constructors, namely the
// default constructor, copy constructor and operator= functions.
//
// This should be used in the private: declarations for a class
// that wants to prevent anyone from instantiating it. This is
// especially useful for classes containing only static methods.
#define DISALLOW_IMPLICIT_CONSTRUCTORS(TypeName) \
  TypeName();                                    \
  DISALLOW_COPY_AND_ASSIGN(TypeName)

}  // namespace cs

#endif  // CS_GLOBALS_H_



