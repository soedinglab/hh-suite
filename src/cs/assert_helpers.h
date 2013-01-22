/*
  Copyright 2012 Andreas Biegert

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

#ifndef CS_ASSERT_HELPERS_H_
#define CS_ASSERT_HELPERS_H_

#include <string.h>

namespace cs {

extern "C" void CS_Fatal(const char* file, int line, const char* format, ...);

// The FATAL, UNREACHABLE and UNIMPLEMENTED macros are useful during
// development, but they should not be relied on in the final product.
#ifndef NDEBUG
#define FATAL(msg)                              \
  CS_Fatal(__FILE__, __LINE__, "%s", (msg))
#define UNIMPLEMENTED()                         \
  CS_Fatal(__FILE__, __LINE__, "unimplemented code")
#define UNREACHABLE()                           \
  CS_Fatal(__FILE__, __LINE__, "unreachable code")
#else
#define FATAL(msg)                              \
  CS_Fatal("", 0, "%s", (msg))
#define UNIMPLEMENTED()                         \
  CS_Fatal("", 0, "unimplemented code")
#define UNREACHABLE() ((void) 0)
#endif


// Used by the ASSERT macro -- should not be called directly.
static inline void AssertHelper(const char* file,
                                int line,
                                const char* source,
                                bool condition) {
  if (!condition)
    CS_Fatal(file, line, "assert(%s) failed", source);
}


// Helper function used by the ASSERT_EQ function when given int
// arguments.  Should not be called directly.
static inline void AssertEqualsHelper(const char* file, int line,
                                      const char* expected_source,
                                      int expected,
                                      const char* value_source,
                                      int value) {
  if (expected != value) {
    CS_Fatal(file, line,
             "ASSERT_EQ(%s, %s) failed\n#   Expected: %i\n#   Found: %i",
             expected_source, value_source, expected, value);
  }
}


// Helper function used by the ASSERT_EQ function when given int
// arguments.  Should not be called directly.
static inline void AssertEqualsHelper(const char* file, int line,
                                      const char* expected_source,
                                      size_t expected,
                                      const char* value_source,
                                      size_t value) {
  if (expected != value) {
    CS_Fatal(file, line,
             "ASSERT_EQ(%s, %s) failed\n#   Expected: %i\n#   Found: %i",
             expected_source, value_source, expected, value);
  }
}


// Helper function used by the ASSERT_NE function when given int
// arguments.  Should not be called directly.
static inline void AssertNonEqualsHelper(const char* file,
                                         int line,
                                         const char* unexpected_source,
                                         int unexpected,
                                         const char* value_source,
                                         int value) {
  if (unexpected == value) {
    CS_Fatal(file, line, "ASSERT_NE(%s, %s) failed\n#   Value: %i",
             unexpected_source, value_source, value);
  }
}


// Helper function used by the ASSERT function when given string
// arguments.  Should not be called directly.
static inline void AssertEqualsHelper(const char* file,
                                      int line,
                                      const char* expected_source,
                                      const char* expected,
                                      const char* value_source,
                                      const char* value) {
  if (strcmp(expected, value) != 0) {
    CS_Fatal(file, line,
             "ASSERT_EQ(%s, %s) failed\n#   Expected: %s\n#   Found: %s",
             expected_source, value_source, expected, value);
  }
}


static inline void AssertNonEqualsHelper(const char* file,
                                         int line,
                                         const char* expected_source,
                                         const char* expected,
                                         const char* value_source,
                                         const char* value) {
  if (expected == value ||
      (expected != NULL && value != NULL && strcmp(expected, value) == 0)) {
    CS_Fatal(file, line, "ASSERT_NE(%s, %s) failed\n#   Value: %s",
             expected_source, value_source, value);
  }
}


// Helper function used by the ASSERT function when given pointer
// arguments.  Should not be called directly.
static inline void AssertEqualsHelper(const char* file,
                                      int line,
                                      const char* expected_source,
                                      void* expected,
                                      const char* value_source,
                                      void* value) {
  if (expected != value) {
    CS_Fatal(file, line,
             "ASSERT_EQ(%s, %s) failed\n#   Expected: %p\n#   Found: %p",
             expected_source, value_source,
             expected, value);
  }
}


static inline void AssertNonEqualsHelper(const char* file,
                                        int line,
                                        const char* expected_source,
                                        void* expected,
                                        const char* value_source,
                                        void* value) {
  if (expected == value) {
    CS_Fatal(file, line, "ASSERT_NE(%s, %s) failed\n#   Value: %p",
             expected_source, value_source, value);
  }
}


// Helper function used by the ASSERT function when given floating
// point arguments.  Should not be called directly.
static inline void AssertEqualsHelper(const char* file,
                                     int line,
                                     const char* expected_source,
                                     double expected,
                                     const char* value_source,
                                     double value) {
  // Force values to 64 bit memory to truncate 80 bit precision on IA32.
  volatile double* exp = new double[1];
  *exp = expected;
  volatile double* val = new double[1];
  *val = value;
  if (*exp != *val) {
    CS_Fatal(file, line,
             "ASSERT_EQ(%s, %s) failed\n#   Expected: %f\n#   Found: %f",
             expected_source, value_source, *exp, *val);
  }
  delete[] exp;
  delete[] val;
}


// Helper function used by the ASSERT function when given floating
// point arguments.  Should not be called directly.
static inline void AssertNonEqualsHelper(const char* file,
                                         int line,
                                         const char* unexpected_source,
                                         double unexpected,
                                         const char* value_source,
                                         double value) {
  // Force values to 64 bit memory to truncate 80 bit precision on IA32.
  volatile double* exp = new double[1];
  *exp = unexpected;
  volatile double* val = new double[1];
  *val = value;
  if (*exp == *val) {
    CS_Fatal(file, line, "ASSERT_NE(%s, %s) failed\n#   Value: %i",
             unexpected_source, value_source, value);
  }
  delete[] exp;
  delete[] val;
}


// The ASSERT macro only generates code in debug builds.
#undef assert
#ifndef NDEBUG
#define assert_result(expr)  AssertHelper(__FILE__, __LINE__, #expr, expr)
#define assert(condition)    AssertHelper(__FILE__, __LINE__, \
                                         #condition, condition)
#define assert_eq(v1, v2)    AssertEqualsHelper(__FILE__, __LINE__, \
                                               #v1, v1, #v2, v2)
#define assert_ne(v1, v2)    AssertNonEqualsHelper(__FILE__, __LINE__, \
                                                  #v1, v1, #v2, v2)
#else
#define assert_result(expr)  (expr)
#define assert(condition)    ((void) 0)
#define assert_eq(v1, v2)    ((void) 0)
#define assert_ne(v1, v2)    ((void) 0)
#endif

}

#endif  // CS_ASSERT_HELPERS_H_
