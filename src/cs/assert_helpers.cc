#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "cs.h"

// Prints a fatal error message and ends the program.
extern "C" void CS_Fatal(const char* file, int line, const char* format, ...) {
  fprintf(stderr, "\n\n#\n# Fatal error in %s, line %d\n# ", file, line);
  va_list arguments;
  va_start(arguments, format);
  vfprintf(stderr, format, arguments);
  va_end(arguments);
  fprintf(stderr, "\n#\n\n");
  abort();
}
