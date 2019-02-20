/*
 * Ffindex
 * written by Andy Hauser <hauser@genzentrum.lmu.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * Ffindex is provided under the Create Commons license "Attribution-ShareAlike
 * 3.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/3.0/
 * 
 * Ffindex is a very simple database for small files. The files are stored
 * concatenated in one big data file, seperated by '\0'. A second file
 * contains a plain text index, giving name, offset and length of of the small
 * files.
 */

#include "ffutil.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

int fferror_print(char *sourcecode_filename, int line, const char *function_name, const char *message)
{
  int myerrno = errno;
  char* errstr = strerror(myerrno);
  fprintf(stderr, "%s:%d %s: %s: %s\n", sourcecode_filename , line, function_name, message, errstr);
  return myerrno;
}


/* remove \n, assumes UNIX line endings! */
char* ffnchomp(char *s, size_t len)
{
  // prevent underflow because of the size_t
  // we want to chomp off the last element
  len = len > 1 ? len - 1 : 0;
  if(s[len] == '\n')
    s[len] = '\0';

  return s;
}

size_t ffcount_lines(const char *filename)
{
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    return 0;
  }

  size_t lines = 0;
  int ch = 0;

  do {
    ch = fgetc(fp);
    if (ch == '\n') {
      lines++;
    }
  } while(ch != EOF);

  if(ch != '\n' && lines != 0) {
    lines++;
  }

  fclose(fp);

  return lines;
}

/* vim: ts=2 sw=2 et
*/
