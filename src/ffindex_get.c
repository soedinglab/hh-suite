/*
 * Written by Andy Hauser.
 */

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1

#include <stdio.h>
#include <limits.h>

#include "ffindex.h"


int main(int argn, char **argv)
{
  if(argn < 3)
  {
    fprintf(stderr, "USAGE: %s data_filename index_filename filename(s)\n", argv[0]);
    return -1;
  }
  char *data_filename  = argv[1];
  char *index_filename = argv[2];

  FILE *data_file  = fopen(data_filename, "r");
  FILE *index_file = fopen(index_filename, "r");

  if( data_file == NULL) { perror(data_filename); return 1; }
  if(index_file == NULL) { perror(index_filename); return 1; }

  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  ffindex_index_t* index = ffindex_index_parse(index_file);
  if(index == NULL)
  {
    perror("no index:");
    return 1;
  }

  for(int i = 3; i < argn; i++)
  {
    char *filename = argv[i];
    FILE *file = ffindex_fopen(data, index, filename);
    char line[LINE_MAX];
    while(fgets(line, LINE_MAX, file) != NULL)
      printf("%s", line); /* XXX Mask nonprintable characters */
  }

  return 0;
}

/* vim: ts=2 sw=2 et
 */
