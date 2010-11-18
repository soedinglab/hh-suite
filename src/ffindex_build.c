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
*/

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "ffindex.h"

void usage(char *program_name)
{
    fprintf(stderr, "USAGE: %s [-as] data_filename index_filename [dirs_to_index ...]\n"
                    "\t-a append\n"
                    "\t-s sort index file\n", program_name);
}

int main(int argn, char **argv)
{
  int append = 0, sort = 0, opt, err = EXIT_SUCCESS;
  while ((opt = getopt(argn, argv, "as")) != -1)
  {
    switch (opt)
    {
      case 'a':
        append = 1;
        break;
      case 's':
        sort = 1;
        break;
      default:
        usage(argv[0]);
        return EXIT_FAILURE;
    }
  }
  if(argn - optind < 2)
  {
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  char *data_filename  = argv[optind++];
  char *index_filename = argv[optind++];
  FILE *data_file, *index_file;

  size_t offset = 0;

  if(append)
  {
    data_file  = fopen(data_filename, "a");
    index_file = fopen(index_filename, "a+");
    if( data_file == NULL) { perror(data_filename); return EXIT_FAILURE; }
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }

    struct stat sb;
    fstat(fileno(data_file), &sb);
    fseek(data_file, sb.st_size, SEEK_SET);
    offset = sb.st_size;

    fstat(fileno(index_file), &sb);
    fseek(index_file, sb.st_size, SEEK_SET);
  }
  else
  {
    data_file  = fopen(data_filename, "w");
    index_file = fopen(index_filename, "w+");
    if( data_file == NULL) { perror(data_filename); return EXIT_FAILURE; }
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
  }

  for(int i = optind; i < argn; i++)
    if(ffindex_insert(data_file, index_file, &offset, argv[i]) < 0)
    {
      perror(argv[i]);
      err = -1;
    }
  fclose(data_file);

  if(sort)
  {
    rewind(index_file);
    ffindex_index_t* index = ffindex_index_parse(index_file);
    fclose(index_file);
    ffindex_sort_index_file(index);
    index_file = fopen(index_filename, "w");
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
    err += ffindex_write(index, index_file);
  }

  return err;
}

/* vim: ts=2 sw=2 et
 */
