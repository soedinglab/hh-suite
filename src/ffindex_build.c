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

#define MAX_FILENAME_LIST_FILES 4096

void usage(char *program_name)
{
    fprintf(stderr, "USAGE: %s [-asv] [-f file]* data_filename index_filename [dirs_to_index ...]\n"
                    "\t-a\tappend\n"
                    "\t-f file\tfile each line containing a filename to index\n"
                    "\t\t\t-f can be specified up to MAX_FILENAME_LIST_FILES times\n"
                    "\t-s\tsort index file\n"
                    "\t-v\tprint version and other info then exit", program_name);
}

int main(int argn, char **argv)
{
  int append = 0, sort = 0, opt, err = EXIT_SUCCESS;
  char* list_filenames[MAX_FILENAME_LIST_FILES];
  int list_filenames_index = 0, version = 0;
  while ((opt = getopt(argn, argv, "asvf:")) != -1)
  {
    switch (opt)
    {
      case 'a':
        append = 1;
        break;
      case 'f':
        list_filenames[list_filenames_index++] = optarg;
        break;
      case 'v':
        version = 1;
      case 's':
        sort = 1;
        break;
      default:
        usage(argv[0]);
        return EXIT_FAILURE;
    }
  }

  if(version == 1)
  {
    printf("%s version %.2f, off_t = %ld bits\n", argv[0], FFINDEX_VERSION, sizeof(off_t) * 8);
    return EXIT_SUCCESS;
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

  /* open index and data file, seek to end if needed */
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

  /* For each list_file insert */
  if(list_filenames_index > 0)
    for(int i = 0; i < list_filenames_index; i++)
    {
      FILE *list_file = fopen(list_filenames[i], "r");
      if( list_file == NULL) { perror(list_filenames[i]); return EXIT_FAILURE; }
      if(ffindex_insert_list_file(data_file, index_file, &offset, list_file) < 0)
      {
        perror(list_filenames[i]);
        err = -1;
      }
    }

  /* For each dir, insert all files into the index */
  for(int i = optind; i < argn; i++)
    if(ffindex_insert_dir(data_file, index_file, &offset, argv[i]) < 0)
    {
      perror(argv[i]);
      err = -1;
    }
  fclose(data_file);

  if(sort)
  {
    rewind(index_file);
    ffindex_index_t* index = ffindex_index_parse(index_file, 0);
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
