/*
 * ffindex_order
 * written by Milot Mirdita <milot@mirdita.de>.
 * Please add your name here if you distribute modified versions.
 * 
 * FFindex is provided under the Create Commons license "Attribution-ShareAlike
 * 4.0", which basically captures the spirit of the Gnu Public License (GPL).
 * 
 * See:
 * http://creativecommons.org/licenses/by-sa/4.0/
 *
 * ffindex_order
 * Reorders the entries in a FFindex data file by the order given by file
 * Each line of the order file must contain a key from the FFindex index.
 * The FFindex data file entries will have the same order as the order file.
*/

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include "ffindex.h"
#include "ffutil.h"


int main(int argc, char **argv)
{
  if(argc < 6)
  {
    fprintf(stderr, "USAGE: %s ORDER_FILENAME DATA_FILENAME INDEX_FILENAME SORTED_DATA_OUT_FILE SORTED_INDEX_OUT_FILE\n"
                    "\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n",
                    argv[0]);
    return -1;
  }
  
  char *order_filename = argv[1];
  
  char *data_filename  = argv[2];
  char *index_filename = argv[3];
  
  char *sorted_data_filename  = argv[4];
  char *sorted_index_filename = argv[5];

  FILE *order_file = fopen(order_filename, "r");

  FILE *data_file  = fopen(data_filename,  "r");
  FILE *index_file = fopen(index_filename, "r");

  FILE *sorted_data_file  = fopen(sorted_data_filename,  "w+");
  FILE *sorted_index_file = fopen(sorted_index_filename, "w+");

  if(order_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], order_filename);  exit(EXIT_FAILURE); }

  if( data_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], data_filename);  exit(EXIT_FAILURE); }
  if(index_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], index_filename);  exit(EXIT_FAILURE); }

  if( sorted_data_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], sorted_data_filename);  exit(EXIT_FAILURE); }
  if(sorted_index_file == NULL) { fferror_print(__FILE__, __LINE__, argv[0], sorted_index_filename);  exit(EXIT_FAILURE); }

  size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  size_t entries = ffcount_lines(index_filename);
  ffindex_index_t* index = ffindex_index_parse(index_file, entries);
  if(index == NULL)  {   
    perror("ffindex_index_parse failed");
    exit(EXIT_FAILURE);
  }


  char message[LINE_MAX];
  char line[LINE_MAX];
  int i = 0;
  size_t offset = 0;
  while (fgets(line, sizeof(line), order_file)) {
    size_t len = strlen(line);
    if (len && (line[len - 1] != '\n')) {
      // line is incomplete
      snprintf(message, LINE_MAX, "Warning: Line %d of order file %s was too long and cut-off.", i, order_filename);
      fferror_print(__FILE__, __LINE__, argv[0], message);
    }

    // remove new line
    char *name = ffnchomp(line, len);
    ffindex_entry_t* entry = ffindex_get_entry_by_name(index, name);

    if (entry != NULL) {
      char* filedata = ffindex_get_data_by_entry(data, entry);
      size_t entryLength = (entry->length == 0 ) ? 0 : entry->length - 1;
      ffindex_insert_memory(sorted_data_file, sorted_index_file, &offset, filedata, entryLength, name);
    }

    i++;
  }
  
  // cleanup
  fclose(sorted_data_file);
  fclose(index_file);
  fclose(data_file);
  fclose(order_file);

  // sort FFindex index
  fclose(sorted_index_file);
  sorted_index_file = fopen(sorted_index_filename, "r+");
  entries = ffcount_lines(index_filename);
  index = ffindex_index_parse(sorted_index_file, entries);
  if(index == NULL)  {   
    perror("ffindex_index_parse failed");
    exit(EXIT_FAILURE);
  }
  fclose(sorted_index_file);
  
  ffindex_sort_index_file(index);
  sorted_index_file = fopen(sorted_index_filename, "w");
  if(sorted_index_file == NULL) {
    perror(sorted_index_filename);
    return EXIT_FAILURE;
  }
  int err = ffindex_write(index, sorted_index_file);
  fclose(sorted_index_file);
 
  return err;
}

/* vim: ts=2 sw=2 et
 */
