/*
 * Written by Andy Hauser.
 */

#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1


#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>

#include "ffindex.h"

/* XXX Use page size? */
#define N 4096

int ffindex_build(FILE *data_file, FILE *index_file, size_t *base_offset, char *input_dir_name)
{
  DIR *dir = opendir(input_dir_name);
  if(dir == NULL)
    return -1;
  size_t input_dir_name_len = strlen(input_dir_name);
  char path[input_dir_name_len + NAME_MAX + 2];
  strncpy(path, input_dir_name, NAME_MAX);
  if(input_dir_name[input_dir_name_len - 1] != '/')
  {
    path[input_dir_name_len] = '/';
    input_dir_name_len += 1;
  }
  size_t offset = *base_offset;
  struct dirent *entry;
  char buffer[N];
  while((entry = readdir(dir)) != NULL)
  {
    if(entry->d_name[0] == '.')
      continue;
    strncpy(path + input_dir_name_len, entry->d_name, NAME_MAX);
    struct stat sb;
    if(stat(path, &sb) == -1)
    {
      perror("stat");
      exit(EXIT_FAILURE);
    }
    if(!S_ISREG(sb.st_mode))
      continue;
    FILE *file = fopen(path, "r");
    if(file == NULL)
      perror(path);

    /* Paste file to data file */
    size_t offset_start = offset;
    size_t read_size;
    while((read_size = fread(buffer, sizeof(char), sizeof(buffer), file)) > 0)
    {
      size_t write_size = fwrite(buffer, sizeof(char), read_size, data_file);
      offset += write_size;
      if(read_size != write_size)
        perror(path); /* XXX handle better */
    }

    /* Seperate by '\0' and make sure at least one byte is written */
    buffer[0] = 0;
    size_t write_size = fwrite(buffer, sizeof(char), 1, data_file);
    offset += 1;

    fprintf(index_file, "%s\t%ld\t%ld\n", entry->d_name, offset_start, offset - offset_start);

    if(ferror(file) != 0 || ferror(data_file) != 0)
    {
      perror(path);
      exit(1);
    }
    fclose(file);
  }
  closedir(dir);
  *base_offset = offset;
  return 0;
}


int ffindex_restore(FILE *data_file, FILE *index_file, char *input_dir_name)
{
  return 0;
}


char* ffindex_mmap_data(FILE *file, size_t* size)
{
  struct stat sb;
  if (fstat(fileno(file), &sb) < 0)
    {
      fprintf(stderr, "ERROR! Could not open file...\n");
      exit(2);
    }
  *size = sb.st_size;
  int fd =  fileno(file);
  if(fd < 0)
    return NULL;
  return (char*)mmap(NULL, *size, PROT_READ, MAP_PRIVATE, fd, 0);
}


static int ffindex_compare_entries_by_name(const void *pentry1, const void *pentry2)
{   
  ffindex_entry_t* entry1 = (ffindex_entry_t*)pentry1;
  ffindex_entry_t* entry2 = (ffindex_entry_t*)pentry2;
  return strncmp(entry1->name, entry2->name, FFINDEX_MAX_ENTRY_NAME_LENTH);
}


ffindex_entry_t* ffindex_bsearch_get_entry(ffindex_index_t *index, char *name)
{
  ffindex_entry_t search;
  strncpy(search.name, name, FFINDEX_MAX_ENTRY_NAME_LENTH);
  return (ffindex_entry_t*) bsearch(&search, index->entries, index->n_entries, sizeof(ffindex_entry_t), ffindex_compare_entries_by_name);
}


ffindex_index_t* ffindex_index_parse(FILE *index_file)
{
  ffindex_index_t *index = (ffindex_index_t*) calloc(1, sizeof(ffindex_index_t));
  if(index == NULL)
  {
    perror("ffindex_index_parse: calloc failed: ");
    exit(EXIT_FAILURE);
  }

  index->file = index_file;
  index->index_data = ffindex_mmap_data(index_file, &(index->index_data_size));
  int i = 0;
  char* d = index->index_data;
  char* end;
  while(d < (index->index_data + index->index_data_size))
  {
    for(int p = 0; *d != '\t'; d++)
      index->entries[i].name[p++] = *d;
    index->entries[i].offset = strtol(d, &end, 10);
    d = end;
    index->entries[i].length  = strtol(d, &end, 10);
    d = end + 1;
    i++;
  }


    /*
  while((n = fscanf(index->file, "%s\t%ld\t%ld\n", index->entries[i].name, &(index->entries[i].offset), &(index->entries[i].length))) == 3)
    i++;

  if(n < 0 && ferror(index->file))
  {
    perror("ffindex_index_parse: ");
    exit(EXIT_FAILURE);
  }
  else if(n > 0)
  {
    fprintf(stderr, "broken index file: wrong numbers of elements in line");
    exit(EXIT_FAILURE);
  }
  */

  index->n_entries = i;

  if(index->n_entries == 0)
    return NULL;

  return index;
}

/*
ffindex_entry_t*  ffindex_linear_get_entry(ffindex_index_t *index, char *name)
{
  return ffindex_get_next_entry_by_name(index, name);
}
*/


char* ffindex_get_filedata(char* data, size_t offset)
{
  return data + offset;
}


FILE* ffindex_fopen(char *data, ffindex_index_t *index, char *filename)
{
  ffindex_entry_t* entry = ffindex_bsearch_get_entry(index, filename);

  if(entry == NULL)
    return NULL;

  char *filedata = ffindex_get_filedata(data, entry->offset);
  return fmemopen(filedata, entry->length, "r");
}


/* vim: ts=2 sw=2 et
*/
