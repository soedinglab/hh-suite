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

#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <limits.h>

#include "fmemopen.h" /* For OS not yet implementing this new standard function */
#include "ffindex.h"

/* XXX Use page size? */
#define FFINDEX_BUFFER_SIZE 4096

#ifdef HH_MAC
size_t strnlen(const char *s, size_t n)
{
  const char *p = (const char *)memchr(s, 0, n);
  return(p ? p-s : n);
}

char* basename(char* path)
{
  char *ptr = strrchr (path, '/');
  return ptr ? ptr + 1 : (char*)path;
}
#endif

/* Insert all file from directory into ffindex */
int ffindex_insert_list_file(FILE *data_file, FILE *index_file, size_t *start_offset, FILE *list_file)
{
  size_t offset = *start_offset;
  char path[PATH_MAX];
  while(fgets(path, PATH_MAX, list_file) != NULL)
  {
    size_t len = strnlen(path, PATH_MAX);
    len -= 1;
    path[len] = '\0'; /* remove \n*/
    ffindex_insert_file(data_file, index_file, &offset, path, basename(path));
  }
  /* update return value */
  *start_offset = offset;
  return 0;
}


/* Insert all file from directory into ffindex */
int ffindex_insert_dir(FILE *data_file, FILE *index_file, size_t *start_offset, char *input_dir_name)
{
  DIR *dir = opendir(input_dir_name);
  if(dir == NULL)
    return -1;
  size_t input_dir_name_len = strnlen(input_dir_name, PATH_MAX);
  char path[PATH_MAX];
  strncpy(path, input_dir_name, NAME_MAX);
  if(input_dir_name[input_dir_name_len - 1] != '/')
  {
    path[input_dir_name_len] = '/';
    input_dir_name_len += 1;
  }
  size_t offset = *start_offset;
  struct dirent *entry;
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
    ffindex_insert_file(data_file, index_file, &offset, path, entry->d_name);
  }
  closedir(dir);

  /* update return value */
  *start_offset = offset;

  return 0;
}


/* Insert one file into ffindex */
int ffindex_insert_file(FILE *data_file, FILE *index_file, size_t *offset, char *path, char *name)
{
    FILE *file = fopen(path, "r");
    if(file == NULL)
      perror(path);

    /* copy and paste file to data file */
    char buffer[FFINDEX_BUFFER_SIZE];
    size_t offset_before = *offset;
    size_t read_size;
    while((read_size = fread(buffer, sizeof(char), sizeof(buffer), file)) > 0)
    {
      size_t write_size = fwrite(buffer, sizeof(char), read_size, data_file);
      *offset += write_size;
      if(read_size != write_size)
        perror(path); /* XXX handle better */
    }

    /* Seperate by '\0' and thus also make sure at least one byte is written */
    buffer[0] = '\0';
    fwrite(buffer, sizeof(char), 1, data_file); /* XXX check for error */
    *offset += 1;

    /* write index entry */
    fprintf(index_file, "%s\t%ld\t%ld\n", name, offset_before, *offset - offset_before);

    if(ferror(file) != 0 || ferror(data_file) != 0)
    {
      perror(path);
      exit(1);
    }
    fclose(file);
    return 0;
}

/* XXX not implemented yet */
int ffindex_restore(FILE *data_file, FILE *index_file, char *input_dir_name)
{
  return -1;
}


char* ffindex_mmap_data(FILE *file, size_t* size)
{
  struct stat sb;
  fstat(fileno(file), &sb);
  *size = sb.st_size;
  int fd =  fileno(file);
  if(fd < 0)
    return (char*)MAP_FAILED;
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
  return (ffindex_entry_t*)bsearch(&search, index->entries, index->n_entries, sizeof(ffindex_entry_t), ffindex_compare_entries_by_name);
}


ffindex_index_t* ffindex_index_parse(FILE *index_file, size_t num_max_entries)
{
  if(num_max_entries == 0)
    num_max_entries = FFINDEX_MAX_INDEX_ENTRIES_DEFAULT;
  size_t nbytes = sizeof(ffindex_index_t) + (sizeof(ffindex_entry_t) * num_max_entries);
  ffindex_index_t *index = (ffindex_index_t *)malloc(nbytes);
  if(index == NULL)
  {
    int myerrno = errno;
    char* errstr = strerror(myerrno);
    fprintf(stderr, "%s:%d ffindex_index_parse: malloc of %ld bytes failed: %s\n", __FILE__, __LINE__, nbytes ,errstr);
    return NULL;
  }

  index->file = index_file;
  index->index_data = ffindex_mmap_data(index_file, &(index->index_data_size));
  index->type = SORTED_FILE; /* Assume a sorted file for now */
  int i = 0;
  char* d = index->index_data;
  char* end;
  /* Faster than scanf per line */
  for(i = 0; d < (index->index_data + index->index_data_size); i++)
  {
    int p;
    for(p = 0; *d != '\t'; d++)
      index->entries[i].name[p++] = *d;
    index->entries[i].name[p] = '\0';
    index->entries[i].offset = strtol(d, &end, 10);
    d = end;
    index->entries[i].length  = strtol(d, &end, 10);
    d = end + 1; /* +1 for newline */
  }

  index->n_entries = i;

  if(index->n_entries == 0)
    return NULL;

  return index;
}


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


void ffindex_sort_index_file(ffindex_index_t *index)
{
  qsort(index->entries, index->n_entries, sizeof(ffindex_entry_t), ffindex_compare_entries_by_name);
}


int ffindex_write(ffindex_index_t* index, FILE* index_file)
{
  for(int i = 0; i < index->n_entries; i++)
  {
    ffindex_entry_t ffindex_entry = index->entries[i];
    fprintf(index_file, "%s\t%ld\t%ld\n", ffindex_entry.name, ffindex_entry.offset, ffindex_entry.length);
  }
  return EXIT_SUCCESS;
}

/* vim: ts=2 sw=2 et
*/
