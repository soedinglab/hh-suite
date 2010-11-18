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

#define FFINDEX_VERSION 0.9
#define FFINDEX_MAX_INDEX_ENTRIES 6000000
#define FFINDEX_MAX_ENTRY_NAME_LENTH 40

enum ffindex_type {PLAIN_FILE, SORTED_FILE, SORTED_ARRAY};

typedef struct ffindex_entry {
  char name[FFINDEX_MAX_ENTRY_NAME_LENTH];
  size_t offset;
  size_t length;
} ffindex_entry_t;

typedef struct ffindex_index {
  enum ffindex_type type;
  char* filename;
  FILE* file;
  char* index_data;
  size_t index_data_size;
  int n_entries;
  ffindex_entry_t entries[FFINDEX_MAX_INDEX_ENTRIES];
} ffindex_index_t;


int ffindex_insert(FILE *data_file, FILE *index_file, size_t *offset, char *input_dir_name);

FILE* ffindex_fopen(char *data, ffindex_index_t *index, char *filename);

char* ffindex_mmap_data(FILE *file, size_t* size);

char* ffindex_get_filedata(char* data, size_t offset);

ffindex_index_t* ffindex_index_parse(FILE *index_file);

ffindex_entry_t* ffindex_bsearch_get_entry(ffindex_index_t *index, char *name);

void ffindex_sort_index_file(ffindex_index_t *index);

int ffindex_write(ffindex_index_t* index, FILE* index_file);
