#define FFINDEX_MAX_INDEX_ENTRIES 6000000
#define FFINDEX_MAX_ENTRY_NAME_LENTH 24

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


int ffindex_build(FILE *data_file, FILE *index_file, size_t *offset, char *input_dir_name);

FILE* ffindex_fopen(char *data, ffindex_index_t *index, char *filename);

char* ffindex_mmap_data(FILE *file, size_t* size);

//int ffindex_get_next_entry_by_name(FILE *index_file, char *entry_name, size_t *offset, size_t *length);

int ffindex_get_entry(FILE *index_file, char *filename, size_t *offset, size_t *length);

char* ffindex_get_filedata(char* data, size_t offset);

ffindex_index_t* ffindex_index_parse(FILE *index_file);

ffindex_entry_t* ffindex_bsearch_get_entry(ffindex_index_t *index, char *name);
