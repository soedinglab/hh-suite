#ifndef FFINDEX_DATABASE_H
#define FFINDEX_DATABASE_H

extern "C" {
#include <ffindex.h>
}

class FFindexDatabase {
public:
    FFindexDatabase(const char* data_filename, const char* index_filename, bool isCompressed);
    virtual ~FFindexDatabase();

    ffindex_index_t* db_index;
    char* db_data;

    char* data_filename;
    const bool isCompressed;

private:
    size_t data_size;
    FILE* db_data_fh;
};

#endif
