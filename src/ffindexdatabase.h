#ifndef FFINDEX_DATABASE_H
#define FFINDEX_DATABASE_H

extern "C" {
#include <ffindex.h>
}

class FFindexDatabase {
public:
    FFindexDatabase(const char* data_filename, const char* index_filename, bool isCompressed);
    virtual ~FFindexDatabase();

    void ensureLinearAccess();

    ffindex_index_t* db_index;
    char* db_data;

    char* data_filename;
    const bool isCompressed;

private:
    size_t data_size;
    FILE* db_data_fh;

    struct compareEntryByOffset {
        bool operator() (const ffindex_entry_t& lhs, const ffindex_entry_t& rhs) const {
            return (lhs.offset < rhs.offset);
        }
    };
};

#endif
