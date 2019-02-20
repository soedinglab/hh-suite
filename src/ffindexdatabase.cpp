#include "ffindexdatabase.h"
#include "hhutil.h"
#include <cstring>

#include <sys/mman.h>

FFindexDatabase::FFindexDatabase(char* data_filename, char* index_filename,
                                 bool isCompressed) {

    this->isCompressed = isCompressed;

    this->data_filename = new char[strlen(data_filename) + 1];
    strcpy(this->data_filename, data_filename);

    db_data_fh = fopen(data_filename, "r");
    if (db_data_fh == NULL) {
        OpenFileError(data_filename, __FILE__, __LINE__, __func__);
    }

    FILE* db_index_fh = fopen(index_filename, "r");
    if (db_index_fh == NULL) {
        OpenFileError(index_filename, __FILE__, __LINE__, __func__);
    }

    size_t ca3m_data_size = CountLinesInFile(index_filename);
    db_index = NULL;
    db_index = ffindex_index_parse(db_index_fh, ca3m_data_size);

    fclose(db_index_fh);

    if (db_index == NULL) {
        HH_LOG(WARNING) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
        HH_LOG(WARNING) << "\tCould not read index file" << index_filename << ". Is the file empty or corrupted?" << std::endl;
//    exit(1);
    }

    db_data = ffindex_mmap_data(db_data_fh, &data_size);
}

FFindexDatabase::~FFindexDatabase() {
    delete[] data_filename;
    munmap(db_data, data_size);
    free(db_index);
    fclose(db_data_fh);
}
