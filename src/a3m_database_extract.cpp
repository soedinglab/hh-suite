/*
 * a3m_reduce_database.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: meiermark
 */


#include "a3m_compress.h"

#include <sstream>

void usage() {
  std::cout << "USAGE: a3m_database_extract -i [ffindex_ca3m_database_prefix] -o [ffindex_a3m_database_prefix] -d [ffindex_sequence_database_prefix] -q [ffindex_header_database_prefix]" << std::endl;
}

int main(int argc, char **argv) {
  bool iflag, dflag, oflag, qflag;
  iflag = dflag = oflag = qflag = false;

  std::string ffindex_header_db_prefix;
  std::string ffindex_sequence_db_prefix;
  std::string ffindex_ca3m_db_prefix;
  std::string ffindex_a3m_db_prefix;

  int c;
  while ((c = getopt(argc, argv, "i:d:o:q:h")) != -1) {
    switch (c) {
      case 'i':
        iflag = 1;
        ffindex_ca3m_db_prefix = optarg;
        break;
      case 'd':
        dflag = 1;
        ffindex_sequence_db_prefix = optarg;
        break;
      case 'o':
        oflag = 1;
        ffindex_a3m_db_prefix = optarg;
        break;
      case 'q':
        qflag = 1;
        ffindex_header_db_prefix = optarg;
        break;
      case 'h':
        usage();
        exit(0);
      case '?':
        if (optopt == 'c')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  if(!iflag || !dflag || !oflag || !qflag) {
    std::cerr << "Missing arguments!" << std::endl;
    usage();
    exit(0);
  }

  //prepare ffindex a3m database
  std::string a3mDataFile = ffindex_a3m_db_prefix+".ffdata";
  std::string a3mIndexFile = ffindex_a3m_db_prefix+".ffindex";

  FILE *a3m_data_fh  = fopen(a3mDataFile.c_str(), "w");
  FILE *a3m_index_fh = fopen(a3mIndexFile.c_str(), "w");

  if (a3m_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex a3m data file! (" << a3mDataFile << ")!" << std::endl;
    exit(1);
  }

  if(a3m_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex a3m index file! (" << a3mIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t a3m_offset = 0;

  //prepare ffindex ca3m database
  std::string ca3mDataFile = ffindex_ca3m_db_prefix+".ffdata";
  std::string ca3mIndexFile = ffindex_ca3m_db_prefix+".ffindex";

  FILE *ca3m_data_fh  = fopen(ca3mDataFile.c_str(), "r");
  FILE *ca3m_index_fh = fopen(ca3mIndexFile.c_str(), "r");

  if (ca3m_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex a3m data file! (" << ca3mDataFile << ")!" << std::endl;
    exit(1);
  }

  if(ca3m_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex a3m index file! (" << ca3mIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t ca3m_offset;
  char* ca3m_data = ffindex_mmap_data(ca3m_data_fh, &ca3m_offset);
  ffindex_index_t* ca3m_index = ffindex_index_parse(ca3m_index_fh, 0);

  if(ca3m_index == NULL) {
    std::cerr << "ERROR: CA3M index (" << ca3mIndexFile << ") could not be loaded!" << std::endl;
    exit(1);
  }

  //prepare ffindex sequence database
  std::string sequenceDataFile = ffindex_sequence_db_prefix+".ffdata";
  std::string sequenceIndexFile = ffindex_sequence_db_prefix+".ffindex";

  FILE *sequence_data_fh  = fopen(sequenceDataFile.c_str(), "r");
  FILE *sequence_index_fh = fopen(sequenceIndexFile.c_str(), "r");

  if (sequence_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex sequence data file! (" << sequenceDataFile << ")!" << std::endl;
    exit(1);
  }

  if(sequence_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex sequence index file! (" << sequenceIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t sequence_data_size;
  char* sequence_data = ffindex_mmap_data(sequence_data_fh, &sequence_data_size);
  ffindex_index_t* sequence_index = ffindex_index_parse(sequence_index_fh, 80000000);

  if(sequence_index == NULL) {
    std::cerr << "ERROR: Sequence index could not be loaded!" << std::endl;
    exit(1);
  }

  //prepare ffindex header database
  std::string headerDataFile = ffindex_header_db_prefix + ".ffdata";
  std::string headerIndexFile = ffindex_header_db_prefix + ".ffindex";

  FILE *header_data_fh = fopen(headerDataFile.c_str(), "r");
  FILE *header_index_fh = fopen(headerIndexFile.c_str(), "r");

  if (header_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex sequence data file! ("
        << headerDataFile << ")!" << std::endl;
    exit(1);
  }

  if (header_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex header index file! ("
        << headerIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t header_data_size;
  char* header_data = ffindex_mmap_data(header_data_fh,
      &header_data_size);
  ffindex_index_t* header_index = ffindex_index_parse(header_index_fh, 1E8);

  if (header_index == NULL) {
    std::cerr << "ERROR: Header index could not be loaded!" << std::endl;
    exit(1);
  }

  //prepare input stream
  size_t ca3m_range_start = 0;
  size_t ca3m_range_end = ca3m_index->n_entries;

  // Foreach entry
  #pragma omp parallel for shared(ca3m_index, ca3m_data, a3m_data_fh, a3m_index_fh, a3m_offset)
  for(size_t entry_index = ca3m_range_start; entry_index < ca3m_range_end; entry_index++)
  {
    ffindex_entry_t* entry = ffindex_get_entry_by_index(ca3m_index, entry_index);
    if(entry == NULL) { perror(entry->name); continue; }

    char* data = ffindex_get_data_by_entry(ca3m_data, entry);

    std::stringstream* out_buffer = new std::stringstream();
    compressed_a3m::extract_a3m(data, entry->length, sequence_index, sequence_data, header_index, header_data, out_buffer);

    std::string out_string = out_buffer->str();

    #pragma omp critical
    {
      ffindex_insert_memory(a3m_data_fh, a3m_index_fh, &a3m_offset, const_cast<char*>(out_string.c_str()), out_string.size(), entry->name);
    }

    delete out_buffer;
  }
  fflush(a3m_index_fh);
  ffsort_index(a3mIndexFile.c_str(), &a3m_index_fh);

  fclose(a3m_data_fh);
  fclose(a3m_index_fh);
}

