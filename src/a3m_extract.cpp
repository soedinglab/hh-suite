/*
 * a3m_extractor.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: meiermark
 */

#include "a3m_compress.h"

#include <iostream>
#include <getopt.h>

void usage() {
  std::cout
      << "USAGE: a3m_extract -i [inputfile|stdin] -o [outputfile|stdout] -d [ffindex_sequence_database_prefix] -q [ffindex_header_database_prefix]"
      << std::endl;
}

int main(int argc, char **argv) {
  bool iflag = false;
  bool dflag = false;
  bool oflag = false;
  bool qflag = false;

  std::string ffindex_sequence_db_prefix;
  std::string ffindex_header_db_prefix;
  std::string output;
  std::string input;

  int c;
  while ((c = getopt(argc, argv, "i:d:o:q:h")) != -1) {
    switch (c) {
      case 'i':
        iflag = 1;
        input = optarg;
        break;
      case 'd':
        dflag = 1;
        ffindex_sequence_db_prefix = optarg;
        break;
      case 'q':
        qflag = 1;
        ffindex_header_db_prefix = optarg;
        break;
      case 'o':
        oflag = optarg;
        output = optarg;
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

  if (!iflag || !dflag || !oflag || !qflag) {
    usage();
    exit(0);
  }

  //prepare ffindex sequence database
  std::string sequenceDataFile = ffindex_sequence_db_prefix + ".ffdata";
  std::string sequenceIndexFile = ffindex_sequence_db_prefix + ".ffindex";

  FILE *sequence_data_fh = fopen(sequenceDataFile.c_str(), "r");
  FILE *sequence_index_fh = fopen(sequenceIndexFile.c_str(), "r");

  if (sequence_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex sequence data file! ("
        << sequenceDataFile << ")!" << std::endl;
    exit(1);
  }

  if (sequence_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex sequence index file! ("
        << sequenceIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t sequence_data_size;
  char* sequence_data = ffindex_mmap_data(sequence_data_fh,
      &sequence_data_size);
  ffindex_index_t* sequence_index = ffindex_index_parse(sequence_index_fh, 1E8);

  if (sequence_index == NULL) {
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
    std::cerr << "Header index could not be loaded!" << std::endl;
    exit(1);
    //TODO: throw error
  }

  //prepare output stream
  std::ostream* out;
  if (output.compare("stdout") != 0) {
    out = new std::ofstream(output.c_str(), std::ios::binary | std::ios::out);
  }
  else {
    out = &std::cout;
  }

  //prepare input stream
  std::istream* in;
  if (input.compare("stdin") != 0) {
    in = new std::ifstream(input.c_str(), std::ios::binary | std::ios::in);
  }
  else {
    in = &std::cin;
  }

  in->seekg(0, in->end);
  size_t stream_size = in->tellg();
  in->seekg(0, in->beg);

  char* buffer_content_c = new char[stream_size];
  in->read(buffer_content_c, stream_size);

  compressed_a3m::extract_a3m(buffer_content_c, stream_size, sequence_index, sequence_data, header_index, header_data, out);
  out->flush();
}

