#include "a3m_compress.h"

#include <iostream>
#include <getopt.h>
#include <sstream>

void usage() {
  std::cout << "USAGE: a3m_reduce -i [inputfile|stdin] -o [outputfile|stdout] -d [ffindex_sequence_database_prefix]" << std::endl;
}

int main(int argc, char **argv) {
  bool iflag, dflag, oflag = false;

  std::string ffindex_sequence_db_prefix;
  std::string output;
  std::string input;

  int c;
  while ((c = getopt(argc, argv, "i:d:o:h")) != -1) {
    switch (c) {
      case 'i':
        iflag = 1;
        input = optarg;
        break;
      case 'd':
        dflag = 1;
        ffindex_sequence_db_prefix = optarg;
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

  if(!iflag || !dflag || !oflag) {
    usage();
    exit(0);
  }

  //prepare ffindex_database
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

  //prepare input stream
  std::istream* in;
  if (input.compare("stdin") != 0) {
    in = new std::ifstream(input.c_str(), std::ios::binary | std::ios::in);
  }
  else {
    in = &std::cin;
  }

  std::stringstream* out_buffer = new std::stringstream();
  int ret = compressed_a3m::compress_a3m(in, sequence_index, sequence_data, out_buffer);

  if(ret) {
    //prepare output
    if (output.compare("stdout") != 0) {
      std::ofstream out(output.c_str(), std::ios::binary | std::ios::out);
      out << out_buffer->str();
      out.close();
    }
    else {
      std::cout << out_buffer->str();
    }
    return 0;
  }
  else {
    std::cerr << "ERROR: Could not compress A3M! ("<< input << ")" << std::endl;
    return 1;
  }
}

