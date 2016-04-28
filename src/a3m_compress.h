/*
 * a3m_reducer.h
 *
 *  Created on: Nov 15, 2013
 *      Author: meiermark
 */

#ifndef A3M_REDUCER_H_
#define A3M_REDUCER_H_

#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <stdint.h>
#include <cstring>
#include <climits>

extern "C" {
  #include <ffindex.h>     // fast index-based database reading
}

namespace compressed_a3m {
  int compress_a3m(std::istream* input, ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data, std::ostream* output);
  int compress_a3m(char* input, size_t input_size, ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data, std::ostream* output);

  int compress_sequence(std::string id, std::string sequence, ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data, std::ostream* output);

  void extract_a3m(char* data, size_t data_size,
      ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data,
      ffindex_index_t* ffindex_header_database_index, char* ffindex_header_data, std::ostream* output);


  unsigned short int get_start_pos(std::string aligned_sequence, char* full_sequence, size_t full_sequence_length);

  void writeU16(std::ostream& file, uint16_t);
  void readU16(char** ptr, uint16_t &result);

  void writeU32(std::ostream& file, uint32_t);
  void readU32(char**ptr, uint32_t &result);
}

std::string &rtrim(std::string &s);

std::string getNameFromHeader(std::string &header);
std::string getShortIdFromHeader(std::string &header);
bool isConsensus(std::string &id);

#endif /* A3M_REDUCER_H_ */
