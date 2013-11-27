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


extern "C" {
  #include <ffindex.h>     // fast index-based database reading
}

namespace compressed_a3m {

void compress_a3m(std::istream* input, ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data, std::ostream* output);
void compress_sequence(std::string id, std::string sequence, ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data, std::ostream* output);
unsigned short int get_start_pos(std::string aligned_sequence, char* full_sequence, size_t full_sequence_length);
}

std::string &rtrim(std::string &s);

std::string getNameFromHeader(std::string &header);
bool isConsensus(std::string &id);

void writeU16(std::ostream& file, uint16_t);


#endif /* A3M_REDUCER_H_ */
