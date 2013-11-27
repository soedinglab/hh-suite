/*
 * a3m_compress.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: meiermark
 */

#include "a3m_compress.h"


void compressed_a3m::compress_a3m(std::istream* input,
    ffindex_index_t* ffindex_sequence_database_index,
    char* ffindex_sequence_database_data, std::ostream* output) {

  bool sequence_flag, consensus_flag = false;
  std::string sequence;
  std::string header;
  std::string id;

  std::string line;
  while (std::getline(*input, line)) {
    //comment - remove comments
    if (line[0] == '#') {
      ;
    }
    //ss_cons - remove ss annotation
    else if (line.substr(0, 8) == ">ss_cons") {
      //TODO: assumption just one line
      getline(*input, line);
      sequence_flag, consensus_flag = false;
    }
    //ss_pred - remove ss annotation
    else if (line.substr(0, 8) == ">ss_pred") {
      //TODO: assumption just one line
      getline(*input, line);
      sequence_flag, consensus_flag = false;
    }
    //header
    else if (line[0] == '>') {
      //process already possible saved sequence/consensus
      if (id.size() != 0) {
        if (consensus_flag) {
          //TODO: assumption header line before other sequences
          output->write(header.c_str(), header.size());
          output->put('\n');
          output->write(sequence.c_str(), sequence.size());
          output->put('\n');
          output->put(';');
        }
        else {
          compressed_a3m::compress_sequence(id, sequence,
              ffindex_sequence_database_index, ffindex_sequence_database_data,
              output);
        }

        sequence.clear();
        header.clear();
        id.clear();
      }

      //get id
      header = line;
      size_t first_ws_index = line.length() - 1;
      for (size_t index = 0; index < line.length(); index++) {
        if (isspace(line[index])) {
          first_ws_index = index;
          break;
        }
      }

      id = getNameFromHeader(header);

      //check if consensus or sequence
      consensus_flag = isConsensus(id);
      sequence_flag = !consensus_flag;
    }
    //sequence
    else if (sequence_flag || consensus_flag) {
      line = rtrim(line);
      sequence += line;
    }
  }

  //process last possible saved sequence/consensus
  if (id.size() != 0) {
    if (consensus_flag) {
      output->write(header.c_str(), header.size());
      output->put('\n');
      output->write(sequence.c_str(), sequence.size());
      output->put('\n');
      //TODO: Warning
    }
    else {
      compressed_a3m::compress_sequence(id, sequence,
          ffindex_sequence_database_index, ffindex_sequence_database_data,
          output);
    }
  }
}

void compressed_a3m::compress_sequence(std::string id,
    std::string aligned_sequence,
    ffindex_index_t* ffindex_sequence_database_index,
    char* ffindex_sequence_database_data, std::ostream* output) {

  ffindex_entry_t* entry = ffindex_get_entry_by_name(
      ffindex_sequence_database_index, const_cast<char*>(id.c_str()));


  if (entry == NULL) {
    //TODO: proper errors
    std::cerr << "WARNING: could not read sequence for " << id
        << " (entry not found)" << std::endl;
    return;
  }

  ffindex_entry_t* entry_zero = ffindex_get_entry_by_index(ffindex_sequence_database_index, 0);
  uint32_t entry_index = entry - entry_zero;

  char* full_sequence = ffindex_get_data_by_entry(
      ffindex_sequence_database_data, entry);
  if (full_sequence == NULL) {
    //TODO: proper errors
    std::cerr << "Warning: could not read sequence for " << id
        << " (data not found)" << std::endl;
    return;
  }

  writeU32(*output, entry_index);

  unsigned short int start_pos = get_start_pos(aligned_sequence, full_sequence,
      entry->length);

  writeU16(*output, start_pos);

  if (start_pos == 0) {
    //TODO: proper errors
    std::cerr << "WARNING: could not match aligned sequence to full sequence! ("
        << id << ")!" << std::endl;
    return;
  }

  //move to first upper case character
  unsigned short int index = 0;
  while ((!isupper(aligned_sequence[index]) || aligned_sequence[index] == '-')
      && index < aligned_sequence.size()) {
    index++;
  }

  unsigned short int first_match_index = index;

  //count gaps
  unsigned short int nr_blocks = 0;
  while (index < aligned_sequence.size()) {
    nr_blocks++;
    while (isupper(aligned_sequence[index]) && index < aligned_sequence.size()) {
      index++;
    }

    char nr_insertions = 0;
    while (islower(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_insertions++;
      index++;
    }

    char nr_gaps = 0;
    while (nr_insertions == 0 && aligned_sequence[index] == '-'
        && index < aligned_sequence.size()) {
      nr_gaps--;
      index++;
    }
  }

  writeU16(*output, nr_blocks);

  index = first_match_index;
  //count gaps
  while (index < aligned_sequence.size()) {
    unsigned short int nr_matches = 0;
    while (isupper(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_matches++;
      index++;
    }

    writeU16(*output, nr_matches);

    char nr_insertions = 0;
    while (islower(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_insertions++;
      index++;
    }

    char nr_gaps = 0;
    while (nr_insertions == 0 && aligned_sequence[index] == '-'
        && index < aligned_sequence.size()) {
      nr_gaps--;
      index++;
    }

    nr_insertions > 0?output->put(nr_insertions):output->put(nr_gaps);
  }
}

//returns pos not index
unsigned short int compressed_a3m::get_start_pos(std::string aligned_sequence,
    char* full_sequence, size_t full_sequence_length) {
  aligned_sequence.erase(
      std::remove(aligned_sequence.begin(), aligned_sequence.end(), '-'),
      aligned_sequence.end());

  for (unsigned short int i = 0; i < full_sequence_length; i++) {
    bool match = true;
    for (unsigned short int j = 0; j < aligned_sequence.size(); j++) {
      if (full_sequence[i+j] != toupper(aligned_sequence[j])) {
        match = false;
        break;
      }
    }

    if (match) {
      return i + 1;
    }
  }

  return 0;
}

// trim from end
std::string& rtrim(std::string &s) {
  s.erase(
      std::find_if(s.rbegin(), s.rend(),
          std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

std::string getNameFromHeader(std::string &header) {
  size_t first_ws_index = header.length();
  for (size_t index = 0; index < header.length(); index++) {
    if (isspace(header[index])) {
      first_ws_index = index;
      break;
    }
  }

  return header.substr(1, first_ws_index - 1);
}

bool isConsensus(std::string &id) {
  return id.substr(id.length() - 10, 10) == "_consensus";
}

void writeU16(std::ostream& file, uint16_t val) {
  unsigned char bytes[2];

  // extract the individual bytes from our value
  bytes[1] = (val) & 0xFF;  // low byte
  bytes[0] = (val >> 8) & 0xFF;  // high byte

  // write those bytes to the file
  file.write((char*) bytes, 2);
}

void readU16(char** ptr, uint16_t &result) {
  unsigned char array[2];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[1] | (array[0] << 8);
}

void writeU32(std::ostream& file, uint32_t val) {
  unsigned char bytes[4];

  // extract the individual bytes from our value
  bytes[3] = (val) & 0xFF;
  bytes[2] = (val >> 8) & 0xFF;
  bytes[1] = (val >> 16) & 0xFF;
  bytes[0] = (val >> 24) & 0xFF;

  // write those bytes to the file
  file.write((char*) bytes, 4);
}

void readU32(char** ptr, uint32_t &result) {
  unsigned char array[4];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;
  array[2] = (unsigned char) (**ptr);
  (*ptr)++;
  array[3] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[3] | (array[2] << 8) << (array[1] << 16) | (array[0] << 24);
}
