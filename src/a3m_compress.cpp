/*
 * a3m_compress.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: meiermark
 */

#include "a3m_compress.h"
#ifdef OPENMP
#include <omp.h>
#endif
int compressed_a3m::compress_a3m(std::istream* input,
    ffindex_index_t* ffindex_sequence_database_index,
    char* ffindex_sequence_database_data, std::ostream* output) {

  bool sequence_flag = false;
  bool consensus_flag = false;
  std::string sequence;
  std::string header;
  std::string id;

  int nr_sequences = 0;
  int nr_consensus = 0;

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
      sequence_flag = consensus_flag = false;
    }
    //ss_pred - remove ss annotation
    else if (line.substr(0, 8) == ">ss_pred") {
      //TODO: assumption just one line
      getline(*input, line);
      sequence_flag = consensus_flag = false;
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
          nr_consensus++;
        }
        else {
          if(compressed_a3m::compress_sequence(id, sequence,
              ffindex_sequence_database_index, ffindex_sequence_database_data,
              output)) {
            nr_sequences++;
          }
        }

        sequence.clear();
        header.clear();
        id.clear();
      }

      header = line;
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
      nr_consensus++;
      //TODO: Warning
    }
    else {
      if(compressed_a3m::compress_sequence(id, sequence,
          ffindex_sequence_database_index, ffindex_sequence_database_data,
          output)) {
        nr_sequences++;
      }
    }
  }

  if(nr_consensus > 1) {
    std::cerr << "More than one consensus sequence found in a3m!" << std::endl;
    return 0;
  }

  if(nr_sequences == 0) {
    std::cerr << "No protein sequences could be compressed in a3m!" << std::endl;
    return 0;
  }

  return 1;
}

int compressed_a3m::compress_a3m(char* input, size_t input_size,
    ffindex_index_t* ffindex_sequence_database_index,
    char* ffindex_sequence_database_data, std::ostream* output) {

  bool sequence_flag = false;
  bool consensus_flag = false;
  std::string sequence;
  std::string header;
  std::string id;

  int nr_sequences = 0;
  int nr_consensus = 0;

  size_t index = 0;
  while (index < input_size) {
    //comment - remove comments
    if (input[index] == '#') {
      if(index == 0) {
        while(input[index] != '\n' && index < input_size) {
          output->put(input[index++]);
        }
        output->put('\n');
      }
      else {
        while(input[index] != '\n' && index < input_size) {
          index++;
        }
      }
    }
    //ss_cons - remove ss annotation
    else if(input[index] == '>' && (strncmp(&input[index], ">ss_pred", 8) == 0 || strncmp(&input[index], ">ss_conf", 8) == 0)) {
      while(index + 1 < input_size && input[index + 1] != '>' && input[index] != '\n') {
        index++;
      }
    }
    //header
    else if (input[index] == '>') {
      //process already possible saved sequence/consensus
      if (id.size() != 0) {
        if (consensus_flag) {
          //TODO: assumption header line before other sequences
          output->write(header.c_str(), header.size());
          output->put('\n');
          output->write(sequence.c_str(), sequence.size());
          output->put('\n');
          output->put(';');
          nr_consensus++;
        }
        else {
          std::string short_id = getShortIdFromHeader(id);
          if(compressed_a3m::compress_sequence(short_id, sequence,
              ffindex_sequence_database_index, ffindex_sequence_database_data,
              output)) {
            nr_sequences++;
          }
        }

        sequence.clear();
        header.clear();
        id.clear();
      }

      size_t start_index = index;
      while(index < input_size && input[index] != '\n') {
        index++;
      }

      //copy line without new line
      header = std::string(&input[start_index], index - start_index);

      id = getNameFromHeader(header);

      //check if consensus or sequence
      consensus_flag = isConsensus(id);
      sequence_flag = !consensus_flag;
    }
    //sequence
    else if (sequence_flag || consensus_flag) {
      size_t start_index = index;
      while(index + 1 < input_size && input[index + 1] != '>' && input[index] != '\0') {
        if(input[index] == '\n') {
          sequence += std::string(&input[start_index], index - start_index);
          start_index = index + 1;
        }
        index++;
      }

      sequence += std::string(&input[start_index], index - start_index);
    }

    index++;
  }

  //process last possible saved sequence/consensus
  if (id.size() != 0) {
    if (consensus_flag) {
      output->write(header.c_str(), header.size());
      output->put('\n');
      output->write(sequence.c_str(), sequence.size());
      output->put('\n');
      nr_consensus++;
      //TODO: Warning
    }
    else {
      std::string short_id = getShortIdFromHeader(id);
      if(compressed_a3m::compress_sequence(short_id, sequence,
          ffindex_sequence_database_index, ffindex_sequence_database_data,
          output)) {
        nr_sequences++;
      }
    }
  }

  if(nr_consensus > 1) {
    std::cerr << "More than one consensus sequence found in a3m!" << std::endl;
    return 0;
  }

  if(nr_sequences == 0) {
    std::cerr << "No protein sequences could be compressed in a3m!" << std::endl;
    return 0;
  }

  return 1;
}


void compressed_a3m::extract_a3m(char* data, size_t data_size,
    ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data,
    ffindex_index_t* ffindex_header_database_index, char* ffindex_header_database_data, std::ostream* output) {

  //read stuff till compressed part
  char last_char = '\0';
  size_t index = 0;
  size_t consensus_length = 0;
  char inConsensus = 0;

  //handle commentary line
  if((*data) == '#') {
    while ((*data) != '\n') {
      output->put((*data));
      last_char = (*data);
      data++;
      index++;
    }

    output->put('\n');
    last_char = '\n';
    data++;
    index++;
  }

  while(!(last_char == '\n' && (*data) == ';') && index < data_size) {
    if((*data) == '\n') {
      inConsensus++;
    }
    else if(inConsensus == 1) {
      consensus_length++;
    }

    output->put((*data));
    last_char = (*data);
    data++;
    index++;
  }

  //get past ';'
  data++;
  index++;

  while(index < data_size - 1) {
    unsigned int entry_index;
    unsigned short int nr_blocks;
    unsigned short int start_pos;

    readU32(&data, entry_index);
    index += 4;

    ffindex_entry_t* sequence_entry = ffindex_get_entry_by_index(ffindex_sequence_database_index, entry_index);
    char* sequence = ffindex_get_data_by_entry(ffindex_sequence_database_data, sequence_entry);

    ffindex_entry_t* header_entry = ffindex_get_entry_by_index(ffindex_header_database_index, entry_index);
    char* header = ffindex_get_data_by_entry(ffindex_header_database_data, header_entry);

    // make sure we always have a valid fasta prefix
    if (header[0] != '>') {
      output->put('>');
    }

    output->write(header, header_entry->length - 1);
    output->put('\n');

    readU16(&data, start_pos);
    index += 2;

    readU16(&data, nr_blocks);
    index += 2;

    size_t actual_pos = start_pos;
    size_t alignment_length = 0;
    for(unsigned short int block_index = 0; block_index < nr_blocks; block_index++) {
      unsigned char nr_matches = (unsigned char)(*data);
      data++;
      index++;

      for(int i = 0; i < nr_matches; i++) {
        output->put(sequence[actual_pos - 1]);
        actual_pos++;
        alignment_length++;
      }

      char nr_insertions_deletions = (*data);
      data++;
      index++;

      if(nr_insertions_deletions > 0) {
        for(int i = 0; i < nr_insertions_deletions; i++) {
          output->put(tolower(sequence[actual_pos - 1]));
          actual_pos++;
        }
      }
      else {
        for(int i = 0; i < -nr_insertions_deletions; i++) {
          output->put('-');
          alignment_length++;
        }
      }
    }

    while(alignment_length < consensus_length) {
      output->put('-');
      alignment_length++;
    }

    output->put('\n');
  }
}

int compressed_a3m::compress_sequence(std::string id,
    std::string aligned_sequence,
    ffindex_index_t* ffindex_sequence_database_index,
    char* ffindex_sequence_database_data, std::ostream* output) {

  ffindex_entry_t* entry = ffindex_get_entry_by_name(
      ffindex_sequence_database_index, const_cast<char*>(id.c_str()));


  if (entry == NULL) {
    //TODO: proper errors
    std::cerr << "WARNING: could not read sequence for " << id
        << " (entry not found)" << std::endl;
    return 0;
  }

  ffindex_entry_t* entry_zero = ffindex_get_entry_by_index(ffindex_sequence_database_index, 0);
  uint32_t entry_index = entry - entry_zero;

  char* full_sequence = ffindex_get_data_by_entry(
      ffindex_sequence_database_data, entry);
  if (full_sequence == NULL) {
    //TODO: proper errors
    std::cerr << "Warning: could not read sequence for " << id
        << " (data not found)" << std::endl;
    return 0;
  }

  unsigned short int start_pos = get_start_pos(aligned_sequence, full_sequence,
      entry->length);

  if (start_pos == 0) {
    //TODO: proper errors
    std::cerr << "WARNING: could not match aligned sequence to full sequence! (" << id << ")!" << std::endl;
    return 0;
  }

  writeU32(*output, entry_index);
  writeU16(*output, start_pos);

  //count blocks
  unsigned short int index = 0;
  unsigned short int nr_blocks = 0;
  while (index < aligned_sequence.size()) {
    int nr_matches = 0;
    while (aligned_sequence[index] != '-' && isupper(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_matches++;
      index++;
    }

    int nr_insertions = 0;
    while (islower(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_insertions++;
      index++;
    }

    int nr_gaps = 0;
    while (nr_insertions == 0 && aligned_sequence[index] == '-'
        && index < aligned_sequence.size()) {
      nr_gaps++;
      index++;
    }

    while(nr_gaps != 0 || nr_insertions != 0 || nr_matches != 0) {
      if(index == aligned_sequence.size() && nr_matches == 0 && nr_insertions == 0) {
        break;
      }

      nr_matches -= std::min(nr_matches, SCHAR_MAX);
      nr_gaps -= std::min(nr_gaps, SCHAR_MAX);
      nr_insertions -= std::min(nr_insertions, SCHAR_MAX);

      nr_blocks++;
    }
  }

  writeU16(*output, nr_blocks);

  index = 0;
  //count gaps
  while (index < aligned_sequence.size()) {
    int nr_matches = 0;
    while (aligned_sequence[index] != '-' && isupper(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_matches++;
      index++;
    }

    int nr_insertions = 0;
    while (islower(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_insertions++;
      index++;
    }

    int nr_gaps = 0;
    while (nr_insertions == 0 && aligned_sequence[index] == '-' && index < aligned_sequence.size()) {
      nr_gaps++;
      index++;
    }

    while(nr_gaps != 0 || nr_insertions != 0 || nr_matches != 0) {
      if(index == aligned_sequence.size() && nr_matches == 0 && nr_insertions == 0) {
        break;
      }

      char print_matches = std::min(nr_matches, SCHAR_MAX);
      char print_gaps = std::min(nr_gaps, SCHAR_MAX);
      char print_insertions = std::min(nr_insertions, SCHAR_MAX);

      nr_matches -= print_matches;
      nr_gaps -= print_gaps;
      nr_insertions -= print_insertions;

      output->put(print_matches);
      print_insertions > 0?output->put(print_insertions):output->put((-1)*print_gaps);
    }
  }

  return 1;
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

std::string getShortIdFromHeader(std::string &id) {
  size_t first_pipe_index = id.length();
  size_t second_pipe_index = id.length();

  for (size_t index = 0; index < id.length(); index++) {
    if (id[index] == '|') {
      if(first_pipe_index == id.length()) {
        first_pipe_index = index;
      }
      else if(second_pipe_index == id.length()) {
        second_pipe_index = index;
      }
    }
  }

  std::string short_id = id;
  if(first_pipe_index != id.length() && second_pipe_index != id.length()) {
    short_id = id.substr(first_pipe_index + 1, second_pipe_index - first_pipe_index - 1);
  }

  return short_id;
}

bool isConsensus(std::string &id) {
  return id.length() > 11 && id.substr(id.length() - 10, 10) == "_consensus";
}

void compressed_a3m::writeU16(std::ostream& file, uint16_t val) {
  unsigned char bytes[2];

  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  bytes[1] = (val >> 8) & 0xFF;  // high byte

  // write those bytes to the file
  file.write((char*) bytes, 2);
}

void compressed_a3m::readU16(char** ptr, uint16_t &result) {
  unsigned char array[2];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8);
}

void compressed_a3m::writeU32(std::ostream& file, uint32_t val) {
  unsigned char bytes[4];

  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;
  bytes[1] = (val >> 8) & 0xFF;
  bytes[2] = (val >> 16) & 0xFF;
  bytes[3] = (val >> 24) & 0xFF;

  // write those bytes to the file
  file.write((char*) bytes, 4);
}

void compressed_a3m::readU32(char** ptr, uint32_t &result) {
  unsigned char array[4];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;
  array[2] = (unsigned char) (**ptr);
  (*ptr)++;
  array[3] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);
}
