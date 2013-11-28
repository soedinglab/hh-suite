/*
 * a3m_compress.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: meiermark
 */

#include "a3m_compress.h"


int compressed_a3m::compress_a3m(std::istream* input,
    ffindex_index_t* ffindex_sequence_database_index,
    char* ffindex_sequence_database_data, std::ostream* output) {

  bool sequence_flag, consensus_flag = false;
  std::string sequence;
  std::string header;
  std::string id;

  int nr_sequences = 0;
  int nr_consensus = 0;

  std::string line;
  while (std::getline(*input, line)) {
    std::cout << "klug" << line << "depp" << std::endl;
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

  bool sequence_flag, consensus_flag = false;
  std::string sequence;
  std::string header;
  std::string id;

  int nr_sequences = 0;
  int nr_consensus = 0;

  size_t index = 0;
  while (index < input_size) {
    //comment - remove comments
    if (input[index] == '#') {
      while(input[index] != '\n' && index < input_size){
        index++;
      }
    }
    //ss_cons - remove ss annotation
    else if(strncmp(&input[index], ">ss_pred", 8) == 0) {
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


void compressed_a3m::extract_a3m(char* data, size_t data_size,
    ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data,
    ffindex_index_t* ffindex_header_database_index, char* ffindex_header_database_data, std::ostream* output) {

  //read stuff till compressed part
  char last_char = '\0';
  size_t index = 0;
  size_t consensus_length = 0;
  char inConsensus = 0;
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

  while(index < data_size) {
    unsigned int entry_index;
    unsigned short int nr_blocks;
    unsigned short int start_pos;
    unsigned short int nr_matches;

    char nr_insertions_deletions;

    readU32(&data, entry_index);
    index += 4;

    ffindex_entry_t* sequence_entry = ffindex_get_entry_by_index(ffindex_sequence_database_index, entry_index);
    char* sequence = ffindex_get_data_by_entry(ffindex_sequence_database_data, sequence_entry);
    //TODO: catch errors

    ffindex_entry_t* header_entry = ffindex_get_entry_by_index(ffindex_header_database_index, entry_index);
    char* header = ffindex_get_data_by_entry(ffindex_header_database_data, header_entry);

    output->write(header, header_entry->length);
    output->put('\n');

    readU16(&data, start_pos);
    index += 2;

    readU16(&data, nr_blocks);
    index += 2;

    size_t actual_pos = start_pos;
    size_t alignment_length = 0;
    for(unsigned short int block_index = 0; block_index < nr_blocks; block_index++) {
      readU16(&data, nr_matches);
      index += 2;

      for(int i = 0; i < nr_matches; i++) {
        output->put(sequence[actual_pos - 1]);
        actual_pos++;
        alignment_length++;
      }

      nr_insertions_deletions = (*data);
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

  writeU32(*output, entry_index);

  unsigned short int start_pos = get_start_pos(aligned_sequence, full_sequence,
      entry->length);

  if (start_pos == 0) {
    //TODO: proper errors
    std::cerr << "WARNING: could not match aligned sequence to full sequence! (" << id << ")!" << std::endl;
    std::cerr << aligned_sequence << std::endl;
    std::cerr << full_sequence << std::endl;
    return 0;
  }

  writeU16(*output, start_pos);


  //move to first upper case character
//  while ((!isupper(aligned_sequence[index]) || aligned_sequence[index] == '-')
//      && index < aligned_sequence.size()) {
//    index++;
//  }


  //count blocks
  unsigned short int index = 0;
  unsigned short int nr_blocks = 0;
  while (index < aligned_sequence.size()) {
    unsigned short int nr_matches = 0;
    while (aligned_sequence[index] != '-' && isupper(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_matches++;
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

    if(!(index == aligned_sequence.size() && nr_matches == 0)) {
      nr_blocks++;
    }
  }

  writeU16(*output, nr_blocks);

  index = 0;
  //count gaps
  while (index < aligned_sequence.size()) {
    unsigned short int nr_matches = 0;
    while (aligned_sequence[index] != '-' && isupper(aligned_sequence[index]) && index < aligned_sequence.size()) {
      nr_matches++;
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

    if(!(index == aligned_sequence.size() && nr_matches == 0)) {
      writeU16(*output, nr_matches);
      nr_insertions > 0?output->put(nr_insertions):output->put(nr_gaps);
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

bool isConsensus(std::string &id) {
  return id.substr(id.length() - 10, 10) == "_consensus";
}

void writeU16(std::ostream& file, uint16_t val) {
  unsigned char bytes[2];

  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;  // low byte
  bytes[1] = (val >> 8) & 0xFF;  // high byte

  // write those bytes to the file
  file.write((char*) bytes, 2);
}

void readU16(char** ptr, uint16_t &result) {
  unsigned char array[2];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8);
}

void writeU32(std::ostream& file, uint32_t val) {
  unsigned char bytes[4];

  // extract the individual bytes from our value
  bytes[0] = (val) & 0xFF;
  bytes[1] = (val >> 8) & 0xFF;
  bytes[2] = (val >> 16) & 0xFF;
  bytes[3] = (val >> 24) & 0xFF;

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

  result = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);
}
