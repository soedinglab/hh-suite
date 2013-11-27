#include "a3m_reducer.h"

int main(int argc, char **argv) {
  bool iflag, dflag, oflag;

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
  std::string sequenceDataFile = ffindex_sequence_db_prefix+"_disorder.ffdata";
  std::string sequenceIndexFile = ffindex_sequence_db_prefix+"_disorder.ffindex";

  FILE *sequence_data_fh  = fopen(sequenceDataFile.c_str(), "r");
  FILE *sequence_index_fh = fopen(sequenceIndexFile.c_str(), "r");

  if (sequence_data_fh == NULL) {
    std::cerr << "Could not open ffindex sequence data file! (" << sequenceDataFile << ")!" << std::endl;
    exit(1);
    //TODO: throw error
  }

  if(sequence_index_fh == NULL) {
    std::cerr << "Could not open ffindex sequence index file! (" << sequenceIndexFile << ")!" << std::endl;
    exit(1);
    //TODO: throw error
  }

  size_t sequence_data_size;
  char* sequence_data = ffindex_mmap_data(sequence_data_fh, &sequence_data_size);
  ffindex_index_t* sequence_index = ffindex_index_parse(sequence_index_fh, 0);

  if(sequence_index == NULL) {
    std::cerr << "Sequence index could not be loaded!" << std::endl;
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
    in = new std::ofstream(input.c_str(), std::ios::binary | std::ios::out);
  }
  else {
    in = &std::cin;
  }

  compressed_a3m::compress_a3m(in, sequence_index, sequence_data, out);
  out->flush();
}

void usage() {
  std::cout << "A3M_Reducer -i [inputfile|stdin] -o [outputfile|stdout] -d [ffindex_sequence_database_prefix]" << std::endl;
}

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
    if (line[0] == "#") {
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
    else if (line[0] == ">") {
      //process already possible saved sequence/consensus
      if (id.size() != 0) {
        if (consensus_flag) {
          //TODO: assumption header line before other sequences
          output << header;
          output << sequence << std::endl;
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

      id = header.substr(1, first_ws_index - 1);

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
      output << header;
      output << sequence << std::endl;
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

  char* full_sequence = ffindex_get_data_by_entry(
      ffindex_sequence_database_data, entry);
  if (full_sequence == NULL) {
    //TODO: proper errors
    std::cerr << "Warning: could not read sequence for " << id
        << " (data not found)" << std::endl;
    return;
  }

  unsigned short int start_pos = get_start_pos(aligned_sequence, full_sequence,
      entry->length);

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

  index = first_match_index;
  //count gaps
  while (index < aligned_sequence.size()) {
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
}

//returns pos not index
unsigned short int compressed_a3m::get_start_pos(std::string aligned_sequence,
    char* full_sequence, size_t full_sequence_length) {
  aligned_sequence.erase(
      std::remove(aligned_sequence.begin(), aligned_sequence.end(), '-'),
      aligned_sequence.end());

  for (short int i = 0; i < full_sequence_length; i++) {
    bool match = true;
    short int j;
    for (j = 0; j < aligned_sequence.size(); j++) {
      if (full_sequence[i] != toupper(aligned_sequence[j])) {
        match = false;
        break;
      }
    }

    if (match) {
      return j + 1;
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
  size_t first_ws_index = header.length() - 1;
  for (size_t index = 0; index < header.length(); index++) {
    if (isspace(line[index])) {
      first_ws_index = index;
      break;
    }
  }

  return header.substr(1, first_ws_index - 1);
}

bool isConsensus(std::string id) {
  return id.substr(id.length() - 11, 10) == "_consensus";
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

