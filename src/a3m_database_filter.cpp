/*
 * a3m_database_filter.cpp
 *
 *  Created on: Feb 27, 2014
 *      Author: meiermark
 */

#include <set>
#include <string.h>
#include <fstream>
#include <sstream>

#include "a3m_compress.h"

extern "C" {
  #include "ffindex.h"
}

void usage() {
  std::cout << "USAGE: a3m_database_filter -i [ffindex_a3m_database_prefix] -o [ffindex_a3m_database_prefix] -s [filter]" << std::endl;
}

int main(int argc, char **argv) {
  bool iflag, sflag, oflag = false;

  std::string set_file;
  std::string ffindex_oa3m_db_prefix;
  std::string ffindex_a3m_db_prefix;

  int c;
  while ((c = getopt(argc, argv, "i:s:o:h")) != -1) {
    switch (c) {
      case 'i':
        iflag = 1;
        ffindex_a3m_db_prefix = optarg;
        break;
      case 's':
        sflag = 1;
        set_file = optarg;
        break;
      case 'o':
        oflag = optarg;
        ffindex_oa3m_db_prefix = optarg;
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

  if(!iflag || !sflag || !oflag) {
    std::cerr << "Missing input!" << std::endl;
    usage();
    exit(1);
  }

  //prepare ffindex a3m output database
  std::string oa3mDataFile = ffindex_oa3m_db_prefix+".ffdata";
  std::string oa3mIndexFile = ffindex_oa3m_db_prefix+".ffindex";

  FILE *oa3m_data_fh  = fopen(oa3mDataFile.c_str(), "w");
  FILE *oa3m_index_fh = fopen(oa3mIndexFile.c_str(), "w");

  if (oa3m_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex ca3m data file! (" << oa3mDataFile << ")!" << std::endl;
    exit(1);
  }

  if(oa3m_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex ca3m index file! (" << oa3mIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t oa3m_offset = 0;

  //prepare ffindex a3m database
  std::string a3mDataFile = ffindex_a3m_db_prefix+".ffdata";
  std::string a3mIndexFile = ffindex_a3m_db_prefix+".ffindex";

  FILE *a3m_data_fh  = fopen(a3mDataFile.c_str(), "r");
  FILE *a3m_index_fh = fopen(a3mIndexFile.c_str(), "r");

  if (a3m_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex a3m data file! (" << a3mDataFile << ")!" << std::endl;
    exit(1);
  }

  if(a3m_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex a3m index file! (" << a3mIndexFile << ")!" << std::endl;
    exit(1);
  }

  size_t a3m_offset;
  char* a3m_data = ffindex_mmap_data(a3m_data_fh, &a3m_offset);
  ffindex_index_t* a3m_index = ffindex_index_parse(a3m_index_fh, 0);

  if(a3m_index == NULL) {
    std::cerr << "ERROR: A3M index could not be loaded!" << std::endl;
    exit(1);
  }

  //prepare filter
  std::set<std::string> filter;
  std::ifstream infile(set_file.c_str());

  std::string line;
  while (std::getline(infile, line)) {
    std::string item = line.substr(0, line.length());
    filter.insert(item);
  }

  infile.close();

  //prepare input stream
  size_t a3m_range_start = 0;
  size_t a3m_range_end = a3m_index->n_entries;

  // Foreach entry
  #pragma omp parallel for shared(a3m_index, a3m_data, oa3m_data_fh, oa3m_index_fh, oa3m_offset)
  for(size_t entry_index = a3m_range_start; entry_index < a3m_range_end; entry_index++)
  {
    //fprintf(stderr, "index %ld\n", entry_index);
    ffindex_entry_t* entry = ffindex_get_entry_by_index(a3m_index, entry_index);
    if(entry == NULL) { perror(entry->name); continue; }

    char* data = ffindex_get_data_by_entry(a3m_data, entry);

    std::stringstream* out_buffer = new std::stringstream();

    size_t nr_sequences = 0;

    for(size_t index = 0; index < entry->length; index++) {
      //write annotation line
      if(data[index] == '#') {
        while(data[index] != '\n' && index < entry->length) {
          out_buffer->put(data[index++]);
        }
        out_buffer->put('\n');
      }
      else if(data[index] == '>') {
        size_t start_index = index;
        while(index < entry->length && data[index] != '\n') {
          index++;
        }

        //copy line without new line
        std::string header = std::string(&data[start_index], index - start_index);
        std::string id = getNameFromHeader(header);
        bool consensus_flag = isConsensus(id);

        std::string short_id = getShortIdFromHeader(header);

        while(index < entry->length - 1 && data[index] != '>') {
          index++;
        }
        if(data[index] == '>' || data[index] == '\0') {
          index--;
        }

        bool passedFilter = false;
        if(filter.find(short_id) != filter.end()) {
          nr_sequences++;
          passedFilter = true;
        }

        if(passedFilter ||
            consensus_flag ||
            id.compare("ss_dssp") == 0 ||
            id.compare("sa_dssp") == 0 ||
            id.compare("ss_pred") == 0 ||
            id.compare("ss_conf") == 0) {
          std::string seq = std::string(&data[start_index], index - start_index);
          out_buffer->write(seq.c_str(), seq.size());
          out_buffer->put('\n');
        }
      }
    }

    if(nr_sequences > 0) {
      std::string out_string = out_buffer->str();
      #pragma omp critical
      {
        ffindex_insert_memory(oa3m_data_fh, oa3m_index_fh, &oa3m_offset, const_cast<char*>(out_string.c_str()), out_string.size(), entry->name);
      }
    }
    else {
      std::cerr << "WARNING: No sequences left for cluster " << entry->name << std::endl;
    }

    delete out_buffer;
  }

  fflush(oa3m_index_fh);
  ffsort_index(oa3mIndexFile.c_str(), oa3m_index_fh);

  fclose(oa3m_index_fh);
  fclose(oa3m_data_fh);
}

