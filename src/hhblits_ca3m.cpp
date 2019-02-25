/*
 * hhblits_mpi.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include "hhblits.h"

#ifdef OPENMP
#include <omp.h>
#endif

struct OutputFFIndex {
    char base[NAMELEN];
    FILE* data_fh;
    FILE* index_fh;
    size_t offset;
    size_t number_entries;
    void (*print)(HHblits&, std::stringstream&);

    void close() {
      char index_filename[NAMELEN];
      snprintf(index_filename, FILENAME_MAX, "%s.ffindex", base);

      fclose(index_fh);
      fclose(data_fh);

      ffsort_index(index_filename);
    }

    void saveOutput(HHblits& hhblits, char* name) {
      std::stringstream out;
      print(hhblits, out);

      std::string tmp = out.str();
      ffindex_insert_memory(data_fh, index_fh, &offset,
          const_cast<char*>(tmp.c_str()), tmp.size(), name);

      fflush(data_fh);
      fflush(index_fh);
      number_entries++;
    }
};


void makeOutputFFIndex(char* par, void (*print)(HHblits&, std::stringstream&),
    std::vector<OutputFFIndex>& outDatabases) {
  if (*par) {
    OutputFFIndex db;

    strcpy(db.base, par);
    db.offset = 0;
    db.print = print;
    db.number_entries = 0;

    char data_filename_out_rank[NAMELEN];
    char index_filename_out_rank[NAMELEN];

    snprintf(data_filename_out_rank, FILENAME_MAX, "%s.ffdata", par);
    snprintf(index_filename_out_rank, FILENAME_MAX, "%s.ffindex", par);

    db.data_fh = fopen(data_filename_out_rank, "w+");
    db.index_fh = fopen(index_filename_out_rank, "w+");

    if (db.data_fh == NULL) {
      HH_LOG(ERROR) << "Could not open datafile " << data_filename_out_rank << "!" << std::endl;
      return;
    }

    if (db.index_fh == NULL) {
      HH_LOG(ERROR) << "Could not open indexfile " << index_filename_out_rank << "!" << std::endl;
      return;
    }

    outDatabases.push_back(db);
  }
}

int main(int argc, const char **argv) {
  Parameters par(argc, argv);
  HHblits::ProcessAllArguments(par);

  char data_filename[NAMELEN];
  char index_filename[NAMELEN];

  strcpy(data_filename, par.infile);
  strcat(data_filename, ".ffdata");

  strcpy(index_filename, par.infile);
  strcat(index_filename, ".ffindex");

  //FILE *data_file = fopen(data_filename, "r");
  //FILE *index_file = fopen(index_filename, "r");

/*
  if (data_file == NULL) {
    HH_LOG(ERROR) << "Input data file " << data_filename << " does not exist!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (index_file == NULL) {
    HH_LOG(ERROR) << "Input index file " << index_filename << " does not exist!" << std::endl;
    exit(EXIT_FAILURE);
  }*/

  //init input ffindex
  /*size_t data_size;
  char *data = ffindex_mmap_data(data_file, &data_size);

  size_t number_input_index_lines = CountLinesInFile(index_filename);
  ffindex_index_t* index = ffindex_index_parse(index_file, number_input_index_lines);
  if (index == NULL) {
    HH_LOG(ERROR) << "Could not parse index from " << index_filename << std::endl;
    exit(EXIT_FAILURE);
  }*/
  
  
  
  
  
  std::string ffindex_header_db_prefix(par.infile); ffindex_header_db_prefix += "_header";
  std::string ffindex_sequence_db_prefix(par.infile); ffindex_sequence_db_prefix += "_sequence";
  std::string ffindex_ca3m_db_prefix(par.infile); ffindex_ca3m_db_prefix += "_ca3m";
  


  //prepare ffindex ca3m database
  std::string ca3mDataFile = ffindex_ca3m_db_prefix + ".ffdata";
  std::string ca3mIndexFile = ffindex_ca3m_db_prefix + ".ffindex";

  FILE *ca3m_data_fh  = fopen(ca3mDataFile.c_str(), "r");
  FILE *ca3m_index_fh = fopen(ca3mIndexFile.c_str(), "r");

  if (ca3m_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex ca3m data file! (" << ca3mDataFile << ")!" << std::endl;
    exit(1);
  }

  if(ca3m_index_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex ca3m index file! (" << ca3mIndexFile << ")!" << std::endl;
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
  std::string sequenceDataFile = ffindex_sequence_db_prefix + ".ffdata";
  std::string sequenceIndexFile = ffindex_sequence_db_prefix + ".ffindex";

  FILE *sequence_data_fh  = fopen(sequenceDataFile.c_str(), "r");
  FILE *sequence_index_fh = fopen(sequenceIndexFile.c_str(), "r");

  if (sequence_data_fh == NULL) {
    std::cerr << "ERROR: Could not open ffindex sequence data file! (" << sequenceDataFile << ")!" << std::endl;
    exit(1);
  } else {
      std::cout<<"Data file: "<< sequenceDataFile <<std::endl;
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


  
  

  std::vector<OutputFFIndex> outputDatabases;
  makeOutputFFIndex(par.outfile, &HHblits::writeHHRFile, outputDatabases);
  makeOutputFFIndex(par.scorefile, &HHblits::writeScoresFile, outputDatabases);
  makeOutputFFIndex(par.pairwisealisfile, &HHblits::writePairwiseAlisFile, outputDatabases);
  makeOutputFFIndex(par.alitabfile, &HHblits::writeAlitabFile, outputDatabases);
  makeOutputFFIndex(par.psifile, &HHblits::writePsiFile, outputDatabases);
  makeOutputFFIndex(par.hhmfile, &HHblits::writeHMMFile, outputDatabases);
  makeOutputFFIndex(par.alnfile, &HHblits::writeA3MFile, outputDatabases);
  makeOutputFFIndex(par.matrices_output_file, &HHblits::writeMatricesFile, outputDatabases);
  makeOutputFFIndex(par.m8file, &HHblits::writeM8, outputDatabases);

  std::vector<HHblitsDatabase*> databases;
  HHblits::prepareDatabases(par, databases);

  //need to save for cleanup in the end
  int threads = par.threads;

#ifdef OPENMP
  omp_set_num_threads(threads);
#endif

  //no openmp parallelization in hhblits methods

  HHblits* hhblits_instances = new HHblits[par.threads];
  par.threads = 1;

  for(int i = 0; i < threads; i++) {
    hhblits_instances[i] = new HHblits(par, databases);
  }
  
  
  //prepare input stream
  size_t ca3m_range_start = 0;
  size_t ca3m_range_end = ca3m_index->n_entries;
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t entry_index = ca3m_range_start; entry_index < ca3m_range_end; entry_index++) {
      
      
      
    
      
  ffindex_entry_t* entry = ffindex_get_entry_by_index(ca3m_index, entry_index);
  char* data = ffindex_get_data_by_entry(ca3m_data, entry);

      /*
      
    
    if(entry == NULL) { perror(entry->name); continue; }

    char* data = ffindex_get_data_by_entry(ca3m_data, entry);

    std::stringstream* out_buffer = new std::stringstream();
    compressed_a3m::extract_a3m(data, entry->length, sequence_index, sequence_data, header_index, header_data, out_buffer);

    std::string out_string = out_buffer->str();
    */

/*
    HMM q(MAXSEQDIS, par.maxres);
    Alignment ali_tmp(MAXSEQ, par.maxres);

    // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
    ali_tmp.ReadCompressed(entry,data,
      sequence_index, seq,
      header_index, header,
      par.mark, par.maxcol);

    // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
    // and store marked sequences in name[k] and seq[k]
    ali_tmp.Compress(par.infile, par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);

    ali_tmp.Shrink();

    // Sort out the nseqdis most dissimilar sequences for display in the output alignments
    ali_tmp.FilterForDisplay(par.max_seqid, par.mark, S, par.coverage, par.qid, par.qsc, par.nseqdis);

    // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
    ali_tmp.N_filtered = ali_tmp.Filter(par.max_seqid, S, par.coverage, par.qid, par.qsc, par.Ndiff);

    if (par.Neff >= 0.999)
    	ali_tmp.FilterNeff(use_global_weights, par.mark, par.cons, par.showcons,
          par.maxres, par.max_seqid, par.coverage, par.Neff, pb, S, Sim);

    // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
    ali_tmp.FrequenciesAndTransitions(q, use_global_weights, par.mark, par.cons,
        par.showcons, par.maxres, pb, Sim);



    input_format = 0;
*/


/*

    if (entry == NULL) {
      HH_LOG(WARNING) << "Could not open entry " << entry_index << " from input ffindex!" << std::endl;
      continue;
    }
*/
    int bin = 0;
#ifdef OPENMP
    bin = omp_get_thread_num();
    omp_set_num_threads(1);
#endif
/*
    FILE* inf = ffindex_fopen_by_entry(data, entry);
    if(inf == NULL) {
      HH_LOG(WARNING) << "Could not open input entry (" << entry->name << ")!" << std::endl;
      continue;
    }
*/
    HH_LOG(INFO) << "Thread " << bin << "\t" << entry->name << std::endl;
    hhblits_instances[bin]->run(entry,data,
      sequence_index, sequence_data,
      header_index, header_data);

    #pragma omp critical
    {
      for (size_t i = 0; i < outputDatabases.size(); i++) {
        outputDatabases[i].saveOutput(*hhblits_instances[bin], entry->name);
      }
    }

    hhblits_instances[bin]->Reset();
  }


  

  for(int i = 0; i < threads; i++) {
    delete hhblits_instances[i];
  }
  delete[] hhblits_instances;

  for (size_t i = 0; i < databases.size(); i++) {
    delete databases[i];
  }
  databases.clear();

  for (size_t i = 0; i < outputDatabases.size(); i++) {
    outputDatabases[i].close();
  }
}

