/*
 * hhalign_app.cpp
 *
 *  Created on: Jun 16, 2014
 *      Author: meiermark
 */

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <omp.h>

extern "C" {
#include <ffindex.h>
}

#include "hhalign_ori.h"
#include "util.h"

struct OutputFFIndex {
    char base[NAMELEN];
    FILE* data_fh;
    FILE* index_fh;
    size_t offset;
    size_t number_entries;
    void (*print)(HHalign&, std::stringstream&);

    void close() {
      fclose(data_fh);
      fclose(index_fh);
    }

    void saveOutput(HHalign& hhalign, char* name) {
      std::stringstream out;
      print(hhalign, out);

      std::string tmp = out.str();
      ffindex_insert_memory(data_fh, index_fh, &offset,
          const_cast<char*>(tmp.c_str()), tmp.size(), name);

      fflush(data_fh);
      fflush(index_fh);
      number_entries++;
    }

    void sort() {
      /* Sort the index entries and write back */
      rewind(index_fh);
      ffindex_index_t* index = ffindex_index_parse(index_fh, number_entries);
      if (index == NULL) {
        //TODO: throw error
      }
      fclose(index_fh);

      ffindex_sort_index_file(index);

      char index_filename[NAMELEN];
      snprintf(index_filename, FILENAME_MAX, "%s.ffindex", base);
      index_fh = fopen(index_filename, "w");

      if (index_fh == NULL) {
        //TODO: throw error
      }

      ffindex_write(index, index_fh);
    }
};


void makeOutputFFIndex(char* par, void (*print)(HHalign&, std::stringstream&),
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
      std::cerr << "could not open datafile " << data_filename_out_rank
          << std::endl;
      return;
    }

    if (db.index_fh == NULL) {
      std::cerr << "could not open indexfile " << index_filename_out_rank
          << std::endl;
      return;
    }

    outDatabases.push_back(db);
  }
}


void readQueriesToTemplates(std::string queries_to_templates_file,
		std::map<std::string, std::vector<std::string>>& map) {
	std::ifstream infile(queries_to_templates_file);

	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#') {
			continue;
		}

		std::vector<std::string> tokens;

		split(line, '\t', tokens);
		if (tokens.size() == 0) {
			continue;
		}

		std::string query = tokens[0];
		tokens.erase(tokens.begin());

		if (map.find(query) != map.end()) {
			(*map.find(query)).second.insert((*map.find(query)).second.end(),
					tokens.begin(), tokens.end());
		} else {
			map[query] = tokens;
		}
	}

	infile.close();
}

int main(int argc, char **argv) {
	Parameters par;
	HHalign::ProcessAllArguments(argc, argv, par);

	char data_filename[NAMELEN];
	char index_filename[NAMELEN];

	strcpy(data_filename, par.infile);
	strcat(data_filename, ".ffdata");

	strcpy(index_filename, par.infile);
	strcat(index_filename, ".ffindex");

	FILE *data_file = fopen(data_filename, "r");
	FILE *index_file = fopen(index_filename, "r");

	if (data_file == NULL) {
		std::cerr << "input data file " << data_filename << " does not exist!"
				<< std::endl;
		exit(EXIT_FAILURE);
	}
	if (index_file == NULL) {
		std::cerr << "input index file " << index_filename << " does not exist!"
				<< std::endl;
		exit(EXIT_FAILURE);
	}

	//init input ffindex
	size_t data_size;
	char *data = ffindex_mmap_data(data_file, &data_size);

	size_t number_input_index_lines = CountLinesInFile(index_filename);
	ffindex_index_t* index = ffindex_index_parse(index_file,
			number_input_index_lines);
	if (index == NULL) {
		std::cerr << "Could not parse index from " << index_filename
				<< std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<OutputFFIndex> outputDatabases;
	makeOutputFFIndex(par.outfile, &HHalign::writeHHRFile, outputDatabases);
//	makeOutputFFIndex(par.scorefile, &HHblits::writeScoresFile, outputDatabases);
//	makeOutputFFIndex(par.pairwisealisfile, &HHblits::writePairwiseAlisFile, outputDatabases);
//	makeOutputFFIndex(par.alitabfile, &HHblits::writeAlitabFile, outputDatabases);
//	makeOutputFFIndex(par.reduced_outfile, &HHblits::writeReducedHHRFile, outputDatabases);
//	makeOutputFFIndex(par.psifile, &HHblits::writePsiFile, outputDatabases);
//	makeOutputFFIndex(par.hhmfile, &HHblits::writeHMMFile, outputDatabases);
//	makeOutputFFIndex(par.alnfile, &HHblits::writeA3MFile, outputDatabases);

	std::vector<HHblitsDatabase*> databases;
	HHblits::prepareDatabases(par, databases);

	std::map<std::string, std::vector<std::string>> queries_to_templates;
	readQueriesToTemplates(std::string(par.queries_to_template_file), queries_to_templates);

	int threads = par.threads;
	omp_set_num_threads(threads);
	par.threads = 1;

	HHalign* hhalign_instances[255];
	for (int i = 0; i < threads; i++) {
		hhalign_instances[i] = new HHalign(par, databases);
	}

//	#pragma omp parallel for
	for (std::map<std::string, std::vector<std::string>>::iterator iterator = queries_to_templates.begin();
			iterator != queries_to_templates.end(); iterator++) {
		std::string query = (*iterator).first;
		std::vector<std::string> templates = (*iterator).second;

		ffindex_entry_t* entry = ffindex_get_entry_by_name(index, const_cast<char*>(query.c_str()));
		if (entry == NULL) {
			std::cerr << "Could not open entry " << query << " from input ffindex!" << std::endl;
			continue;
		}

		int bin = omp_get_thread_num();
		omp_set_num_threads(1);

		hhalign_instances[bin]->Reset();

		FILE* inf = ffindex_fopen_by_entry(data, entry);
		hhalign_instances[bin]->run(inf, entry->name, templates);
		fclose(inf);

		#pragma omp critical
		{
			for (size_t i = 0; i < outputDatabases.size(); i++) {
				outputDatabases[i].saveOutput(*hhalign_instances[bin], entry->name);
			}
		}
	}

	fclose(data_file);
	fclose(index_file);

	for (size_t i = 0; i < outputDatabases.size(); i++) {
		outputDatabases[i].sort();
		outputDatabases[i].close();
	}
}
