/*
 * hhblits_mpi.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#include "hhsearch.h"
#include "hhalign.h"

#ifdef OPENMP
#include <omp.h>
#endif

struct OutputFFIndex {
    char base[NAMELEN];
    FILE *data_fh;
    FILE *index_fh;
    size_t offset;
    size_t number_entries;

    void (*print)(HHblits &, std::stringstream &);

    void close() {
        char index_filename[NAMELEN];
        snprintf(index_filename, FILENAME_MAX, "%s.ffindex", base);

        fclose(index_fh);
        fclose(data_fh);

        ffsort_index(index_filename);
    }

    void saveOutput(HHblits &hhblits, char *name) {
        std::stringstream out;
        print(hhblits, out);

        std::string tmp = out.str();
        ffindex_insert_memory(data_fh, index_fh, &offset, const_cast<char *>(tmp.c_str()), tmp.size(), name);

        fflush(data_fh);
        fflush(index_fh);
        number_entries++;
    }
};


void makeOutputFFIndex(char *par, void (*print)(HHblits &, std::stringstream &), std::vector<OutputFFIndex> &outDatabases) {
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
#ifdef HHSEARCH
    HHsearch::ProcessAllArguments(par);
#elif HHALIGN
    HHalign::ProcessAllArguments(par);
#else
    HHblits::ProcessAllArguments(par);
#endif

    std::string data_filename(par.infile);
    data_filename.append(".ffdata");

    std::string index_filename(par.infile);
    index_filename.append(".ffindex");

    FFindexDatabase reader(data_filename.c_str(), index_filename.c_str(), false);
    reader.ensureLinearAccess();

    std::vector<HHblitsDatabase*> databases;
#ifdef HHSEARCH
    HHsearch::prepareDatabases(par, databases);
#elif HHALIGN
#else
    HHblits::prepareDatabases(par, databases);
#endif

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

    // only parallelize over queries, not per query
    int threads = par.threads;
    par.threads = 1;

#pragma omp parallel num_threads(threads)
    {
#ifdef HHSEARCH
        HHblits app(par, databases);
#elif HHALIGN
        HHalign app(par);
#else
        HHblits app(par, databases);
#endif

        int bin = 0;
#ifdef OPENMP
        bin = omp_get_thread_num();
        omp_set_num_threads(1);
#endif

#pragma omp for schedule(dynamic, 1)
        for (size_t entry_index = 0; entry_index < reader.db_index->n_entries; entry_index++) {
            ffindex_entry_t *entry = ffindex_get_entry_by_index(reader.db_index, entry_index);
            if (entry == NULL) {
                HH_LOG(WARNING) << "Could not open entry " << entry_index << " from input ffindex!" << std::endl;
                continue;
            }

            FILE *inf = ffindex_fopen_by_entry(reader.db_data, entry);
            if (inf == NULL) {
                HH_LOG(WARNING) << "Could not open input entry (" << entry->name << ")!" << std::endl;
                continue;
            }

            HH_LOG(INFO) << "Thread " << bin << "\t" << entry->name << std::endl;
            app.run(inf, entry->name);

#pragma omp critical
            {
                for (size_t i = 0; i < outputDatabases.size(); ++i) {
                    outputDatabases[i].saveOutput(app, entry->name);
                }
            }

            app.Reset();
        }
    }

    for (size_t i = 0; i < outputDatabases.size(); ++i) {
        outputDatabases[i].close();
    }

    for (size_t i = 0; i < databases.size(); ++i) {
        delete databases[i];
    }
    databases.clear();
}

