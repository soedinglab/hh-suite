/*
 * hhblits_mpi.cpp
 *
 *  Created on: Apr 1, 2014
 *      Author: Markus Meier (markus.meier@mpibpc.mpg.de)
 */

#include "hhsearch.h"
#include "hhalign.h"

#include <mpi.h>

extern "C" {
#include <ffindex.h>
#include <mpq/mpq.h>
}

#ifdef OPENMP
#include <omp.h>
#endif

struct OutputFFIndex {
    char base[NAMELEN];
    FILE *data_fh;
    FILE *index_fh;
    size_t offset;

    void (*print)(HHblits &, std::stringstream &);

    void close() {
        fclose(data_fh);
        fclose(index_fh);
    }

    void saveOutput(HHblits &hhblits, char *name) {
        std::stringstream out;
        print(hhblits, out);

        std::string tmp = out.str();
        ffindex_insert_memory(data_fh, index_fh, &offset,
                              const_cast<char *>(tmp.c_str()), tmp.size(), name);

        fflush(data_fh);
        fflush(index_fh);
    }
};

void makeOutputFFIndex(char *par, const int mpi_rank,
                       void (*print)(HHblits &, std::stringstream &),
                       std::vector<OutputFFIndex> &outDatabases) {
    if (*par) {
        OutputFFIndex db;

        strcpy(db.base, par);
        db.offset = 0;
        db.print = print;

        char data_filename_out_rank[NAMELEN];
        char index_filename_out_rank[NAMELEN];

        snprintf(data_filename_out_rank, FILENAME_MAX, "%s.ffdata.%d", par,
                 mpi_rank);
        snprintf(index_filename_out_rank, FILENAME_MAX, "%s.ffindex.%d", par,
                 mpi_rank);

        db.data_fh = fopen(data_filename_out_rank, "w+");
        db.index_fh = fopen(index_filename_out_rank, "w+");

        if (db.data_fh == NULL) {
            HH_LOG(WARNING) << "Could not open datafile " << data_filename_out_rank << std::endl;
            return;
        }

        if (db.index_fh == NULL) {
            HH_LOG(WARNING) << "Could not open indexfile " << index_filename_out_rank << std::endl;
            return;
        }

        outDatabases.push_back(db);
    }
}

void merge_splits(const char *prefix) {
    if (*prefix) {
        char data_filename[FILENAME_MAX];
        char index_filename[FILENAME_MAX];

        snprintf(data_filename, FILENAME_MAX, "%s.ffdata", prefix);
        snprintf(index_filename, FILENAME_MAX, "%s.ffindex", prefix);

        ffmerge_splits(data_filename, index_filename, 1, MPQ_size - 1, true);
    }
}

struct HHblits_MPQ_Wrapper {
    char *data;
    ffindex_index_t *index;
    HHblits *hhblits;
    std::vector<OutputFFIndex> *outputDatabases;

    HHblits_MPQ_Wrapper(char *data, ffindex_index_t *index, HHblits &hhblits,
                        std::vector<OutputFFIndex> &outputDatabases) {
        this->data = data;
        this->index = index;
        this->hhblits = &hhblits;
        this->outputDatabases = &outputDatabases;
    }

    void Payload(const size_t start, const size_t end) {
        // Foreach entry in the input file
        for (size_t entry_index = start; entry_index < end; entry_index++) {
            ffindex_entry_t *entry = ffindex_get_entry_by_index(index, entry_index);
            if (entry == NULL) {
                continue;
            }

            hhblits->Reset();

            FILE *inf = ffindex_fopen_by_entry(data, entry);
            hhblits->run(inf, entry->name);
            fclose(inf);

            for (size_t i = 0; i < outputDatabases->size(); i++) {
                outputDatabases->operator[](i).saveOutput(*hhblits, entry->name);
            }
        }
    }
};

void static payload(void *env, const size_t start, const size_t end) {
    HHblits_MPQ_Wrapper *hhblits_wrapper = (HHblits_MPQ_Wrapper *) env;
    hhblits_wrapper->Payload(start, end);
}

int main(int argc, char **argv) {
    Parameters par;
#ifdef HHSEARCH
    HHsearch::ProcessAllArguments(argc, argv, par);
#elif HHALIGN
    HHalign::ProcessAllArguments(argc, argv, par);
#else
    HHblits::ProcessAllArguments(argc, argv, par);
#endif

    // hhblits_mpi will be parallelized with openmpi, no other parallelization
    par.threads = 1;
#ifdef OPENMP
    omp_set_num_threads(1);
#endif

    std::string data_filename(par.infile);
    data_filename.append(".ffdata");

    std::string index_filename(par.infile);
    index_filename.append(".ffindex");

    FFindexDatabase reader(data_filename.c_str(), index_filename.c_str(), false);
    reader.ensureLinearAccess();

    int mpq_status = MPQ_Init(argc, argv, reader.db_index->n_entries);

    if (mpq_status == MPQ_SUCCESS) {
        if (MPQ_rank == MPQ_MASTER) {
            MPQ_Master(1);
        } else {
            std::vector<OutputFFIndex> outputDatabases;
            makeOutputFFIndex(par.outfile, MPQ_rank, &HHblits::writeHHRFile,
                              outputDatabases);
            makeOutputFFIndex(par.scorefile, MPQ_rank, &HHblits::writeScoresFile,
                              outputDatabases);
            makeOutputFFIndex(par.pairwisealisfile, MPQ_rank,
                              &HHblits::writePairwiseAlisFile, outputDatabases);
            makeOutputFFIndex(par.alitabfile, MPQ_rank, &HHblits::writeAlitabFile,
                              outputDatabases);
            makeOutputFFIndex(par.psifile, MPQ_rank, &HHblits::writePsiFile,
                              outputDatabases);
            makeOutputFFIndex(par.hhmfile, MPQ_rank, &HHblits::writeHMMFile,
                              outputDatabases);
            makeOutputFFIndex(par.alnfile, MPQ_rank, &HHblits::writeA3MFile,
                              outputDatabases);
            makeOutputFFIndex(par.matrices_output_file, MPQ_rank, &HHblits::writeMatricesFile,
                              outputDatabases);
            makeOutputFFIndex(par.m8file, MPQ_rank, &HHblits::writeM8,
                              outputDatabases);

            std::vector<HHblitsDatabase*> databases;
#ifdef HHSEARCH
            HHsearch::prepareDatabases(par, databases);
            HHblits app(par, databases);
#elif HHalign
            HHalign app(par);
#else
            HHblits::prepareDatabases(par, databases);
            HHblits app(par, databases);
#endif

            HHblits_MPQ_Wrapper wrapper(reader.db_data, reader.db_index, app, outputDatabases);
            MPQ_Worker(payload, &wrapper);

            for (size_t i = 0; i < outputDatabases.size(); i++) {
                outputDatabases[i].close();
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (MPQ_rank == MPQ_MASTER) {
            merge_splits(par.outfile);
            merge_splits(par.scorefile);
            merge_splits(par.pairwisealisfile);
            merge_splits(par.alitabfile);
            merge_splits(par.psifile);
            merge_splits(par.hhmfile);
            merge_splits(par.alnfile);
            merge_splits(par.matrices_output_file);
            merge_splits(par.m8file);
        }
    } else {
        if (mpq_status == MPQ_ERROR_NO_WORKERS) {
            fprintf(stderr, "MPQ_Init: Needs at least one worker process.\n");
            exit(EXIT_FAILURE);
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

