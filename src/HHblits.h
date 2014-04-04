/*
 * HHblits.h
 *
 *  Created on: Apr 1, 2014
 *      Author: meiermark
 */

#ifndef HHBLITS_H_
#define HHBLITS_H_

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <sstream>
#include <vector>
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror(), strerror(errno)
#include <cassert>
#include <stdexcept>
#include <map>
#include <omp.h>

#include <emmintrin.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::vector;
using std::pair;

#ifdef PTHREAD
#include <pthread.h>  // POSIX pthread functions and data structures
#endif

#include <sys/time.h>

extern "C" {
#include <ffindex.h>
}

#include "cs.h"
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"
#include "abstract_state_matrix.h"

#include "hhdecl.h"
#include "list.h"
#include "hash.h"
#include "util.h"
#include "hhutil.h"

#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfullalignment.h" // class FullAlignment
#include "hhhitlist.h"   // class Hit

#include "hhmatrices.h"
#include "hhfunc.h"      // some functions common to hh programs

const int MAXBINS = 384; // maximum number of bins (positions in thread queue)
const char HHBLITS_REFERENCE[] =
		"Remmert M., Biegert A., Hauser A., and Soding J.\nHHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.\nNat. Methods 9:173-175 (2011)\n";

class HHblits {
public:
	HHblits(int argc, char** argv);
	virtual ~HHblits();

	void Init(int argc, char** argv);
	void Reset();

	void SetDatabase(char* db_base);

    //writer for non-mpi version
    void writeHHRFile(char* hhrFile);
    void writeAlisFile(char* basename);
    void writeScoresFile(char* scoresFile);
    void writePairwiseAlisFile(char* pairwieseAlisFile, char outformat);
    void writeAlitabFile(char* alitabFile);
    void writeReducedHHRFile(char* reduceHHRFile);
    void writePsiFile(char* psiFile);
    void writeHMMFile(char* HMMFile);
    void writeA3MFile(char* A3MFile);

    //output writer for mpi version
    std::map<int, Alignment>& getAlis();
    std::stringstream* writeHHRFile();
    std::stringstream* writeScoresFile();
    std::stringstream* writePairwiseAlisFile(char outformat);
    std::stringstream* writeAlitabFile();
    std::stringstream* writeReducedHHRFile();
    std::stringstream* writePsiFile();
    std::stringstream* writeHMMFile();
    std::stringstream* writeA3MFile();

	void run(FILE* query_fh, char* query_path);

private:
	const int SHORT_BIAS = 32768;
	const int NUMCOLSTATES = cs::AS219::kSize;

	FILE* db_data_file;
	unsigned char* db_data;

	// number of sequences in prefilter database file
	size_t num_dbs = 0;
	// array containing all sequence names in prefilter db file
	char** dbnames;
	// pointer to first letter of next sequence in db_data
	unsigned char** first;
	// length of next sequence
	int* length;

	unsigned char* qc;     // extended column state query profile as char
	int W;                 //
	unsigned short* qw;    // extended column state query profile as short int
	int Ww;

	cs::ContextLibrary<cs::AA> *cs_lib;

	// HHblits variables
	int v1 = v;                               // verbose mode
	int num_rounds = 2;                   // number of iterations
	bool last_round = false;                // set to true in last iteration
	bool already_seen_filter = true;   // Perform filtering of already seen HHMs
	bool realign_old_hits = false; // Realign old hits in last round or use previous alignments

	char input_format = 0; // Set to 1 if input in HMMER format (has already pseudocounts)

	float neffmax = 10;                     // Break if Neff > Neffmax

	char config_file[NAMELEN];
	char infile[NAMELEN];
	char query_hhmfile[NAMELEN];             // -qhmm output file
	bool alitab_scop = false;        // Write only SCOP alignments in alitabfile

	//database filenames
	char db_base[NAMELEN];
	char dbcs_base[NAMELEN];
	char dbcs_index_filename[NAMELEN];
	char dbcs_data_filename[NAMELEN];

	char dbhhm_base[NAMELEN];
	char dbhhm_index_filename[NAMELEN];
	char dbhhm_data_filename[NAMELEN];

	char dba3m_base[NAMELEN];
	char dba3m_index_filename[NAMELEN];
	char dba3m_data_filename[NAMELEN];

	//compressed a3m stuff
	bool use_compressed_a3m;
	char dbca3m_base[NAMELEN];
	char dbca3m_index_filename[NAMELEN];
	char dbca3m_data_filename[NAMELEN];

	char dbuniprot_base[NAMELEN];
	char dbuniprot_header_index_filename[NAMELEN];
	char dbuniprot_header_data_filename[NAMELEN];
	char dbuniprot_sequence_index_filename[NAMELEN];
	char dbuniprot_sequence_data_filename[NAMELEN];

	// Needed for fast index reading
	size_t data_size;
	FILE *dba3m_data_file;
	FILE *dba3m_index_file;
	ffindex_index_t* dba3m_index = NULL;
	char* dba3m_data;

	FILE *dbhhm_data_file;
	FILE *dbhhm_index_file;
	ffindex_index_t* dbhhm_index = NULL;
	char* dbhhm_data;

	// compressed a3m stuff
	size_t ca3m_data_offset;
	FILE* dbca3m_data_file;
	FILE* dbca3m_index_file;
	ffindex_index_t* dbca3m_index = NULL;
	char* dbca3m_data;

	size_t uniprot_header_data_offset;
	FILE* dbuniprot_header_data_file;
	FILE* dbuniprot_header_index_file;
	ffindex_index_t* dbuniprot_header_index = NULL;
	char* dbuniprot_header_data;

	size_t uniprot_sequence_data_offset;
	FILE* dbuniprot_sequence_data_file;
	FILE* dbuniprot_sequence_index_file;
	ffindex_index_t* dbuniprot_sequence_index = NULL;
	char* dbuniprot_sequence_data;

	char** dbfiles_new;
	char** dbfiles_old;
	int ndb_new = 0;
	int ndb_old = 0;

	// HHsearch variables
	int threads = 2; // number of compute pthreads during Viterbi and Realign (apart from the main thread which reads from db file) and # OpenMP threads; 0:no multithreading
    int omp_threads = 2;                    // number of OpenMP threads to start
	int bins; // number of bins; jobs gets allocated to a FREE bin were they are waiting for execution by a thread

	// Create query HMM with maximum of par.maxres match states
	HMM* q = NULL;
	 // Create query HMM with maximum of par.maxres match states (needed for prefiltering)
	HMM* q_tmp = NULL;

	Alignment Qali; // output A3M generated by merging A3M alignments for significant hits to the query alignment
	Alignment Qali_allseqs; // output A3M alignment with no sequence filtered out (only active with -all option)

	HMM* t[MAXBINS]; // Each bin has a template HMM allocated that was read from the database file
	Hit* hit[MAXBINS]; // Each bin has an object of type Hit allocated with a separate dynamic programming matrix (memory!!)
	int* format; // format[bin] = 0 if in HHsearch format => add pcs; format[bin] = 1 if in HMMER format => no pcs

	HitList hitlist; // list of hits with one Hit object for each pairwise comparison done
	HitList reducedHitlist;
	int N_searched;           // Number of HMMs searched
	std::map<int, Alignment> alis;

	void help(char all = 0);
	void ProcessArguments(int argc, char** argv);

	void ReadQueryA3MFile();
	void getTemplateA3M(char* entry_name, long& ftellpos, Alignment& tali);
	void getTemplateHMMFromA3M(char* entry_name, char use_global_weights, long& ftellpos, int& format, HMM* t);
	void getTemplateHMM(char* entry_name, char use_global_weights, long& ftellpos, int& format, HMM* t);
	void DoViterbiSearch(char *dbfiles[], int ndb, Hash<Hit>* previous_hits, bool alignByWorker = true);
	void ViterbiSearch(char *dbfiles[], int ndb, Hash<Hit>* previous_hits, int db_size);

	void RescoreWithViterbiKeepAlignment(int db_size, Hash<Hit>* previous_hits);
	void perform_realign(char *dbfiles[], int ndb, Hash<char>* premerged_hits);

    void wiggleQSC(HitList& hitlist, int n_redundancy, Alignment& Qali,
            char inputformat, float* qsc, size_t nqsc,
            HitList& reducedFinalHitList);
    void wiggleQSC(int n_redundancy, float* qsc, size_t nqsc, HitList& reducedFinalHitList);
	void reduceRedundancyOfHitList(int n_redundancy, int query_length,
			HitList& hitlist, HitList& reducedHitList);
	void recalculateAlignmentsForDifferentQSC(HitList& hitlist, Alignment& Qali,
			char inputformat, float* qsc, size_t nqsc,
			HitList& recalculated_hitlist);
	void uniqueHitlist(HitList& hitlist);

	//Worker
	void AlignByWorker(int bin);
	void PerformViterbiByWorker(int bin, Hash<Hit>* previous_hits);
	void RealignByWorker(int bin);

	//Prefilter
	int ungapped_sse_score(const unsigned char* query_profile,
			const int query_length, const unsigned char* db_sequence,
			const int dbseq_length, const unsigned char score_offset, __m128i* workspace);

			int swStripedByte(unsigned char *querySeq,
			int queryLength,
			unsigned char *dbSeq,
			int dbLength,
			unsigned short gapOpen,
			unsigned short gapExtend,
			__m128i *pvHLoad,
			__m128i *pvHStore,
			__m128i *pvE,
			unsigned short bias);

			void init_no_prefiltering();
			void checkCSFormat(size_t nr_checks);
			void init_prefilter();
			void prefilter_db(Hash<Hit>* previous_hits);
			void stripe_query_profile();
};

#endif /* HHBLITS_H_ */
