/*
 * hhprefilter.h
 *
 *  Created on: Apr 17, 2014
 *      Author: meiermark
 */

#ifndef HHPREFILTER_H_
#define HHPREFILTER_H_

#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif
namespace hh {
class Prefilter;
}

#include "hhhmm.h"
#include "hash.h"
#include "hhhit.h"
#include "hhdatabase.h"
#include "simd.h"

namespace hh {

const int SHORT_BIAS = 32768;
const int NUMCOLSTATES = cs::AS219::kSize;

//////////////////////////////////////////////////////////////////////////////////////////
//     This file contains code adapted from Michael Farrar
//     (http://sites.google.com/site/farrarmichael/smith-waterman). His code is marked.
//     The copy right of his code is shown below:

//     Copyright 2006, by Michael Farrar.  All rights reserved. The SWSSE2
//     program and documentation may not be sold or incorporated into a
//     commercial product, in whole or in part, without written consent of
//     Michael Farrar.
//
//     For further information regarding permission for use or reproduction,
//     please contact Michael Farrar at:
//
//         farrar.michael@gmail.com
//
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
//     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
//     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
//     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

//     Reference:
//     Farrar M. Striped Smith-Waterman speeds database searches
//     six times over other SIMD implementations. Bioinformatics. 2007, 23:156-61.
//
//     Michael Farrar died unexpectedly in December 2010.
//     Many thanks posthumously for your great code!
//     Johannes


class Prefilter {

public:
	Prefilter(const char* cs_library, FFindexDatabase* cs219_database);
	virtual ~Prefilter();

	static void init_no_prefiltering(FFindexDatabase* cs219_database, std::vector<std::pair<int, std::string> >& prefiltered_entries);
	static void init_selected(FFindexDatabase* cs219_database, std::vector<std::string> templates, std::vector<std::pair<int, std::string> >& prefiltered_entries);

	void prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits,
			const int threads, const int prefilter_gap_open, const int prefilter_gap_extend,
			const int prefilter_score_offset, const int prefilter_bit_factor, const double prefilter_evalue_thresh,
			const double prefilter_evalue_coarse_thresh, const int preprefilter_smax_thresh,
            const int min_prefilter_hits, const int maxnumdb, const float R[20][20],
			std::vector<std::pair<int, std::string> >& new_prefilter_hits, std::vector<std::pair<int, std::string> >& old_prefilter_hits);

private:
	cs::ContextLibrary<cs::AA> *cs_lib;

	// number of sequences in prefilter database file
	size_t num_dbs;

	// array containing all sequence names in prefilter db file
	char** dbnames;

	// pointer to first letter of next sequence in db_data
	unsigned char** first;

	// length of next sequence
	int* length;

	// extended column state query profile as char
//	unsigned char* qc;
//	int W;

	void init_prefilter(FFindexDatabase* cs219_database);

	int ungapped_sse_score(const unsigned char* query_profile,
		const int query_length, const unsigned char* db_sequence,
		const int dbseq_length, const unsigned char score_offset, simd_int* workspace);

	int swStripedByte(unsigned char *querySeq,
		int queryLength,
		unsigned char *dbSeq,
		int dbLength,
		unsigned short gapOpen,
		unsigned short gapExtend,
		simd_int *pvHLoad,
		simd_int *pvHStore,
		simd_int *pvE,
		unsigned short bias);

	void checkCSFormat(size_t nr_checks);
	void stripe_query_profile(HMM* q_tmp, const int prefilter_score_offset, const int prefilter_bit_factor, const int W, unsigned char* qc);
};

} /* namespace hh */

#endif /* HHPREFILTER_H_ */
