/*
 * hhprefilter.h
 *
 *  Created on: Apr 17, 2014
 *      Author: meiermark
 */
#ifndef HHPREFILTER_H_
#define HHPREFILTER_H_

#include <sstream>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

#include "hhhmm.h"
#include "hash.h"
#include "hhhit.h"
#include "simd.h"
#include "ffindexdatabase.h"

//////////////////////////////////////////////////////////////////////////////////////////
//   The function swStripedByte contains code adapted from Mengyao Zhao
//   The MIT License
//   Copyright (c) 2012-2015 Boston College.
//   Permission is hereby granted, free of charge, to any person obtaining
//   a copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to
//   permit persons to whom the Software is furnished to do so, subject to
//   the following conditions:
//   The above copyright notice and this permission notice shall be
//   included in all copies or substantial portions of the Software.
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.


// The 2-clause BSD License
//   Copyright 2006 Michael Farrar.  
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//   
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//   
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//   
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


class Prefilter {
public:
	Prefilter(const std::string& cs_library, FFindexDatabase* cs219_database);
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

#endif /* HHPREFILTER_H_ */
