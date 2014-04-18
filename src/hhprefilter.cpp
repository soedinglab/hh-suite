/*
 * hhprefilter.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: meiermark
 */

#include "hhprefilter.h"

namespace hh {

#define SWAP(tmp, arg1, arg2) tmp = arg1; arg1 = arg2; arg2 = tmp;

//TODO: read cs_library once and give pointer to it

Prefilter::Prefilter(const char* cs_library, FFindexDatabase* cs219_database) {
	// Prepare column state lib (context size =1 )
	FILE* fin = fopen(cs_library, "r");
	if (!fin)
		OpenFileError(cs_library, __FILE__, __LINE__, __func__);

	cs_lib = new cs::ContextLibrary<cs::AA>(fin);
	fclose(fin);

	cs::TransformToLin(*cs_lib);

	init_prefilter(cs219_database);
}

Prefilter::~Prefilter() {
	if (use_prefilter) {
		free(length);
		free(first);

		for (size_t n = 0; n < num_dbs; n++)
			delete[] dbnames[n];
		free(dbnames);
	}

	delete cs_lib;
}

int Prefilter::swStripedByte(unsigned char *querySeq, int queryLength,
		unsigned char *dbSeq, int dbLength, unsigned short gapOpen,
		unsigned short gapExtend, __m128i *pvHLoad, __m128i *pvHStore,
		__m128i *pvE, unsigned short bias) {
	int i, j;
	int score;

	int cmp;
	int iter = (queryLength + 15) / 16;

	__m128i *pv;

	__m128i vE, vF, vH;

	__m128i vMaxScore;
	__m128i vBias;
	__m128i vGapOpen;
	__m128i vGapExtend;

	__m128i vTemp;
	__m128i vZero;

	__m128i *pvScore;

	__m128i *pvQueryProf = (__m128i *) querySeq;

	/* Load the bias to all elements of a constant */
	vBias = _mm_set1_epi8(bias);

	/* Load gap opening penalty to all elements of a constant */
	vGapOpen = _mm_set1_epi8(gapOpen);

	/* Load gap extension penalty to all elements of a constant */
	vGapExtend = _mm_set1_epi8(gapExtend);

	vMaxScore = _mm_setzero_si128();
	vZero = _mm_setzero_si128();

	/* Zero out the storage vector */
	for (i = 0; i < iter; ++i) {
		_mm_store_si128(pvE + i, vMaxScore);
		_mm_store_si128(pvHStore + i, vMaxScore);
	}

	for (i = 0; i < dbLength; ++i) {
		/* fetch first data asap. */
		pvScore = pvQueryProf + dbSeq[i] * iter;

		/* zero out F. */
		vF = _mm_setzero_si128();

		/* load the next h value */
		vH = _mm_load_si128(pvHStore + iter - 1);
		vH = _mm_slli_si128(vH, 1);

		pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		for (j = 0; j < iter; ++j) {
			/* load values of vF and vH from previous row (one unit up) */
			vE = _mm_load_si128(pvE + j);

			/* add score to vH */
			vH = _mm_adds_epu8(vH, *(pvScore++));
			vH = _mm_subs_epu8(vH, vBias);

			/* Update highest score encountered this far */
			vMaxScore = _mm_max_epu8(vMaxScore, vH);

			/* get max from vH, vE and vF */
			vH = _mm_max_epu8(vH, vE);
			vH = _mm_max_epu8(vH, vF);

			/* save vH values */
			_mm_store_si128(pvHStore + j, vH);

			/* update vE value */
			vH = _mm_subs_epu8(vH, vGapOpen);
			vE = _mm_subs_epu8(vE, vGapExtend);
			vE = _mm_max_epu8(vE, vH);

			/* update vF value */
			vF = _mm_subs_epu8(vF, vGapExtend);
			vF = _mm_max_epu8(vF, vH);

			/* save vE values */
			_mm_store_si128(pvE + j, vE);

			/* load the next h value */
			vH = _mm_load_si128(pvHLoad + j);
		}

		/* reset pointers to the start of the saved data */
		j = 0;
		vH = _mm_load_si128(pvHStore);

		/*  the computed vF value is for the given column.  since */
		/*  we are at the end, we need to shift the vF value over */
		/*  to the next column. */
		vF = _mm_slli_si128(vF, 1);
		vTemp = _mm_subs_epu8(vH, vGapOpen);
		vTemp = _mm_subs_epu8(vF, vTemp);
		vTemp = _mm_cmpeq_epi8(vTemp, vZero);
		cmp = _mm_movemask_epi8(vTemp);

		while (cmp != 0xffff) {
			vE = _mm_load_si128(pvE + j);

			vH = _mm_max_epu8(vH, vF);

			/* save vH values */
			_mm_store_si128(pvHStore + j, vH);

			/*  update vE incase the new vH value would change it */
			vH = _mm_subs_epu8(vH, vGapOpen);
			vE = _mm_max_epu8(vE, vH);
			_mm_store_si128(pvE + j, vE);

			/* update vF value */
			vF = _mm_subs_epu8(vF, vGapExtend);

			++j;
			if (j >= iter) {
				j = 0;
				vF = _mm_slli_si128(vF, 1);
			}

			vH = _mm_load_si128(pvHStore + j);

			vTemp = _mm_subs_epu8(vH, vGapOpen);
			vTemp = _mm_subs_epu8(vF, vTemp);
			vTemp = _mm_cmpeq_epi8(vTemp, vZero);
			cmp = _mm_movemask_epi8(vTemp);
		}
	}

	/* find largest score in the vMaxScore vector */
	vTemp = _mm_srli_si128(vMaxScore, 8);
	vMaxScore = _mm_max_epu8(vMaxScore, vTemp);
	vTemp = _mm_srli_si128(vMaxScore, 4);
	vMaxScore = _mm_max_epu8(vMaxScore, vTemp);
	vTemp = _mm_srli_si128(vMaxScore, 2);
	vMaxScore = _mm_max_epu8(vMaxScore, vTemp);
	vTemp = _mm_srli_si128(vMaxScore, 1);
	vMaxScore = _mm_max_epu8(vMaxScore, vTemp);

	/* store in temporary variable */
	score = _mm_extract_epi16(vMaxScore, 0);
	score = score & 0x00ff;

	/* return largest score */
	return score;
}

// d = i-j+LT-1 is index of diagonal
int Prefilter::ungapped_sse_score(const unsigned char* query_profile,
		const int query_length, const unsigned char* db_sequence,
		const int dbseq_length, const unsigned char score_offset, __m128i* workspace)
{
	int i; // position in query bands (0,..,W-1)
	int j;// position in db sequence (0,..,dbseq_length-1)
	int W = (query_length + 15) / 16;// width of bands in query and score matrix = hochgerundetes LQ/16

	__m128i *p;
	__m128i S;// 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15)
	__m128i Smax = _mm_setzero_si128();
	__m128i Soffset;// all scores in query profile are shifted up by Soffset to obtain pos values
	__m128i *s_prev, *s_curr;// pointers to Score(i-1,j-1) and Score(i,j), resp.
	__m128i *qji;// query profile score in row j (for residue x_j)
	__m128i *s_prev_it, *s_curr_it;
	__m128i *query_profile_it = (__m128i *) query_profile;
	__m128i Zero = _mm_setzero_si128();

	// Load the score offset to all 16 unsigned byte elements of Soffset
	Soffset = _mm_set1_epi8(score_offset);

	// Initialize  workspace to zero
	for (i=0, p=workspace; i < 2*W; ++i)
	_mm_store_si128(p++, Zero);

	s_curr = workspace;
	s_prev = workspace + W;

	for (j=0; j<dbseq_length; ++j)// loop over db sequence positions
	{
		// Get address of query scores for row j
		qji = query_profile_it + db_sequence[j]*W;

		// Load the next S value
		S = _mm_load_si128(s_curr + W - 1);
		S = _mm_slli_si128(S, 1);

		// Swap s_prev and s_curr, smax_prev and smax_curr
		SWAP(p,s_prev,s_curr);

		s_curr_it = s_curr;
		s_prev_it = s_prev;

		for (i=0; i<W; ++i)// loop over query band positions
		{
			// Saturated addition and subtraction to score S(i,j)
			S = _mm_adds_epu8(S, *(qji++));// S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset)
			S = _mm_subs_epu8(S, Soffset);// S(i,j) = max(0, S(i,j) - Soffset)
			_mm_store_si128(s_curr_it++, S);// store S to s_curr[i]
			Smax = _mm_max_epu8(Smax, S);// Smax(i,j) = max(Smax(i,j), S(i,j))

			// Load the next S and Smax values
			S = _mm_load_si128(s_prev_it++);
		}
	}

	/* find largest score in the Smax vector */
	S = _mm_srli_si128 (Smax, 8);
	Smax = _mm_max_epu8 (Smax, S);
	S = _mm_srli_si128 (Smax, 4);
	Smax = _mm_max_epu8 (Smax, S);
	S = _mm_srli_si128 (Smax, 2);
	Smax = _mm_max_epu8 (Smax, S);
	S = _mm_srli_si128 (Smax, 1);
	Smax = _mm_max_epu8 (Smax, S);

	/* store in temporary variable */
	int score = _mm_extract_epi16 (Smax, 0);
	score = score & 0x00ff;

	/* return largest score */
	return score;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Pull out all names from prefilter db file and copy into dbfiles_new for full HMM-HMM comparison
///////////////////////////////////////////////////////////////////////////////////////////////////
void Prefilter::init_no_prefiltering(FFindexDatabase* cs219_database, std::vector<std::string>& prefiltered_entries) {
	ffindex_index_t* db_index = cs219_database->db_index;

	for (size_t n = 0; n < db_index->n_entries; n++) {
		ffindex_entry_t* entry = ffindex_get_entry_by_index(db_index, n);
		prefiltered_entries.push_back(std::string(entry->name));
	}

	if (v >= 2)
		std::cout << "Searching " << prefiltered_entries.size() << " database HHMs without prefiltering" << std::endl;
}


//////////////////////////////////////////////////////////////
// Reading in column state sequences for prefiltering
//////////////////////////////////////////////////////////////
void Prefilter::init_prefilter(FFindexDatabase* cs219_database) {
	// Set up variables for prefiltering
	num_dbs = cs219_database->db_index->n_entries;
	first = (unsigned char**) memalign(16, num_dbs * sizeof(unsigned char*));
	length = (int*) memalign(16, num_dbs * sizeof(int));
	dbnames = (char**) memalign(16, num_dbs * sizeof(char*));
	for (size_t n = 0; n < num_dbs; n++) {
		ffindex_entry_t* entry = ffindex_get_entry_by_index(cs219_database->db_index, n);
		first[n] = (unsigned char*) ffindex_get_data_by_entry(cs219_database->db_data, entry);
		length[n] = entry->length - 1;
		dbnames[n] = new char[strlen(entry->name) + 1];
		strcpy(dbnames[n], entry->name);
	}

	//check if cs219 format is new binary format
	checkCSFormat(5);

	if (v >= 2) {
		printf("Searching %zu column state sequences.\n", num_dbs);
	}
}


void Prefilter::checkCSFormat(size_t nr_checks) {
	for (size_t n = 0; n < std::min(nr_checks, num_dbs); n++) {
		if (first[n][0] == '>') {
			nr_checks--;
		}
	}

	if (nr_checks == 0) {
		std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
				<< __func__ << ":" << std::endl;
		std::cerr << "\tYour cs database is in an old format!" << std::endl;
		std::cerr << "\tThis format is no longer supportet!" << std::endl;
		std::cerr << "\tCorrespond to the user manual!" << std::endl;
		exit(1);
	}
}


////////////////////////////////////////////////////////////////////////
// Prepare query profile for prefitering
////////////////////////////////////////////////////////////////////////
void Prefilter::stripe_query_profile(HMM* q_tmp, const int prefilter_score_offset, const int prefilter_bit_factor) {
	int LQ = q_tmp->L;
	int a, h, i, j, k;

	// Build query profile with 219 column states
	float** query_profile = new float*[LQ + 1];
	for (i = 0; i < LQ + 1; ++i)
		query_profile[i] = (float*) memalign(16, NUMCOLSTATES * sizeof(float),
				"the query profile during prefiltering");

	const cs::ContextLibrary<cs::AA>& lib = *cs_lib;

	// log S(i,k) = log ( SUM_a p(i,a) * p(k,a) / f(a) )   k: column state, i: pos in ali, a: amino acid
	for (i = 0; i < LQ; ++i)
		for (k = 0; k < NUMCOLSTATES; ++k) {
			float sum = 0;
			for (a = 0; a < 20; ++a)
				sum += (q_tmp->p[i][a] * lib[k].probs[0][a]) / q_tmp->pav[a];
			query_profile[i + 1][k] = sum;
		}

	/////////////////////////////////////////
	// Stripe query profile with chars
	qc = (unsigned char*) memalign(16,
			(NUMCOLSTATES + 1) * (LQ + 15) * sizeof(unsigned char),
			"the striped query profile during prefiltering"); // query profile (states + 1 because of ANY char)
	W = (LQ + 15) / 16;   // band width = hochgerundetes LQ/16

	for (a = 0; a < NUMCOLSTATES; ++a) {
		h = a * W * 16;
		for (i = 0; i < W; ++i) {
			j = i;
			for (k = 0; k < 16; ++k) {
				if (j >= LQ)
					qc[h] = (unsigned char) prefilter_score_offset;
				else {
					float dummy = flog2(query_profile[j + 1][a])
							* prefilter_bit_factor
							+ prefilter_score_offset + 0.5;
					// if (dummy>255.0) qc[h] = 255;
					// else if (dummy<0) qc[h] = 0;
					// else qc[h] = (unsigned char) dummy;  // 1/3 bits & make scores >=0 everywhere
					qc[h] = (unsigned char) fmax(0.0, fmin(255.0, dummy));
				}
				++h;
				j += W;
			}
		}
	}

	// Add extra ANY-state (220'th state)
	h = NUMCOLSTATES * W * 16;
	for (i = 0; i < W; ++i) {
		j = i;
		for (k = 0; k < 16; ++k) {
			if (j >= LQ)
				qc[h] = (unsigned char) prefilter_score_offset;
			else
				qc[h] = (unsigned char) (prefilter_score_offset - 1);
			h++;
			j += W;
		}
	}

	//////////////////////////////////////////////+
	// Stripe query profile with shorts
	unsigned short* qw = (unsigned short*) memalign(16,
			(NUMCOLSTATES + 1) * (LQ + 7) * sizeof(unsigned short),
			"the striped 2B query profile during prefiltering"); // query profile (states + 1 because of ANY char)
	int Ww = (LQ + 7) / 8;

	/////////////////////////////////////////
	// Stripe query profile
	for (a = 0; a < NUMCOLSTATES; ++a) {
		h = a * Ww * 8;
		for (i = 0; i < Ww; ++i) {
			j = i;
			for (k = 0; k < 8; ++k) {
				if (j >= LQ)
					qw[h] = 0;
				else {
					float dummy = flog2(query_profile[j + 1][a])
							* prefilter_bit_factor;
					qw[h] = (unsigned short) dummy; // 1/3 bits & make scores >=0 everywhere
				}
				++h;
				j += Ww;
			}
		}
	}

	// Add extra ANY-state
	h = NUMCOLSTATES * Ww * 8;
	for (i = 0; i < Ww; ++i) {
		j = i;
		for (k = 0; k < 8; ++k) {
			if (j >= LQ)
				qw[h] = 0;
			else
				qw[h] = (unsigned short) -1;
			h++;
			j += W;
		}
	}

	free(qw);

	for (i = 0; i < LQ + 1; ++i)
		free(query_profile[i]);
	delete[] query_profile;
}
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Main prefilter function
////////////////////////////////////////////////////////////////////////
void Prefilter::prefilter_db(HMM* q_tmp, Hash<Hit>* previous_hits,
		const int threads, const int prefilter_gap_open, const int prefilter_gap_extend,
		const int prefilter_score_offset, const int prefilter_bit_factor, const double prefilter_evalue_thresh,
		const double prefilter_evalue_coarse_thresh, const int preprefilter_smax_thresh, const int min_prefilter_hits, const float R[20][20],
		std::vector<std::string>& new_prefilter_hits, std::vector<std::string>& old_prefilter_hits) {

	Hash<char>* doubled = new Hash<char>;
	doubled->New(16381, 0);

	stripe_query_profile(q_tmp, prefilter_score_offset, prefilter_bit_factor);

	__m128i ** workspace = new __m128i *[threads];

	int score;
	double evalue;
	std::vector < std::pair<double, int> > first_prefilter;
	std::vector < std::pair<double, int> > hits;

	int thread_id = 0;
	int count_dbs = 0;
	int gap_init = prefilter_gap_open + prefilter_gap_extend;
	int gap_extend = prefilter_gap_extend;
	int LQ = q_tmp->L;
	const float log_qlen = flog2(LQ);
	const double factor = (double) num_dbs * LQ;

	for (int i = 0; i < threads; i++)
		workspace[i] = (__m128i *) memalign(16, 3 * (LQ + 15) * sizeof(char),
				"the dynamic programming workspace during prefiltering");

#pragma omp parallel for schedule(static) private(score, thread_id)
	// Loop over all database sequences
	for (size_t n = 0; n < num_dbs; n++) {
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif

		// Perform search step
		score = ungapped_sse_score(qc, LQ, first[n], length[n],
				prefilter_score_offset, workspace[thread_id]);

		score = score
				- (int) (prefilter_bit_factor
						* (log_qlen + flog2(length[n])));

#pragma omp critical
		first_prefilter.push_back(std::pair<double, int>(score, n));

		if (v >= 2 && !(n % 100000)) {
			std::cout << ".";
			std::cout.flush();
		}
	}

	//filter after calculation of ungapped sse score to include at least min_prefilter_hits
	std::vector<std::pair<double, int> >::iterator it;

	sort(first_prefilter.begin(), first_prefilter.end());
	std::reverse(first_prefilter.begin(), first_prefilter.end());

	std::vector<std::pair<double, int> >::iterator first_prefilter_begin_erase =
			first_prefilter.end();
	std::vector<std::pair<double, int> >::iterator first_prefilter_end_erase =
			first_prefilter.end();
	count_dbs = 0;
	for (it = first_prefilter.begin(); it < first_prefilter.end(); it++) {
		if (count_dbs >= min_prefilter_hits
				&& (*it).first < preprefilter_smax_thresh) {
			first_prefilter_begin_erase = it;
			break;
		} else {
			count_dbs++;
		}
	}

	first_prefilter.erase(first_prefilter_begin_erase,
			first_prefilter_end_erase);

	if (v >= 2) {
		printf(
				"\nHMMs passed 1st prefilter (gapless profile-profile alignment)  : %6i\n",
				count_dbs);
	}

#pragma omp parallel for schedule(static) private(evalue, score, thread_id)
	// Loop over all database sequences
//  for (int n = 0; n < count_dbs; n++) {
	for (it = first_prefilter.begin(); it < first_prefilter.end(); it++) {
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif

		int n = (*it).second;

		// Perform search step
		score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend,
				workspace[thread_id], workspace[thread_id] + W,
				workspace[thread_id] + 2 * W, prefilter_score_offset);

		evalue = factor * length[n] * fpow2(-score / prefilter_bit_factor);

		if (evalue < prefilter_evalue_coarse_thresh) {
#pragma omp critical
			hits.push_back(std::pair<double, int>(evalue, n));
		}
	}

	//filter after calculation of evalues to include at least min_prefilter_hits
	sort(hits.begin(), hits.end());

	std::vector<std::pair<double, int> >::iterator second_prefilter_begin_erase =
			hits.end();
	std::vector<std::pair<double, int> >::iterator second_prefilter_end_erase =
			hits.end();
	count_dbs = 0;
	for (it = hits.begin(); it < hits.end(); it++) {
		if (count_dbs >= min_prefilter_hits
				&& (*it).first > prefilter_evalue_thresh) {
			second_prefilter_begin_erase = it;
			break;
		} else {
			count_dbs++;
		}
	}

	hits.erase(second_prefilter_begin_erase, second_prefilter_end_erase);

	count_dbs = 0;

	for (it = hits.begin(); it < hits.end(); it++) {
		// Add hit to dbfiles
		char name[NAMELEN];
		strcpy(name, dbnames[(*it).second]);

		char db_name[NAMELEN];
		strcpy(db_name, name);

		if (!doubled->Contains(db_name)) {
			doubled->Add(db_name);

			// check, if DB was searched in previous rounds
			strcat(name, "__1");  // irep=1
			if (previous_hits->Contains(name)) {
				old_prefilter_hits.push_back(std::string(db_name));
			} else {
				new_prefilter_hits.push_back(std::string(db_name));
			}
		}
	}

	// Free memory
	free(qc);
	for (int i = 0; i < threads; i++)
		free(workspace[i]);
	delete[] workspace;
	if (doubled)
		delete doubled;
}

} /* namespace hh */
