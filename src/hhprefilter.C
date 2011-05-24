// hhprefilter.C
//
// some function for prefiltering in HHblits

#define SWAP(tmp, arg1, arg2) tmp = arg1; arg1 = arg2; arg2 = tmp;

struct ali_pos {
  int q_start;
  int q_stop;
  int t_start;
  int t_stop;
  float evalue;
};

// Needed variables

#define SHORT_BIAS 32768

int LDB = 0;              // number of characters of input cs-database
int num_dbs = 0;
Hash<char>* doubled;

int pos;
int block_count;
char actual_hit[NAMELEN];
int* block;

char** dbnames;
int LQ;
unsigned char* qc;
unsigned char* X;
unsigned char** first;
int* length;
int W;
unsigned short* qw;
int Ww;

// fast Smith-Waterman with SSE
int swStripedByte(unsigned char   *querySeq,
		  int              queryLength,
		  unsigned char   *dbSeq,
		  int              dbLength,
		  unsigned short   gapOpen,
		  unsigned short   gapExtend,
		  __m128i         *pvHLoad,
		  __m128i         *pvHStore,
		  __m128i         *pvE,
		  unsigned short   bias)
{
  int     i, j;
  int     score;
  
  int     cmp;
  int     iter = (queryLength + 15) / 16;
  
  __m128i *pv;
  
  __m128i vE, vF, vH;
  
  __m128i vMaxScore;
  __m128i vBias;
  __m128i vGapOpen;
  __m128i vGapExtend;
  
  __m128i vTemp;
  __m128i vZero;
  
  __m128i *pvScore;

  __m128i *pvQueryProf = (__m128i*) querySeq;
  
  /* Load the bias to all elements of a constant */
  vBias = _mm_set1_epi8(bias);
  
  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi8(gapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi8(gapExtend);
  
  vMaxScore = _mm_setzero_si128();
  vZero = _mm_setzero_si128();
  
  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i)
    {
      _mm_store_si128 (pvE + i, vMaxScore);
      _mm_store_si128 (pvHStore + i, vMaxScore);
    }

  for (i = 0; i < dbLength; ++i)
    {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;
      
      /* zero out F. */
      vF = _mm_setzero_si128();
      
      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 1);
      
      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;
      
      for (j = 0; j < iter; ++j)
        {
	  /* load values of vF and vH from previous row (one unit up) */
	  vE = _mm_load_si128 (pvE + j);
	  
	  /* add score to vH */
	  vH = _mm_adds_epu8 (vH, *(pvScore++));
	  vH = _mm_subs_epu8 (vH, vBias);
	  
	  /* Update highest score encountered this far */
	  vMaxScore = _mm_max_epu8 (vMaxScore, vH);

	  /* get max from vH, vE and vF */
	  vH = _mm_max_epu8 (vH, vE);
	  vH = _mm_max_epu8 (vH, vF);

	  /* save vH values */
	  _mm_store_si128 (pvHStore + j, vH);
	  
	  /* update vE value */
	  vH = _mm_subs_epu8 (vH, vGapOpen);
	  vE = _mm_subs_epu8 (vE, vGapExtend);
	  vE = _mm_max_epu8 (vE, vH);
	  
	  /* update vF value */
	  vF = _mm_subs_epu8 (vF, vGapExtend);
	  vF = _mm_max_epu8 (vF, vH);
	  
	  /* save vE values */
	  _mm_store_si128 (pvE + j, vE);
	  
	  /* load the next h value */
	  vH = _mm_load_si128 (pvHLoad + j);
        }
      
      /* reset pointers to the start of the saved data */
      j = 0;
      vH = _mm_load_si128 (pvHStore);
      
      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */
      vF = _mm_slli_si128 (vF, 1);
      vTemp = _mm_subs_epu8 (vH, vGapOpen);
      vTemp = _mm_subs_epu8 (vF, vTemp);
      vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
      cmp  = _mm_movemask_epi8 (vTemp);
      
      while (cmp != 0xffff) 
        {
	  vE = _mm_load_si128 (pvE + j);
	  
	  vH = _mm_max_epu8 (vH, vF);
	  
	  /* save vH values */
	  _mm_store_si128 (pvHStore + j, vH);
	  
	  /*  update vE incase the new vH value would change it */
	  vH = _mm_subs_epu8 (vH, vGapOpen);
	  vE = _mm_max_epu8 (vE, vH);
	  _mm_store_si128 (pvE + j, vE);
	  
	  /* update vF value */
	  vF = _mm_subs_epu8 (vF, vGapExtend);
	  
	  ++j;
	  if (j >= iter)
            {
	      j = 0;
	      vF = _mm_slli_si128 (vF, 1);
            }
	  
	  vH = _mm_load_si128 (pvHStore + j);
	  
	  vTemp = _mm_subs_epu8 (vH, vGapOpen);
	  vTemp = _mm_subs_epu8 (vF, vTemp);
	  vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
	  cmp  = _mm_movemask_epi8 (vTemp);
        }
    }
  
  /* find largest score in the vMaxScore vector */
  vTemp = _mm_srli_si128 (vMaxScore, 8);
  vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
  vTemp = _mm_srli_si128 (vMaxScore, 4);
  vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
  vTemp = _mm_srli_si128 (vMaxScore, 2);
  vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
  vTemp = _mm_srli_si128 (vMaxScore, 1);
  vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
  
  /* store in temporary variable */
  score = _mm_extract_epi16 (vMaxScore, 0);
  score = score & 0x00ff;

  /* return largest score */
  return score;
}

int swStripedWord_backtrace(int              queryLength,
			    unsigned char   *dbSeq,
			    int              dbLength,
			    unsigned short   gapOpen,
			    unsigned short   gapExtend,
			    __m128i         *pvQueryProf,
			    __m128i         *pvHLoad,
			    __m128i         *pvHStore,
			    __m128i         *pvE,
			    ali_pos          *res)
{

  int     i, j;
  int     cmp = 0;
  int     iter = (queryLength + 7) / 8;
  int     LQ = iter*8;
    
  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore;
  __m128i vGapOpen;
  __m128i vGapExtend;

  __m128i vMin;
  __m128i vTemp;

  __m128i *pvScore;

  __m128i Bt;
  __m128i Bttmp;
  __m128i v0x01 = _mm_set1_epi16(0x0001);
  __m128i v0x02 = _mm_set1_epi16(0x0010);
  __m128i v0x03 = _mm_set1_epi16(0x0100);
  __m128i v0x04 = _mm_set1_epi16(0x1000);

  __m128i* Hmatrix = (__m128i*)memalign(16,dbLength*LQ*sizeof(unsigned short));   // 2GB für 35000*35000 (titin)
  __m128i* Btmatrix = (__m128i*)memalign(16,dbLength*LQ*sizeof(unsigned short));
    
  __m128i* Hmatrix_it = Hmatrix;
  __m128i* Btmatrix_it = Btmatrix;

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi16(gapOpen);
    
  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi16(gapExtend);

  /*  load vMaxScore with the zeros.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short. */
  vMaxScore = _mm_set1_epi16(-1);
  vMaxScore = _mm_slli_epi16 (vMaxScore, 15);

  vMin = _mm_shuffle_epi32 (vMaxScore, 0);
  vMin = _mm_srli_si128 (vMin, 14);

  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i)
    {
      _mm_store_si128 (pvE + i, vMaxScore);
      _mm_store_si128 (pvHStore + i, vMaxScore);
    }

  for (i = 0; i < dbLength; ++i)
    {
      Hmatrix_it = Hmatrix + iter*i;
      Btmatrix_it = Btmatrix + iter*i;

      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;
      
      /* zero out F. */
      vF = _mm_set1_epi16(-1);
      vF = _mm_slli_epi16 (vF, 15);
      
      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 2);
      vH = _mm_or_si128 (vH, vMin);
      
      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;
      
      for (j = 0; j < iter; ++j)
        {
	  /* load values of vF and vH from previous row (one unit up) */
	  vE = _mm_load_si128 (pvE + j);
	  
	  /* add score to vH */
	  vH = _mm_adds_epi16 (vH, *pvScore++);
	  
	  /* Update highest score encountered this far */
	  vMaxScore = _mm_max_epi16 (vMaxScore, vH);
	  
	  /* get max from vH, vE and vF */
	  vH = _mm_max_epi16 (vH, vE);
	  vH = _mm_max_epi16 (vH, vF);

	  // Set backtrace register
	  Bt = _mm_setzero_si128();
	  Bttmp = _mm_cmpeq_epi16(vH, vE);
	  Bt = _mm_and_si128(Bttmp, v0x01);
	  Bttmp = _mm_cmpeq_epi16(vH, vF);
	  Bttmp = _mm_and_si128(Bttmp, v0x02);
	  Bt = _mm_or_si128(Bt,Bttmp);
	  
	  /* save vH values */
	  _mm_store_si128 (pvHStore + j, vH);
	  _mm_store_si128 (Hmatrix_it + j, vH);

	  /* update vE value */
	  vH = _mm_subs_epi16 (vH, vGapOpen);
	  vE = _mm_subs_epi16 (vE, vGapExtend);
	  vE = _mm_max_epi16 (vE, vH);
	  
	  /* update vF value */
	  vF = _mm_subs_epi16 (vF, vGapExtend);
	  vF = _mm_max_epi16 (vF, vH);

	  // Set backtrace register for gap matrices E and F
	  Bttmp = _mm_cmpeq_epi16(vH, vE);
	  Bttmp = _mm_and_si128(Bttmp, v0x03);
	  Bt = _mm_or_si128(Bt,Bttmp);
	  Bttmp = _mm_cmpeq_epi16(vH, vF);
	  Bttmp = _mm_and_si128(Bttmp, v0x04);
	  Bt = _mm_or_si128(Bt,Bttmp);
	  
	  /* save Bt values */
	  _mm_store_si128 (Btmatrix_it + j, Bt);

	  /* save vE values */
	  _mm_store_si128 (pvE + j, vE);

	  /* load the next h value */
	  vH = _mm_load_si128 (pvHLoad + j);
        }
      
      /* reset pointers to the start of the saved data */
      j = 0;
      vH = _mm_load_si128 (pvHStore + j);
      
      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */
      vF = _mm_slli_si128 (vF, 2);
      vF = _mm_or_si128 (vF, vMin);
      vTemp = _mm_subs_epi16 (vH, vGapOpen);
      vTemp = _mm_cmpgt_epi16 (vF, vTemp);
      cmp  = _mm_movemask_epi8 (vTemp);
      while (cmp != 0x0000) 
        {
	  vE = _mm_load_si128 (pvE + j);
	  
	  vH = _mm_max_epi16 (vH, vF);

	  // Set backtrace register
	  Bt = _mm_setzero_si128();
	  Bttmp = _mm_cmpeq_epi16(vH, vE);
	  Bt = _mm_and_si128(Bttmp, v0x01);
	  Bttmp = _mm_cmpeq_epi16(vH, vF);
	  Bttmp = _mm_and_si128(Bttmp, v0x02);
	  Bt = _mm_or_si128(Bt,Bttmp);

	  /* save vH values */
	  _mm_store_si128 (pvHStore + j, vH);
	  _mm_store_si128 (Hmatrix_it + j, vH);

	  /*  update vE incase the new vH value would change it */
	  vH = _mm_subs_epi16 (vH, vGapOpen);
	  vE = _mm_max_epi16 (vE, vH);

	  _mm_store_si128 (pvE + j, vE);

	  /* update vF value */
	  vF = _mm_subs_epi16 (vF, vGapExtend);

	  // Set backtrace register for gap matrices E and F
	  Bttmp = _mm_cmpeq_epi16(vH, vE);
	  Bttmp = _mm_and_si128(Bttmp, v0x03);
	  Bt = _mm_or_si128(Bt,Bttmp);
	  Bttmp = _mm_cmplt_epi16(vF, vH);
	  Bttmp = _mm_and_si128(Bttmp, v0x04);
	  Bt = _mm_or_si128(Bt,Bttmp);
	  
	  /* save Bt values */
	  _mm_store_si128 (Btmatrix_it + j, Bt);

	  ++j;
	  if (j >= iter)
            {
	      j = 0;
	      vF = _mm_slli_si128 (vF, 2);
	      vF = _mm_or_si128 (vF, vMin);
            }
	  
	  vH = _mm_load_si128 (pvHStore + j);
	  
	  vTemp = _mm_subs_epi16 (vH, vGapOpen);
	  vTemp = _mm_cmpgt_epi16 (vF, vTemp);
	  cmp  = _mm_movemask_epi8 (vTemp);
        }
    }
    
  // Extract all alignments with score > threshold
  __m128i vMax;
  __m128i Dt;
  __m128i Dtmax;
  __m128i Dq;
  __m128i Dqmax;
  __m128i Tmp;
  __m128i Tmp2;
  __m128i One = _mm_set1_epi16(1);
  __m128i vIter = _mm_setr_epi16(0,iter,2*iter,3*iter,4*iter,5*iter,6*iter,7*iter);
  
  int num = 0;
  int crossout_thresh = (int)(par.block_shading_space/1.3);
  double factor = (double)par.dbsize * (double)queryLength * (double)dbLength;

  int b,c,d,k,l;
  int q_pos, t_pos, pos;
  int d1, d2;
  int t_start, q_start, t_end, q_end;
  int kstart, lstart, h;

  while (num < 10) {

    // Find maximum score with position
    vMax = _mm_set1_epi16(-SHORT_BIAS);
    Dtmax = _mm_setzero_si128();
    Dqmax = _mm_setzero_si128();

    Hmatrix_it = Hmatrix;
      
    for (i = 0; i < dbLength; ++i)
      {
	//printf("DB-pos: %3i\n",i);
	Dq = vIter;
	Dt = _mm_set1_epi16(i);
	  
	for (j = 0; j < iter; ++j)
	  {
	    vH = _mm_load_si128(Hmatrix_it++);
	      
	    Tmp = _mm_cmpgt_epi16(vH, vMax);
	    vMax = _mm_max_epi16(vH, vMax);
	      
	    Tmp2 = _mm_and_si128(Tmp,Dt);
	    Dtmax = _mm_max_epi16(Dtmax, Tmp2);
	      
	    Dqmax = _mm_andnot_si128(Tmp,Dqmax);  // clear boxes, where new maximum is found
	    Tmp2 = _mm_and_si128(Tmp,Dq);
	    Dqmax = _mm_max_epi16(Dqmax, Tmp2);
	      
	    Dq = _mm_adds_epi16(Dq, One);
	  }
      }
      
    int max_score = _mm_extract_epi16(vMax,0) + SHORT_BIAS;
    int max_pos_q = _mm_extract_epi16(Dqmax,0);
    int max_pos_t = _mm_extract_epi16(Dtmax,0);

//     for (b = 1; b<8; ++b) {
//       if ((_mm_extract_epi16(vMax,b) + SHORT_BIAS) > max_score) {
// 	max_score = _mm_extract_epi16(vMax,b) + SHORT_BIAS;
// 	max_pos_q = _mm_extract_epi16(Dqmax,b);
// 	max_pos_t = _mm_extract_epi16(Dtmax,b);
//       }
//     }

    // MUST use this way because of Mac!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((_mm_extract_epi16(vMax,1) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,1) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,1);
      max_pos_t = _mm_extract_epi16(Dtmax,1);
    }
    if ((_mm_extract_epi16(vMax,2) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,2) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,2);
      max_pos_t = _mm_extract_epi16(Dtmax,2);
    }
    if ((_mm_extract_epi16(vMax,3) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,3) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,3);
      max_pos_t = _mm_extract_epi16(Dtmax,3);
    }
    if ((_mm_extract_epi16(vMax,4) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,4) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,4);
      max_pos_t = _mm_extract_epi16(Dtmax,4);
    }
    if ((_mm_extract_epi16(vMax,5) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,5) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,5);
      max_pos_t = _mm_extract_epi16(Dtmax,5);
    }
    if ((_mm_extract_epi16(vMax,6) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,6) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,6);
      max_pos_t = _mm_extract_epi16(Dtmax,6);
    }
    if ((_mm_extract_epi16(vMax,7) + SHORT_BIAS) > max_score) {
      max_score = _mm_extract_epi16(vMax,7) + SHORT_BIAS;
      max_pos_q = _mm_extract_epi16(Dqmax,7);
      max_pos_t = _mm_extract_epi16(Dtmax,7);
    }


    double evalue = factor * fpow2(-max_score/par.prefilter_bit_factor);
      
    if (num == 0 || evalue < par.prefilter_evalue_thresh) {
	
      // Backtrace alignment
      short *Btmatrix_sit = (short*) Btmatrix;
      short *Hmatrix_sit = (short*) Hmatrix;
      
      q_pos = max_pos_q;
      t_pos = max_pos_t;
      
      pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
      while (q_pos != 0 && t_pos != 0 && Hmatrix_sit[pos] != -SHORT_BIAS) {
	if ((Btmatrix_sit[pos]&0x0001) != 0) {
	  --t_pos;
	  pos -= LQ;
	  while (q_pos != 0 && t_pos != 0 && Hmatrix_sit[pos] != -SHORT_BIAS) {
	    if ((Btmatrix_sit[pos]&0x0100) == 0) {
	      --t_pos;
	      pos -= LQ;
	    } else {
	      --q_pos;
	      --t_pos;
	      break;
	    }
	  }
	} else if ((Btmatrix_sit[pos]&0x0010) != 0) {
	  --q_pos;
	  pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
	  while (q_pos != 0 && t_pos != 0 && Hmatrix_sit[pos] != -SHORT_BIAS) {
	    if ((Btmatrix_sit[pos]&0x1000) == 0) {
	      --q_pos;
	      pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
	    } else {
	      --q_pos;
	      --t_pos;
	      break;
	    }
	  }
	} else {
	  --q_pos;
	  --t_pos;
	} 
	pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
      }

      // Add 1 to each position because HHsearch counts sequences from position 1
      res[num].q_start = q_pos + 1;
      res[num].t_start = t_pos + 1;
      res[num].q_stop = imin(queryLength,max_pos_q + 1);
      res[num].t_stop = max_pos_t + 1;
      res[num++].evalue = evalue;

      // Cross out cells in H-matrix
      d1 = imin(t_pos-q_pos,max_pos_t-max_pos_q) - crossout_thresh;
      d2 = imax(t_pos-q_pos,max_pos_t-max_pos_q) + crossout_thresh;

      t_start = imax(0, t_pos - crossout_thresh);
      q_start = imax(0, q_pos - crossout_thresh);
      t_end = imin(dbLength, max_pos_t + crossout_thresh);
      q_end = imin(LQ, max_pos_q + crossout_thresh);

      kstart = (q_start % iter) * 8;
      lstart = q_start / iter;
      h = t_start*LQ;

      for (b = t_start; b < t_end; ++b) {
	k = kstart;
	l = lstart + h;
	d = b - q_end;
	for (c = b-q_start; c > d; --c) {
	  if (d1 < c && c < d2) {
	    Hmatrix_sit[k + l] = -SHORT_BIAS;
	  }
	  k+=8;
	  if (k == LQ) {
	    k = 0;
	    ++l;
	  }
	}
	h += LQ;
      }
    } else {
      break;
    }
  }
    
  free(Hmatrix);
  free(Btmatrix);

  return num;
}


// d = i-j+LT-1 is index of diagonal
int ungapped_sse_score(const unsigned char* query_profile,
		       const int            query_length,
		       const unsigned char* db_sequence,
		       const int            dbseq_length,
		       const unsigned char  score_offset,
		       __m128i*             workspace)
{
  int i; // position in query bands (0,..,W-1)
  int j; // position in db sequence (0,..,dbseq_length-1)
  int W = (query_length + 15) / 16; // width of bands in query and score matrix = hochgerundetes LQ/16
  
  __m128i *p;
  __m128i S;                  // 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15)
  __m128i Smax = _mm_setzero_si128();
  __m128i Soffset;            // all scores in query profile are shifted up by Soffset to obtain pos values 
  __m128i *s_prev, *s_curr;       // pointers to Score(i-1,j-1) and Score(i,j), resp.
  __m128i *qji;               // query profile score in row j (for residue x_j)
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

  for (j=0; j<dbseq_length; ++j) // loop over db sequence positions
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

      for (i=0; i<W; ++i) // loop over query band positions
        {
	  // Saturated addition and subtraction to score S(i,j)
	  S = _mm_adds_epu8(S, *(qji++));      // S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset) 
	  S = _mm_subs_epu8(S, Soffset);       // S(i,j) = max(0, S(i,j) - Soffset)
	  _mm_store_si128(s_curr_it++, S);       // store S to s_curr[i]
	  Smax = _mm_max_epu8(Smax, S);        // Smax(i,j) = max(Smax(i,j), S(i,j))
	  	  
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

void init_no_prefiltering()
{
  // Get DBsize
  char tmp_file[NAMELEN];
  strcpy(tmp_file, db);
  strcat(tmp_file, ".sizes");
  FILE* fin = fopen(tmp_file, "r");
  if (!fin)
    {
      FILE *stream;
      // Get DB-size
      command = "cat " + (string)db + " |grep \">\" |wc -l";
      stream = popen(command.c_str(), "r");
      ptr=fgets(line, LINELEN, stream);
      par.dbsize = strint(ptr);
      pclose(stream);
    } 
  else 
    {
      ptr=fgets(line, LINELEN, fin);
      par.dbsize = strint(ptr);
      fclose(fin);
    }
  
  if (par.dbsize > MAXNUMDB_NO_PREFILTER)
    {cerr<<endl<<"Error in "<<program_name<<": Without prefiltering, the max. number of database HHMs is "<<MAXNUMDB_NO_PREFILTER<<" (actual: "<<par.dbsize<<")\n"; exit(4);}

  char word[NAMELEN];
  FILE* dbf = NULL;
  dbf = fopen(db,"rb");
  if (!dbf) OpenFileError(db);

  while(fgetline(line,LINELEN,dbf)) // read HMM files in pal file
    {
      if (line[0]=='>')
	{
	  strwrd(word,line+1);

	  // Add hit to dbfiles
	  char db_name[NAMELEN];
	  
	  if (!strncmp(word,"cl|",3))   // kClust formatted database (NR20, NR30, UNIPROT20)
	    {
	      substr(db_name,word,3,11);
	      strcat(db_name,".hhm");
	    }
	  else                              // other database
	    {
	      strcpy(db_name,word);
	      strtr(db_name,"|", "_");
	      strcat(db_name,".");
	      strcat(db_name,db_ext);
	    }
	  
	  dbfiles_new[ndb_new]=new(char[strlen(db_name)+1]);
	  strcpy(dbfiles_new[ndb_new],db_name);
	  ndb_new++;
	}
    }
  fclose(dbf);

  cerr<<"Init without prefiltering! "<<ndb_new<<" HHMs in database!\n";

}
      

void init_prefilter()
{
  // Get Prefilter Pvalue (Evalue / Par.Dbsize)
  char tmp_file[NAMELEN];
  strcpy(tmp_file, db);
  strcat(tmp_file, ".sizes");
  FILE* fin = fopen(tmp_file, "r");
  if (!fin)
    {
      FILE *stream;
      // Get DB-size
      command = "cat " + (string)db + " |grep \"^>\" |wc -l";
      stream = popen(command.c_str(), "r");
      ptr=fgets(line, LINELEN, stream);
      par.dbsize = strint(ptr);
      pclose(stream);
      // Get DB-length
      command = "cat " + (string)db + " |grep -v \"^>\" |wc -c";
      stream = popen(command.c_str(), "r");
      ptr=fgets(line, LINELEN, stream);
      LDB = strint(ptr);
      pclose(stream);
    } 
  else 
    {
      ptr=fgets(line, LINELEN, fin);
      par.dbsize = strint(ptr);
      LDB = strint(ptr);
      fclose(fin);
    }

  if (par.dbsize == 0 || LDB == 0)
    {cerr<<endl<<"Error! Could not determine DB-size of prefilter db ("<<db<<")\n"; exit(4);}
	    
  par.hhblits_prefilter_logpval=-log(par.prefilter_evalue_thresh / (float)par.dbsize);

  if (v>1) printf("\nFull database size: %6i\n",par.dbsize);

  X = (unsigned char*)memalign(16,LDB*sizeof(unsigned char));                 // database string (concatenate all DB-seqs)
  first = (unsigned char**)memalign(16,(2*par.dbsize)*sizeof(unsigned char*));    // first characters of db sequences
  length = (int*)memalign(16,(2*par.dbsize)*sizeof(int));                         // lengths of db sequences
  dbnames = new char*[par.dbsize*2];                                              // names of db sequences

  /////////////////////////////////////////
  // Read in database
  num_dbs = 0;
  int len = 0;
  int pos = 0;
  int c;
  char word[NAMELEN];
  FILE* dbf = NULL;
  dbf = fopen(db,"rb");
  if (!dbf) OpenFileError(db);

  while(fgetline(line,LINELEN,dbf)) // read prefilter database
    {
      if (line[0] == '>')  // Header
	{
	  if (len > 0)           // if it is not the first sequence
	    length[num_dbs++] = len;
	  len = 0;
	      
	  strwrd(word,line+1);

	  dbnames[num_dbs]=new(char[strlen(word)+1]);
	  strcpy(dbnames[num_dbs],word);
	      
	  first[num_dbs] = X + pos;
	}
      else
	{
	  char* linep = line;
	  while (*linep!=0)
	    {
	      c = *linep < 0 ? *linep + 256 : *linep;
	      if (cs::AS219::kValidChar[c])
	      	{
	      	  X[pos++]=(unsigned char)(cs::AS219::kCharToInt[c]);
	      	  len++;
	      	}
	      else
	      	cerr<<endl<<"WARNING: invalid symbol \'"<<*linep<<"\' of "<<db<<"\n";

	      linep++;
	    }
	  
	}
    }

  if (len > 0)
    length[num_dbs++] = len;
      
  fclose(dbf);
}

void stripe_query_profile()
{
  LQ=q_tmp->L;
  float** query_profile = NULL;
  int a,h,i,j,k;

  // Add Pseudocounts
  if (!*par.clusterfile) {
    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
    q_tmp->PreparePseudocounts();
    // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    q_tmp->AddAminoAcidPseudocounts(par.pcm,par.pca,par.pcb,1);
  } else {
    // Add context specific pseudocounts (now always used, because clusterfile is necessary)
    q_tmp->AddContextSpecificPseudocounts(5,par.pre_pca,par.pre_pcb,1);
  }
      
  q_tmp->CalculateAminoAcidBackground();

  // Build query profile with 62 column states
  query_profile = new float*[LQ+1];
  for (i=0; i<LQ+1; ++i) 
    query_profile[i]=(float*) memalign(16,cs::AS219::kSize*sizeof(float));
      
  const cs::ContextLibrary<cs::AA>& lib = *cs_lib;

  // log (S(i,c)) = log ( SUM_a p(i,a) * p(c,a) / f(a) )   c: column state, i: pos in ali, a: amino acid
  for (i=0; i<LQ; ++i)
    for (k=0; k<(int)cs::AS219::kSize; ++k)
      {
	float sum = 0;
	for (a=0; a<20; ++a)
	  sum += ((q_tmp->p[i][a] * lib[k].probs[0][a]) / q_tmp->pav[a]);
	query_profile[i+1][k] = sum;
      }
      
  /////////////////////////////////////////
  // Stripe query profile with chars
  qc = (unsigned char*)memalign(16,(par.prefilter_states+1)*(LQ+15)*sizeof(unsigned char));   // query profile (states + 1 because of ANY char)
  W = (LQ+15) / 16;   // band width = hochgerundetes LQ/16
  
  for (a=0; a < par.prefilter_states; ++a)  
    {
      h = a*W*16;
      for (i=0; i < W; ++i) 
  	{
  	  j = i;
  	  for (k = 0; k < 16; ++k)
  	    {
  	      if (j >= LQ)
  		qc[h]=(unsigned char) par.prefilter_score_offset;
  	      else
  		{
  		  float dummy = flog2(query_profile[j+1][a])*par.prefilter_bit_factor + par.prefilter_score_offset + 0.5;
  		  if (dummy>255.99) dummy = 255.5;
  		  if (dummy<0) dummy = 0.0;
  		  qc[h] = (unsigned char) dummy;  // 1/3 bits & make scores >=0 everywhere
  		}
  	      ++h;
  	      j+=W;
  	    }
  	}
    }

  // Add extra ANY-state
  h = par.prefilter_states*W*16;
  for (i=0; i < W; ++i) 
    {
      j = i;
      for (k = 0; k < 16; ++k)
	{
	  if (j >= LQ)
	    qc[h]=(unsigned char) par.prefilter_score_offset;
	  else
	    qc[h]=(unsigned char) (par.prefilter_score_offset - 1);
	  h++;
	  j+=W;
	}
    }
  
  //////////////////////////////////////////////+
  // Stripe query profile with shorts
  qw = (unsigned short*)memalign(16,(par.prefilter_states+1)*(LQ+7)*sizeof(unsigned short));   // query profile (states + 1 because of ANY char)
  Ww = (LQ+7) / 8;
  
  /////////////////////////////////////////
  // Stripe query profile
  for (a=0; a < par.prefilter_states; ++a)  
    {
      h = a*Ww*8;
      for (i=0; i < Ww; ++i) 
	{
	  j = i;
	  for (k = 0; k < 8; ++k)
	    {
	      if (j >= LQ)
		qw[h] = 0;
	      else
		{
		  float dummy = flog2(query_profile[j+1][a])*par.prefilter_bit_factor;
		  qw[h] = (unsigned short) dummy;  // 1/3 bits & make scores >=0 everywhere
		}
	      ++h;
	      j+=Ww;
	    }
	}
    }
      
  // Add extra ANY-state
  h = par.prefilter_states*Ww*8;
  for (i=0; i < Ww; ++i) 
    {
      j = i;
      for (k = 0; k < 8; ++k)
	{
	  if (j >= LQ)
	    qw[h] = 0;
	  else
	    qw[h] = (unsigned short) -1;
	  h++;
	  j+=W;
	}
    }

  for (i=0; i<LQ+1; ++i) 
    free(query_profile[i]);
  delete[] query_profile;
}

void prefilter_with_SW_evalue_preprefilter_backtrace()
{
  stripe_query_profile();
  
  int* prefiltered_hits = new int[par.dbsize+1];
  int* backtrace_hits = new int[MAXNUMDB+1];

  __m128i** workspace = new(__m128i*[cpu]);

  int score;
  double evalue;
  vector<pair<double, int> > hits;

  int thread_id = 0;
  int count_dbs = 0;
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;
  const float log_qlen = flog2(LQ);
  const double factor = (double)par.dbsize * LQ;

  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
#pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#endif
  
      // Perform search step
      score = ungapped_sse_score(qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace[thread_id]);
  
      score = score - par.prefilter_bit_factor * (log_qlen + flog2(length[n]));
    
      if (score > par.preprefilter_smax_thresh)
	{
#pragma omp critical
	  prefiltered_hits[count_dbs++] = n;
	}
    }
  if (v>=2)
    {
      printf("Hits passed prefilter 1 (gapless profile-profile alignment): %6i\n", count_dbs);
      //printf("%6i hits through preprefilter!\n", count_dbs);
    }
  if (print_elapsed) ElapsedTimeSinceLastCall("(ungapped preprefilter)");
  
#pragma omp parallel for schedule(static) private(evalue, score, thread_id)
  for (int n = 0; n < count_dbs; n++)     // Loop over all database sequences
    {
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[prefiltered_hits[n]], length[prefiltered_hits[n]], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      evalue = factor * length[prefiltered_hits[n]] * fpow2(-score/par.prefilter_bit_factor);
 
      if (evalue < par.prefilter_evalue_thresh)
	{
#pragma omp critical
	  hits.push_back(pair<double,int>(evalue, prefiltered_hits[n]));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<double, int> >::iterator it;
  count_dbs = 0;
  
  for ( it=hits.begin() ; it < hits.end(); it++ )
    {
      backtrace_hits[count_dbs++] = (*it).second;

      // Add hit to dbfiles
      char tmp_name[NAMELEN];
      char db_name[NAMELEN];
      ptr=strwrd(tmp_name,dbnames[(*it).second]);
      
      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30, UNIPROT20)
	{
	  substr(db_name,tmp_name,3,11);
	  strcat(db_name,".hhm");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".");
	  strcat(db_name,db_ext);
	}
      
      if (! doubled->Contains(db_name))
	{
	  doubled->Add(db_name);
	  // check, if DB was searched in previous rounds 
	  strcat(tmp_name,"__1");  // irep=1
	  if (previous_hits->Contains(tmp_name))
	    {
	      dbfiles_old[ndb_old]=new(char[strlen(db_name)+1]);
	      strcpy(dbfiles_old[ndb_old],db_name);
	      ndb_old++;
	    }
	  else 
	    {
	      dbfiles_new[ndb_new]=new(char[strlen(db_name)+1]);
	      strcpy(dbfiles_new[ndb_new],db_name);
	      ndb_new++;
	    }
	}

      if (count_dbs >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i, hits = %6i)\n",MAXNUMDB, (int)hits.size());
	  break;
	}
    }

  if (print_elapsed) ElapsedTimeSinceLastCall("(SW prefilter)");

  if (block_filter)
    {
      // Run SW with backtrace
      for (int i = 0; i < cpu; i++) {
	free(workspace[i]);
	workspace[i] = (__m128i*)memalign(16,3*(LQ+7)*sizeof(short));
      }
      __m128i *qw_it = (__m128i*) qw;
      
#pragma omp parallel for schedule(static) private(block, block_count, thread_id)
      for (int n = 0; n < count_dbs; n++)     // Loop over all database sequences
	{
#ifdef _OPENMP
	  thread_id = omp_get_thread_num();
#endif
            
	  // Perform backtrace, if one of the profiles has length > 2*par.block_shading_space
	  if (LQ > 2*par.block_shading_space || length[backtrace_hits[n]] > 2*par.block_shading_space)
	    {
	      ali_pos *res = new ali_pos[10];
	      
	      // Perform search step
	      int num_res = swStripedWord_backtrace(LQ, first[backtrace_hits[n]], length[backtrace_hits[n]], gap_init, gap_extend, qw_it, workspace[thread_id], workspace[thread_id] + Ww, workspace[thread_id] + 2*Ww, res);
	      
	      if (num_res > 0) 
		{
		  char tmp_name[NAMELEN];
		  ptr=strwrd(tmp_name,dbnames[backtrace_hits[n]]);
		  block = new(int[400]);
		  block_count = 0;
		  for (int a = 0; a < num_res; a++) 
		    {
		      if (block_count >= 400) { continue; }
		      // Get block of HSP
		      block[block_count++]=res[a].q_start;
		      block[block_count++]=res[a].q_stop;
		      block[block_count++]=res[a].t_start;
		      block[block_count++]=res[a].t_stop;
		    }
#pragma omp critical
		  {
		    par.block_shading->Add(tmp_name,block);
		    par.block_shading_counter->Add(tmp_name,block_count);
		  }
		}
	      delete[] res;
	    }
	}
      if (print_elapsed) ElapsedTimeSinceLastCall("(SW backtrace prefilter)");

      free(qw);
    }

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
  delete[] prefiltered_hits;
  delete[] backtrace_hits;
}

// Main prefilter function
void prefilter_db()
{
  doubled = new(Hash<char>);
  doubled->New(16381,0);
  for (int idb=0; idb<ndb_new; idb++) delete[](dbfiles_new[idb]);
  for (int idb=0; idb<ndb_old; idb++) delete[](dbfiles_old[idb]);
  ndb_new = ndb_old = 0;
  block_count=0;
  //block = new(int[400]);
  block = NULL;
  strcpy(actual_hit,"");

  par.block_shading->Reset();
  while (!par.block_shading->End())
    delete[] (par.block_shading->ReadNext()); 
  par.block_shading->New(16381,NULL);
  par.block_shading_counter->New(16381,NULL);
  
  prefilter_with_SW_evalue_preprefilter_backtrace();

  if(doubled) delete doubled;
}
