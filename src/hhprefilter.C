// hhprefilter.C
//
// some function for prefiltering in HHblits

#define SWAP(tmp, arg1, arg2) tmp = arg1; arg1 = arg2; arg2 = tmp;

struct triple {
  int first;
  int second;
  int third;
};

struct ali_pos {
  int q_start;
  int q_stop;
  int t_start;
  int t_stop;
  float evalue;
};

#define SHORT_BIAS 32768

int dbsize = 0;
int LDB = 0;
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

// int
// swStripedByte(int              queryLength,
//               unsigned char   *dbSeq,
//               int              dbLength,
//               unsigned short   gapOpen,
//               unsigned short   gapExtend,
//               __m128i         *pvQueryProf,
//               __m128i         *pvHLoad,
//               __m128i         *pvHStore,
//               __m128i         *pvE,
//               unsigned short   bias)
// {
//     int     i, j;
//     int     score;

//     int     dup;
//     int     cmp;
//     int     iter = (queryLength + 15) / 16;
    
//     __m128i *pv;

//     __m128i vE, vF, vH;

//     __m128i vMaxScore;
//     __m128i vBias;
//     __m128i vGapOpen;
//     __m128i vGapExtend;

//     __m128i vTemp;
//     __m128i vZero;

//     __m128i *pvScore;

//     /* Load the bias to all elements of a constant */
//     dup    = (bias << 8) | (bias & 0x00ff);
//     vBias = _mm_insert_epi16 (vBias, dup, 0);
//     vBias = _mm_shufflelo_epi16 (vBias, 0);
//     vBias = _mm_shuffle_epi32 (vBias, 0);

//     /* Load gap opening penalty to all elements of a constant */
//     dup    = (gapOpen << 8) | (gapOpen & 0x00ff);
//     vGapOpen = _mm_insert_epi16 (vGapOpen, dup, 0);
//     vGapOpen = _mm_shufflelo_epi16 (vGapOpen, 0);
//     vGapOpen = _mm_shuffle_epi32 (vGapOpen, 0);

//     /* Load gap extension penalty to all elements of a constant */
//     dup    = (gapExtend << 8) | (gapExtend & 0x00ff);
//     vGapExtend = _mm_insert_epi16 (vGapExtend, dup, 0);
//     vGapExtend = _mm_shufflelo_epi16 (vGapExtend, 0);
//     vGapExtend = _mm_shuffle_epi32 (vGapExtend, 0);

//     vMaxScore = _mm_xor_si128 (vMaxScore, vMaxScore);

//     vZero = _mm_xor_si128 (vZero, vZero);

//     /* Zero out the storage vector */
//     for (i = 0; i < iter; i++)
//     {
//         _mm_store_si128 (pvE + i, vMaxScore);
//         _mm_store_si128 (pvHStore + i, vMaxScore);
//     }

//     for (i = 0; i < dbLength; ++i)
//     {
//         /* fetch first data asap. */
//         pvScore = pvQueryProf + dbSeq[i] * iter;

//         /* zero out F. */
//         vF = _mm_xor_si128 (vF, vF);

//         /* load the next h value */
//         vH = _mm_load_si128 (pvHStore + iter - 1);
//         vH = _mm_slli_si128 (vH, 1);

//         pv = pvHLoad;
//         pvHLoad = pvHStore;
//         pvHStore = pv;

//         for (j = 0; j < iter; j++)
//         {
//             /* load values of vF and vH from previous row (one unit up) */
//             vE = _mm_load_si128 (pvE + j);

//             /* add score to vH */
//             vH = _mm_adds_epu8 (vH, *(pvScore + j));
//             vH = _mm_subs_epu8 (vH, vBias);

//             /* Update highest score encountered this far */
//             vMaxScore = _mm_max_epu8 (vMaxScore, vH);

//             /* get max from vH, vE and vF */
//             vH = _mm_max_epu8 (vH, vE);
//             vH = _mm_max_epu8 (vH, vF);

//             /* save vH values */
//             _mm_store_si128 (pvHStore + j, vH);

//             /* update vE value */
//             vH = _mm_subs_epu8 (vH, vGapOpen);
//             vE = _mm_subs_epu8 (vE, vGapExtend);
//             vE = _mm_max_epu8 (vE, vH);

//             /* update vF value */
//             vF = _mm_subs_epu8 (vF, vGapExtend);
//             vF = _mm_max_epu8 (vF, vH);

//             /* save vE values */
//             _mm_store_si128 (pvE + j, vE);

//             /* load the next h value */
//             vH = _mm_load_si128 (pvHLoad + j);
//         }

//         /* reset pointers to the start of the saved data */
//         j = 0;
//         vH = _mm_load_si128 (pvHStore + j);

//         /*  the computed vF value is for the given column.  since */
//         /*  we are at the end, we need to shift the vF value over */
//         /*  to the next column. */
//         vF = _mm_slli_si128 (vF, 1);
//         vTemp = _mm_subs_epu8 (vH, vGapOpen);
//         vTemp = _mm_subs_epu8 (vF, vTemp);
//         vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
//         cmp  = _mm_movemask_epi8 (vTemp);

//         while (cmp != 0xffff) 
//         {
//             vE = _mm_load_si128 (pvE + j);

//             vH = _mm_max_epu8 (vH, vF);

//             /* save vH values */
//             _mm_store_si128 (pvHStore + j, vH);

//             /*  update vE incase the new vH value would change it */
//             vH = _mm_subs_epu8 (vH, vGapOpen);
//             vE = _mm_max_epu8 (vE, vH);
//             _mm_store_si128 (pvE + j, vE);

//             /* update vF value */
//             vF = _mm_subs_epu8 (vF, vGapExtend);

//             j++;
//             if (j >= iter)
//             {
//                 j = 0;
//                 vF = _mm_slli_si128 (vF, 1);
//             }

//             vH = _mm_load_si128 (pvHStore + j);

//             vTemp = _mm_subs_epu8 (vH, vGapOpen);
//             vTemp = _mm_subs_epu8 (vF, vTemp);
//             vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
//             cmp  = _mm_movemask_epi8 (vTemp);
//         }
//     }

//     /* find largest score in the vMaxScore vector */
//     vTemp = _mm_srli_si128 (vMaxScore, 8);
//     vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
//     vTemp = _mm_srli_si128 (vMaxScore, 4);
//     vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
//     vTemp = _mm_srli_si128 (vMaxScore, 2);
//     vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
//     vTemp = _mm_srli_si128 (vMaxScore, 1);
//     vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);

//     /* store in temporary variable */
//     score = _mm_extract_epi16 (vMaxScore, 0);
//     score = score & 0x00ff;

//     /*  check if we might have overflowed */
//     if (score + bias >= 255)
//     {
//         score = 255;
//     }

//     /* return largest score */
//     return score;
// }

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

  //int *matrix = new int[dbLength*iter*16];
  
  for (i = 0; i < dbLength; ++i)
    {
      //printf("\nDB-Char: %c", i2aa(dbSeq[i]));
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
	  
	  //printf("\nVH  ");
	  //for (int b = 0; b<16; b++) // loop over bytes in XMM registers
	    //matrix[i*iter*16+j*16+b] = (int)*(((unsigned char*) &vH) + b);
	    //printf("%2i|",(int)*(((unsigned char*) &vH) + b));

	  /* Update highest score encountered this far */
	  vMaxScore = _mm_max_epu8 (vMaxScore, vH);
	  // if (i == 0) {
	  //   printf("\nVMaxScore  ");
	  //   for (int b = 0; b<16; b++) // loop over bytes in XMM registers
	  //     printf(" %3i |",(int)*(((unsigned char*) &vMaxScore) + b));
	  // }

	  /* get max from vH, vE and vF */
	  vH = _mm_max_epu8 (vH, vE);
	  vH = _mm_max_epu8 (vH, vF);
	  // if (i == 0) {
	  //   printf("\nvH after max with vE & vF  ");
	  //   for (int b = 0; b<16; b++) // loop over bytes in XMM registers
	  //     printf(" %3i |",(int)*(((unsigned char*) &vH) + b));
	  // }

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
      // printf("\nVMaxScore  ");
      // for (int b = 0; b<16; b++) // loop over bytes in XMM registers
      //  	printf(" %3i |",(int)*(((unsigned char*) &vMaxScore) + b));
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

  // Print vH matrix
  // for (i = 0; i < dbLength; ++i)
  //   {
  //     printf("\nDB-Char: %c  ", i2aa(dbSeq[i]));
  //     for (int b = 0; b<16; b++)
  // 	{
  // 	  for (j = 0; j < iter; j++)
  // 	    {
  // 	      printf("%2i|", matrix[i*iter*16+j*16+b]);
  // 	    }
  // 	}
  //   }
    
  /* return largest score */
  return score;
}

// int swStripedWord_backtrace(int              queryLength,
// 			    unsigned char   *dbSeq,
// 			    int              dbLength,
// 			    unsigned short   gapOpen,
// 			    unsigned short   gapExtend,
// 			    __m128i         *pvQueryProf,
// 			    __m128i         *pvHLoad,
// 			    __m128i         *pvHStore,
// 			    __m128i         *pvE,
// 			    ali_pos          *res)
// {

//     int     i, j;
//     // int     score;

//     int     cmp = 0;
//     int     iter = (queryLength + 7) / 8;
//     int     LQ = iter*8;
    
//     __m128i *pv;

//     __m128i vE, vF, vH;

//     __m128i vMaxScore;
//     __m128i vGapOpen;
//     __m128i vGapExtend;

//     __m128i vMin;
//     __m128i vTemp;

//     __m128i *pvScore;

//     __m128i Bt;
//     __m128i Bttmp;
//     __m128i v0x01 = _mm_set1_epi16(0x0001);
//     __m128i v0x02 = _mm_set1_epi16(0x0010);
//     __m128i v0x03 = _mm_set1_epi16(0x0100);
//     __m128i v0x04 = _mm_set1_epi16(0x1000);

//     __m128i* Hmatrix = (__m128i*)memalign(16,dbLength*LQ*sizeof(unsigned short));
//     __m128i* Btmatrix = (__m128i*)memalign(16,dbLength*LQ*sizeof(unsigned short));
    
//     __m128i* Hmatrix_it = Hmatrix;
//     __m128i* Btmatrix_it = Btmatrix;

//     /* Load gap opening penalty to all elements of a constant */
//     vGapOpen = _mm_set1_epi16(gapOpen);
    
//     /* Load gap extension penalty to all elements of a constant */
//     vGapExtend = _mm_set1_epi16(gapExtend);

//     /*  load vMaxScore with the zeros.  since we are using signed */
//     /*  math, we will bias the maxscore to -32768 so we have the */
//     /*  full range of the short. */
//     vMaxScore = _mm_set1_epi16(-1);
//     vMaxScore = _mm_slli_epi16 (vMaxScore, 15);

//     vMin = _mm_shuffle_epi32 (vMaxScore, 0);
//     vMin = _mm_srli_si128 (vMin, 14);

//     /* Zero out the storage vector */
//     for (i = 0; i < iter; ++i)
//     {
//         _mm_store_si128 (pvE + i, vMaxScore);
//         _mm_store_si128 (pvHStore + i, vMaxScore);
//     }

//     for (i = 0; i < dbLength; ++i)
//     {
//       Hmatrix_it = Hmatrix + iter*i;
//       Btmatrix_it = Btmatrix + iter*i;

//       /* fetch first data asap. */
//       pvScore = pvQueryProf + dbSeq[i] * iter;
      
//       /* zero out F. */
//       vF = _mm_set1_epi16(-1);
//       vF = _mm_slli_epi16 (vF, 15);
      
//       /* load the next h value */
//       vH = _mm_load_si128 (pvHStore + iter - 1);
//       vH = _mm_slli_si128 (vH, 2);
//       vH = _mm_or_si128 (vH, vMin);
      
//       pv = pvHLoad;
//       pvHLoad = pvHStore;
//       pvHStore = pv;
      
//       for (j = 0; j < iter; ++j)
//         {
// 	  /* load values of vF and vH from previous row (one unit up) */
// 	  vE = _mm_load_si128 (pvE + j);
	  
// 	  /* add score to vH */
// 	  vH = _mm_adds_epi16 (vH, *pvScore++);
	  
// 	  /* Update highest score encountered this far */
// 	  vMaxScore = _mm_max_epi16 (vMaxScore, vH);
	  
// 	  /* get max from vH, vE and vF */
// 	  vH = _mm_max_epi16 (vH, vE);
// 	  vH = _mm_max_epi16 (vH, vF);

// 	  // Set backtrace register
// 	  Bt = _mm_setzero_si128();
// 	  Bttmp = _mm_cmpeq_epi16(vH, vE);
// 	  Bt = _mm_and_si128(Bttmp, v0x01);
// 	  Bttmp = _mm_cmpeq_epi16(vH, vF);
// 	  Bttmp = _mm_and_si128(Bttmp, v0x02);
// 	  Bt = _mm_or_si128(Bt,Bttmp);
	  
// 	  /* save vH values */
// 	  _mm_store_si128 (pvHStore + j, vH);
// 	  _mm_store_si128 (Hmatrix_it + j, vH);

// 	  /* update vE value */
// 	  vH = _mm_subs_epi16 (vH, vGapOpen);
// 	  vE = _mm_subs_epi16 (vE, vGapExtend);
// 	  vE = _mm_max_epi16 (vE, vH);
	  
// 	  /* update vF value */
// 	  vF = _mm_subs_epi16 (vF, vGapExtend);
// 	  vF = _mm_max_epi16 (vF, vH);

// 	  // Set backtrace register for gap matrices E and F
// 	  Bttmp = _mm_cmpeq_epi16(vH, vE);
// 	  Bttmp = _mm_and_si128(Bttmp, v0x03);
// 	  Bt = _mm_or_si128(Bt,Bttmp);
// 	  Bttmp = _mm_cmpeq_epi16(vH, vF);
// 	  Bttmp = _mm_and_si128(Bttmp, v0x04);
// 	  Bt = _mm_or_si128(Bt,Bttmp);
	  
// 	  /* save Bt values */
// 	  _mm_store_si128 (Btmatrix_it + j, Bt);

// 	  /* save vE values */
// 	  _mm_store_si128 (pvE + j, vE);

// 	  /* load the next h value */
// 	  vH = _mm_load_si128 (pvHLoad + j);
//         }
      
//       /* reset pointers to the start of the saved data */
//       j = 0;
//       vH = _mm_load_si128 (pvHStore + j);
      
//       /*  the computed vF value is for the given column.  since */
//       /*  we are at the end, we need to shift the vF value over */
//       /*  to the next column. */
//       vF = _mm_slli_si128 (vF, 2);
//       vF = _mm_or_si128 (vF, vMin);
//       vTemp = _mm_subs_epi16 (vH, vGapOpen);
//       vTemp = _mm_cmpgt_epi16 (vF, vTemp);
//       cmp  = _mm_movemask_epi8 (vTemp);
//       while (cmp != 0x0000) 
//         {
// 	  vE = _mm_load_si128 (pvE + j);
	  
// 	  vH = _mm_max_epi16 (vH, vF);

// 	  // Set backtrace register
// 	  Bt = _mm_setzero_si128();
// 	  Bttmp = _mm_cmpeq_epi16(vH, vE);
// 	  Bt = _mm_and_si128(Bttmp, v0x01);
// 	  Bttmp = _mm_cmpeq_epi16(vH, vF);
// 	  Bttmp = _mm_and_si128(Bttmp, v0x02);
// 	  Bt = _mm_or_si128(Bt,Bttmp);

// 	  /* save vH values */
// 	  _mm_store_si128 (pvHStore + j, vH);
// 	  _mm_store_si128 (Hmatrix_it + j, vH);

// 	  /*  update vE incase the new vH value would change it */
// 	  vH = _mm_subs_epi16 (vH, vGapOpen);
// 	  vE = _mm_max_epi16 (vE, vH);

// 	  _mm_store_si128 (pvE + j, vE);

// 	  /* update vF value */
// 	  vF = _mm_subs_epi16 (vF, vGapExtend);

// 	  // Set backtrace register for gap matrices E and F
// 	  Bttmp = _mm_cmpeq_epi16(vH, vE);
// 	  Bttmp = _mm_and_si128(Bttmp, v0x03);
// 	  Bt = _mm_or_si128(Bt,Bttmp);
// 	  Bttmp = _mm_cmplt_epi16(vF, vH);
// 	  Bttmp = _mm_and_si128(Bttmp, v0x04);
// 	  Bt = _mm_or_si128(Bt,Bttmp);
	  
// 	  /* save Bt values */
// 	  _mm_store_si128 (Btmatrix_it + j, Bt);

// 	  ++j;
// 	  if (j >= iter)
//             {
// 	      j = 0;
// 	      vF = _mm_slli_si128 (vF, 2);
// 	      vF = _mm_or_si128 (vF, vMin);
//             }
	  
// 	  vH = _mm_load_si128 (pvHStore + j);
	  
// 	  vTemp = _mm_subs_epi16 (vH, vGapOpen);
// 	  vTemp = _mm_cmpgt_epi16 (vF, vTemp);
// 	  cmp  = _mm_movemask_epi8 (vTemp);
//         }
//     }
    
//     // printf("Query-pos:");
//     // for (int b = 0; b < queryLength; b++) {
//     //   printf("%5i ", b);
//     // }
//     // printf("\n");

//     // short *tmp = (short*) Btmatrix;
//     // for (int b = 0; b < dbLength; b++) {
//     //   printf("DB: %3i   : ", b);
//     //   for (int c = 0; c < 8; c++) {
//     // 	for (int d = 0; d < iter; d++) {
//     // 	  int i = b*iter*8 + d * 8 + c;
//     // 	  if ((tmp[i]&0x1000) != 0) { printf("1"); } else { printf("0"); }
//     // 	  if ((tmp[i]&0x0100) != 0) { printf("1"); } else { printf("0"); }
//     // 	  if ((tmp[i]&0x0010) != 0) { printf("1"); } else { printf("0"); }
//     // 	  if ((tmp[i]&0x0001) != 0) { printf("1"); } else { printf("0"); }
//     // 	  printf(" |", tmp[i]);
//     // 	}
//     //   }
//     //   printf("\n");
//     // }
//     // printf("\n");

//     // /* find largest score in the vMaxScore vector */
//     // vTemp = _mm_srli_si128 (vMaxScore, 8);
//     // vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
//     // vTemp = _mm_srli_si128 (vMaxScore, 4);
//     // vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
//     // vTemp = _mm_srli_si128 (vMaxScore, 2);
//     // vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);

//     // /* store in temporary variable */
//     // score = (short) _mm_extract_epi16 (vMaxScore, 0);

//     // printf("Largest Score: %4i!\n", score + SHORT_BIAS);

//     // Extract all alignments with score > threshold
//     __m128i vMax;
//     __m128i Dt;
//     __m128i Dtmax;
//     __m128i Dq;
//     __m128i Dqmax;
//     __m128i Tmp;
//     __m128i Tmp2;
//     __m128i One = _mm_set1_epi16(1);
//     __m128i vIter = _mm_setr_epi16(0,iter,2*iter,3*iter,4*iter,5*iter,6*iter,7*iter);
  
//     int num = 0;
//     int crossout_thresh = (int)(par.block_shading_space/1.3);
//     double factor = (double)dbsize * (double)queryLength * (double)dbLength;

//     int b,c,d,k,l;
//     int q_pos, t_pos, pos;
//     int d1, d2;
//     int t_start, q_start, t_end, q_end;
//     int kstart, lstart, h;

//     while (num < 10) {

//       // printf("Query-pos: ");
//       // for (b = 0; b < queryLength; b++) {
//       // 	printf("%3i ", b);
//       // }
//       // printf("\n");
//       // tmp = (short*) Hmatrix;
//       // for (b = 0; b < dbLength; b++) {
//       // 	printf("DB: %3i  : ", b);
//       // 	for (c = 0; c < 8; c++) {
//       // 	  for (d = 0; d < iter; d++) {
//       // 	  int i = b*iter*8 + d * 8 + c;
//       // 	  printf("%3i|", tmp[i]+SHORT_BIAS);
//       // 	  }
//       // 	}
//       // 	printf("\n");
//       // }
//       // printf("\n");
      
//       // Find maximum score with position
//       vMax = _mm_set1_epi16(-SHORT_BIAS);
//       Dtmax = _mm_setzero_si128();
//       Dqmax = _mm_setzero_si128();

//       Hmatrix_it = Hmatrix;
      
//       for (i = 0; i < dbLength; ++i)
// 	{
// 	  //printf("DB-pos: %3i\n",i);
// 	  Dq = vIter;
// 	  Dt = _mm_set1_epi16(i);
	  
// 	  for (j = 0; j < iter; ++j)
// 	    {
// 	      vH = _mm_load_si128(Hmatrix_it++);
	      
// 	      Tmp = _mm_cmpgt_epi16(vH, vMax);
// 	      vMax = _mm_max_epi16(vH, vMax);
	      
// 	      Tmp2 = _mm_and_si128(Tmp,Dt);
// 	      Dtmax = _mm_max_epi16(Dtmax, Tmp2);
	      
// 	      Dqmax = _mm_andnot_si128(Tmp,Dqmax);  // clear boxes, where new maximum is found
// 	      Tmp2 = _mm_and_si128(Tmp,Dq);
// 	      Dqmax = _mm_max_epi16(Dqmax, Tmp2);
	      
// 	      Dq = _mm_adds_epi16(Dq, One);
// 	    }
// 	}
      
//       // printf("Max and pos:\n");
//       // for (int b = 0; b<8; b++) // loop over bytes in XMM registers
//       // 	printf("Score: %3i  (q:%3i t:%3i) \n",*(((short*) &vMax) + b) + SHORT_BIAS,*(((short*) &Dqmax) + b),*(((short*) &Dtmax) + b));
//       // printf("\n");
      
//       int max_score = _mm_extract_epi16(vMax,0) + SHORT_BIAS;
//       int max_pos_q = _mm_extract_epi16(Dqmax,0);
//       int max_pos_t = _mm_extract_epi16(Dtmax,0);

//       for (b = 1; b<8; ++b) {
//       	if ((_mm_extract_epi16(vMax,b) + SHORT_BIAS) > max_score) {
//       	  max_score = _mm_extract_epi16(vMax,b) + SHORT_BIAS;
// 	  max_pos_q = _mm_extract_epi16(Dqmax,b);
// 	  max_pos_t = _mm_extract_epi16(Dtmax,b);
//       	}
//       }

//       double evalue = factor * fpow2(-max_score/par.prefilter_bit_factor);
      
//       //printf("Score: %5i  Evalue: %8.4g   (num: %2i)\n", max_score, evalue, num);

//       if (num == 0 || evalue < par.prefilter_smax_thresh) {
	
// 	// Backtrace alignment
// 	short *Btmatrix_sit = (short*) Btmatrix;
// 	short *Hmatrix_sit = (short*) Hmatrix;

// 	q_pos = max_pos_q;
// 	t_pos = max_pos_t;

// 	pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
// 	while (q_pos != 0 && t_pos != 0 && Hmatrix_sit[pos] != -SHORT_BIAS) {
// 	  if ((Btmatrix_sit[pos]&0x0001) != 0) {
// 	    --t_pos;
// 	    pos -= LQ;
// 	    while (q_pos != 0 && t_pos != 0 && Hmatrix_sit[pos] != -SHORT_BIAS) {
// 	      if ((Btmatrix_sit[pos]&0x0100) == 0) {
// 		--t_pos;
// 		pos -= LQ;
// 	      } else {
// 		--q_pos;
// 		--t_pos;
// 		break;
// 	      }
// 	    }
// 	  } else if ((Btmatrix_sit[pos]&0x0010) != 0) {
// 	    --q_pos;
// 	    pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
// 	    while (q_pos != 0 && t_pos != 0 && Hmatrix_sit[pos] != -SHORT_BIAS) {
// 	      if ((Btmatrix_sit[pos]&0x1000) == 0) {
// 		--q_pos;
// 		pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
// 	      } else {
// 		--q_pos;
// 		--t_pos;
// 		break;
// 	      }
// 	    }
// 	  } else {
// 	    --q_pos;
// 	    --t_pos;
// 	  } 
// 	  pos = t_pos * LQ + (q_pos%iter) * 8 + (q_pos/iter);
// 	}

// 	// Add 1 to each position because HHsearch counts sequences from position 1
// 	res[num].q_start = q_pos + 1;
// 	res[num].t_start = t_pos + 1;
// 	res[num].q_stop = imin(queryLength,max_pos_q + 1);
// 	res[num].t_stop = max_pos_t + 1;
// 	res[num++].evalue = evalue;

// 	// Cross out cells in H-matrix
// 	d1 = imin(t_pos-q_pos,max_pos_t-max_pos_q) - crossout_thresh;
// 	d2 = imax(t_pos-q_pos,max_pos_t-max_pos_q) + crossout_thresh;

// 	t_start = imax(0, t_pos - crossout_thresh);
// 	q_start = imax(0, q_pos - crossout_thresh);
// 	t_end = imin(dbLength, max_pos_t + crossout_thresh);
// 	q_end = imin(LQ, max_pos_q + crossout_thresh);

// 	kstart = (q_start % iter) * 8;
// 	lstart = q_start / iter;
// 	h = t_start*LQ;

// 	for (b = t_start; b < t_end; ++b) {
// 	  k = kstart;
// 	  l = lstart + h;
// 	  d = b - q_end;
// 	  for (c = b-q_start; c > d; --c) {
// 	    if (d1 < c && c < d2) {
// 	      Hmatrix_sit[k + l] = -SHORT_BIAS;
// 	    }
// 	    k+=8;
// 	    if (k == LQ) {
// 	      k = 0;
// 	      ++l;
// 	    }
// 	  }
// 	  h += LQ;
// 	}
// 	// for (int b = t_start; b < t_end; b++) {
// 	//   int h = b*iter*8;
// 	//   for (int c = q_start; c < q_end; c++) {
// 	//     if ((b-c) > d1 && (b-c) < d2) {
// 	//       Hmatrix_sit[h + (c%iter)*8 + (c/iter)] = -SHORT_BIAS;
// 	//     }
// 	//   }
// 	// }
//       } else {
// 	break;
//       }
//     }
    
//     free(Hmatrix);
//     free(Btmatrix);

//     return num;
// }

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

    __m128i* Hmatrix = (__m128i*)memalign(16,dbLength*LQ*sizeof(unsigned short));
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
    double factor = (double)dbsize * (double)queryLength * (double)dbLength;

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

      for (b = 1; b<8; ++b) {
      	if ((_mm_extract_epi16(vMax,b) + SHORT_BIAS) > max_score) {
      	  max_score = _mm_extract_epi16(vMax,b) + SHORT_BIAS;
	  max_pos_q = _mm_extract_epi16(Dqmax,b);
	  max_pos_t = _mm_extract_epi16(Dtmax,b);
      	}
      }

      double evalue = factor * fpow2(-max_score/par.prefilter_bit_factor);
      
      if (num == 0 || evalue < par.prefilter_smax_thresh) {
	
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


////////////////////////////////////////////////////////////////////////////////////////////////////////
// d = i-j+LT-1 is index of diagonal
void ungapped_sse2(unsigned char*       query_profile,
                   const int            query_length,
                   const unsigned char* db_sequence,
                   const int            dbseq_length,
                   const unsigned char  score_offset,
                   __m128i*             workspace,
	 	   unsigned char*       smax)
{
  int i; // position in query bands (0,..,W-1)
  int j; // position in db sequence (0,..,dbseq_length-1)
  int W = (query_length + 15) / 16; // width of bands in query and score matrix = hochgerundetes LQ/16
  
  __m128i *p;
  __m128i S;                       // 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15)
  __m128i Smax;                    // 16 unsigned bytes holding Smax(b*W+i,j) (b=0,..,15)
  __m128i Soffset;                 // all scores in query profile are shifted up by Soffset to obtain pos values 
  __m128i *s_prev, *s_curr;        // pointers to Score(i-1,j-1) and Score(i,j), resp.
  __m128i *smax_prev, *smax_curr;  // pointers to Score_max(i-1,j-1) and Score_max(i,j), resp.
  __m128i *qji;                    // query profile score in row j (for residue x_j)
  __m128i *s_prev_it, *s_curr_it;       
  __m128i *smax_prev_it, *smax_curr_it; 
  __m128i *query_profile_it = (__m128i *) query_profile;
  __m128i Zero = _mm_setzero_si128();

  unsigned char* smax_store = smax + (16*W) + dbseq_length - 2;

  // Load the score offset to all 16 unsigned byte elements of Soffset
  Soffset = _mm_set1_epi8(score_offset);
  
  // Initialize  workspace to zero 
  for (i=0, p=workspace; i < 4*W; ++i) 
    _mm_store_si128(p++, Zero);
  
  s_curr = workspace;
  s_prev = workspace + W;
  smax_curr = workspace + 2*W;
  smax_prev = workspace + 3*W;

  for (j=0; j<dbseq_length; ++j) // loop over db sequence positions
    {
      // Get address of query scores for row j 
      qji = query_profile_it + db_sequence[j]*W;

      // Load the next S value
      S = _mm_load_si128(s_curr + W - 1);
      S = _mm_slli_si128(S, 1);

      // Load the next Smax value
      Smax = _mm_load_si128(smax_curr + W - 1);
      Smax = _mm_slli_si128(Smax, 1);
            
      // Swap s_prev and s_curr, smax_prev and smax_curr
      SWAP(p,s_prev,s_curr);
      SWAP(p,smax_prev,smax_curr);
      
      s_curr_it = s_curr;
      s_prev_it = s_prev;
      smax_curr_it = smax_curr;
      smax_prev_it = smax_prev;

      for (i=0; i<W; ++i) // loop over query band positions
        {
	  // Saturated addition and subtraction to score S(i,j)
	  S = _mm_adds_epu8(S, *(qji++));        // S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset) 
	  S = _mm_subs_epu8(S, Soffset);         // S(i,j) = max(0, S(i,j) - Soffset)
	  _mm_store_si128(s_curr_it++, S);       // store S to s_curr[i]

	  // Update highest scores in diagonal encountered so far
	  Smax = _mm_max_epu8(Smax, S);          // Smax(i,j) = max(Smax(i,j), S(i,j))
	  _mm_store_si128(smax_curr_it++, Smax); // store Smax to Smax[i]
	  
	  // Load the next S and Smax values 
	  S = _mm_load_si128(s_prev_it++);    
	  Smax = _mm_load_si128(smax_prev_it++); 
        }
      
      // Store highest byte of smax_curr at pos smax_store
      *(smax_store--) = *(((char*) (smax_curr+W))-1); 
    }
  
  // Store Smax values for negative diagonals 
  char* smax_it = (char*) smax_curr;
  for (i=0; i<W; i++)
    {
      smax_store = smax+i;
      for (int b = 0; b<16; b++) // loop over bytes in XMM registers
	{
	  *smax_store = *(smax_it++); 
	  smax_store+=W;
	}
    }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// d = i-j+LT-1 is index of diagonal
int ungapped_sse_score(const unsigned char*       query_profile,
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

int identify_region(const int LQ, const int LDB, unsigned char* smax, triple* res)
{
  __m128i G = _mm_set1_epi16(par.prefilter_gap_open);
  __m128i H = _mm_set1_epi16(-par.prefilter_gap_extend);
  __m128i S0;
  __m128i* smax_tmp;
  __m128i sum = _mm_setzero_si128();
  
  short* smax_stripes = (short*)memalign(16,(LQ+LDB+15)*sizeof(short));
  int W = (LQ+LDB-1+7)/8;
  
  // Store Smax in stripes
  int a = 0;
  for (int i=0; i < W; i++) 
    {
      int j = i;
      for (int k = 0; k < 8; k++)
	{
	  if (j >= (LQ+LDB-1))
	    smax_stripes[a]=0;
	  else
	    smax_stripes[a]=(short) smax[j];
	  a++;
	  j+=W;
	}
    }

  // Calculate average smax
  smax_tmp = (__m128i *) smax_stripes;

  for (int d=0; d<W; d++)
    sum = _mm_adds_epu16(sum, *(smax_tmp++));

  sum = _mm_hadds_epi16(sum,sum);
  sum = _mm_hadds_epi16(sum,sum);
  sum = _mm_hadds_epi16(sum,sum);

  S0 = _mm_set1_epi16(*((short*) &sum) / (LQ+LDB-1));
  S0 = _mm_adds_epi16(S0,G);

  // Subtract S0 from Smax and store max(-h, smax - smax_avg - g)
  smax_tmp = (__m128i *) smax_stripes;
  for (int d=0; d<W; d++)
    {
      sum = _mm_subs_epi16(*smax_tmp, S0);
      sum = _mm_max_epi16(sum, H);
      _mm_store_si128(smax_tmp++,sum);
    }

  // Calculate R(d) in all segments parallel
  __m128i R = _mm_setzero_si128();
  __m128i R_prev = _mm_setzero_si128();
  __m128i R_max = _mm_setzero_si128();  
  __m128i D_start = _mm_setzero_si128();
  __m128i D_stop = _mm_setzero_si128();
  __m128i D;
  __m128i D_tmp = _mm_setzero_si128();
  __m128i Tmp = _mm_setzero_si128();;
  __m128i One = _mm_set1_epi16(1);
  uint64_t stop_when_0[] = {1,1};
  __m128i* stop_when_0_ref = (__m128i *) stop_when_0;

  int max_it = 6;

  for (int it = 0; it <= max_it; it++) 
    {
      D = _mm_setr_epi16(0,W,2*W,3*W,4*W,5*W,6*W,7*W);
      D = _mm_subs_epi16(D, One);
      smax_tmp = (__m128i *) smax_stripes;
      for (int d = 0; d < W; d++) 
	{
	  D = _mm_adds_epi16(D, One);
	  
	  // Whenever R[b]==0, set D_tmp[b] = pos
	  Tmp = _mm_cmpeq_epi16(R, _mm_setzero_si128());
	  Tmp = _mm_and_si128(Tmp,D);
	  D_tmp = _mm_max_epi16(D_tmp,Tmp);

	  // Calculate new R(d) = max{0, R(d-1)-h, R(d-1) + Smax(d) - <Smax> - g}
	  Tmp = _mm_load_si128(smax_tmp++);
	  R = _mm_adds_epi16(R, Tmp);
	  R = _mm_max_epi16(R, _mm_setzero_si128());

	  // Whenever R[b] > R_max[b], set D_stop[b] = pos
	  Tmp = _mm_cmpgt_epi16(R, R_max);
	  Tmp = _mm_and_si128(Tmp,D);
	  D_stop = _mm_max_epi16(D_stop, Tmp);
	  Tmp = _mm_cmpgt_epi16(R, R_max);
	  Tmp = _mm_and_si128(Tmp,D_tmp);
	  D_start = _mm_max_epi16(D_start, Tmp);
	  
	  R_max = _mm_max_epi16(R_max, R);
	}

      Tmp = _mm_cmpgt_epi16(R,R_prev);
      _mm_store_si128(stop_when_0_ref, Tmp);
      if (it == max_it || (stop_when_0[0] == 0 && stop_when_0[1] == 0))
	break;

      // Shift register
      R_prev = R;
      R = _mm_slli_si128(R, 2);
      D_start = _mm_slli_si128(D_start, 2);
      D_tmp = _mm_setzero_si128();
      R_max = _mm_setzero_si128();

    }

  // Return all start- and end-positions, where R_max is above threshold
  int num = 0;
  for (int i = 0; i < 8; i++) {
    if ((int)*(((short*) &R_max)+i) > par.prefilter_rmax_thresh) {
      res[num].first = (int)*(((short*) &D_start)+i);
      res[num].second = (int)*(((short*) &D_stop)+i);
      res[num++].third = (int)*(((short*) &R_max)+i);    
    }
  }

  free(smax_stripes);

  return num;
}

void init_prefilter()
{
   // Get Prefilter Pvalue (Evalue / DBsize)
  FILE *stream;
  command = (string)blast + "/fastacmd -d " + (string)db + " -I";
  stream = popen(command.c_str(), "r");
  while (fgets(line, LINELEN, stream))
    {
      char tmp_number[NAMELEN];
      char tmp_name[NAMELEN];
      ptr=strscn(line); 
      ptr=strwrd(tmp_number,ptr);
      ptr=strwrd(tmp_name,ptr);
      if (!strcmp(tmp_name,"sequences;")) {
	strtrd(tmp_number,",");
	dbsize = atoi(tmp_number);
	ptr=strwrd(tmp_number,ptr);
	strtrd(tmp_number,",");
	LDB = atoi(tmp_number);
	break;
      }
    }
  pclose(stream);
  if (dbsize == 0 || LDB == 0)
    {cerr<<endl<<"Error! Could not determine DB-size of prefilter db ("<<db<<")\n"; exit(4);}
  par.hhblits_prefilter_logpval=-log(e_psi / (float)dbsize);

  if (!(!strcmp(pre_mode,"csblast") || !strcmp(pre_mode,"blast")))
    {
      printf("Prefilter DB with %6i sequences with %s ...\n",dbsize,pre_mode);

      X = (unsigned char*)memalign(16,LDB*sizeof(unsigned char));                             // database string (concatenate all DB-seqs
      first = (unsigned char**)memalign(16,(2*dbsize)*sizeof(unsigned char*));                // first characters of db sequences
      length = (int*)memalign(16,(2*dbsize)*sizeof(int));                                     // lengths of db sequences
      dbnames = new char*[dbsize*2];

      /////////////////////////////////////////
      // Read in database
      num_dbs = 0;
      int len = 0;
      int pos = 0;
      int blocknum = 0;
      char word[NAMELEN];
      FILE* dbf = NULL;
      dbf = fopen(db,"rb");
      if (!dbf) OpenFileError(db);
      while(fgetline(line,LINELEN,dbf)) // read HMM files in pal file
	{
	  if (line[0]=='>')
	    {
	      if (len > 0)           // if it is not the first sequence
		length[num_dbs++] = len;
	      len = 0;
	      blocknum = 0;
	      
	      strwrd(word,line+1);
	      dbnames[num_dbs]=new(char[strlen(word)+5]);
	      strcpy(dbnames[num_dbs],word);
	      strcat(dbnames[num_dbs]," 0");
	      
	      first[num_dbs] = X + pos;
	    }
	  else
	    {
	      int h = 0;
	      while (h<LINELEN && line[h]>'\0')
		{
		  if (aa2i(line[h])>=0) // ignore white-space characters ' ', \t and \n (aa2i()==-1)
		    {
		      X[pos++]=(unsigned char)(aa2i(line[h])); //  AS62::kCharToInt[line[h]]
		      len++;
		    }
		  else 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" of "<<db<<"\n";
		  h++;
		  if (len == par.prefilter_lmax && (!strcmp(pre_mode,"only_prefilt_ungapped_SW") || !strcmp(pre_mode,"ungapped_SW") || !strcmp(pre_mode,"only_prefilt_ungapped_SW_no_region") || !strcmp(pre_mode,"ungapped_SW_no_region")))
		    {
		      blocknum++;
		      length[num_dbs++] = len;
		      len = par.prefilter_db_overlap;
		      first[num_dbs] = X + pos - par.prefilter_db_overlap;
		      dbnames[num_dbs]=new(char[strlen(word)+10]);
		      stringstream ss_tmp;
		      ss_tmp << word << " " << (blocknum*(par.prefilter_lmax-par.prefilter_db_overlap));
		      strcpy(dbnames[num_dbs],ss_tmp.str().c_str());
		    }
		}
	    }
	}
      if (len > 0)
	length[num_dbs++] = len;
      
      fclose(dbf);
      
    }
}

void init_sse_prefiltering()
{
  LQ=q.L;
  qc = (unsigned char*)memalign(16,(par.prefilter_states+1)*(LQ+15)*sizeof(unsigned char));   // query profile
  W = (LQ+15) / 16;   // band width = hochgerundetes LQ/16
  
  // Add Pseudocounts
  if (!*par.clusterfile) { //compute context-specific pseudocounts?
    // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
    q_tmp.PreparePseudocounts();
    // Add amino acid pseudocounts to query: p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    q_tmp.AddAminoAcidPseudocounts(2,1.5,2,1);
  } else {
    // Add context specific pseudocounts
    q_tmp.AddContextSpecificPseudocounts(2,1.5,2,1);
  }
  
  q_tmp.CalculateAminoAcidBackground();

  // Divide by NullModel
  q_tmp.IncludeNullModelInHMM(q_tmp,q_tmp,1);
  //q_tmp.IncludeNullModelInHMM(q_tmp,q_tmp,0);
        
  // printf("\n\nQuery profile:\n        ");
  // for (int j=1; j <= LQ; j++)
  //   printf("%4i  ", j);
  // printf("\n");
  // for (int a=0; a < par.prefilter_states; a++)  
  //   {
  //     printf("a=%3i   ",a);
  //     for (int j=1; j <= LQ; j++)
  // 	{
  // 	  printf("%5.2f|",flog2(q_tmp.p[j][a])*2.0);
  // 	}
  //     printf("\n");
  //   }

  // printf("\n\nQuery profile:\n      ");
  // for (int j=1; j <= LQ; j++)
  //   printf(" %c|", q_tmp.seq[0][j]);
  // printf("\n");
  // for (int a=0; a < par.prefilter_states; a++)  
  //   {
  //     printf("a=%c   ",i2aa(a));
  //     for (int j=1; j <= LQ; j++)
  // 	{
  // 	  float dummy = flog2(q_tmp.p[j][a])*par.prefilter_bit_factor + par.prefilter_score_offset+0.5;
  // 	  if (dummy>255.99) dummy = 255.5;
  // 	  if (dummy<0) dummy = 0.0;
  // 	  unsigned char c = (unsigned char) dummy;  // 1/3 bits & make scores >=0 everywhere
  // 	  printf("%2i|",(int)c);
  // 	}
  //     printf("\n");
  //   }

  // /////////////////////////////////////////
  // // Stripe query profile
  for (int a=0; a < par.prefilter_states; a++)  
    {
      int h = a*W*16;
      for (int i=0; i < W; i++) 
  	{
  	  int j = i;
  	  for (int k = 0; k < 16; k++)
  	    {
  	      if (j >= LQ)
  		qc[h]=(unsigned char) par.prefilter_score_offset;
  	      else
  		{
  		  float dummy = flog2(q_tmp.p[j+1][a])*par.prefilter_bit_factor + par.prefilter_score_offset + 0.5;
  		  if (dummy>255.99) dummy = 255.5;
  		  if (dummy<0) dummy = 0.0;
  		  qc[h] = (unsigned char) dummy;  // 1/3 bits & make scores >=0 everywhere
  		}
  	      h++;
  	      j+=W;
  	    }
  	}
   }

  // For K=20 alphabet add extra X-state
  int h = par.prefilter_states*W*16;
  for (int i=0; i < 16*W; i++)
    {
      if (i >= LQ)
  	qc[h]=(unsigned char) par.prefilter_score_offset;
      else
  	qc[h]=(unsigned char) (par.prefilter_score_offset - 1);
      h++;
    }

  if (!strcmp(pre_mode,"SW_evalue_preprefilter_backtrace"))
    {
      qw = (unsigned short*)memalign(16,(par.prefilter_states+1)*(LQ+7)*sizeof(unsigned short));   // query profile
      Ww = (LQ+7) / 8;

      // /////////////////////////////////////////
      // // Stripe query profile
      for (int a=0; a < par.prefilter_states; a++)  
	{
	  int h = a*Ww*8;
	  for (int i=0; i < Ww; i++) 
	    {
	      int j = i;
	      for (int k = 0; k < 8; k++)
		{
		  if (j >= LQ)
		    qw[h] = 0;
		  else
		    {
		      float dummy = flog2(q_tmp.p[j+1][a])*par.prefilter_bit_factor;
		      qw[h] = (unsigned short) dummy;  // 1/3 bits & make scores >=0 everywhere
		    }
		  h++;
		  j+=Ww;
		}
	    }
	}
      
      // For K=20 alphabet add extra X-state
      int h = par.prefilter_states*Ww*8;
      for (int i=0; i < 8*Ww; i++)
	{
	  if (i >= LQ)
	    qw[h]=0;
	  else
	    qw[h]= (unsigned short) -1;
	  h++;
	}
    }
  // int matrix[500];
  // FILE* inf = NULL;
  // inf = fopen("/net/cluster/user/michael/tmp/sse_test/farrar/multi-threaded/blosum62.mat","rb");
  // if (!inf) OpenFileError("/net/cluster/user/michael/tmp/sse_test/farrar/multi-threaded/blosum62.mat");
  // while(fgetline(line,LINELEN,inf))
  //   {
  //     if (line[0] != '#')
  // 	break;
  //   }
  // int num = 0;
  // while(fgetline(line,LINELEN,inf))
  //   {
  //     ptr = line;
  //     for (int i = 0; i < 20; i++)
  // 	{
  // 	  matrix[num*20 + i] = strint(ptr);
  // 	}
  //     num++;
  //   }
  // fclose(inf);

  // printf("Matrix:\n");
  // for (int i = 0; i < 20; i++) 
  //   printf("%2i|", i);
  // printf("\n");
  // for (int i = 0; i < 20; i++) 
  //   printf("---");
  // printf("\n");
  // for (int i = 0; i < 20 * 20; i++) 
  //   {
  //     if ((i % 20) == 0)
  // 	printf("\n");
  //     printf("%2i|", matrix[i]);
  //   }


  // /////////////////////////////////////////
  // // Stripe query profile
  // for (int a=0; a < par.prefilter_states; a++)  
  //   {
  //     int h = a*W*16;
  //     for (int i=0; i < W; i++) 
  // 	{
  // 	  int j = i;
  // 	  for (int k = 0; k < 16; k++)
  // 	    {
  // 	      if (j >= LQ)
  // 		//qc[h]=0;
  // 		qc[h]=par.prefilter_score_offset;
  // 	      else
  // 		{
  // 		  int dummy = matrix[aa2i(q_tmp.seq[0][j+1])*20+a] + par.prefilter_score_offset;
  // 		  qc[h] = (unsigned char) dummy;
  // 		}
  // 	      h++;
  // 	      j+=W;
  // 	    }
  // 	}
  //   }

  // // For K=20 alphabet add extra X-state
  // int h = par.prefilter_states*W*16;
  // for (int i=0; i < 16*W; i++)
  //   {
  //     if (i >= LQ)
  // 	qc[h]=(unsigned char) par.prefilter_score_offset;
  //     else
  // 	qc[h]=(unsigned char) (par.prefilter_score_offset - 1);
  //     h++;
  //   }

  // printf("\n\nStriped Query profile:\n      ");
  // for (int j = 0; j < W; ++j) {
  //   for (int k = j; k < W*16; k += W) {
  //     if (k < LQ)
  // 	printf(" %c|", q_tmp.seq[0][k+1]);
  //     else
  // 	printf(" -|");
  //   }
  // }
  // printf("\n");
  // for (int a=0; a < par.prefilter_states; a++)  
  //   {
  //     printf("a=%c   ",i2aa(a));
  //     for (int j=0; j < W*16; j++)
  // 	{
  // 	  printf("%2i|",qc[a*W*16+j]);
  // 	}
  //     printf("\n");
  //   }
}

void prefilter_with_ungapped_SW()
{
  strcpy(par.block_shading_mode,"SSE");
  
  init_sse_prefiltering();

  unsigned char* smax = (unsigned char*)memalign(16,(LQ+par.prefilter_lmax+15)*sizeof(unsigned char));
  for (int a = 0; a < (LQ+par.prefilter_lmax+15); a++)
    smax[a] = 0;
  __m128i* workspace = (__m128i*)memalign(16,4*(LQ+15)*sizeof(char));
  
  triple* res = new triple[8];
  int num_res;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];
  int block_start;
  
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      // Perform search step
      ungapped_sse2 (qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace, smax);
      
      // identify region with SSE 2 when one smax[d] > smax_thresh
      bool check = false;
      for(int d=0; d<LQ+length[n]-1; d++)
	if (smax[d] > par.prefilter_smax_thresh)
	  {
	    check = true;
	    break;
	  }
      if (check) {
	num_res = identify_region(LQ, length[n], smax, res);
	if (num_res > 0) {
	  
	  ptr=strwrd(tmp_name,dbnames[n]);
	  block_start = strint(ptr);
	  
	  if (block_filter)
	    {
	      // Write block for template
	      if (strcmp(actual_hit,"") && strcmp(actual_hit,tmp_name))   // New template
		{
		  // printf("Add to block shading   key: %s    data:",actual_hit);
		  // for (int i = 0; i < block_count; i++)
		  //   printf(" %i,",block[i]);
		  // printf("\n");
		  par.block_shading->Add(actual_hit,block);
		  par.block_shading_counter->Add(actual_hit,block_count);
		  block = new(int[400]);
		  block_count = 0;
		}
	      if (block_count >= 400) { continue; }
	      strcpy(actual_hit,tmp_name);
	      
	      // Extract only best hit for each region
	      int beg = res[0].first;
	      int end = res[0].second;
	      int rmax = res[0].third;
	      for (int i = 1; i < num_res; i++) 
		{
		  if (res[i].first == beg)
		    {
		      if (res[i].third > rmax)
			{
			  rmax = res[i].third;
			  end = res[i].second;
			}
		    }
		  else
		    {
		      //printf("%4i-%4i (i-j: %4i - %4i), %4i |", beg, end, beg-length[n]+1, end-length[n]+1, rmax);
		      block[block_count++]=block_start;
		      block[block_count++]=beg-length[n]+1;  // d1
		      block[block_count++]=end-length[n]+1;  // d2
		      beg = res[i].first;
		      end = res[i].second;
		      rmax = res[i].third;
		    }
		}
	      //printf("%4i-%4i (i-j: %4i - %4i), %4i |", beg, end, beg-length[n]+1, end-length[n]+1, rmax);
	      block[block_count++]=block_start;
	      block[block_count++]=beg-length[n]+1;  // d1
	      block[block_count++]=end-length[n]+1;  // d2
	    }
	  
	  if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	    {
	      substr(tmp,tmp_name,3,11);
	      substr(db_name,tmp,0,1);
	      strcat(db_name,"/");
	      strcat(db_name,tmp);
	      strcat(db_name,".db");
	    }
	  else                              // other database
	    {
	      strcpy(db_name,tmp_name);
	      strtr(db_name,"|", "_");
	      strcat(db_name,".hhm");
	    }
	  
	  if (! doubled->Contains(db_name))
	    {
	      doubled->Add(db_name);
	      // check, if DB was searched in previous rounds 
	      strcat(tmp_name,"__1");  // irep=1
	      if (previous_hits->Contains(tmp_name))
		{
		  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
		  strcpy(dbfiles_old[ndb_old],dbhhm);
		  strcat(dbfiles_old[ndb_old],"/");
		  strcat(dbfiles_old[ndb_old],db_name);
		  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
		  ndb_old++;
		}
	      else 
		{
		  dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
		  strcpy(dbfiles_new[ndb_new],dbhhm);
		  strcat(dbfiles_new[ndb_new],"/");
		  strcat(dbfiles_new[ndb_new],db_name);
		  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
		  ndb_new++;
		}
	    }
	}
      }
    }
  
  if (block_filter && strcmp(actual_hit,""))   // New template
    {
      // printf("Add to block shading   key: %s    data:",actual_hit);
      // for (int i = 0; i < block_count; i++)
      //   printf(" %i,",block[i]);
      // printf("\n");
      par.block_shading->Add(actual_hit,block);
      par.block_shading_counter->Add(actual_hit,block_count);
    }
  
  // Free memory
  free(qc);
  free(smax);
  free(workspace);
  delete[](res);
}

void prefilter_with_ungapped_SW_only()
{
  init_sse_prefiltering();

  unsigned char* smax = (unsigned char*)memalign(16,(LQ+par.prefilter_lmax+15)*sizeof(unsigned char));
  for (int a = 0; a < (LQ+par.prefilter_lmax+15); a++)
    smax[a] = 0;
  __m128i* workspace = (__m128i*)memalign(16,4*(LQ+15)*sizeof(char));
  
  triple* res = new triple[8];
  int num_res;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];

  char old_db_name[NAMELEN];
  strcpy(old_db_name, "");
  int score = 0;

  FILE* outf=NULL;
  outf=fopen(par.outfile,"w");
  fprintf(outf,"        Name        Score\n");
  
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      // Perform search step
      ungapped_sse2 (qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace, smax);
      
      // identify region with SSE 2 when one smax[d] > smax_thresh
      bool check = false;
      for(int d=0; d<LQ+length[n]-1; d++)
	if (smax[d] > par.prefilter_smax_thresh)
	  {
	    check = true;
	    break;
	  }
      if (check) {
	num_res = identify_region(LQ, length[n], smax, res);
	if (num_res > 0) {

	  ptr=strwrd(tmp_name,dbnames[n]);
	  
	  if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	    substr(db_name,tmp_name,3,11);
	  else                              // other database
	    strcpy(db_name,tmp_name);

	  if (strcmp(db_name,old_db_name) && strcmp(old_db_name,""))
	    {
	      fprintf(outf,"%-20s   %4i\n", old_db_name, score);
	      score = 0;
	    }
	  
	  strcpy(old_db_name,db_name);

	  for (int i = 0; i < num_res; i++)
	    {
	      if (res[i].third > score)
		score = res[i].third;
	    }
	  
	}
      }
    }

  if (strcmp(db_name,old_db_name) && strcmp(old_db_name,""))
    fprintf(outf,"%-20s   %4i\n", old_db_name, score);
	  
  fclose(outf);
  printf("Scorefile of prefiltering-hits written to %s!\n",par.outfile);
  
  // Free memory
  free(qc);
  free(smax);
  free(workspace);
  delete[](res);
}

void prefilter_with_ungapped_SW_no_region()
{
  init_sse_prefiltering();
  
  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);
  
  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,2*(LQ+15)*sizeof(char));
  
  int score;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<int, string> > hits;

  #pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif
  
      // Perform search step
      score = ungapped_sse_score (qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace[thread_id]);
      
      if (score > par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<int,string>(score, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());
  
  vector<pair<int, string> >::reverse_iterator it;
  
  for ( it=hits.rbegin() ; it < hits.rend(); it++ )
    {
      //printf("Hit %20s with score %4i\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());
      
      ptr=strwrd(tmp_name,tmp);
      
      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      if (! doubled->Contains(db_name))
	{
	  doubled->Add(db_name);

	  // check, if DB was searched in previous rounds 
	  strcat(tmp_name,"__1");  // irep=1
	  if (previous_hits->Contains(tmp_name))
	    {
	      dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	      strcpy(dbfiles_old[ndb_old],dbhhm);
	      strcat(dbfiles_old[ndb_old],"/");
	      strcat(dbfiles_old[ndb_old],db_name);
	      if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	      ndb_old++;
	    }
	  else 
	    {
	      dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	      strcpy(dbfiles_new[ndb_new],dbhhm);
	      strcat(dbfiles_new[ndb_new],"/");
	      strcat(dbfiles_new[ndb_new],db_name);
	      if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	      ndb_new++;
	    }
	}

      if (ndb_new >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
	  break;
	  //exit(4);
	}
    }

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}

void prefilter_with_SW_evalue_preprefilter()
{
  int* prefiltered_hits = new int[1000000];

  init_sse_prefiltering();
  
  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);
  
  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int score;
  int thread_id = 0;
  int count_dbs = 0;

  float log_qlen = flog2(LQ);

  #pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif
  
      // Perform search step
      score = ungapped_sse_score (qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace[thread_id]);
  
      score = score - par.prefilter_bit_factor * (log_qlen + flog2(length[n]));
    
      if (score > par.preprefilter_smax_thresh)
	{
          #pragma omp critical
	  prefiltered_hits[count_dbs++] = n;
	}
    }
  printf("%6i hits through preprefilter!\n", count_dbs);
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  double evalue;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<double, string> > hits;

  const double factor = (double)dbsize * LQ;

  #pragma omp parallel for schedule(static) private(evalue, score, thread_id)
  for (int n = 0; n < count_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[prefiltered_hits[n]], length[prefiltered_hits[n]], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      evalue = factor * length[prefiltered_hits[n]] * fpow2(-score/par.prefilter_bit_factor);
 
      if (evalue < par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<double,string>(evalue, string(dbnames[prefiltered_hits[n]])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<double, string> >::iterator it;
  
  for ( it=hits.begin() ; it < hits.end(); it++ )
    {
      //printf("Hit %20s with score %6.2g\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      // check, if DB was searched in previous rounds 
      strcat(tmp_name,"__1");  // irep=1
      if (previous_hits->Contains(tmp_name))
	{
	  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_old[ndb_old],dbhhm);
	  strcat(dbfiles_old[ndb_old],"/");
	  strcat(dbfiles_old[ndb_old],db_name);
	  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	  ndb_old++;
	}
      else 
	{
	 dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_new[ndb_new],dbhhm);
	  strcat(dbfiles_new[ndb_new],"/");
	  strcat(dbfiles_new[ndb_new],db_name);
	  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	  ndb_new++;
	}
      
      if (ndb_new >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
	  break;
	  //exit(4);
	}
    }

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
  delete[] prefiltered_hits;
}

void prefilter_with_SW_evalue_preprefilter_backtrace()
{
  int* prefiltered_hits = new int[1000000];
  int* backtrace_hits = new int[MAXNUMDB+1];

  init_sse_prefiltering();
  
  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);
  
  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int score;
  int thread_id = 0;
  int count_dbs = 0;

  float log_qlen = flog2(LQ);

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
  printf("%6i hits through preprefilter!\n", count_dbs);
  if (print_elapsed) ElapsedTimeSinceLastCall("(ungapped preprefilter)");
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  double evalue;

  const double factor = (double)dbsize * LQ;

  vector<pair<double, int> > hits;

  #pragma omp parallel for schedule(static) private(evalue, score, thread_id)
  for (int n = 0; n < count_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[prefiltered_hits[n]], length[prefiltered_hits[n]], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      evalue = factor * length[prefiltered_hits[n]] * fpow2(-score/par.prefilter_bit_factor);
 
      if (evalue < par.prefilter_smax_thresh)
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
      char tmp[NAMELEN];
      ptr=strwrd(tmp_name,dbnames[(*it).second]);
      
      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      if (! doubled->Contains(db_name))
	{
	  doubled->Add(db_name);
	  // check, if DB was searched in previous rounds 
	  strcat(tmp_name,"__1");  // irep=1
	  if (previous_hits->Contains(tmp_name))
	    {
	      dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	      strcpy(dbfiles_old[ndb_old],dbhhm);
	      strcat(dbfiles_old[ndb_old],"/");
	      strcat(dbfiles_old[ndb_old],db_name);
	      if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	      ndb_old++;
	    }
	  else 
	    {
	      dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	      strcpy(dbfiles_new[ndb_new],dbhhm);
	      strcat(dbfiles_new[ndb_new],"/");
	      strcat(dbfiles_new[ndb_new],db_name);
	      if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	      ndb_new++;
	    }
	}


      if (count_dbs >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
	  break;
	}
    }

  if (print_elapsed) ElapsedTimeSinceLastCall("(SW prefilter)");

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
      //if (1)
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
		// printf("Add to block shading   key: %s    data:",tmp_name);
		// for (int i = 0; i < block_count; i++)
		//   printf(" %i,",block[i]);
		// printf("\n");
		par.block_shading->Add(tmp_name,block);
		par.block_shading_counter->Add(tmp_name,block_count);
	      }
	    }
	  delete[] res;
	}
    }
  if (print_elapsed) ElapsedTimeSinceLastCall("(SW backtrace prefilter)");

  // Free memory
  free(qc);
  free(qw);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
  delete[] prefiltered_hits;
  delete[] backtrace_hits;
}

// void prefilter_with_ungapped_SW_no_region()
// {
//   int num_dbs = init_sse_prefiltering();
  
//   __m128i* workspace = (__m128i*)memalign(16,2*(LQ+15)*sizeof(char));
  
//   char db_name[NAMELEN];
//   char tmp_name[NAMELEN];
//   char tmp[NAMELEN];
  
//   char old_db_name[NAMELEN];
//   strcpy(old_db_name, "");
//   int score = 0;
  
//   for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
//     {
//       // Perform search step
//       int tmp_score = ungapped_sse_score (qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace);
      
//       if (tmp_score > par.prefilter_smax_thresh)
// 	{
// 	  ptr=strwrd(tmp_name,dbnames[n]);
	  
// 	  if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
// 	    {
// 	      substr(tmp,tmp_name,3,11);
// 	      substr(db_name,tmp,0,1);
// 	      strcat(db_name,"/");
// 	      strcat(db_name,tmp);
// 	      strcat(db_name,".db");
// 	    }
// 	  else                              // other database
// 	    {
// 	      strcpy(db_name,tmp_name);
// 	      strtr(db_name,"|", "_");
// 	      strcat(db_name,".hhm");
// 	    }
	  
// 	  if (strcmp(db_name,old_db_name) && strcmp(old_db_name,""))
// 	    {
// 	      // check, if DB was searched in previous rounds 
// 	      strcat(tmp_name,"__1");  // irep=1
// 	      if (previous_hits->Contains(tmp_name))
// 		{
// 		  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(old_db_name)+2]);
// 		  strcpy(dbfiles_old[ndb_old],dbhhm);
// 		  strcat(dbfiles_old[ndb_old],"/");
// 		  strcat(dbfiles_old[ndb_old],old_db_name);
// 		  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
// 		  ndb_old++;
// 		}
// 	      else 
// 		{
// 		  dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(old_db_name)+2]);
// 		  strcpy(dbfiles_new[ndb_new],dbhhm);
// 		  strcat(dbfiles_new[ndb_new],"/");
// 		  strcat(dbfiles_new[ndb_new],old_db_name);
// 		  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
// 		  ndb_new++;
// 		}
// 	      score = 0;
// 	    }
// 	  if (ndb_new >= MAXNUMDB) 
// 	    {
// 	      printf("\nERROR! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
// 	      exit(4);
// 	    }
	  
// 	  strcpy(old_db_name,db_name);
	  
// 	  if (tmp_score > score)
// 	    score = tmp_score;
// 	}
//     }
  
//   if (strcmp(db_name,old_db_name) && strcmp(old_db_name,""))
//     {
//       // check, if DB was searched in previous rounds 
//       strcat(tmp_name,"__1");  // irep=1
//       if (previous_hits->Contains(tmp_name))
// 	{
// 	  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(old_db_name)+2]);
// 	  strcpy(dbfiles_old[ndb_old],dbhhm);
// 	  strcat(dbfiles_old[ndb_old],"/");
// 	  strcat(dbfiles_old[ndb_old],old_db_name);
// 	  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
// 	  ndb_old++;
// 	}
//       else 
// 	{
// 	  dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(old_db_name)+2]);
// 	  strcpy(dbfiles_new[ndb_new],dbhhm);
// 	  strcat(dbfiles_new[ndb_new],"/");
// 	  strcat(dbfiles_new[ndb_new],old_db_name);
// 	  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
// 	  ndb_new++;
// 	}
//     }
  
//   // Free memory
//   free(qc);
//   free(X);
//   free(length);
//   free(first);
//   free(workspace);
//   for (int n = 0; n < num_dbs; n++)
//     delete[](dbnames[n]);
//   delete[](dbnames);
// }

void prefilter_with_ungapped_SW_no_region_only()
{
  init_sse_prefiltering();

  __m128i** workspace = new(__m128i*[cpu]);
  
  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,2*(LQ+15)*sizeof(char));
  
  int score;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<int, string> > hits;

  #pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif
  
      // Perform search step
      score = ungapped_sse_score (qc, LQ, first[n], length[n], par.prefilter_score_offset, workspace[thread_id]);
      
      #pragma omp critical
      hits.push_back(pair<int,string>(score, string(dbnames[n])));
    }

  sort(hits.begin(), hits.end());
  
  vector<pair<int, string> >::reverse_iterator it;
  
  FILE* outf=NULL;
  outf=fopen(par.outfile,"w");
  fprintf(outf,"        Name        Score\n");

  for ( it=hits.rbegin() ; it < hits.rend(); it++ )
    {
      //printf("Hit %20s with score %4i\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());
      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	substr(db_name,tmp_name,3,11);
      else                              // other database
	strcpy(db_name,tmp_name);
      
      if (! doubled->Contains(db_name))
	{
	  doubled->Add(db_name);
	  fprintf(outf,"%-20s   %4i\n", db_name, (*it).first);
	}
    }

  fclose(outf);
  printf("Scorefile of prefiltering-hits written to %s!\n",par.outfile);
  
  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}


void prefilter_with_BLAST()
{
  if (v>=2) printf("Pre-filtering with %s ...\n",(strcmp(pre_mode,"csblast"))?"PSI-BLAST":"CS-BLAST");
  stringstream ss;
  if (!is_regular_file(tmp_psifile.c_str())) 
    {
      if (strcmp(pre_mode,"csblast")) 
	ss << blast << "/blastpgp -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 2> /dev/null";
      else
	ss << csblast << "/csblast -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -D " << csblast_db << " --blast-path " << blast << " --no-penalty 2> /dev/null";
    } 
  else
    {
      if (strcmp(pre_mode,"csblast")) 
	ss << blast << "/blastpgp -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -B " << tmp_psifile << " 2> /dev/null";
      else
	ss << csblast << "/csblast -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -B " << tmp_psifile << " -D " << csblast_db << " --blast-path " << blast << " --no-penalty 2> /dev/null";
    }
  
  command = ss.str();
  if (v>=4) cout << "Command: " << command << "\n";
  
  // and extract prefiltered HHMs
  FILE *stream = popen(command.c_str(), "r");
  while (fgets(line, LINELEN, stream))
    {
      ptr=strscn(line);
      ptr=strcut(ptr);
      if (ptr==NULL)
	{cerr<<endl<<"Error in "<<program_name<<": Cannot parse result of pre-filtering!\n"; exit(6);}
      char tmp_name[NAMELEN];
      char db_name[NAMELEN];
      char tmp[NAMELEN];
      ptr=strwrd(tmp_name,ptr);
      
      if (block_filter)
	{
	  // Write block for template
	  if (strcmp(actual_hit,"") && strcmp(actual_hit,tmp_name))   // New template
	    {
	      //printf("Add to block shading   key: %s    data:",actual_hit);
	      // for (int i = 0; i < block_count; i++)
	      //   printf(" %i,",block[i]);
	      // printf("\n");
	      par.block_shading->Add(actual_hit,block);
	      par.block_shading_counter->Add(actual_hit,block_count);
	      block = new(int[400]);
	      block_count = 0;
	    }
	  if (block_count >= 400) { continue; }
	  strcpy(actual_hit,tmp_name);
	  // Get block of HSP
	  ptr=strwrd(tmp,ptr); // sequence identity
	  pos=strint(ptr); // ali length
	  pos=strint(ptr); // mismatches
	  pos=strint(ptr); // gap openings
	  pos=strint(ptr); // query start
	  block[block_count++]=pos;
	  pos=strint(ptr); // query end
	  block[block_count++]=pos;
	  pos=strint(ptr); // subject start
	  block[block_count++]=pos;
	  pos=strint(ptr); // subject end
	  block[block_count++]=pos;	    
	}
      
      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      if (! doubled->Contains(db_name))
	{
	  doubled->Add(db_name);
	  // check, if DB was searched in previous rounds 
	  strcat(tmp_name,"__1");  // irep=1
	  if (previous_hits->Contains(tmp_name))
	    {
	      dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	      strcpy(dbfiles_old[ndb_old],dbhhm);
	      strcat(dbfiles_old[ndb_old],"/");
	      strcat(dbfiles_old[ndb_old],db_name);
	      if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	      ndb_old++;
	    }
	  else 
	    {
	      dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	      strcpy(dbfiles_new[ndb_new],dbhhm);
	      strcat(dbfiles_new[ndb_new],"/");
	      strcat(dbfiles_new[ndb_new],db_name);
	      if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	      ndb_new++;
	    }
	}
    }
  pclose(stream);
  if (block_filter && strcmp(actual_hit,""))   // New template
    {
      // printf("Add to block shading   key: %s    data:",actual_hit);
      // for (int i = 0; i < block_count; i++)
      //   printf(" %i,",block[i]);
      // printf("\n");
      par.block_shading->Add(actual_hit,block);
      par.block_shading_counter->Add(actual_hit,block_count);
    }
}

void prefilter_with_BLAST_only()
{
  if (v>=2) printf("Pre-filtering with %s ...\n",(strcmp(pre_mode,"csblast"))?"PSI-BLAST":"CS-BLAST");
  stringstream ss;
  if (!is_regular_file(tmp_psifile.c_str())) 
    {
      if (strcmp(pre_mode,"csblast")) 
	ss << blast << "/blastpgp -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 2> /dev/null";
      else
	ss << csblast << "/csblast -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -D " << csblast_db << " --blast-path " << blast << " --no-penalty 2> /dev/null";
    } 
  else
    {
      if (strcmp(pre_mode,"csblast")) 
	ss << blast << "/blastpgp -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -B " << tmp_psifile << " 2> /dev/null";
      else
	ss << csblast << "/csblast -d " << db << " -a " << cpu << " -b " << N_PSI << " -e " << e_psi << " -i " << tmp_infile << " -m 8 -B " << tmp_psifile << " -D " << csblast_db << " --blast-path " << blast << " --no-penalty 2> /dev/null";
    }
  
  command = ss.str();
  if (v>=4) cout << "Command: " << command << "\n";
  
  // and extract prefiltered HHMs
  FILE* outf=NULL;
  outf=fopen(par.outfile,"w");
  fprintf(outf,"        Name        Score\n");

  float evalue;
  float score;

  FILE *stream = popen(command.c_str(), "r");
  while (fgets(line, LINELEN, stream))
    {
      ptr=strscn(line);
      ptr=strcut(ptr);
      if (ptr==NULL)
	{cerr<<endl<<"Error in "<<program_name<<": Cannot parse result of pre-filtering!\n"; exit(6);}
      char tmp_name[NAMELEN];
      char db_name[NAMELEN];
      char tmp[NAMELEN];
      ptr=strwrd(tmp_name,ptr);
      
      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	//substr(db_name,tmp_name,3,11);
	continue;
      else                              // other database
	strcpy(db_name,tmp_name);
	  
      if (! doubled->Contains(db_name))
	{
	  doubled->Add(db_name);      

	  ptr=strwrd(tmp,ptr); // sequence identity
	  pos=strint(ptr); // ali length
	  pos=strint(ptr); // mismatches
	  pos=strint(ptr); // gap openings
	  pos=strint(ptr); // query start
	  pos=strint(ptr); // query end
	  pos=strint(ptr); // subject start
	  pos=strint(ptr); // subject end
	  ptr=strwrd(tmp,ptr); // sequence identity
	  evalue = atof(tmp);
	  ptr=strwrd(tmp,ptr); // sequence identity
	  score = atof(tmp);
	  
	  fprintf(outf,"%-20s   %6.2f   %6.2g\n", db_name, score, evalue);
	  //printf("Hit %-20s   has score %6.2f  (E-value: %6.2g)\n", db_name, score, evalue);
	}
    }
  pclose(stream);
  fclose(outf);
  printf("Scorefile of prefiltering-hits written to %s!\n",par.outfile);
}

void prefilter_with_SW()
{
  init_sse_prefiltering();

  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  int score;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<int, string> > hits;

#pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      //score = score / flog2(length[n]);
 
      if (score > par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<int,string>(score, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<int, string> >::reverse_iterator it;
  
  for ( it=hits.rbegin() ; it < hits.rend(); it++ )
    {
      //printf("Hit %20s with score %4i\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      // check, if DB was searched in previous rounds 
      strcat(tmp_name,"__1");  // irep=1
      if (previous_hits->Contains(tmp_name))
	{
	  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_old[ndb_old],dbhhm);
	  strcat(dbfiles_old[ndb_old],"/");
	  strcat(dbfiles_old[ndb_old],db_name);
	  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	  ndb_old++;
	}
      else 
	{
	 dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_new[ndb_new],dbhhm);
	  strcat(dbfiles_new[ndb_new],"/");
	  strcat(dbfiles_new[ndb_new],db_name);
	  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	  ndb_new++;
	}
      
      if (ndb_new >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
	  break;
	  //exit(4);
	}
    }

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}

void prefilter_with_SW_score()
{
  init_sse_prefiltering();

  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  int score;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  int log_qlen = flog2(LQ);

  vector<pair<int, string> > hits;

#pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      score = score - par.prefilter_bit_factor * (log_qlen + flog2(length[n]));
 
      if (score > par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<int,string>(score, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<int, string> >::reverse_iterator it;
  
  for ( it=hits.rbegin() ; it < hits.rend(); it++ )
    {
      //printf("Hit %20s with score %4i\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      // check, if DB was searched in previous rounds 
      strcat(tmp_name,"__1");  // irep=1
      if (previous_hits->Contains(tmp_name))
	{
	  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_old[ndb_old],dbhhm);
	  strcat(dbfiles_old[ndb_old],"/");
	  strcat(dbfiles_old[ndb_old],db_name);
	  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	  ndb_old++;
	}
      else 
	{
	 dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_new[ndb_new],dbhhm);
	  strcat(dbfiles_new[ndb_new],"/");
	  strcat(dbfiles_new[ndb_new],db_name);
	  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	  ndb_new++;
	}
      
      if (ndb_new >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
	  break;
	  //exit(4);
	}
    }

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}

void prefilter_with_SW_evalue()
{
  init_sse_prefiltering();

  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  int score;
  double evalue;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<double, string> > hits;

  const double factor = (double)dbsize * LQ;

#pragma omp parallel for schedule(static) private(evalue, score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      evalue = factor * length[n] * fpow2(-score/par.prefilter_bit_factor);
 
      if (evalue < par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<double,string>(evalue, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<double, string> >::iterator it;
  
  for ( it=hits.begin() ; it < hits.end(); it++ )
    {
      //printf("Hit %20s with score %6.2g\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	{
	  substr(tmp,tmp_name,3,11);
	  substr(db_name,tmp,0,1);
	  strcat(db_name,"/");
	  strcat(db_name,tmp);
	  strcat(db_name,".db");
	}
      else                              // other database
	{
	  strcpy(db_name,tmp_name);
	  strtr(db_name,"|", "_");
	  strcat(db_name,".hhm");
	}
      
      // check, if DB was searched in previous rounds 
      strcat(tmp_name,"__1");  // irep=1
      if (previous_hits->Contains(tmp_name))
	{
	  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_old[ndb_old],dbhhm);
	  strcat(dbfiles_old[ndb_old],"/");
	  strcat(dbfiles_old[ndb_old],db_name);
	  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
	  ndb_old++;
	}
      else 
	{
	 dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
	  strcpy(dbfiles_new[ndb_new],dbhhm);
	  strcat(dbfiles_new[ndb_new],"/");
	  strcat(dbfiles_new[ndb_new],db_name);
	  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
	  ndb_new++;
	}
      
      if (ndb_new >= MAXNUMDB) 
	{
	  printf("\nWARNING! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
	  break;
	  //exit(4);
	}
    }

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}

// void prefilter_with_SW()
// {
//   block_filter = false;
  
//   int num_dbs = init_sse_prefiltering();

//   __m128i* workspace = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
//   int score;
//   char db_name[NAMELEN];
//   char tmp_name[NAMELEN];
//   char tmp[NAMELEN];

//   for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
//     {
//       // Perform search step
//       //score = swStripedByte(qc, LQ, first[n], length[n], -par.prefilter_gap_open, -par.prefilter_gap_extend, workspace, workspace + W, workspace + 2*W, par.prefilter_score_offset);
//       score = swStripedByte(qc, LQ, first[n], length[n], par.prefilter_gap_open + par.prefilter_gap_extend, par.prefilter_gap_extend, workspace, workspace + W, workspace + 2*W, par.prefilter_score_offset);
      
//       if (score > par.prefilter_smax_thresh)
// 	{
// 	  ptr=strwrd(tmp_name,dbnames[n]);
	  
// 	  if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
// 	    {
// 	      substr(tmp,tmp_name,3,11);
// 	      substr(db_name,tmp,0,1);
// 	      strcat(db_name,"/");
// 	      strcat(db_name,tmp);
// 	      strcat(db_name,".db");
// 	    }
// 	  else                              // other database
// 	    {
// 	      strcpy(db_name,tmp_name);
// 	      strtr(db_name,"|", "_");
// 	      strcat(db_name,".hhm");
// 	    }
	  
// 	  if (! doubled->Contains(db_name))
// 	    {
// 	      doubled->Add(db_name);
// 	      // check, if DB was searched in previous rounds 
// 	      strcat(tmp_name,"__1");  // irep=1
// 	      if (previous_hits->Contains(tmp_name))
// 		{
// 		  dbfiles_old[ndb_old]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
// 		  strcpy(dbfiles_old[ndb_old],dbhhm);
// 		  strcat(dbfiles_old[ndb_old],"/");
// 		  strcat(dbfiles_old[ndb_old],db_name);
// 		  if (ndb_old<5 && ndb_old>0 && access(dbfiles_old[ndb_old],R_OK)) OpenFileError(dbfiles_old[ndb_old]); // file not readable?
// 		  ndb_old++;
// 		}
// 	      else 
// 		{
// 		  dbfiles_new[ndb_new]=new(char[strlen(dbhhm)+strlen(db_name)+2]);
// 		  strcpy(dbfiles_new[ndb_new],dbhhm);
// 		  strcat(dbfiles_new[ndb_new],"/");
// 		  strcat(dbfiles_new[ndb_new],db_name);
// 		  if (ndb_new<5 && ndb_new>0 && access(dbfiles_new[ndb_new],R_OK)) OpenFileError(dbfiles_new[ndb_new]); // file not readable?
// 		  ndb_new++;
// 		}
// 	    }
// 	  if (ndb_new >= MAXNUMDB) 
// 	    {
// 	      printf("\nERROR! To many hits through prefilter! (MAXNUM = %6i)\n",MAXNUMDB);
// 	      exit(4);
// 	    }
// 	}
//     }

//   // Free memory
//   free(qc);
//   free(X);
//   free(length);
//   free(first);
//   free(workspace);
//   for (int n = 0; n < num_dbs; n++)
//     delete[](dbnames[n]);
//   delete[](dbnames);
// }

void prefilter_with_SW_only()
{
  init_sse_prefiltering();

  __m128i** workspace = new(__m128i*[cpu]);

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  int score;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<int, string> > hits;

  #pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
      
      if (score > par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<int,string>(score, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<int, string> >::reverse_iterator it;
  
  FILE* outf=NULL;
  outf=fopen(par.outfile,"w");
  fprintf(outf,"        Name        Score\n");

  for ( it=hits.rbegin() ; it < hits.rend(); it++ )
    {
      //printf("Hit %20s with score %4i\n", ((*it).second).c_str(), (*it).first);

      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	//substr(db_name,tmp_name,3,11);
	continue;
      else                              // other database
	strcpy(db_name,tmp_name);
	  
      fprintf(outf,"%-20s   %4i\n", db_name, (*it).first);
    }
  fclose(outf);
  printf("Scorefile of prefiltering-hits written to %s!\n",par.outfile);

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}

void prefilter_with_SW_score_only()
{
  init_sse_prefiltering();

  __m128i** workspace = new(__m128i*[cpu]);

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  int score;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  int log_qlen = flog2(LQ);

  vector<pair<int, string> > hits;

  #pragma omp parallel for schedule(static) private(score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);

      score = score - par.prefilter_bit_factor * (log_qlen + flog2(length[n]));
      
      if (score > par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<int,string>(score, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<int, string> >::reverse_iterator it;
  
  FILE* outf=NULL;
  outf=fopen(par.outfile,"w");
  fprintf(outf,"        Name        Score\n");

  for ( it=hits.rbegin() ; it < hits.rend(); it++ )
    {
      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	//substr(db_name,tmp_name,3,11);
	continue;
      else                              // other database
	strcpy(db_name,tmp_name);
	  
      fprintf(outf,"%-20s   %4i\n", db_name, (*it).first);
    }
  fclose(outf);
  printf("Scorefile of prefiltering-hits written to %s!\n",par.outfile);

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
}

void prefilter_with_SW_evalue_only()
{
  init_sse_prefiltering();

  if (print_elapsed) ElapsedTimeSinceLastCall("(init prefiltering)");

  __m128i** workspace = new(__m128i*[cpu]);

  for (int i = 0; i < cpu; i++)
    workspace[i] = (__m128i*)memalign(16,3*(LQ+15)*sizeof(char));
  
  int gap_init = par.prefilter_gap_open + par.prefilter_gap_extend;
  int gap_extend = par.prefilter_gap_extend;

  int score;
  double evalue;
  int thread_id = 0;
  char db_name[NAMELEN];
  char tmp_name[NAMELEN];
  char tmp[NAMELEN];

  vector<pair<double, string> > hits;

  const double factor = (double)dbsize * LQ;

#pragma omp parallel for schedule(static) private(evalue, score, thread_id)
  for (int n = 0; n < num_dbs; n++)     // Loop over all database sequences
    {
      #ifdef _OPENMP
      thread_id = omp_get_thread_num();
      #endif

      // Perform search step
      score = swStripedByte(qc, LQ, first[n], length[n], gap_init, gap_extend, workspace[thread_id], workspace[thread_id] + W, workspace[thread_id] + 2*W, par.prefilter_score_offset);
     
      evalue = factor * length[n] * fpow2(-score/par.prefilter_bit_factor);
 
      if (evalue < par.prefilter_smax_thresh)
	{
          #pragma omp critical
	  hits.push_back(pair<double,string>(evalue, string(dbnames[n])));
	}
    }

  sort(hits.begin(), hits.end());

  vector<pair<double, string> >::iterator it;
  
  FILE* outf=NULL;
  outf=fopen(par.outfile,"w");
  fprintf(outf,"        Name        Score\n");

  for ( it=hits.begin() ; it < hits.end(); it++ )
    {
      strcpy(tmp, ((*it).second).c_str());

      ptr=strwrd(tmp_name,tmp);

      if (!strncmp(tmp_name,"cl|",3))   // kClust formatted database (NR20, NR30)
	//substr(db_name,tmp_name,3,11);
	continue;
      else                              // other database
	strcpy(db_name,tmp_name);
	  
      fprintf(outf,"%-20s   %6.2g\n", db_name, (*it).first);
    }
  fclose(outf);
  printf("Scorefile of prefiltering-hits written to %s!\n",par.outfile);

  // Free memory
  free(qc);
  for (int i = 0; i < cpu; i++)
    free(workspace[i]);
  delete[] workspace;
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
  block = new(int[400]);
  strcpy(actual_hit,"");

  par.block_shading->Reset();
  while (!par.block_shading->End())
    delete[] (par.block_shading->ReadNext()); 
  par.block_shading->New(16381,NULL);
  par.block_shading_counter->New(16381,NULL);
  
  if (!strcmp(pre_mode,"ungapped_SW"))  // SSE case
    prefilter_with_ungapped_SW();
  else if (!strcmp(pre_mode,"only_prefilt_ungapped_SW"))  // Smith-Waterman case
    prefilter_with_ungapped_SW_only();
  else if (!strcmp(pre_mode,"ungapped_SW_no_region"))  // Smith-Waterman case
    prefilter_with_ungapped_SW_no_region();
  else if (!strcmp(pre_mode,"only_prefilt_ungapped_SW_no_region"))  // Smith-Waterman case
    prefilter_with_ungapped_SW_no_region_only();
  else if (!strcmp(pre_mode,"SW"))  // Smith-Waterman case
    prefilter_with_SW();
  else if (!strcmp(pre_mode,"SW_score"))  // Smith-Waterman case
    prefilter_with_SW_score();
  else if (!strcmp(pre_mode,"SW_evalue"))  // Smith-Waterman case
    prefilter_with_SW_evalue();
  else if (!strcmp(pre_mode,"SW_evalue_preprefilter"))  // Smith-Waterman case
    prefilter_with_SW_evalue_preprefilter();
  else if (!strcmp(pre_mode,"SW_evalue_preprefilter_backtrace"))  // Smith-Waterman case
    {
      if (block_filter) 
	prefilter_with_SW_evalue_preprefilter_backtrace();
      else
	prefilter_with_SW_evalue_preprefilter();
    }
  else if (!strcmp(pre_mode,"only_prefilt_SW"))  // Smith-Waterman case
    prefilter_with_SW_only();
  else if (!strcmp(pre_mode,"only_prefilt_SW_score"))  // Smith-Waterman case
    prefilter_with_SW_score_only();
  else if (!strcmp(pre_mode,"only_prefilt_SW_evalue"))  // Smith-Waterman case
    prefilter_with_SW_evalue_only();
  else if (!strcmp(pre_mode,"combi"))  // Smith-Waterman case
    {
      prefilter_with_SW();
      prefilter_with_BLAST();
    }
  else if (!strcmp(pre_mode,"only_prefilt_csblast"))  // Smith-Waterman case
    prefilter_with_BLAST_only();
  else
    prefilter_with_BLAST();

  if(doubled) delete doubled;
}
