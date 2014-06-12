//
//  hhhmm4.C
//  GitHHSuite
//
//  Created by Martin Steinegger on 24.06.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#ifndef HMMSIMD_c
#define HMMSIMD_c
#include "hhhmmsimd.h"

#include <algorithm>
#include <float.h>    // FLT_MIN
#include "simd.h"
#include "hhutil.h"


/////////////////////////////////////////////////////////////////////////////////////
//// Class HMM
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
HMMSimd::HMMSimd(int maxres)
{
    
    p = new float*[maxres];         // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
    for (int i=0; i<maxres; i++) p[i]=(float*) malloc_simd_float(VEC_SIZE*NAA*sizeof(float));  // align memory on 16b boundaries for SSE2
    tr = malloc_simd_float(VEC_SIZE*maxres*sizeof(float)*7);
//    tr_m2i = malloc_simd_float(VEC_SIZE*maxres*sizeof(float));
//    tr_m2d = malloc_simd_float(VEC_SIZE*maxres*sizeof(float));
//    tr_d2m = malloc_simd_float(VEC_SIZE*maxres*sizeof(float));
//    tr_d2d = malloc_simd_float(VEC_SIZE*maxres*sizeof(float));
//    tr_i2m = malloc_simd_float(VEC_SIZE*maxres*sizeof(float));
//    tr_i2i = malloc_simd_float(VEC_SIZE*maxres*sizeof(float));
    this->seqarr = new HMM*[VEC_SIZE];
    this->L = 0;
    lengths = malloc_simd_int(VEC_SIZE * sizeof(int));
    this->maxres = maxres;
}



/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
HMMSimd::~HMMSimd()
{
    //M2M M2I M2D I2M I2I D2M D2D
    free(tr);
    free(lengths);
    for (int i=0; i<maxres; i++) if (p[i])  free(p[i]);  else break;
    delete[] p;
    delete[] seqarr;
}

void HMMSimd::MapOneHMM(HMM *seq){
    std::vector<HMM *> tmp;
    for (unsigned int i = 0; i < VEC_SIZE; i++) {
        tmp.push_back(seq);
    }
    MapHMMVector(tmp);
}

HMM* HMMSimd::GetHMM(int elem){
    return this->seqarr[elem];
}

// Maps four HMMs to a HMM4 
void HMMSimd::MapHMMVector(std::vector<HMM *> hmms){
    
    if(hmms.size()>VEC_SIZE){
        std::cerr<<"Error: More than expected HMMs should be mapped. Please report this bug to developers\n";
        exit(3);
    }
    this->L = 0;

    int * lenghts_ptr = (int *)this->lengths;


    for (unsigned int i = 0; i < hmms.size(); i++) {
        HMM * hmm = hmms[i];
        this->L = std::max(hmm->L,this->L);

        // store all lengths of the HMMs in a vector
        lenghts_ptr[i] = hmm->L;

        this->seqarr[i] = hmm;
    }


    float * tr_scalar = (float *)this->tr;
    
    for(unsigned int seq_i = 0; seq_i < hmms.size(); seq_i++){
        HMM * curr = seqarr[seq_i];
        for(int i = 0; i <= curr->L; i++){
            const unsigned int start_pos = i * VEC_SIZE * 7;
//            const simd_float t_m2m = simdf32_load((float *)&t->tr_m2m[j-1]);
//            const simd_float t_d2m = simdf32_load((float *)&t->tr_d2m[j-1]);
//            const simd_float t_i2m = simdf32_load((float *)&t->tr_i2m[j-1]);
//            const simd_float t_m2d = simdf32_load((float *)&t->tr_m2d[j-1]);
//            const simd_float t_d2d = simdf32_load((float *)&t->tr_d2d[j-1]);
//            const simd_float t_m2i = simdf32_load((float *)&t->tr_m2i[j]);
//            const simd_float t_i2i = simdf32_load((float *)&t->tr_i2i[j]);
//          cache line optimized order   
            tr_scalar[start_pos + 0 * VEC_SIZE + seq_i] = curr->tr[i][I2I];
            tr_scalar[start_pos + 1 * VEC_SIZE + seq_i] = curr->tr[i][M2I];
            tr_scalar[start_pos + 2 * VEC_SIZE + seq_i] = curr->tr[i][M2M];
            tr_scalar[start_pos + 3 * VEC_SIZE + seq_i] = curr->tr[i][M2D];
            tr_scalar[start_pos + 4 * VEC_SIZE + seq_i] = curr->tr[i][D2M];
            tr_scalar[start_pos + 5 * VEC_SIZE + seq_i] = curr->tr[i][D2D];
            tr_scalar[start_pos + 6 * VEC_SIZE + seq_i] = curr->tr[i][I2M];
            for(int aa_i=0; aa_i < NAA;aa_i++){
                p[i][aa_i*VEC_SIZE+seq_i]=curr->p[i][aa_i];
            }
        }
        for(int i = curr->L+1; i < this->L+1; i++){
            const unsigned int start_pos = i * VEC_SIZE * 7;
            tr_scalar[start_pos+0*VEC_SIZE+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+1*VEC_SIZE+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+2*VEC_SIZE+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+3*VEC_SIZE+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+4*VEC_SIZE+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+5*VEC_SIZE+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+6*VEC_SIZE+seq_i] = -FLT_MAX;

            for(int aa_i=0; aa_i < NAA;aa_i++){
                p[i][aa_i*VEC_SIZE+seq_i]=0;
            }
        }
    }
}

#endif
