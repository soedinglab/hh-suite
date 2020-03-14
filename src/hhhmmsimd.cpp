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
#include "log.h"


/////////////////////////////////////////////////////////////////////////////////////
//// Class HMM
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
HMMSimd::HMMSimd(int maxres)
{
    
    p = new float*[maxres];         // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
    // align memory on 16b boundaries for SSE2
    p_data = (float *) malloc_simd_float(maxres * VECSIZE_FLOAT * NAA * sizeof(float));
    size_t currPos = 0;
    for (int i=0; i<maxres; i++){
        p[i] = p_data + currPos;
        currPos += (VECSIZE_FLOAT * NAA);
    }
    tr = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float)*7);
//    tr_m2i = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float));
//    tr_m2d = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float));
//    tr_d2m = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float));
//    tr_d2d = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float));
//    tr_i2m = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float));
//    tr_i2i = malloc_simd_float(VECSIZE_FLOAT*maxres*sizeof(float));
    this->seqarr = new HMM*[VECSIZE_FLOAT];
    this->L = 0;
    lengths = malloc_simd_int(VECSIZE_FLOAT * sizeof(int));
    for(unsigned int i = 0; i < VECSIZE_FLOAT; i++){
       ((int*)lengths)[i] = 0;
    }
    this->dssp_index = (unsigned char *) malloc_simd_int(maxres * VECSIZE_FLOAT * sizeof(unsigned char));
    this->pred_index = (unsigned char *) malloc_simd_int(maxres * VECSIZE_FLOAT * sizeof(unsigned char));

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
    free(dssp_index);
    free(pred_index);
    free(p_data);
    delete[] p;
    delete[] seqarr;
}

void HMMSimd::MapOneHMM(HMM *seq){
    std::vector<HMM *> tmp;
    for (unsigned int i = 0; i < VECSIZE_FLOAT; i++) {
        tmp.push_back(seq);
    }
    MapHMMVector(tmp);
}

HMM* HMMSimd::GetHMM(int elem){
    return this->seqarr[elem];
}

// Maps four HMMs to a HMM4 
void HMMSimd::MapHMMVector(std::vector<HMM *> hmms){
    
    if(hmms.size() > VECSIZE_FLOAT){
      HH_LOG(ERROR) << "More than expected HMMs should be mapped. Please report this bug to developers\n";
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
            const unsigned int start_pos = i * VECSIZE_FLOAT * 7;
//            const simd_float t_m2m = simdf32_load((float *)&t->tr_m2m[j-1]);
//            const simd_float t_d2m = simdf32_load((float *)&t->tr_d2m[j-1]);
//            const simd_float t_i2m = simdf32_load((float *)&t->tr_i2m[j-1]);
//            const simd_float t_m2d = simdf32_load((float *)&t->tr_m2d[j-1]);
//            const simd_float t_d2d = simdf32_load((float *)&t->tr_d2d[j-1]);
//            const simd_float t_m2i = simdf32_load((float *)&t->tr_m2i[j]);
//            const simd_float t_i2i = simdf32_load((float *)&t->tr_i2i[j]);
//          cache line optimized order   
            tr_scalar[start_pos + 0 * VECSIZE_FLOAT + seq_i] = curr->tr[i][I2I];
            tr_scalar[start_pos + 1 * VECSIZE_FLOAT + seq_i] = curr->tr[i][M2I];
            tr_scalar[start_pos + 2 * VECSIZE_FLOAT + seq_i] = curr->tr[i][M2M];
            tr_scalar[start_pos + 3 * VECSIZE_FLOAT + seq_i] = curr->tr[i][M2D];
            tr_scalar[start_pos + 4 * VECSIZE_FLOAT + seq_i] = curr->tr[i][D2M];
            tr_scalar[start_pos + 5 * VECSIZE_FLOAT + seq_i] = curr->tr[i][D2D];
            tr_scalar[start_pos + 6 * VECSIZE_FLOAT + seq_i] = curr->tr[i][I2M];
            for(int aa_i=0; aa_i < NAA;aa_i++){
                p[i][aa_i*VECSIZE_FLOAT+seq_i]=curr->p[i][aa_i];
            }
            if (i > 0){
                pred_index[(i - 1) * VECSIZE_FLOAT + seq_i] = (unsigned char) curr->ss_pred[i] * MAXCF + curr->ss_conf[i];
                dssp_index[(i - 1) * VECSIZE_FLOAT + seq_i] = (unsigned char) curr->ss_dssp[i];
            }
        }
        for(int i = curr->L+1; i < this->L+1; i++){
            const unsigned int start_pos = i * VECSIZE_FLOAT * 7;
            tr_scalar[start_pos+0*VECSIZE_FLOAT+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+1*VECSIZE_FLOAT+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+2*VECSIZE_FLOAT+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+3*VECSIZE_FLOAT+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+4*VECSIZE_FLOAT+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+5*VECSIZE_FLOAT+seq_i] = -FLT_MAX;
            tr_scalar[start_pos+6*VECSIZE_FLOAT+seq_i] = -FLT_MAX;

            for(int aa_i=0; aa_i < NAA;aa_i++){
                p[i][aa_i*VECSIZE_FLOAT+seq_i] = 0;
            }
            pred_index[(i - 1) * VECSIZE_FLOAT + seq_i] = 0;
            dssp_index[(i - 1) * VECSIZE_FLOAT + seq_i] = 0;
        }
    }
    for(unsigned int seq_i = hmms.size(); seq_i < VECSIZE_FLOAT; seq_i++){
        for(int i = 1; i < this->L+1; i++){
            pred_index[(i - 1) * VECSIZE_FLOAT + seq_i] = 0;
            dssp_index[(i - 1) * VECSIZE_FLOAT + seq_i] = 0;
        }
    }
}

#endif
