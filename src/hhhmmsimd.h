// hhhmm.h

#ifndef HMMSIMD_h
#define HMMSIMD_h
#include "hhhmm.h"
#include "simd.h"

#include <emmintrin.h>


class HMMSimd
{
public:
    
    HMMSimd(int maxres);
    ~HMMSimd();
    HMMSimd& operator=(HMMSimd&);
    float** p;                // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
    int L;                    // length of HMM = number of match states; set in declaration of HMM object
    simd_int * lengths;					// length of HMMs respectively
    //M2M M2I M2D I2M I2I D2M D2D
    simd_float * tr;               // tr_m2m[i] = log2 of transition probabilities M2M
    // contains index lookup for SS score dssp
    unsigned char *dssp_index;
    // contains index lookup for SS score prediction
    unsigned char *pred_index;


//    simd_float * tr_m2m;               // tr_m2m[i] = log2 of transition probabilities M2M
//    simd_float * tr_m2i;               // tr_m2i[i] = log2 of transition probabilities M2I
//    simd_float * tr_m2d;               // tr_m2i[i] = log2 of transition probabilities M2d
//    simd_float * tr_d2m;               // tr_d2m[i] = log2 of transition probabilities D2M
//    simd_float * tr_d2d;               // tr_d2d[i] = log2 of transition probabilities D2D
//    simd_float * tr_i2m;               // tr_i2m[i] = log2 of transition probabilities I2M
//    simd_float * tr_i2i;               // tr_i2i[i] = log2 of transition probabilities I2I 

    // Maps HMMs to a HMMSimd
    void MapHMMVector(std::vector<HMM *> hmms);
    // Maps one HMM to a HMMSimd 
    void MapOneHMM(HMM *seq);
    // returns pointer to HMM
    HMM* GetHMM(int elem);
private:
    // pointer to sequences
    HMM ** seqarr;
    int maxres;

    // contains profile data (aligned in memory)
    float *p_data;
};



#endif
