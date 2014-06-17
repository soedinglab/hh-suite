//
//  hhviterbifourmatrix.c
//  GitHHSuite
//
//  Created by Martin Steinegger on 19.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef HHVITERBIMATRIX_c
#define HHVITERBIMATRIX_c
#include "hhviterbimatrix.h"



ViterbiMatrix::ViterbiMatrix(){
    this->bCO_MI_DG_IM_GD_MM_vec=NULL;
    this->cellOff = false;
}


ViterbiMatrix::~ViterbiMatrix(){
}



/////////////////////////////////////////////////////////////////////////////////////
//// Allocate memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void ViterbiMatrix::AllocateBacktraceMatrix(int Nq, int Nt)
{
    int i;
    this->bCO_MI_DG_IM_GD_MM_vec=new(unsigned char*[Nq]);
    
    for (i=0; i<Nq; ++i)
    {
        this->bCO_MI_DG_IM_GD_MM_vec[i]=(unsigned char *)malloc_simd_float(VEC_SIZE*Nt*sizeof(unsigned char));
        memset(this->bCO_MI_DG_IM_GD_MM_vec[i], 0, VEC_SIZE*Nt*sizeof(unsigned char));
        if (!this->bCO_MI_DG_IM_GD_MM_vec[i])
        {
            fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
            fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
            fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
            exit(3);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
//// Delete memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void ViterbiMatrix::DeleteBacktraceMatrix(int Nq)
{
    int i;
    for (i=0; i<Nq; ++i) {
        free(this->bCO_MI_DG_IM_GD_MM_vec[i]);
    }
    delete[] this->bCO_MI_DG_IM_GD_MM_vec;
    this->bCO_MI_DG_IM_GD_MM_vec=NULL;
}

#endif
