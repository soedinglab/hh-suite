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
    this->max_query_length = 0;
    this->max_template_length = 0;
}


ViterbiMatrix::~ViterbiMatrix(){
  DeleteBacktraceMatrix();
}



/////////////////////////////////////////////////////////////////////////////////////
//// Allocate memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void ViterbiMatrix::AllocateBacktraceMatrix(int Nq, int Nt)
{
  Nq += 1;
  Nt += 1;

    if(Nq > max_query_length || Nt > max_template_length) {
      DeleteBacktraceMatrix();
    }
    else {
      return;
    }

    this->bCO_MI_DG_IM_GD_MM_vec = new unsigned char*[Nq];
    
    for (int i=0; i<Nq; ++i) {
        this->bCO_MI_DG_IM_GD_MM_vec[i]=(unsigned char *)malloc_simd_float(VEC_SIZE*Nt*sizeof(unsigned char));
        memset(this->bCO_MI_DG_IM_GD_MM_vec[i], 0, VEC_SIZE*Nt*sizeof(unsigned char));
        if (!this->bCO_MI_DG_IM_GD_MM_vec[i]) {
            fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
            fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
            fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
            exit(3);
        }
    }

    max_query_length = Nq;
    max_template_length = Nt;
}



/////////////////////////////////////////////////////////////////////////////////////
//// Delete memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void ViterbiMatrix::DeleteBacktraceMatrix() {
  if(max_query_length == 0) {
    return;
  }

  for (int i=0; i < max_query_length; ++i) {
      free(this->bCO_MI_DG_IM_GD_MM_vec[i]);
  }
  delete[] this->bCO_MI_DG_IM_GD_MM_vec;
  this->bCO_MI_DG_IM_GD_MM_vec = NULL;

  max_query_length = 0;
  max_template_length = 0;
}

#endif
