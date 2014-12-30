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

    max_query_length = ICEIL(Nq+2,VECSIZE_FLOAT);
    max_template_length = ICEIL(((Nt+2) *VEC_SIZE),VECSIZE_FLOAT);
    // Allocate posterior prob matrix (matrix rows are padded to make them aligned to multiples of ALIGN_FLOAT)
    bCO_MI_DG_IM_GD_MM_vec = malloc_matrix<unsigned char>(max_query_length, max_template_length);
    if (!bCO_MI_DG_IM_GD_MM_vec)
        MemoryError("m_probabilities", "hhviterbimatrix.cpp", 55, "ViterbiMatrix::allocateMatrix");

}



/////////////////////////////////////////////////////////////////////////////////////
//// Delete memory for dynamic programming matrix
/////////////////////////////////////////////////////////////////////////////////////
void ViterbiMatrix::DeleteBacktraceMatrix() {
    if(max_query_length == 0) {
        return;
    }

    free(this->bCO_MI_DG_IM_GD_MM_vec);
    this->bCO_MI_DG_IM_GD_MM_vec = NULL;

    max_query_length = 0;
    max_template_length = 0;
}

#endif
