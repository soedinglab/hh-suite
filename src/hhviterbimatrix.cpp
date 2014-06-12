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


/* a=target variable, b=bit number to act upon 0-n */
#define BIT_SET(a,b) ((a) |= (1<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1<<(b)))
#define BIT_CHECK(a,b) ((a) & (1<<(b)))

/* x=target variable, y=mask */
#define BITMASK_SET(x,y) ((x) |= (y))
#define BITMASK_CLEAR(x,y) ((x) &= (~(y)))
#define BITMASK_FLIP(x,y) ((x) ^= (y))
#define BITMASK_CHECK(x,y) ((x) & (y))

ViterbiMatrix::ViterbiMatrix(){
    this->bCO_MI_DG_IM_GD_MM_vec=NULL;
    this->cellOff = false;
}


ViterbiMatrix::~ViterbiMatrix(){
}

//TODO: inline
bool ViterbiMatrix::ViterbiMatrix::hasCellOff(){
    return this->cellOff;
}

//TODO: inline
unsigned char * ViterbiMatrix::ViterbiMatrix::getRow(int row){
    return this->bCO_MI_DG_IM_GD_MM_vec[row];
}

//TODO: inline
void ViterbiMatrix::setCellOff(bool value){
    this->cellOff = value;
}

//TODO: inline
void ViterbiMatrix::setCellOff(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem],7);
        this->setCellOff(true);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem], 7);
    }
}

inline void ViterbiMatrix::setMatIns(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem],6);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem], 6);
    }
}

inline void ViterbiMatrix::setDelGap(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem],5);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem], 5);
    }
}

inline void ViterbiMatrix::setInsMat(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem],4);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem], 4);
    }
}

inline void ViterbiMatrix::setGapDel(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem],3);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem], 3);
    }
}

inline void ViterbiMatrix::setMatMat(int row,int col,int elem,unsigned char value){
    //0xF8 11111000
//    this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem]&=(0xF8 ^ value);
		unsigned char c = 0xF8;
		this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem] =
				(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem] & c) | value;
}

inline bool ViterbiMatrix::getCellOff(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 128);
}

//TODO: inline
bool ViterbiMatrix::getMatIns(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 64);
}

//TODO: inline
bool ViterbiMatrix::getDelGap(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 32);
}

//TODO: inline
bool ViterbiMatrix::getInsMat(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 16);
}

//TODO: inline
bool ViterbiMatrix::getGapDel(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 8);
}

//TODO: inline
int  ViterbiMatrix::getMatMat(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VEC_SIZE)+elem];
    return (int) (vCO_MI_DG_GD_MM & 7);
}

inline void ViterbiMatrix::printCellOff(int row_size,int col_size,int elem){
    for(int row = 0; row < row_size; row++){
        for(int col = 0; col < col_size;col++){
            std::cout << getCellOff(row,col,elem) << " ";
        }
        std::cout << std::endl;
    }
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
        this->bCO_MI_DG_IM_GD_MM_vec[i]=(unsigned char *)malloc_simd_float(VEC_SIZE*Nt*sizeof(unsigned char));;
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
        delete[] this->bCO_MI_DG_IM_GD_MM_vec[i];
    }
    delete[] this->bCO_MI_DG_IM_GD_MM_vec;
    this->bCO_MI_DG_IM_GD_MM_vec=NULL;
}

#endif
