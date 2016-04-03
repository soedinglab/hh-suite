//
//  hhviterbimatrix.h
//
//  Created by Martin Steinegger on 19.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef HHVITERBIMATRIX_h
#define HHVITERBIMATRIX_h

#include "hhutil.h"
#include "simd.h"
#include "hhhmmsimd.h"

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



class ViterbiMatrix {
public:  

    // Constructor (only set pointers to NULL)
    ViterbiMatrix();
    ~ViterbiMatrix();
    const static char STOP=0;
    const static char MM=2;
    const static char GD=3;
    const static char IM=4;
    const static char DG=5;
    const static char MI=6;
    const static bool GD_MM=true;
    const static bool GD_OTHER=false;
    const static bool IM_MM=true;
    const static bool IM_OTHER=false;
    const static bool DG_MM=true;
    const static bool DG_OTHER=false;
    const static bool MI_MM=true;
    const static bool MI_OTHER=false;


    unsigned char * getRow(int row); 
    void AllocateBacktraceMatrix(int Nq, int Nt);
    void DeleteBacktraceMatrix();

    bool getCellOff(int row,int col,int elem); 
    bool getMatIns(int row,int col,int elem); 
    bool getGapDel(int row,int col,int elem); 
    bool getInsMat(int row,int col,int elem); 
    bool getDelGap(int row,int col,int elem); 
    int  getMatMat(int row,int col,int elem); 
    
    void setCellOff(int row,int col,int elem,bool value);
    void setMatIns(int row,int col,int elem,bool value);
    void setDelGap(int row,int col,int elem,bool value);
    void setInsMat(int row,int col,int elem,bool value);
    void setGapDel(int row,int col,int elem,bool value);
    void setMatMat(int row,int col,int elem,unsigned char value);
    
    bool hasCellOff();
    void setCellOff(bool value);

    void printCellOff(int row_size,int col_size,int elem);

private:
    //CO	MI	DG	IM	GD	MM
    //1     1	1	1	1	1	1	1
    unsigned char ** bCO_MI_DG_IM_GD_MM_vec;
    // flag to indecated if cellOff is activ or not
    bool cellOff;

    int max_query_length;
    int max_template_length;

};

#include "hhviterbimatrix-inl.h"


#endif
