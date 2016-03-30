/*
 * hhbiterbimatrix-inl.h
 *
 *  Created on: Jun 16, 2014
 *      Author: meiermark
 */

#ifndef HHBITERBIMATRIX_INL_H_
#define HHBITERBIMATRIX_INL_H_


inline bool ViterbiMatrix::ViterbiMatrix::hasCellOff(){
    return this->cellOff;
}


inline unsigned char * ViterbiMatrix::ViterbiMatrix::getRow(int row){
    return this->bCO_MI_DG_IM_GD_MM_vec[row];
}


inline void ViterbiMatrix::setCellOff(bool value){
    this->cellOff = value;
}


inline void ViterbiMatrix::setCellOff(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem],7);
        this->setCellOff(true);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem], 7);
    }
}

inline void ViterbiMatrix::setMatIns(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem],6);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem], 6);
    }
}

inline void ViterbiMatrix::setDelGap(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem],5);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem], 5);
    }
}

inline void ViterbiMatrix::setInsMat(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem],4);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem], 4);
    }
}

inline void ViterbiMatrix::setGapDel(int row,int col,int elem,bool value){
    if(value){
        BIT_SET(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem],3);
    }else{
        BIT_CLEAR(this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem], 3);
    }
}


inline void ViterbiMatrix::setMatMat(int row,int col,int elem,unsigned char value){
    //0xF8 11111000
//    this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem]&=(0xF8 ^ value);
        unsigned char c = 0xF8;
        this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem] =
                (this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem] & c) | value;
}

inline bool ViterbiMatrix::getCellOff(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 128);
}


inline bool ViterbiMatrix::getMatIns(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 64);
}


inline bool ViterbiMatrix::getDelGap(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 32);
}


inline bool ViterbiMatrix::getInsMat(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 16);
}


inline bool ViterbiMatrix::getGapDel(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem];
    return (bool) (vCO_MI_DG_GD_MM & 8);
}


inline int  ViterbiMatrix::getMatMat(int row,int col,int elem){
    unsigned char vCO_MI_DG_GD_MM=this->bCO_MI_DG_IM_GD_MM_vec[row][(col*VECSIZE_FLOAT)+elem];
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

#endif /* HHBITERBIMATRIX_INL_H_ */
