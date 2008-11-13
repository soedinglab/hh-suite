/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "amino_acid.h"

AminoAcid::AminoAcid() : NAA(20), NAA_ANY(21), ANY(20), GAP(21), ENDGAP(22) {
    i2aa  = new char[NAA_ANY];
    al2i  = new int[NAA];
    al2aa = new char[NAA];
    
    for(int i=0; i<NAA_ANY; ++i) i2aa[i]  = __I2AA[i];
    for(int i=0; i<NAA; ++i)     al2i[i]  = __AL2I[i];
    for(int i=0; i<NAA; ++i)     al2aa[i] = __AL2AA[i];

    const size_t char_size = static_cast<int>(pow(2, 8*sizeof(char)));
    aa2i = new int[char_size];
    for(size_t i=0; i<char_size; ++i) aa2i[i]=-1;   
    for(int i=0; i<NAA; ++i) aa2i[i2aa[i]]=i;
    aa2i['X']=ANY;
    aa2i['J']=ANY;
    aa2i['O']=ANY;
    aa2i['U']=aa2i['C'];
    aa2i['B']=aa2i['D'];
    aa2i['Z']=aa2i['E'];
}

AminoAcid::~AminoAcid(){
    delete [] i2aa;
    delete [] aa2i;
    delete [] al2i;
    delete [] al2aa;
}
