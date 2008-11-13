#ifndef AB_AMINO_ACID_H
#define AB_AMINO_ACID_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//Amino acid class that encapsulates amino acid meta data

#include <cstring> //size_t
#include <cctype>  //toupper

static const char __I2AA[]  = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X'};
static const int  __AL2I[]  = { 0 , 4 , 3 , 6 , 13, 7 , 8 , 9 , 11, 10, 12, 2 , 14, 5 , 1 , 15, 16, 19, 17, 18};
static const char __AL2AA[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

class AminoAcid{
public:

    AminoAcid();
    ~AminoAcid();

    inline const char * const get_i2aa() const { return i2aa; }
    const int  * const get_aa2i() const { return aa2i; }
    inline const char * const get_al2aa() const { return al2aa; }
    const int  * const get_al2i() const { return al2i; }
    const int get_naa() const { return NAA; }
    const int get_naa_any() const { return NAA_ANY; }
    const int get_any() const { return ANY; }
    const int get_gap() const { return GAP; }
    const int get_endgap() const { return ENDGAP; }
    const int naa() const { return NAA; }
    const int naa_any() const { return NAA_ANY; }
    const int any() const { return ANY; }
    const int gap() const { return GAP; }
    const int endgap() const { return ENDGAP; }
    inline const bool is_any(char aa) const { return aa2i[toupper(aa)]==ANY; }
    inline const bool is_any(int i) const { return i==ANY; }

private:
    const int NAA;
    const int NAA_ANY;
    const int ANY;
    const int GAP;
    const int ENDGAP;

    char *i2aa;
    int *aa2i;
    int *al2i;
    char *al2aa;
};

#endif
