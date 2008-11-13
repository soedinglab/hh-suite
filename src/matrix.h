#ifndef AB_MATRIX_H
#define AB_MATRIX_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//Amino acid substitution matrix class
//Exact values for substitution scores are computed from the raw BLOSUM data
//Export BLOSUM_MATRICES='/dir/' to point to the directory where the blosumXX.out files reside
//The files are available at ftp://ftp.ncbi.nih.gov/repository/blocks/unix/blosum/blosum.tar.Z
//The dummy 'static_blosum62' can be alternatively used

#include <sstream>      //stringstream
#include <fstream>      //file streams
#include <iostream>     //ostream, istream
#include <cstdio>       //printf, fprintf
#include <cstdlib>      //getenv

#include "dummy_matrix.h"
#include "my_exception.h"

class Matrix{
public:
    enum mtype{
        blosum30  = 0,
        blosum35  = 1,
        blosum40  = 2,
        blosum45  = 3,
        blosum50  = 4,
        blosum55  = 5,
        blosum60  = 6,
        blosum62  = 7,
        blosum65  = 8,
        blosum70  = 9,
        blosum75  = 10,
        blosum80  = 11,
        blosum85  = 12,
        blosum90  = 13,
        blosum95  = 14,
        blosum100 = 15,
        static_blosum62   = 16
    };
    Matrix(const mtype t=Matrix::static_blosum62) throw (std::exception);
    ~Matrix();

    //rounds values in matrix to integer values
    void round_bit_scores();

    //pointers to data
    const float** const get_matrix() const;
    const float** const get_prob_matrix() const;
    const float* const  get_p_back() const;
    const float** const  get_p_cond() const;
    const char * const get_i2aa() const;
    const int  * const get_aa2i() const;
    const size_t get_aa_dim() const;

    //prints float values of matrix
    std::ostream& print_p_back(std::ostream&);
    //sets matrix to original blosum bit-scores
    void reset();

    //prints float values of matrix
    std::ostream& print(std::ostream&);

    //prints float values of matrix
    std::ostream& print_frequencies(std::ostream&);
    //prints rounded bit-scores of matrix[][]
    std::ostream& print_int(std::ostream&);

    //set matrix[][] with arbitrary values
    void read_matrix(const float src[][21]);
    std::ostream& print_debug(std::ostream&);

private:
    const size_t AMINOACID_DIM;
    float **matrix;
    float **original;
    float *p_background;
    float **p_cond;
    float scale;
    void _read_blosum_matrix(const char*) throw (std::exception);
    void _copy_dummy();
    void _copy();
    const std::string _get_fn(const mtype);
    char *i2aa;
    int *aa2i;
};

#endif
