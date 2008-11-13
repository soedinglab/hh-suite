#ifndef AB_SEQUENCE_H
#define AB_SEQUENCE_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//A container class for protein sequences

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "my_exception.h"
#include "amino_acid.h"
#include "matrix.h"

//2MB
#define READ_HBUFFER_SIZE 2097152
//1MB
#define READ_SBUFFER_SIZE 1048576

class Sequence
{
public:
    enum format
    {
        fasta     = 0,
        clustal   = 1, //not implemented
        stockholm = 2  //not implemented
    };

public:
    //constructors for numeric transformation
    Sequence( char* header,
              char* sequence,
              const int* const aa2int,
              const char* const int2aa,
              size_t ADIM,
              size_t index=0
        ) throw (std::exception);

    Sequence( const char* const header,
              const char* const sequence,
              AminoAcid *aa,
                size_t index=0
        ) throw (std::exception);

    Sequence( const char* fn,
              const format&,
              AminoAcid *aa,
              size_t index=0
        ) throw (std::exception);

    //create random seq of length 'length' with respect to background amino acid frequencies
    Sequence( const size_t length,
              AminoAcid *aa,
              Matrix *m,
              size_t index=0
        ) throw (std::exception);


    virtual ~Sequence();

    //inline functions HAVE TO reside in header files to be accessible from other cpp files
    inline const char * const get_sequence() const { return seq; }
    inline const char * const get_header()   const { return header; }
    inline char *sequence()                        { return seq; }
    inline const size_t length() const             { return len; }
    inline const size_t index()  const             { return idx; }
    static size_t get_length( const char* const ptr, const size_t l);
    void set_index(const size_t);

    std::ostream& print_debug(std::ostream&);
    std::ostream& write(std::ostream&, const size_t width=60);
    std::ostream& write_sequence(std::ostream&, const size_t width=60);

    void replace_J();
private:
    const size_t NAA;
    const int*  const aa2i;
    const char* const i2aa;

    char *header;
    char *seq;
    size_t len;
    size_t hlen;
    size_t idx;

    void read_fasta(const char*) throw (std::exception);
    void init( const char* const head,
               const size_t hl,
               const char* const sequence,
               const size_t sl) throw (std::exception);
};

#endif

