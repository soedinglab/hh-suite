#ifndef AB_PROFILE_H
#define AB_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//A class for sequence profiles

#include <string>
#include <fstream>
#include <iostream>
#include <limits>   // infinity
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cfloat>

#include "amino_acid.h"
#include "matrix.h"
#include "sequence.h"
#include "my_exception.h"

class Profile
{
public:
    //constructor for dummy profile object
    Profile() throw (std::exception);
    //constructor for empty profile object
    Profile( AminoAcid* aa, Matrix* m ) throw (std::exception);
    // constructor for profile of given length with probabilities initialized to zero
    Profile( size_t l, AminoAcid* aa, Matrix* m ) throw (std::exception);
    //constructor for profile object from serialized profile object
    Profile( std::istream& in,
             AminoAcid* aa,
             Matrix* m ) throw (std::exception);
    //constructor for profile object from sequence object
    Profile( const Sequence& sequence,
             AminoAcid* aa,
             Matrix* m ) throw (std::exception);
    // copy constructor
    Profile(const Profile& rhs);
    // destructor
    virtual ~Profile();
    // assignment operator
    virtual Profile& operator=(const Profile& rhs);

    inline float** get_profile() const { return p; }
    inline float** get_log_profile() const { return logp; }
    inline const size_t get_length() const { return len; }
    inline const size_t length() const { return len; }
    inline const size_t get_center() const { return center; }
    inline const size_t get_index()  const { return index; }
    inline const std::string get_name()  const { return name; }
    inline const std::string get_file()  const { return file; }
    inline const std::string get_date()  const { return date; }
    inline const char* const get_sequence()  const { return seq; }
    inline float* const get_neff_per_column() const { return neffi; }
    inline const float get_neff()  const { return neff; }

    void set_index(const size_t);
    void set_name(const std::string);
    void set_file(const std::string);

    virtual std::ostream& write(std::ostream& out);
    virtual std::ostream& write_hmm(std::ostream& out);
    virtual std::ostream& write_binary(std::ostream& out);
    virtual inline std::ostream& print_debug(std::ostream& out) { return  write(out); }

    void reset() throw (std::exception);
    void log2lin() throw (std::exception);
    void lin2log() throw (std::exception);
    void add_pseudocounts(float pca=1.0, float pcb=1.5);
    void normalize() throw (std::exception);
    void calculate_number_of_effective_sequences() throw (std::exception);

protected:
    static const int SCALE_FAC = 1000;
    const int NAA;
    const int NAA_ANY;
    const int ANY;
    const int GAP;
    const int* aa2i;
    const char* i2aa;
    const int* al2i;
    const char* al2aa;
    const float** p_cond;
    const float* p_background;

    float **p;     // matrix with amino acid probabilities
    float **logp;  // matrix with amino acid probabilities in log space
    size_t len;    // number of columns in the profile
    size_t index;  // index of profile window
    size_t center; // index of central profile column

    std::string name;
    std::string file;
    std::string date;
    char *seq;
    float *neffi;  //number of effective sequences per column
    float neff;    //number of effective sequences in profile

    void init() throw (std::exception);
    void read(std::istream& in) throw (std::exception);
    virtual void read_header(std::istream& in) throw (std::exception);
    virtual void read_columns(std::istream& in) throw (std::exception);

    virtual std::ostream& write_header(std::ostream& out);
    virtual std::ostream& write_columns(std::ostream& out);
    virtual std::ostream& write_hmm_header(std::ostream& out);
    virtual std::ostream& write_hmm_columns(std::ostream& out);

    virtual void copy_internals(const Profile& rhs) throw (std::exception);
    virtual void free_memory();
};

#endif

