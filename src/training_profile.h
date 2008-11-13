#ifndef AB_TRAINING_PROFILE_H
#define AB_TRAINING_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//A class for training profiles for EM-clustering

#include <string>
#include <fstream>
#include <iostream>
#include <limits>   // infinity
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdio>

#include "profile.h"
#include "my_exception.h"
#include "amino_acid.h"

class TrainingProfile : public Profile {
public:
    //constructor for training profile object from serialized profile object
    TrainingProfile( std::istream& in,
                     AminoAcid* aa,
                     Matrix* m ) throw (std::exception);
   // copy constructor
    TrainingProfile(const TrainingProfile& rhs);
    // destructor
    virtual ~TrainingProfile();
    // assignment operator
    virtual TrainingProfile& operator=(const TrainingProfile& rhs);

//    inline float* const get_neff_per_column() const { return neffi; }
//    inline const float get_neff()  const { return neff; }
    inline const float get_phi()  const { return phi; }
    inline void set_phi(float p) { phi=p; }
    void lgamma_correction();

protected:
    float phi;     //how well is this training profile described by the current cluster profiles

    virtual void copy_internals(const TrainingProfile& rhs) throw (std::exception);

//    virtual void read_header(std::istream& in) throw (std::exception);
//    virtual void read_columns(std::istream& in) throw (std::exception);

//    virtual std::ostream& write_header(std::ostream& out);
//    virtual std::ostream& write_columns(std::ostream& out, size_t colwidth=8);
};

#endif
