/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "training_profile.h"

// constructor for training profile object from input stream
TrainingProfile::TrainingProfile( std::istream& in,
                                  AminoAcid* aa,
                                  Matrix* m) throw (std::exception) :
    Profile::Profile(aa, m),
    phi(0.0) {
    read(in);
}

TrainingProfile::~TrainingProfile(){
}

TrainingProfile::TrainingProfile(const TrainingProfile& rhs) { copy_internals(rhs); }

TrainingProfile& TrainingProfile::operator=(const TrainingProfile& rhs) {
    // identity test, if self-assignment, do nothing
    if (this == &rhs) return *this;
    copy_internals(rhs);
    return *this;
}

void TrainingProfile::copy_internals(const TrainingProfile& rhs) throw (std::exception)  {
    Profile::copy_internals(rhs);
    phi = rhs.phi;
}

void TrainingProfile::lgamma_correction() {
    for(size_t i=0; i<len; ++i)
        for(int a=0; a<NAA; ++a)
            logp[i][a] = (fast_log_gamma(neffi[i]*p[i][a]+1.0)-fast_log_gamma(neffi[i]+1)/static_cast<float>(NAA)) / (neffi[i]*p[i][a]);
}
