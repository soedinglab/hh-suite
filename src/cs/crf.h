/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_CRF_H_
#define CS_CRF_H_

#include "count_profile.h"
#include "context_profile-inl.h"
#include "context_library-inl.h"
#include "crf_state.h"
#include "pseudocounts-inl.h"
#include "sequence.h"

namespace cs {

// Forward declarations
template<class Abc>
class Crf;

// Strategy class for initializing a CRF
template<class Abc>
class CrfInit {
  public:
    CrfInit() {}
    virtual ~CrfInit() {}
    virtual void operator() (Crf<Abc>& crf) = 0;
};


// A container class for CRF states to represent the most common
// sequence motifs in a training database of proteins/DNA sequences.
template<class Abc>
class Crf {
  public:
    typedef CrfState<Abc>* StateIter;
    typedef const CrfState<Abc>* ConstStateIter;

    // Constructs an empty CRF of given dimenions.
    Crf(size_t size, size_t wlen);

    // Constructs a CRF from serialized data read from input stream.
    explicit Crf(FILE* fin);

    // Constructs CRF with a specific init-strategy encapsulated by an
    // initializer.
    Crf(size_t size, size_t wlen, CrfInit<Abc>& init);

    // Constructs CRF using a context library.
    Crf(const ContextLibrary<Abc>& lib, double weight_center = 1.6, double weight_decay = 0.85)
              : wlen_(lib.wlen()), states_(lib.size(), CrfState<Abc>(lib.wlen())) {
        for (size_t i = 0; i < lib.size(); ++i) 
            states_[i] = CrfState<Abc>(lib[i], weight_center, weight_decay);                         
    } 

    // Deallocates states in profile vector
    virtual ~Crf() {}

    // Returns the number of profiles in the fully assembled CRF
    size_t size() const { return states_.size(); }

    // Returns the number of columns in each context profile.
    size_t wlen() const { return wlen_; }

    // Returns total number of weights in this CRF. Note that context weights
    // and pseudocount weights of letter ANY are not accounted for since these are
    // held fix at zero anyway.
    size_t nweights() const { return size() * (1 + (wlen() + 1) * Abc::kSize); }

    // Returns index of central profile column.
    size_t center() const { return (wlen_ - 1) / 2; }

    // Accessor methods for state i, where i is from interval [0,size].
    CrfState<Abc>& operator[](size_t i) { return states_[i]; }
    const CrfState<Abc>& operator[](size_t i) const { return states_[i]; }

    // Initializes profile at index 'idx' with given profile.
    void SetState(size_t idx, const CrfState<Abc>& s);

    // Returns an iterator to a list of pointers to profiles.
    StateIter begin() { return &states_[0]; }

    // Returns an iterator pointing past the end of pointers to profiles.
    StateIter end() { return &states_[0] + states_.size(); }

    // Returns a const iterator over pointers of profiles.
    ConstStateIter begin() const { return &states_[0]; }

    // Returns a const iterator pointing past the end of pointers to profiles.
    ConstStateIter end() const { return &states_[0] + states_.size(); }

    // Writes the CRF in serialization format to output stream.
    void Write(FILE* fout) const;

  private:
    // Initializes the library from serialized data read from stream.
    void Read(FILE* fin);

    size_t wlen_;                           // size of context window.
    Vector< CrfState<Abc> > states_;  // states ordered by index.
};  // Crf


// Prints the library in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const Crf<Abc>& crf) {
    out << "CRF" << std::endl;
    out << "size:\t" << crf.size() << std::endl;
    out << "wlen:\t" << crf.wlen() << std::endl;
    for (size_t k = 0; k < crf.size(); ++k) out << crf[k];
    return out;
}

}  // namespace cs

#endif  // CS_CRF_H_
