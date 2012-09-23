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

using std::vector;

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


// Strategy for initializing CRF by sampling from training set, optionally
// adding pseudocounts.
template<class Abc, class TrainingPair>
class SamplingCrfInit : public CrfInit<Abc> {
  public:
    typedef std::vector<TrainingPair> TrainingSet;

    SamplingCrfInit(const TrainingSet& trainset,
                    const Pseudocounts<Abc>& pc,
                    Admix& admix,
                    const SubstitutionMatrix<Abc>& sm,
                    unsigned int seed = 0, 
                    double weight_center = 1.6,
                    double weight_decay = 0.85)
            : trainset_(trainset),
              pc_(pc),
              admix_(admix),
              sm_(sm),
              seed_(seed),
              ran_(seed),
              weight_center_(weight_center),
              weight_decay_(weight_decay) {}

    virtual ~SamplingCrfInit() {}

    virtual void operator() (Crf<Abc>& crf);

  private:
    const TrainingSet& trainset_;
    const Pseudocounts<Abc>& pc_;
    Admix& admix_;
    const SubstitutionMatrix<Abc>& sm_;
    const unsigned int seed_;
    Ran ran_;
    const double weight_center_;
    const double weight_decay_;
};  // SamplingCrfInit


// FIXME: why doesn't this compile?!
// Comparison function for sorting context profiles by prior in descending order
// template<class Abc>
// struct ContextProfilePriorCompare :
//       public std::binary_function<ContextProfile<Abc>, ContextProfile<Abc>, bool> {
//   bool operator() (const ContextProfile<Abc>& lh,
//                    const ContextProfile<Abc>& rh) const {
//     return rh.prior < lh.prior;
//   }
// };

// template<class Abc>
// bool ContextProfilePriorCompare(const ContextProfile<Abc>& lhs,
//                                 const ContextProfile<Abc>& rhs) {
//   return rhs.prior < lhs.prior;
// }

// Strategy that uses CRF states from a CRF to initialize
// CRF states.
template<class Abc>
class CrfBasedCrfInit : public CrfInit<Abc> {
  public:
    CrfBasedCrfInit(const Crf<Abc>& crf)
            : crf_states_(crf.begin(), crf.end()) {}

    virtual ~CrfBasedCrfInit() {}

    virtual void operator() (Crf<Abc>& crf) {
      if (crf_states_.size() < crf.size())
          throw Exception("Too few CRF states for CRF initialization!");
      vector<const CrfState<Abc>* > crf_states;
      for (size_t i = 0; i < crf_states_.size(); ++i)
        crf_states.push_back(&crf_states_[i]);
      if (crf_states_.size() > crf.size()) {
        // Use CRF states with the highest probability
        sort(crf_states.begin(), crf_states.end(), CompareCrfStates);
      }
      for (size_t k = 0; k < crf.size(); ++k)
        crf.SetState(k, *crf_states[k]);
    }

  private:

    static bool CompareCrfStates(const CrfState<Abc>* p, const CrfState<Abc>* q) {
      return p->bias_weight > q->bias_weight;
    }
        
    const vector<CrfState<Abc> > crf_states_;
};  // class CrfBasedCrfInit

// Strategy that uses context profiles from a profile library to initialize
// CRF states.
template<class Abc>
class LibraryBasedCrfInit : public CrfInit<Abc> {
  public:
    LibraryBasedCrfInit(const ContextLibrary<Abc>& lib,
                        double weight_center = 1.6,
                        double weight_decay  = 0.85,
                        unsigned int seed = 0)
            : profiles_(lib.begin(), lib.end()),
              weight_center_(weight_center), weight_decay_(weight_decay), seed_(seed) {}

    virtual ~LibraryBasedCrfInit() {}

    virtual void operator() (Crf<Abc>& crf) {
      if (profiles_.size() < crf.size())
          throw Exception("Too few profiles in context library for CRF initialization!");
      vector<const ContextProfile<Abc>* > crf_profiles;
      for (size_t i = 0; i < profiles_.size(); ++i)
        crf_profiles.push_back(&profiles_[i]);
      if (profiles_.size() > crf.size()) {
        // Use context profiles with the highest probability
        sort(crf_profiles.begin(), crf_profiles.end(), CompareProfiles);
      }
      for (size_t k = 0; k < crf.size(); ++k)
        crf.SetState(k, CrfState<Abc>(*crf_profiles[k], weight_center_, weight_decay_));
    }

  private:

    static bool CompareProfiles(const ContextProfile<Abc>* p, const ContextProfile<Abc>* q) {
      return p->prior > q->prior;
    }
        
    const vector<ContextProfile<Abc> > profiles_;
    const SubstitutionMatrix<Abc>* sm_;
    const double weight_center_;
    const double weight_decay_;
    const unsigned int seed_;
};  // class LibraryBasedCrfInit


// Strategy that initializes CRF weights by sammpling from gaussian distribution
template<class Abc>
class GaussianCrfInit : public CrfInit<Abc> {
  public:
    GaussianCrfInit(double sigma,
                    const SubstitutionMatrix<Abc>& sm,
                    unsigned int seed = 0)
            : sigma_(sigma), sm_(sm), seed_(seed), gaussian_(0, sigma, seed) {}

    virtual ~GaussianCrfInit() {}

    virtual void operator() (Crf<Abc>& crf);

  protected:
    double sigma_;
    const SubstitutionMatrix<Abc>& sm_;
    unsigned int seed_;
    Gaussian gaussian_;

};  // class GaussianCrfInit

}  // namespace cs

#endif  // CS_CRF_H_
