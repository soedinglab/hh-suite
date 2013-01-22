/*
  Copyright 2009 Andreas Biegert

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

#ifndef CS_CONTEXT_LIBRARY_H_
#define CS_CONTEXT_LIBRARY_H_

#include "abstract_state_matrix.h"
#include "co_emission.h"
#include "context_profile.h"
#include "pseudocounts-inl.h"

namespace cs {

// Forward declarations
template<class Abc>
class ContextLibrary;

template<class Abc>
class Emission;

template<class Abc>
class Crf;

template<class Abc>
class CrfState;

// Strategy class for initializing a context library
template<class Abc>
class LibraryInit {
 public:
  LibraryInit() {}
  virtual ~LibraryInit() {}
  virtual void operator() (ContextLibrary<Abc>& lib) const = 0;
};

// A container of K context profiles representing the most common
// sequence motifs in a training database of proteins/DNA sequences.
template<class Abc>
class ContextLibrary {
 public:
  typedef ContextProfile<Abc>* ProfileIter;
  typedef const ContextProfile<Abc>* ConstProfileIter;

  // Constructs an empty profile library of given dimenions.
  ContextLibrary(size_t size, size_t wlen);

  // Constructs a profile library from serialized data read from input stream.
  explicit ContextLibrary(FILE* fin);

  // Constructs profile library with a specific init-strategy encapsulated by an
  // initializer.
  ContextLibrary(size_t size, size_t wlen, const LibraryInit<Abc>& init);

  // Nothing to do here
  virtual ~ContextLibrary() {}

  // Returns the number of profiles in the fully assembled profile library
  size_t size() const { return profiles_.size(); }

  // Returns the number of columns in each context profile.
  size_t wlen() const { return wlen_; }

  // Returns index of central profile column.
  size_t center() const { return (wlen_ - 1) / 2; }

  // Accessor methods for state i, where i is from interval [0,size].
  ContextProfile<Abc>& operator[](size_t i) { return profiles_[i]; }
  const ContextProfile<Abc>& operator[](size_t i) const { return profiles_[i]; }

  // Initializes profile at index 'k' with given profile.
  void SetProfile(size_t k, const ContextProfile<Abc>& p);

  // Returns an iterator pointing to beginning of profiles.
  ProfileIter begin() { return &profiles_[0]; }

  // Returns an iterator pointing past the end of profiles.
  ProfileIter end() { return &profiles_[0] + profiles_.size(); }

  // Returns a const iterator pointing to beginning of profiles.
  ConstProfileIter begin() const { return &profiles_[0]; }

  // Returns a const iterator pointing past the end of profiles.
  ConstProfileIter end() const { return &profiles_[0] + profiles_.size(); }

  // Writes the profile library in serialization format to output stream.
  void Write(FILE* fout) const;

  // Sorts context states by relative entropy of central column and assigns
  // new state indices according to this new ordering
  void SortByEntropy();

 private:
   // Initializes the library from serialized data read from stream.
  void Read(FILE* fin);

  size_t wlen_;                            // size of context window.
  Vector<ContextProfile<Abc> > profiles_;  // context profiles ordered by index.
};  // ContextLibrary


// Prints the library in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const ContextLibrary<Abc>& lib) {
  out << "ContextLibrary" << std::endl;
  out << "size:\t" << lib.size() << std::endl;
  out << "wlen:\t" << lib.wlen() << std::endl;
  for (size_t k = 0; k < lib.size(); ++k) out << lib[k];
  return out;
}

// Transforms probabilites in context profiles to log-space and sets 'is_log' flag.
template<class Abc>
void TransformToLog(ContextLibrary<Abc>& lib);

// Transforms probabilites in context profiles to lin-space and sets 'is_log' flag.
template<class Abc>
void TransformToLin(ContextLibrary<Abc>& lib);

// Calculates posterior probs for a context library and sequence window X_i
// centered at index 'i' and writes them to array 'pp'. Caller is responsible for
// making sure that 'pp' has sufficient length. Return value is log sum of all
// individual emission terms. The third template parameter specifies the central
// position of the context window. Note: For PO-HMMs this is not an ordinary size_t
// but a vertex descriptor.
template<class Abc, class ContextInput, class CenterPos>
double CalculatePosteriorProbs(const ContextLibrary<Abc>& lib,
                               const Emission<Abc>& emission,
                               const ContextInput& input,
                               CenterPos i,
                               double* pp);

// Strategy for initializing library by sampling from training set of count
// profiles, optionally adding pseudocounts.
template<class Abc>
class SamplingLibraryInit : public LibraryInit<Abc> {
 public:
  typedef std::vector< CountProfile<Abc> > TrainingSet;

  SamplingLibraryInit(const TrainingSet& trainset,
                      const Pseudocounts<Abc>& pc,
                      const Admix& admix,
                      unsigned int seed = 0)
      : trainset_(trainset),
        pc_(pc),
        admix_(admix),
        seed_(seed) {}

  virtual ~SamplingLibraryInit() {}

  virtual void operator() (ContextLibrary<Abc>& lib) const;

 private:
  const TrainingSet& trainset_;
  const Pseudocounts<Abc>& pc_;
  const Admix& admix_;
  const unsigned int seed_;
};  // SamplingLibraryInit


// Strategy that initializes profile probs by sammpling from gaussian distribution
// with mean at background frequencies.
template<class Abc>
class GaussianLibraryInit : public LibraryInit<Abc> {
 public:
  GaussianLibraryInit(double sigma,
                      const SubstitutionMatrix<Abc>& sm,
                      unsigned int seed = 0)
      : sigma_(sigma), sm_(sm), seed_(seed) {}

  virtual ~GaussianLibraryInit() {}

  virtual void operator() (ContextLibrary<Abc>& lib) const;

 protected:
  double sigma_;
  const SubstitutionMatrix<Abc>& sm_;
  unsigned int seed_;
};  // class GaussianLibraryInit

// Uses a CRF for initialization
template<class Abc>
class CrfBasedLibraryInit : public LibraryInit<Abc> {
 public:
  CrfBasedLibraryInit(
      const Crf<Abc>& crf, 
      double wcenter = 1.6,
      double wdecay  = 0.85,
      double neff    = 1.0) :
    crf_(crf), wcenter_(wcenter), wdecay_(wdecay), neff_(neff) {}

  virtual ~CrfBasedLibraryInit() {}

  virtual void operator() (ContextLibrary<Abc>& lib) const;

 protected:
  const Crf<Abc>& crf_;
  const double wcenter_;
  const double wdecay_;
  const double neff_;
};  // class CrfBasedLibraryInit

// Translate a sequence or count profile into an abstract state sequence.
template<class AS, class Abc, class CountsInput>
Sequence<AS> TranslateIntoStateSequence(const CountsInput& input,
                                        const ContextLibrary<Abc>& lib,
                                        const Emission<Abc>& emission);

// Translate a sequence or count profile into an abstract state profile given
// the learned context-to-state mutation probabilities in 'matrix'.
template<class AS, class Abc, class CountsInput>
Profile<AS> TranslateIntoStateProfile(const CountsInput& input,
                                      const ContextLibrary<Abc>& lib,
                                      const Emission<Abc>& emission,
                                      const AbstractStateMatrix<AS>& matrix);

// Learns a color-space SOM from a full-blown context-library.
template<class Abc>
void LearnContextMap(const ContextLibrary<Abc>& lib,
                     ContextLibrary<Abc>& som,
                     const CoEmission<Abc>& co_emission,
                     int nsteps,          // number of learning steps
                     double sigma = 5.0,  // initial neighborhood gaussian sigma
                     double alpha = 0.1,  // initial learning rate
                     double tau1 = 0.0,   // timescale parameter for sigma
                     double tau2 = 0.0,   // timescale parameter for alpha
                     unsigned int seed = 0);

// Assigns each context profile in given context lib an RBG color based learned SOM
template<class Abc>
void AssignContextColors(ContextLibrary<Abc>& lib,
                         const ContextLibrary<Abc>& som,
                         const CoEmission<Abc>& co_emission,
                         double color_offset = 0.2);

// Assigns each context profile a unique name based on its position in learned SOM
template<class Abc>
void AssignContextNames(ContextLibrary<Abc>& lib,
                        const ContextLibrary<Abc>& som,
                        const CoEmission<Abc>& co_emission);

}  // namespace cs

#endif  // CS_CONTEXT_LIBRARY_H_
