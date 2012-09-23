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

#ifndef CS_CRF_STATE_H_
#define CS_CRF_STATE_H_

#include "profile_column.h"
#include "profile-inl.h"
#include "context_profile-inl.h"
#include "substitution_matrix.h"

namespace cs {

template<class Abc>
struct CrfState {
    // Default construction
    CrfState() : bias_weight(0.0) {};

    // Constructs a CRF state with 'len' columns
    explicit CrfState(size_t len)
            : bias_weight(0.0), context_weights(len) {
        assert(len & 1);
    }

    // Construction from serialized profile read from input stream.
    explicit CrfState(FILE* fin) { Read(fin); }

    // Constructs a CRF from a profile of probabilities.
    CrfState(double prior, Profile<Abc> prof, ProfileColumn<Abc> col,
             double weight_center = 1.0, double weight_decay = 1.0) {
      Init(prior, prof, col, weight_center, weight_decay);
    }

    // Constructs a CRF from a ContextProfile.
    CrfState(ContextProfile<Abc> p, double weight_center = 1.6, double weight_decay = 0.85) {
      if (p.is_log) TransformToLin(p);
      Init(p.prior, p.probs, ProfileColumn<Abc>(p.probs[(p.length() - 1) / 2]), 
          weight_center, weight_decay);
    }

    // Initializes the CRF with context weights based on the values in profile 'prof'
    // and pseudocount weights based on values in profile column 'col' with
    // prior probability 'prior'. We assume that all arguments are in lin space. The 
    // column weights are defined by wcenter and wdecay.
    void Init(double prior, Profile<Abc> prof, ProfileColumn<Abc> col,
             double weight_center = 1.0, double weight_decay = 1.0) {
        assert(prof.length() & 1);
        context_weights = Profile<Abc>(prof.length());

        Normalize(prof, 1.0);
        Normalize(col, 1.0);

        bias_weight = log(MAX(prior, DBL_MIN));

        const size_t c = (length() - 1) / 2;
        double weights[length()];
        for (size_t i = 1; i <= c; ++i) {
          double weight = weight_center * pow(weight_decay, i);
          weights[c - i] = weight;
          weights[c + i] = weight;
        }     
        weights[c] = weight_center;
        for (size_t j = 0; j < length(); ++j) {
          for (size_t a = 0; a < Abc::kSize; ++a)
            context_weights[j][a] = weights[j] * log(MAX(prof[j][a], DBL_MIN));
          context_weights[j][Abc::kAny] = 0.0;
        }

        for (size_t a = 0; a < Abc::kSize; ++a)
          pc_weights[a] = log(MAX(col[a], DBL_MIN));
        pc_weights[Abc::kAny] = 0.0;
        UpdatePseudocounts(*this);
    }

    // Initializes count profile with a serialized profile read from stream.
    void Read(FILE* fin);

    // Initializes count profile with a serialized profile read from stream.
    void Write(FILE* fin) const;

    // Returns number of context weights columns.
    inline size_t length() const { return context_weights.length(); }

    // Compares two CRF states.
    bool operator< (const CrfState<Abc>& other) const {
      return bias_weight < other.bias_weight;
    }


    std::string name;               // name of this state
    double bias_weight;             // bias weight lamda_k of this state
    Profile<Abc> context_weights;   // context weights lamda_k(j,a)
    ProfileColumn<Abc> pc_weights;  // unnormalized logs of pseudocounts
    ProfileColumn<Abc> pc;          // predicted pseudocounts at central column
};

// Prints CRF state weights in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const CrfState<Abc>& crf);

// Updates pseudocount emission probs in given CRF state based on pc_weights.
template<class Abc>
void UpdatePseudocounts(CrfState<Abc>& state);

// Calculates context score between a CRF state and a sequence window
template<class Abc>
double ContextScore(const Profile<Abc>& context_weights,
                    const Sequence<Abc>& seq,
                    size_t idx,
                    size_t center);

// Calculates context score between a CRF state and a count profile window
template<class Abc>
double ContextScore(const Profile<Abc>& context_weights,
                    const CountProfile<Abc>& cp,
                    size_t idx,
                    size_t center);

}  // namespace cs

#endif  // CS_CRF_STATE_H_
