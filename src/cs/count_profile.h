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

#ifndef CS_COUNT_PROFILE_H_
#define CS_COUNT_PROFILE_H_

#include "alignment.h"
#include "sequence.h"
#include "substitution_matrix-inl.h"
#include "profile-inl.h"

namespace cs {

template<class Abc>
struct CountProfile {
    // Constructs a count profile with 'len' columns
    explicit CountProfile(size_t len = 0) : counts(len, 0.0), neff(len, 1.0) {}

    // Construction from serialized profile read from input stream.
    explicit CountProfile(FILE* fin) { Read(fin); }

    // Constructs a profile representation of given sequence.
    explicit CountProfile(const Sequence<Abc>& seq)
            : counts(seq.length(), 0.0f), neff(seq.length(), 1.0f) {
        for (size_t i = 0; i < seq.length(); ++i) counts[i][seq[i]] = 1.0f;
    }

    // Constructs a count profile from a profile and sets effective number
    // of sequences to one. Usefull for constructing a count profile from
    // a pseudocount factory.
    CountProfile(const Profile<Abc>& p) : counts(p), neff(p.length(), 1.0) {}

    // Construction from alignment with specified sequence weighting method
    CountProfile(const Alignment<Abc>& ali, bool pos_weights = true, bool neff_sum_pairs = false);

    // Creates profile by copying subprofile starting at index 'idx' for 'len' cols
    CountProfile(const CountProfile& other, size_t idx, size_t len)
            : counts(len), neff(len) {
        for (size_t i = 0; i < len; ++i) {
            neff[i] = other.neff[idx + i];
            for (size_t a = 0; a < Abc::kSizeAny; ++a)
                counts[i][a] = other.counts[idx + i][a];
        }
    }

    // Initializes count profile with a serialized profile read from stream.
    void Read(FILE* fin);

    // Writes serialized count profile to stream.
    void Write(FILE* fout) const;

    // Returns number of columns.
    size_t length() const { return counts.length(); }


    std::string name;               // optional name descriptor
    Profile<Abc> counts;            // absolute counts of alphabet letters
    Vector<double> neff;            // effective number of sequences at column i
};

// Returns the average number of effective sequences in given count profile
template<class Abc>
double Neff(const CountProfile<Abc>& cp);

// Builds and returns a consensus sequence of the given count profile by
// calculating at each position the alphabet character that deviates most strongly
// from its background probability.
template<class Abc>
std::string ConsensusSequence(const CountProfile<Abc>& cp,
                              const SubstitutionMatrix<Abc>& sm);

// Builds and returns a conservation string for given count profile that
// indicates conservation of residues by uppercase, lowercase, and '~'
template<class Abc>
std::string ConservationSequence(const CountProfile<Abc>& cp,
                                 const SubstitutionMatrix<Abc>& sm);

// Prints counts and neff in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const CountProfile<Abc>& cp);


}  // namespace cs

#endif  // CS_COUNT_PROFILE_H_
