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

#ifndef CS_CONTEXT_LIBRARY_INL_H_
#define CS_CONTEXT_LIBRARY_INL_H_

#include "context_profile-inl.h"
#include "pseudocounts-inl.h"
#include "context_library.h"
#include "ran.h"

namespace cs {

template<class Abc>
ContextLibrary<Abc>::ContextLibrary(size_t size, size_t wlen)
    : wlen_(wlen), profiles_(size, ContextProfile<Abc>(wlen)) {}

template<class Abc>
ContextLibrary<Abc>::ContextLibrary(FILE* fin)
    : wlen_(0), profiles_() {
  Read(fin);
}

template<class Abc>
ContextLibrary<Abc>::ContextLibrary(size_t size, size_t wlen,
                                    const LibraryInit<Abc>& init)
    : wlen_(wlen), profiles_(size, ContextProfile<Abc>(wlen)) {
  init(*this);
}

template<class Abc>
void ContextLibrary<Abc>::Read(FILE* fin) {
  // Parse and check header information
  if (!StreamStartsWith(fin, "ContextLibrary"))
      throw Exception("Stream does not start with class id 'ContextLibrary'!");

  char buffer[KB];
  size_t size = 0;
  if (cs::fgetline(buffer, KB, fin))
    size = ReadInt(buffer, "SIZE", "Unable to parse context library 'SIZE'!");
  if (cs::fgetline(buffer, KB, fin))
    wlen_ = ReadInt(buffer, "LENG", "Unable to parse context library 'LENG'!");

  // Read context profiles
  profiles_.Resize(size);
  size_t k = 0;
  for (; k < size && !feof(fin); ++k)
    profiles_[k] = ContextProfile<Abc>(fin);
  if (k != size)
    throw Exception("Serialized context library should have %i profiles but"
                    "actually has %i!", size, k);
}

template<class Abc>
void ContextLibrary<Abc>::Write(FILE* fout) const {
  // Write header
  fputs("ContextLibrary\n", fout);
  fprintf(fout, "SIZE\t%d\n", static_cast<int>(size()));
  fprintf(fout, "LENG\t%d\n", static_cast<int>(wlen()));
  // Serialize profiles
  for (size_t k = 0; k < profiles_.size(); ++k) profiles_[k].Write(fout);
}

// Transforms probabilites in context profiles to log-space and sets 'is_log' flag.
template<class Abc>
inline void TransformToLog(ContextLibrary<Abc>& lib) {
  for (size_t k = 0; k < lib.size(); ++k) TransformToLog(lib[k]);
}

// Transforms probabilites in context profiles to lin-space and sets 'is_log' flag.
template<class Abc>
inline void TransformToLin(ContextLibrary<Abc>& lib) {
  for (size_t k = 0; k < lib.size(); ++k) TransformToLin(lib[k]);
}

template<class Abc, class CountsInput, class CenterPos>
double CalculatePosteriorProbs(const ContextLibrary<Abc>& lib,
                               const Emission<Abc>& emission,
                               const CountsInput& input,
                               CenterPos i,
                               double* pp) {
  // Calculate posterior probability ppi[k] of state k given context window
  // around position 'i'
  double max = -FLT_MAX;
  for (size_t k = 0; k < lib.size(); ++k) {
    assert(lib[k].is_log);
    pp[k] = lib[k].prior + emission(lib[k].probs, input, i);
    if (pp[k] > max)
      max = pp[k];  // needed for log-sum-exp trick
  }

  // Log-sum-exp trick begins here
  double sum = 0.0;
  for (size_t k = 0; k < lib.size(); ++k)
    sum += exp(pp[k] - max);
  double tmp = max + log(sum);
  for (size_t k = 0; k < lib.size(); ++k)
    pp[k] = exp(pp[k] - tmp);

  return tmp;
}

template<class AS, class Abc, class CountsInput>
Sequence<AS> TranslateIntoStateSequence(const CountsInput& input,
                                        const ContextLibrary<Abc>& lib,
                                        const Emission<Abc>& emission) {
  Sequence<AS> as_seq(input.length());
  Vector<double> pp(AS::kSize);

  for (size_t i = 0; i < input.length(); ++i) {
    // Calculate posterior probabilities given sequence window around 'i'
    double max = -FLT_MAX;
    size_t k_max = 0;
    for (size_t k = 0; k < lib.size(); ++k) {
      assert(lib[k].is_log);
      pp[k] = lib[k].prior + emission(lib[k].probs, input, i);
      k_max = (pp[k] > max ) ? k : k_max;
      max   = (pp[k] > max ) ? pp[k] : max;
    }
    as_seq[i] = k_max;
  }
  return as_seq;
}


}  // namespace cs

#endif  // CS_CONTEXT_LIBRARY_INL_H_
