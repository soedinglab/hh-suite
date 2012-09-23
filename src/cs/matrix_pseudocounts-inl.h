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

#ifndef CS_MATRIX_PSEUDOCOUNTS_INL_H_
#define CS_MATRIX_PSEUDOCOUNTS_INL_H_

#include "matrix_pseudocounts.h"

#include "count_profile-inl.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "substitution_matrix-inl.h"

namespace cs {

template<class Abc>
void MatrixPseudocounts<Abc>::AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const {
    assert_eq(seq.length(), p.length());
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to sequence ...";

    for(size_t i = 0; i < p.length(); ++i) {
        for(size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = m_.r(a, seq[i]);
    }
}

template<class Abc>
void MatrixPseudocounts<Abc>::AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const {
    assert_eq(cp.counts.length(), p.length());
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to profile ...";

    for(size_t i = 0; i < cp.counts.length(); ++i) {
        for(size_t a = 0; a < Abc::kSize; ++a) {
            double sum = 0.0;
            for(size_t b = 0; b < Abc::kSize; ++b)
                sum += m_.r(a,b) * cp.counts[i][b] / cp.neff[i];
            p[i][a] = sum;
        }
    }
}

}  // namespace cs

#endif  // CS_MATRIX_PSEUDOCOUNTS_INL_H_
