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

#ifndef CS_MATRIX_PSEUDOCOUNTS_H_
#define CS_MATRIX_PSEUDOCOUNTS_H_

#include "count_profile.h"
#include "profile-inl.h"
#include "pseudocounts-inl.h"
#include "sequence.h"
#include "substitution_matrix.h"

namespace cs {

// Substitution matrix pseudocounts factory.
template<class Abc>
class MatrixPseudocounts : public Pseudocounts<Abc> {
  public:
    MatrixPseudocounts(const SubstitutionMatrix<Abc>& m) : m_(m) {};

    virtual ~MatrixPseudocounts() {}

  private:
    virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const;

    virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const;

    // Substitution matrix with conditional probabilities for pseudocounts.
    const SubstitutionMatrix<Abc>& m_;

    DISALLOW_COPY_AND_ASSIGN(MatrixPseudocounts);
};  // MatrixPseudocounts

}  // namespace cs

#endif  // CS_MATRIX_PSEUDOCOUNTS_H_
