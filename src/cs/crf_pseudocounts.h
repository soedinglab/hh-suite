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

#ifndef CS_CRF_PSEUDOCOUNTS_H_
#define CS_CRF_PSEUDOCOUNTS_H_

#include "count_profile-inl.h"
#include "emission.h"
#include "profile-inl.h"
#include "pseudocounts-inl.h"
#include "sequence-inl.h"
#include "crf-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from a context CRF
template<class Abc>
class CrfPseudocounts : public Pseudocounts<Abc> {
 public:
  CrfPseudocounts(const Crf<Abc>& lib);

  virtual ~CrfPseudocounts() {}

  virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const;

  virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const;

 private:
  // CRF with context weights and pseudocount emission weights.
  const Crf<Abc>& crf_;

  DISALLOW_COPY_AND_ASSIGN(CrfPseudocounts);
};  // CrfPseudocounts

}  // namespace cs

#endif  // CS_CRF_PSEUDOCOUNTS_H_
