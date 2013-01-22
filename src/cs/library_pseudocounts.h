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

#ifndef CS_CONTEXT_LIB_PSEUDOCOUNTS_H_
#define CS_CONTEXT_LIB_PSEUDOCOUNTS_H_

#include "count_profile-inl.h"
#include "emission.h"
#include "profile-inl.h"
#include "pseudocounts-inl.h"
#include "sequence-inl.h"
#include "context_library-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from a library
// of context profiles.
template<class Abc>
class LibraryPseudocounts : public Pseudocounts<Abc> {
  public:
    LibraryPseudocounts(const ContextLibrary<Abc>& lib, double weight_center = 1.6, double weight_decay = 0.85);

    virtual ~LibraryPseudocounts() {}

    virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const;

    virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const;

  private:
    // Profile library with context profiles.
    const ContextLibrary<Abc>& lib_;
    // Needed to compute emission probabilities of context profiles.
    const Emission<Abc> emission_;

    DISALLOW_COPY_AND_ASSIGN(LibraryPseudocounts);
};  // LibraryPseudocounts

}  // namespace cs

#endif  // CS_CONTEXT_LIB_PSEUDOCOUNTS_H_
