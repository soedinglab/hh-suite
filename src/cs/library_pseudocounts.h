// Copyright 2009, Andreas Biegert

#ifndef CS_CONTEXT_LIB_PSEUDOCOUNTS_H_
#define CS_CONTEXT_LIB_PSEUDOCOUNTS_H_

#include "count_profile-inl.h"
#include "emission.h"
#include "profile-inl.h"
#include "pseudocounts.h"
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

    virtual void AddToSequence(const Sequence<Abc>& seq, const Admix& pca, Profile<Abc>& p) const;

    virtual void AddToProfile(const CountProfile<Abc>& cp, const Admix& pca, Profile<Abc>& p) const;

  private:
    // Profile library with context profiles.
    const ContextLibrary<Abc>& lib_;
    // Needed to compute emission probabilities of context profiles.
    const Emission<Abc> emission_;

    DISALLOW_COPY_AND_ASSIGN(LibraryPseudocounts);
};  // LibraryPseudocounts

}  // namespace cs

#endif  // CS_CONTEXT_LIB_PSEUDOCOUNTS_H_
