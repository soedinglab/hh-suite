// Copyright 2009, Andreas Biegert

#ifndef CS_CONTEXT_LIB_PSEUDOCOUNTS_INL_H_
#define CS_CONTEXT_LIB_PSEUDOCOUNTS_INL_H_

#include "library_pseudocounts.h"

namespace cs {

template<class Abc>
LibraryPseudocounts<Abc>::LibraryPseudocounts(const ContextLibrary<Abc>& lib,
                                              double weight_center,
                                              double weight_decay)
        : lib_(lib), emission_(lib.wlen(), weight_center, weight_decay) {}

template<class Abc>
void LibraryPseudocounts<Abc>::AddToSequence(const Sequence<Abc>& seq,
                                             const Admix& pca,
                                             Profile<Abc>& p) const {
    assert_eq(seq.length(), p.length());
    LOG(INFO) << "Adding library pseudocounts to sequence ...";

    const double tau = pca(1.0);  // effective number of sequences is one
    Matrix<double> pp(seq.length(), lib_.size(), 0.0);  // posterior probabilities
    Vector<double> pc(Abc::kSize);                      // pseudocount vector P(a|X_i)

    // Calculate and add pseudocounts for each sequence window X_i separately
    for (size_t i = 0; i < seq.length(); ++i) {
        double* ppi = &pp[i][0];
        // Calculate posterior probability of state k given sequence window around 'i'
        CalculatePosteriorProbs(lib_, emission_, seq, i, ppi);
        // Calculate pseudocount vector P(a|X_i)
        Assign(pc, 0.0);
        for (size_t k = 0; k < lib_.size(); ++k)
            for(size_t a = 0; a < Abc::kSize; ++a)
                pc[a] += ppi[k] * lib_[k].pc[a];
        // FIXME: is normalization here really needed?
        Normalize(&pc[0], Abc::kSize);
        // Add pseudocounts to sequence
        for(size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = (1.0 - tau) * (seq[i] == a ? 1.0 : 0.0) + tau * pc[a];
    }
}

template<class Abc>
void LibraryPseudocounts<Abc>::AddToProfile(const CountProfile<Abc>& cp,
                                            const Admix& pca,
                                            Profile<Abc>& p) const {
    assert_eq(cp.counts.length(), p.length());
    LOG(INFO) << "Adding library pseudocounts to profile ...";

    Matrix<double> pp(cp.counts.length(), lib_.size(), 0.0);  // posterior probs
    Vector<double> pc(Abc::kSize, 0.0);                       // pseudocount vector P(a|X_i)

    // Calculate and add pseudocounts for each sequence window X_i separately
    for (size_t i = 0; i < cp.counts.length(); ++i) {
        double* ppi = &pp[i][0];
        // Calculate posterior probability of state k given sequence window around 'i'
        CalculatePosteriorProbs(lib_, emission_, cp, i, ppi);
        // Calculate pseudocount vector P(a|X_i)
        Assign(pc, 0.0);
        for (size_t k = 0; k < lib_.size(); ++k)
            for(size_t a = 0; a < Abc::kSize; ++a)
                pc[a] += ppi[k] * lib_[k].pc[a];
        // FIXME: is normalization here really needed?
        Normalize(&pc[0], Abc::kSize);
        // Add pseudocounts to profile
        double tau = pca(cp.neff[i]);
        for(size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = (1.0 - tau) * cp.counts[i][a] / cp.neff[i] + tau * pc[a];
    }
}

}  // namespace cs

#endif  // CS_CONTEXT_LIB_PSEUDOCOUNTS_INL_H_
