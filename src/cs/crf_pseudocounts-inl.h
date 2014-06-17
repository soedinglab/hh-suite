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

#ifndef CS_CRF_PSEUDOCOUNTS_INL_H_
#define CS_CRF_PSEUDOCOUNTS_INL_H_

#include "crf_pseudocounts.h"

namespace cs {

template<class Abc>
CrfPseudocounts<Abc>::CrfPseudocounts(const Crf<Abc>& crf) : crf_(crf) {}

template<class Abc>
void CrfPseudocounts<Abc>::AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const {
  assert_eq(seq.length(), p.length());
  LOG(INFO) << "Adding CRF pseudocounts to sequence ...";

  const size_t center = crf_.center();
  Matrix<double> pp(seq.length(), crf_.size(), 0.0);  // posterior probabilities
  int len = static_cast<int>(seq.length());

  // Calculate and add pseudocounts for each sequence window X_i separately
#pragma omp parallel for schedule(static)
  for (int i = 0; i < len; ++i) {
    double* ppi = &pp[i][0];

    // Calculate posterior probability ppi[k] of state k given sequence window
    // around position 'i'
    double max = -DBL_MAX;
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = crf_[k].bias_weight +
        ContextScore(crf_[k].context_weights, seq, i, center);
      if (ppi[k] > max)
        max = ppi[k];  // needed for log-sum-exp trick
    }

    // Log-sum-exp trick begins here
    double sum = 0.0;
    for (size_t k = 0; k < crf_.size(); ++k)
      sum += exp(ppi[k] - max);
    double tmp = max + log(sum);
    // Calculate pseudocount vector P(a|X_i)
    double* pc = p[i];
    for (size_t a = 0; a < Abc::kSize; ++a) pc[a] = 0.0;
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = exp(ppi[k] - tmp);
      for(size_t a = 0; a < Abc::kSize; ++a)
        pc[a] += ppi[k] * crf_[k].pc[a];
    }
    Normalize(&pc[0], Abc::kSize);
  }
}

template<class Abc>
void CrfPseudocounts<Abc>::AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const {
  assert_eq(cp.counts.length(), p.length());
  LOG(INFO) << "Adding library pseudocounts to profile ...";

  const size_t center = crf_.center();
  Matrix<double> pp(cp.counts.length(), crf_.size(), 0.0);  // posterior probs
  int len = static_cast<int>(cp.length());

  // Calculate and add pseudocounts for each sequence window X_i separately
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < len; ++i) {
    double* ppi = &pp[i][0];

    // Calculate posterior probability ppi[k] of state k given sequence window
    // around position 'i'
    double max = -DBL_MAX;
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = crf_[k].bias_weight +
        ContextScore(crf_[k].context_weights, cp, i, center);
      if (ppi[k] > max)
        max = ppi[k];  // needed for log-sum-exp trick
    }

     // Log-sum-exp trick begins here
    double sum = 0.0;
    for (size_t k = 0; k < crf_.size(); ++k)
      sum += exp(ppi[k] - max);
    double tmp = max + log(sum);
    // Calculate pseudocount vector P(a|X_i)
    double* pc = p[i];
    for (size_t a = 0; a < Abc::kSize; ++a) pc[a] = 0.0;
    for (size_t k = 0; k < crf_.size(); ++k) {
      ppi[k] = exp(ppi[k] - tmp);
      for(size_t a = 0; a < Abc::kSize; ++a)
        pc[a] += ppi[k] * crf_[k].pc[a];
    }
    Normalize(&pc[0], Abc::kSize);  
  }
}

}  // namespace cs

#endif  // CS_LIBRARY_PSEUDOCOUNTS_INL_H_
