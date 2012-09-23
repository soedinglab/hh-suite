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

#ifndef CS_EMISSION_H_
#define CS_EMISSION_H_

#include "count_profile-inl.h"
#include "profile-inl.h"
#include "profile_column.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Functor for calculating multinomial emission probabilities for context profiles.
template<class Abc>
class Emission {

  public:

    // Constructs an emission object with positional window weights.
    Emission(size_t wlen,
	     double w_center = 1.6,
	     double w_decay = 0.85,
	     const SubstitutionMatrix<Abc>* sm = NULL)
	    : center_((wlen - 1) / 2), weights_(wlen), logp_(0.0) {
	assert(wlen & 1);

	weights_[center_] = w_center;
	for (size_t i = 1; i <= center_; ++i) {
	    double weight = w_center * pow(w_decay, i);
	    weights_[center_ - i] = weight;
	    weights_[center_ + i] = weight;
	}

	if (sm) {
	    for (size_t a = 0; a < Abc::kSize; ++a)
		logp_[a] = log(sm->p(a));
	}
    }

    // Calculates the sum of positional weights.
    float GetSumWeights() const {
	double sum = 0.0;
	for (size_t i = 0; i < weights_.size(); ++i) sum += weights_[i];
	return sum;
    }

    // Calculates the log of the probability that profile 'p' emits the sequence
    // window centered at index 'idx' in 'seq'. Note that the normalization factor
    // that is usualy used in multinomial distributions is left out since it
    // cancels out anyway.
    double operator() (const Profile<Abc>& p,
		       const Sequence<Abc>& seq,
		       size_t idx) const {
	assert(p.length() & 1);
	assert_eq(weights_.size(), p.length());

	const size_t beg = MAX(0, static_cast<int>(idx - center_));
	const size_t end = MIN(seq.length(), idx + center_ + 1);
	double rv = 0.0;
	for(size_t i = beg, j = beg - idx + center_; i < end; ++i, ++j) {
	    rv += weights_[j] * (p[j][seq[i]] - logp_[seq[i]]);
	}
	return rv;
    }

    // Calculates the log of the probability that profile 'p' emits the counts in
    // a window centered at index 'idx' in 'c'. Note that the normalization factor
    // that is usualy used in multinomial distributions is left out since it
    // cancels out anyway.
    double operator() (const Profile<Abc>& p,
		       const CountProfile<Abc>& cp,
		       size_t idx) const {
	assert(p.length() & 1);
	assert_eq(weights_.size(), p.length());

	const size_t beg = MAX(0, static_cast<int>(idx - center_));
	const size_t end = MIN(cp.counts.length(), idx + center_ + 1);
	double sum, rv = 0.0;
	for(size_t i = beg, j = beg - idx + center_; i < end; ++i, ++j) {
	    sum = 0.0;
	    for (size_t a = 0; a < Abc::kSize; ++a)
		sum += cp.counts[i][a] * (p[j][a] - logp_[a]);
	    rv += weights_[j] * sum;
	}
	return rv;
    }

  private:
    size_t center_;            // index of central column in context window
    Vector<double> weights_;   // positional window weights
    ProfileColumn<Abc> logp_;  // log of background frequencies
};  // class Emission

}  // namespace cs

#endif  // CS_EMISSION_H_
