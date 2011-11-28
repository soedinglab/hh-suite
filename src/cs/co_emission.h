// Copyright 2009, Andreas Biegert

#ifndef CS_CO_EMISSION_H_
#define CS_CO_EMISSION_H_

#include "profile-inl.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Functor for calculating co-emission probabilities for two profiles.
template<class Abc>
class CoEmission {
public:

  // Constructs an emission object with positional window weights.
  CoEmission(size_t wlen,
             const SubstitutionMatrix<Abc>& sm,
             double w_center = 1.6,
             double w_decay = 0.85)
      : center_((wlen - 1) / 2), weights_(wlen), sm_(sm) {
    assert(wlen & 1);

    weights_[center_] = w_center;
    for (size_t i = 1; i <= center_; ++i) {
      double weight = w_center * pow(w_decay, i);
      weights_[center_ - i] = weight;
      weights_[center_ + i] = weight;
    }
  }

  // Calculates the sum of positional weights.
  float GetSumWeights() const {
    double sum = 0.0;
    for (size_t i = 0; i < weights_.size(); ++i) sum += weights_[i];
    return sum;
  }

  // Calculates weighted sum of log co-emissions probs over all profile columns.
  // Probabilities in both profiles are assumed to be in lin-space.
  double operator() (const Profile<Abc>& p,
                     const Profile<Abc>& q) const {
    assert(p.length() & 1);
    assert_eq(weights_.size(), p.length());
    assert_eq(p.length(), q.length());

    double rv = 0.0;
    for (size_t j = 0; j < weights_.size(); ++j) {
      double prob_j = 0.0;
      for (size_t a = 0; a < Abc::kSize; ++a)
        prob_j += (p[j][a] * q[j][a]) / sm_.p(a);
      rv += weights_[j] * log(prob_j);
    }

    return rv;
  }

 private:
  size_t center_;                      // index of central column in context window
  Vector<double> weights_;             // positional window weights
  const SubstitutionMatrix<Abc>& sm_;  // neeeded for background frequencies

  DISALLOW_COPY_AND_ASSIGN(CoEmission);
};  // class CoEmission

}  // namespace cs

#endif  // CS_CO_EMISSION_H_
