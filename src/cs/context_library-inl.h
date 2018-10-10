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

#include "abstract_state_matrix-inl.h"
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
void ContextLibrary<Abc>::SetProfile(size_t k, const ContextProfile<Abc>& p) {
  assert_eq(wlen(), p.probs.length());
  assert(k < size());
  profiles_[k] = p;
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

template<class Abc>
void ContextLibrary<Abc>::SortByEntropy() {
  typedef std::pair<double, int> EntropyIndexPair;
  std::vector<EntropyIndexPair> pairs;

  for (size_t k = 0; k < profiles_.size(); ++k)
    pairs.push_back(std::make_pair(Entropy(profiles_[k].pc), k));
  sort(pairs.begin(), pairs.end());

  Vector<ContextProfile<Abc> > profiles_sorted(profiles_.size());
  for (size_t k = 0; k < pairs.size(); ++k)
    profiles_sorted[k] = profiles_[pairs[k].second];
  profiles_ = profiles_sorted;
}


template<class Abc>
void SamplingLibraryInit<Abc>::operator() (ContextLibrary<Abc>& lib) const {
  LOG(DEBUG) << "Initializing context library with by sampling " << lib.size()
             << " profile windows from training set ...";

  assert(trainset_.size() >= lib.size());
  Ran ran(seed_);
  Vector<bool> used(trainset_.size(), false);

  size_t k = 0;
  while (k < lib.size()) {
    size_t r = ran(lib.size());
    assert(r < trainset_.size());
    assert_eq(lib.wlen(), trainset_[r].counts.length());

    if (!used[r]) {
      ContextProfile<Abc> p(pc_.AddTo(trainset_[r], admix_));
      LOG(DEBUG) << p;
      p.prior = 1.0 / lib.size();
      lib.SetProfile(k, p);

      used[r] = true;
      ++k;
    }
  }

  LOG(DEBUG) << lib;
}


template<class Abc>
void GaussianLibraryInit<Abc>::operator() (ContextLibrary<Abc>& lib) const {
  Gaussian gauss(0, sigma_, seed_);

  for (size_t k = 0; k < lib.size(); ++k) {
    ContextProfile<Abc> cp(lib.wlen());
    cp.is_log = true;
    cp.prior = log(1.0 / lib.size());
    for (size_t j = 0; j < lib.wlen(); ++j) {
      for (size_t a = 0; a < Abc::kSize; ++a)
        cp.probs[j][a] = log(sm_.p(a)) + gauss();
      cp.probs[j][Abc::kAny] = 0.0;
    }
    TransformToLin(cp);
    Normalize(cp.probs, 1.0);
    lib.SetProfile(k, cp);
  }
}

template<class Abc>
void CrfBasedLibraryInit<Abc>::operator() (ContextLibrary<Abc>& lib) const {
  if (crf_.size() < lib.size())
    throw Exception("Too few context profiles for CRF initialization!");

  // Precompute column weights
  size_t c = crf_.center();
  double cw[crf_.wlen()];
  cw[c] = wcenter_;
  for (size_t j = 1; j <= c; ++j) {
    cw[c - j] = cw[c - j + 1] * wdecay_;
    cw[c + j] = cw[c - j];
  }   

  // Initialize all context profiles
  double prior[lib.size()];
  for (size_t k = 0; k < lib.size(); ++k) {
    const CrfState<Abc>& state = crf_[k];
    ContextProfile<Abc> cp(crf_.wlen());
    cp.is_log = true;

    prior[k] = 0.0;
    for (size_t i = 0; i < crf_.wlen(); ++i) {
      // Compute delta_k(i)
      double col[Abc::kSize];
      double max = -DBL_MAX;
      for (size_t a = 0; a < Abc::kSize; ++a) {
        col[a] = state.context_weights[i][a] / cw[i];
        if (col[a] > max) max = col[a];
      }
      double delta = 0.0;
      for (size_t a = 0; a < Abc::kSize; ++a) 
        delta += exp(col[a] - max);
      delta = -cw[i] * (max + log(delta));

      // log(p_k(i,a))
      for (size_t a = 0; a < Abc::kSize; ++a)
        cp.probs[i][a] = (state.context_weights[i][a] + delta) / cw[i];
      cp.probs[i][Abc::kAny] = 0.0;

      // needed for computing the prior
      prior[k] += delta;
    }
    // Set the pc column
    for (size_t a = 0; a < Abc::kSizeAny; ++a)
      cp.pc[a] = cp.probs[c][a];

    lib.SetProfile(k, cp);
  }

  // Compute log(prior)
  double max = -DBL_MAX;
  for (size_t k = 0; k < lib.size(); ++k) {
    prior[k] = crf_[k].bias_weight - neff_ * prior[k];
    if (prior[k] > max) max = prior[k];
  }
  double delta = 0.0;
  for (size_t k = 0; k < lib.size(); ++k) 
    delta += exp(prior[k] - max);
  delta = max + log(delta);
  for (size_t k = 0; k < lib.size(); ++k) {
    ContextProfile<Abc>& cp = lib[k];
    cp.prior = prior[k] - delta;
    TransformToLin(cp);
    Normalize(cp.probs, 1.0);
  }
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

template<class AS, class Abc, class CountsInput>
Profile<AS> TranslateIntoStateProfile(const CountsInput& input,
                                      const ContextLibrary<Abc>& lib,
                                      const Emission<Abc>& emission,
                                      const AbstractStateMatrix<AS>& matrix) {
  Profile<AS> asp(input.length(), 0.0);
  Vector<double> pp(lib.size());

  // For each i calculate posteriors and translate them into abstract state probs
  for (size_t i = 0; i < input.length(); ++i) {
    CalculatePosteriorProbs(lib, emission, input, i, &pp[0]);
    for (size_t a = 0; a < AS::kSize; ++a)
      for (size_t k = 0; k < lib.size(); ++k)
        asp[i][a] += pp[k] * matrix.r(a,k);
  }

  return asp;
}


// Learns a SOM from a full-blown context-library.
template<class Abc>
void LearnContextMap(const ContextLibrary<Abc>& lib,
                     ContextLibrary<Abc>& som,
                     const CoEmission<Abc>& co_emission,
                     int nsteps,    // number of learning steps
                     double sigma,  // initial neighborhood gaussian sigma
                     double alpha,  // initial learning rate
                     double tau1,   // timescale parameter for sigma
                     double tau2,   // timescale parameter for alpha
                     unsigned int seed) {
  assert_eq(lib.wlen(), som.wlen());
  const int num_colors = iround(pow(som.size(), 1.0 / 3));  // colors per channel
  if (tau1 == 0) tau1 = nsteps;
  if (tau2 == 0) tau2 = nsteps;
  Ran ran(seed);  // for picking input profiles

  // Perform 'nsteps' training iterations
  for (int n = 0; n < nsteps; ++n) {
    // Pick an input profile at random
    const ContextProfile<Abc>& x = lib[ran(lib.size())];
    // Find best matching SOM profile
    int rmax = 0, gmax = 0, bmax = 0;
    double pmax = -FLT_MAX;
    for (int r = 0; r < num_colors; ++r) {
        for (int g = 0; g < num_colors; ++g) {
          for (int b = 0; b < num_colors; ++b) {
            size_t k = r * num_colors * num_colors + g * num_colors + b;
            assert(!x.is_log);
            assert(!som[k].is_log);
            double p = co_emission(x.probs, som[k].probs);
            if (p > pmax) { rmax = r; gmax = g; bmax = b; pmax = p; }
          }
        }
    }
    // Set global learning rate 'alpha_n' and neighborhood width 'sigma_n'
    double alpha_n = alpha * exp(-n / tau2);
    double sigma_n = sigma * exp(-n / tau1);
    LOG(INFO) << strprintf("alpha[%3zu]=%8.5f  sigma[%3zu]=%8.5f",
                            n, alpha_n, n, sigma_n);
    // Update SOM profiles in neighborhood of best matching SOM profile
    for (int r = 0; r < num_colors; ++r) {
      for (int g = 0; g < num_colors; ++g) {
        for (int b = 0; b < num_colors; ++b) {
          size_t k = r * num_colors * num_colors + g * num_colors + b;
          ContextProfile<Abc>& cp = som[k];
          double sqr_dist = SQR(rmax - r) + SQR(gmax - g) + SQR(bmax - b);
          double h = exp(-sqr_dist / (2 * SQR(sigma_n)));
          LOG(INFO) << strprintf("n=%zu  dist=%8.5f  k=%zu  h=%8.5f",
                                  n, sqrt(sqr_dist), k, h);
          cp.prior += alpha_n * h * (x.prior - cp.prior);
          for (size_t i = 0; i < cp.probs.length(); ++i)
            for (size_t a = 0; a < Abc::kSize; ++a)
              cp.probs[i][a] += alpha_n * h * (x.probs[i][a] - cp.probs[i][a]);
        }
      }
    }
  }
  // Set pc probs in SOM profiles
  for (size_t k = 0; k < som.size(); ++k)
    for (size_t a = 0; a < Abc::kSize; ++a)
      som[k].pc[a] = som[k].probs[som.center()][a];
}

// Maps each context profile in 'lib' to an RBG color based color-space SOM 'som'
template<class Abc>
void AssignContextColors(ContextLibrary<Abc>& lib,
                         const ContextLibrary<Abc>& som,
                         const CoEmission<Abc>& co_emission,
                         double color_offset) {
  const int num_colors = iround(pow(som.size(), 1.0 / 3));  // colors per channel

  // Each context profile receives the color of the best matching SOM profile.
  for (size_t k = 0; k < lib.size(); ++k) {
    int rmax = 0, gmax = 0, bmax = 0;
    double pmax = -FLT_MAX;
    for (int r = 0; r < num_colors; ++r) {
      for (int g = 0; g < num_colors; ++g) {
        for (int b = 0; b < num_colors; ++b) {
          size_t l = r * num_colors * num_colors + g * num_colors + b;
          assert(!lib[k].is_log);
          assert(!som[l].is_log);
          double p = co_emission(lib[k].probs, som[l].probs);
          if (p > pmax) { rmax = r; gmax = g; bmax = b; pmax = p; }
        }
      }
    }
    // Set color value of profile 'k'
    Color color(color_offset + (1.0 - color_offset) * (rmax + 1) / num_colors,
                color_offset + (1.0 - color_offset) * (gmax + 1) / num_colors,
                color_offset + (1.0 - color_offset) * (bmax + 1) / num_colors);
    lib[k].color = color;
  }
}

// Assigns each context profile a unique name based on its position in learned SOM
template<class Abc>
void AssignContextNames(ContextLibrary<Abc>& lib,
                        const ContextLibrary<Abc>& som,
                        const CoEmission<Abc>& co_emission) {
  const char kVowels[] = { 'A', 'E', 'I', 'O', 'U', 'Y'};
  const char kConsonants[] = { 'B', 'C', 'D', 'F', 'G', 'H' };
  const int num_colors = iround(pow(som.size(), 1.0 / 3));  // colors per channel
  Vector<int> suffixes(som.size(), 0);

  // Each context profile receives the color of the best matching SOM profile.
  for (size_t k = 0; k < lib.size(); ++k) {
    int rmax = 0, gmax = 0, bmax = 0;
    double pmax = -FLT_MAX;
    for (int r = 0; r < num_colors; ++r) {
      for (int g = 0; g < num_colors; ++g) {
        for (int b = 0; b < num_colors; ++b) {
          size_t l = r * num_colors * num_colors + g * num_colors + b;
          assert(!lib[k].is_log);
          assert(!som[l].is_log);
          double p = co_emission(lib[k].probs, som[l].probs);
          if (p > pmax) { rmax = r; gmax = g; bmax = b; pmax = p; }
        }
      }
    }
    // Increment suffix counter
    size_t l = rmax * num_colors * num_colors + gmax * num_colors + bmax;
    suffixes[l]++;

    // Build name from position in lattice as "vowel + consonant + vowel + suffix"
    std::string name;
    name += kVowels[num_colors - rmax - 1];
    name += kConsonants[num_colors - gmax - 1];
    name += kVowels[num_colors - bmax - 1];
    name += strprintf("%d", suffixes[l]);

    lib[k].name = name;
  }
}

}  // namespace cs

#endif  // CS_CONTEXT_LIBRARY_INL_H_
