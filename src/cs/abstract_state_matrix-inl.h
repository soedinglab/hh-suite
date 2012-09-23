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

#ifndef CS_ABSTRACT_STATE_MATRIX_INL_H_
#define CS_ABSTRACT_STATE_MATRIX_INL_H_

#include "abstract_state_matrix.h"

namespace cs {

template<class AS>
AbstractStateMatrix<AS>::AbstractStateMatrix(std::string matrixfile)
    : num_contexts_(0), py_(AS::kSizeAny, 0.0) {
  // Read matrix file
  char buffer[MB];
  FILE* fin = fopen(matrixfile.c_str(), "r");
  if (!fin)
    throw Exception("Can't read abstract state matrix '%s'!", matrixfile.c_str());

  // Determine number of contexts
  fgetline(buffer, MB, fin);  // skip alphabet description line
  while (fgetline(buffer, KB, fin))
    if (strscn(buffer)) num_contexts_++;

  // Resize and assign matrices and vectors
  q_.Assign(num_contexts_, AS::kSizeAny, 0.0);
  s_.Assign(num_contexts_, AS::kSizeAny, 0.0);
  r_.Assign(AS::kSizeAny, num_contexts_, 0.0);
  px_.Assign(num_contexts_, 0.0);

  rewind(fin);
  const char* ptr = NULL;
  size_t k = 0;
  fgetline(buffer, MB, fin);  // skip alphabet description line
  while (fgetline(buffer, KB, fin)) {
    ptr = buffer;
    for (size_t s = 0; s < AS::kSize; ++s)
      q_[k][s] = strtof(ptr);
    k++;
  }
  fclose(fin);

  // Init remaining members based on target freqs
  Init();
}

template<class AS>
void AbstractStateMatrix<AS>::Init() {
  // Check transition probability matrix, renormalize P
  double sumab = 0.0;
  for (size_t a = 0; a < num_contexts_; a++)
    for (size_t b = 0; b < AS::kSize; ++b) sumab += q_[a][b];
  for (size_t a = 0; a < num_contexts_; ++a)
    for (size_t b = 0; b < AS::kSize; ++b) q_[a][b] /= sumab;

  // Calculate background frequencies
  for (size_t a = 0; a < num_contexts_; ++a) {
    px_[a] = 0.0;
    for (size_t b = 0; b < AS::kSize; ++b) px_[a] += q_[a][b];
  }
  Normalize(&px_[0], num_contexts_);
  for (size_t b = 0; b < AS::kSize; ++b) {
    py_[b] = 0.0;
    for (size_t a = 0; a < num_contexts_; ++a) py_[b] += q_[a][b];
  }
  Normalize(&py_[0], AS::kSize);

  // Precompute conditional probability matrix R "abstract state b given context a"
  for (size_t a = 0; a < num_contexts_; ++a)
    for (size_t b = 0; b < AS::kSize; ++b)
      r_[b][a] = q_[a][b] / px_[a]; // P(b|a)
  for (size_t a = 0; a < num_contexts_; ++a)
    r_[AS::kAny][a] = 1.0; // P(X|a) = 1.0

  // Calculate scoring matrix as S[a][b] = (1 / lambda) * log(P(a,b) / (P(a)*P(b)))
  const double lambda = log(2) / 2.0;
  for (size_t a = 0; a < num_contexts_; ++a)
    for (size_t b = 0; b < AS::kSize; ++b)
      s_[a][b] = (1.0 / lambda) * log(q_[a][b] / (px_[a] * py_[b]));

  LOG(DEBUG1) << *this;
}

template<class AS>
void AbstractStateMatrix<AS>::Print(std::ostream& out) const {
  out << "Background frequencies on query side (in %):\n";
  for (size_t a = 0; a < num_contexts_; ++a)
    out << a << "\t";
  out << std::endl;
  for (size_t a = 0; a < num_contexts_; ++a)
    out << strprintf("%-.4f\t", 100.0 * px_[a]);

  out << "\nBackground frequencies on database side (in %):\n";
  for (size_t a = 0; a < AS::kSize; ++a)
    out << AS::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < AS::kSize; ++a)
    out << strprintf("%-.2f\t", 100.0 * py_[a]);

  out << "\nSubstitution matrix log( q(k,s) / p(k)*p(s) ) (scaled by 1/lambda):\n";
  for (size_t b = 0; b < AS::kSize; ++b)
    out << AS::kIntToChar[b] << "\t";
  out << std::endl;
  for (size_t a = 0; a < num_contexts_; ++a) {
    for (size_t b = 0; b < AS::kSize; ++b)
      out << strprintf("%+5.2f\t", s_[a][b]);
    out << std::endl;
  }

  out << "Target frequency matrix q(k,s) (in %):\n";
  for (size_t b = 0; b < AS::kSize; ++b)
    out << AS::kIntToChar[b] << "\t";
  out << std::endl;
  for (size_t a = 0; a < num_contexts_; ++a) {
    for (size_t b = 0; b < AS::kSize; ++b)
      out << strprintf("%-.2f\t", 100.0f * q_[a][b]);
    out << std::endl;
  }

  out << "Matrix of conditional probabilities P(s|k) = P(k,s)/p(k) (in %):\n";
  for (size_t b = 0; b < AS::kSize; ++b)
    out << AS::kIntToChar[b] << "\t";
  out << std::endl;
  for (size_t a = 0; a < num_contexts_; ++a) {
    for (size_t b = 0; b < AS::kSize; ++b)
      out << strprintf("%-.1f\t", 100.0f * r_[b][a]);
    out << std::endl;
  }
}

}  // namespace cs

#endif  // CS_ABSTRACT_STATE_MATRIX_INL_H_
