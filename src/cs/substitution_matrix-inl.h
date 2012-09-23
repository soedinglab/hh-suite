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

#ifndef CS_SUBSTITUTION_MATRIX_INL_H_
#define CS_SUBSTITUTION_MATRIX_INL_H_

#include "substitution_matrix.h"

namespace cs {

template<class Abc>
SubstitutionMatrix<Abc>::SubstitutionMatrix(double l)
    : q_(Abc::kSizeAny, Abc::kSizeAny, 0.0),
      s_(Abc::kSizeAny, Abc::kSizeAny, 0.0),
      rx_(Abc::kSizeAny, Abc::kSizeAny, 0.0),
      ry_(Abc::kSizeAny, Abc::kSizeAny, 0.0),
      px_(Abc::kSizeAny, 0.0),
      py_(Abc::kSizeAny, 0.0),
      lambda_(l) {}

template<class Abc>
SubstitutionMatrix<Abc>::~SubstitutionMatrix() {}

template<class Abc>
void SubstitutionMatrix<Abc>::InitFromTargetFreqs() {
  // Check transition probability matrix, renormalize P
  double sumab = 0.0;
  for (size_t a = 0; a < Abc::kSize; a++)
    for (size_t b = 0; b < Abc::kSize; ++b) sumab += q_[a][b];
  for (size_t a = 0; a < Abc::kSize; ++a)
    for (size_t b = 0; b < Abc::kSize; ++b) q_[a][b] /= sumab;

  // Calculate background frequencies
  for (size_t a = 0; a < Abc::kSize; ++a) {
    px_[a] = 0.0;
    for (size_t b = 0; b < Abc::kSize; ++b) px_[a] += q_[a][b];
  }
  Normalize(&px_[0], Abc::kSize);
  for (size_t b = 0; b < Abc::kSize; ++b) {
    py_[b] = 0.0;
    for (size_t a = 0; a < Abc::kSize; ++a) py_[b] += q_[a][b];
  }
  Normalize(&py_[0], Abc::kSize);

  // Precompute conditional probability matrix Px
  for (size_t a = 0; a < Abc::kSize; ++a)
    for (size_t b = 0; b < Abc::kSize; ++b)
      rx_[b][a] = q_[a][b] / px_[a]; // Px(b|a)
  for (size_t b = 0; b < Abc::kSize; ++b)
    rx_[b][Abc::kAny] = py_[b]; // Px(b|X) = Py(b)
  for (size_t a = 0; a < Abc::kSize; ++a)
    rx_[Abc::kAny][a] = 1.0; // Px(X|a) = 1.0

  // Precompute conditional probability matrix Py
  for (size_t a = 0; a < Abc::kSize; ++a)
    for (size_t b = 0; b < Abc::kSize; ++b)
      ry_[a][b] = q_[a][b] / py_[b]; // Py(a|b)
  for (size_t a = 0; a < Abc::kSize; ++a)
    ry_[a][Abc::kAny] = px_[a]; // Py(a|X) = Px(a)
  for (size_t b = 0; b < Abc::kSize; ++b)
    ry_[Abc::kAny][b] = 1.0; // Py(X|b) = 1.0

  // Calculate scoring matrix as S[a][b] = (1 / lambda) * log(P(a,b) / (P(a)*P(b)))
  for (size_t a = 0; a < Abc::kSize; ++a)
    for (size_t b = 0; b < Abc::kSize; ++b)
      s_[a][b] = (1.0 / lambda_) * log(q_[a][b] / (px_[a] * py_[b]));

  LOG(DEBUG1) << *this;
}

// template<class Abc>
// void SubstitutionMatrix<Abc>::InitFromScoresAndBackgroundFreqs() {
//   // Calculate target frequencies
//   for (size_t a = 0; a < Abc::kSize; ++a)
//     for (size_t b = 0; b < Abc::kSize; ++b)
//       p_[a][b] = pow(2.0, s_[a][b]) * f_[a] * f_[b];

//   // Check transition probability matrix, renormalize P
//   double sumab = 0.0f;
//   for (size_t a = 0; a < Abc::kSize; a++)
//     for (size_t b = 0; b < Abc::kSize; ++b) sumab += p_[a][b];
//   for (size_t a = 0; a < Abc::kSize; ++a)
//     for (size_t b = 0; b < Abc::kSize; ++b) p_[a][b] /= sumab;

//   // Precompute matrix R for amino acid pseudocounts:
//   for (size_t a = 0; a < Abc::kSize; ++a)
//     for (size_t b = 0; b < Abc::kSize; ++b)
//       r_[a][b] = p_[a][b] / f_[b]; // R[a][b] = P(a|b)

//   LOG(DEBUG1) << *this;
// }

template<class Abc>
void SubstitutionMatrix<Abc>::Print(std::ostream& out) const {
  out << "Background frequencies on query side (in %):\n";
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << Abc::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << strprintf("%-.1f\t", 100.0f * px_[a]);

  out << "\nBackground frequencies on database side (in %):\n";
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << Abc::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << strprintf("%-.1f\t", 100.0f * py_[a]);

  out << "\nSubstitution matrix log( q(a,b) / p(a)*p(b) ) (scaled by 1/lambda):\n";
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << Abc::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < Abc::kSize; ++a) {
    for (size_t b = 0; b < Abc::kSize; ++b)
      out << strprintf("%+5.1f\t", s_[a][b]);
    out << std::endl;
  }

  out << "Target frequency matrix q(a,b) (in %):\n";
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << Abc::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < Abc::kSize; ++a) {
    for (size_t b = 0; b < Abc::kSize; ++b)
      out << strprintf("%-.2f\t", 100.0f * q_[a][b]);
    out << std::endl;
  }

  out << "Matrix of conditional probs on query side P(a|b) = P(a,b)/p(b) (in %):\n";
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << Abc::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < Abc::kSize; ++a) {
    for (size_t b = 0; b < Abc::kSize; ++b)
      out << strprintf("%-.1f\t", 100.0f * rx_[a][b]);
    out << std::endl;
  }

  out << "Matrix of conditional probs on database side P(a|b) = P(a,b)/p(b) (in %):\n";
  for (size_t a = 0; a < Abc::kSize; ++a)
    out << Abc::kIntToChar[a] << "\t";
  out << std::endl;
  for (size_t a = 0; a < Abc::kSize; ++a) {
    for (size_t b = 0; b < Abc::kSize; ++b)
      out << strprintf("%-.1f\t", 100.0f * ry_[a][b]);
    out << std::endl;
  }
}

template<class Abc>
double Neff(const SubstitutionMatrix<Abc>& sm) {
  double neff = 0.0;
  for (size_t b = 0; b < Abc::kSize; ++b) {
    double tmp = 0.0;
    for (size_t a = 0; a < Abc::kSize; ++a)
      if (sm.r(a,b) > FLT_MIN) tmp -= sm.r(a,b) * log2(sm.r(a,b));
    neff += sm.p(b) * pow(2.0, tmp);
  }
  return neff;
}

}  // namespace cs

#endif  // CS_SUBSTITUTION_MATRIX_INL_H_
