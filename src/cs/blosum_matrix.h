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

#ifndef CS_BLOSUM_MATRIX_H_
#define CS_BLOSUM_MATRIX_H_

#include "substitution_matrix-inl.h"

namespace cs {

enum BlosumType {
  BLOSUM45 = 0,
  BLOSUM62 = 1,
  BLOSUM80 = 2
};

// BLOSUM family of substitution matrices for  class for substitution matrix
// classes.
class BlosumMatrix : public SubstitutionMatrix<AA> {
 public:
  BlosumMatrix(BlosumType matrix = BLOSUM62);
  virtual ~BlosumMatrix() {}

 private:
  // Initializes the matrix from target frequencies in raw data array.
  void Init(const float* blosum_xx);
};

// Converts a BLOSUM matrix string to a BLOSUM matrix type.
BlosumType BlosumTypeFromString(const std::string& s);

}  // namespace cs

#endif  // CS_BLOSUM_MATRIX_H_
