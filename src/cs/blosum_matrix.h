// Copyright 2009, Andreas Biegert

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
