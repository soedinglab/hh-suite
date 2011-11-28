// Copyright 2009, Andreas Biegert

#ifndef CS_SUBSTITUTION_MATRIX_H_
#define CS_SUBSTITUTION_MATRIX_H_

namespace cs {

// Abstract base class for substitution matrix classes.
template<class Abc>
class SubstitutionMatrix {
 public:
  SubstitutionMatrix(double l = 1.0);
  virtual ~SubstitutionMatrix() = 0;

  // Return substitution score s(a,b).
  double s(size_t a, size_t b) const { return s_[a][b]; }

  // Returns target frequency q(a,b).
  double q(size_t a, size_t b) const { return q_[a][b]; }

  // Returns conditional probability P(a|b) "a given b" for symmetric matrices.
  double r(size_t a, size_t b) const { return rx_[a][b]; }

  // Returns conditional probability P(a|b) "a given b" on query side
  double rx(size_t a, size_t b) const { return rx_[a][b]; }

  // Returns conditional probability P(a|b) "a given b" on database side
  double ry(size_t a, size_t b) const { return ry_[a][b]; }

  // Returns equilibrium frequency of alphbet letter a for symmetric matrices.
  double p(size_t a) const { return px_[a]; }

  // Returns equilibrium frequency of alphbet letter a on query side.
  double px(size_t a) const { return px_[a]; }

  // Returns equilibrium frequency of alphbet letter a on database side.
  double py(size_t a) const { return py_[a]; }

  // Prints the substitution matrix in human readable format to stream.
  friend std::ostream& operator<< (std::ostream& out,
                                   const SubstitutionMatrix& m) {
    m.Print(out);
    return out;
  }

 protected:
  // Initializes the other matrix data members from target frequencies.
  void InitFromTargetFreqs();

  // Target frequency matrix P(a,b).
  Matrix<double> q_;
  // Substitution matrix S(a,b).
  Matrix<double> s_;
  // Conditional probability matrix on query side.
  Matrix<double> rx_;
  // Conditional probability matrix on database side.
  Matrix<double> ry_;
  // Background frequencies of alphabet on query side.
  Vector<double> px_;
  // Background frequencies of alphabet on database side.
  Vector<double> py_;
  // Scaling factor lambda
  double lambda_;

 private:
  // Prints the substitution matrix in human-readable format to output stream.
  virtual void Print(std::ostream& out) const;
};

// Returns the average number of effective sequences for given substitution matrix
template<class Abc>
double Neff(const SubstitutionMatrix<Abc>& sm);

}  // namespace cs

#endif  // CS_SUBSTITUTION_MATRIX_H_
