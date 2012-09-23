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

#ifndef CS_ABSTRACT_STATE_MATRIX_H_
#define CS_ABSTRACT_STATE_MATRIX_H_

namespace cs {

// A substitution matrix between contexts on the query side and abstract states
// on the database side.
template<class AS>
class AbstractStateMatrix {
 public:
  AbstractStateMatrix(std::string matrixfile);

  // Return substitution score s(a,b) of context k and abstract state s.
  double s(size_t k, size_t s) const { return s_[k][s]; }

  // Returns target frequency q(k,s) of context k and abstract state s.
  double q(size_t k, size_t s) const { return q_[k][s]; }

  // Returns conditional probability P(s|k) "s given k".
  double r(size_t s, size_t k) const { return r_[s][k]; }

  // Returns equilibrium frequency of context k on query side.
  double px(size_t k) const { return px_[k]; }

  // Returns equilibrium frequency of abstract state s on database side.
  double py(size_t s) const { return py_[s]; }

  // Returns number of context profiles on query side.
  size_t num_contexts() const { return num_contexts_; }

  // Prints the substitution matrix in human readable format to stream.
  friend std::ostream& operator<< (std::ostream& out,
                                   const AbstractStateMatrix& m) {
    m.Print(out);
    return out;
  }

 protected:
  // Initializes the other matrix data members from target frequencies.
  void Init();

  // Number of contexts on query side
  size_t num_contexts_;
  // Target frequency matrix P(k,s).
  Matrix<double> q_;
  // Substitution matrix S(k,s).
  Matrix<double> s_;
  // Conditional probability matrix "s given k".
  Matrix<double> r_;
  // Background frequencies of contexts on query side.
  Vector<double> px_;
  // Background frequencies of abstract states on database side.
  Vector<double> py_;

 private:
  // Prints the substitution matrix in human-readable format to output stream.
  virtual void Print(std::ostream& out) const;
};

}  // namespace cs

#endif  // CS_ABSTRACT_STATE_MATRIX_H_
