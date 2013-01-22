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

#ifndef CS_MATRIX_H_
#define CS_MATRIX_H_

#include "stride_iter.h"

#include <valarray>
#include <numeric>
#include <algorithm>

using std::valarray;

template<class T>
class matrix {
 public:
  // public typedefs
  typedef T value_type;
  typedef matrix self;
  typedef value_type* iterator;
  typedef const value_type* const_iterator;
  typedef value_type* row_type;
  typedef stride_iter<value_type*> col_type;
  typedef const value_type* const_row_type;
  typedef stride_iter<const value_type*> const_col_type;

  // constructors
  matrix() : num_rows_(0), num_cols_(0), m_() {}
  matrix(int r, int c) : num_rows_(r), num_cols_(c), m_(r * c) {}
  matrix(const self& x)
      : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x.m_) {}

  matrix(int r, int c, const T& val)
      : num_rows_(r), num_cols_(c), m_(r * c) {
    for (int i = 0; i < r * c; ++i) m_[i] = val;
  }

  template<typename New_T>
  explicit matrix(const valarray<New_T>& x)
      : num_rows_(x.size()), num_cols_(1), m_(x.size()) {
    for (int i =0 ; i < x.size(); ++i)
      m_[i] = x[i];
  }

  // allow construction from matricies of other types
  template<typename New_T>
  explicit matrix(const matrix<New_T>& x)
      : num_rows_(x.num_rows_), num_cols_(x.num_cols_), m_(x.size()) {
    copy(x.begin(), x.end(), m_.begin());
  }

  // public functions
  int num_rows() const { return num_rows_; }
  int num_cols() const { return num_cols_; }
  int size() const { return num_rows_ * num_cols_; }

  // element access
  row_type row_begin(int n) { return &m_[n * num_cols()]; }
  row_type row_end(int n) { return row_begin(n) + num_cols(); }
  col_type col_begin(int n) { return col_type(&m_[n], num_cols()); }
  col_type col_end(int n) { return col_begin(n) + num_cols(); }
  const_row_type row_begin(int n) const { return &m_[n * num_cols()]; }
  const_row_type row_end(int n) const { return row_begin(n) + num_cols(); }
  const_col_type col_begin(int n) const {
    return const_col_type(&m_[n], num_cols());
  }
  const_col_type col_end(int n) const { return col_begin(n) + num_cols(); }
  iterator begin() { return &m_[0]; }
  iterator end() { return begin() + size(); }
  const_iterator begin() const { return &m_[0]; }
  const_iterator end() const { return begin() + size(); }
  void resize(int r, int c, T v = T()) {
    m_.resize(r * c, v);
    num_rows_ = r;
    num_cols_ = c;
  }

  // operators
  self& operator=(const self& x) {
    m_.resize(x.size());
    m_ = x.m_;
    num_rows_ = x.num_rows_;
    num_cols_ = x.num_cols_;
    return *this;
  }
  self& operator=(value_type x) { m_ = x; return *this; }
  row_type operator[](int n) { return row_begin(n); }
  const_row_type operator[](int n) const { return row_begin(n); }
  self& operator+=(const self& x) { m_ += x.m_; return *this; }
  self& operator-=(const self& x) { m_ -= x.m_; return *this; }
  self& operator+=(value_type x) { m_ += x; return *this; }
  self& operator-=(value_type x) { m_ -= x; return *this; }
  self& operator*=(value_type x) { m_ *= x; return *this; }
  self& operator/=(value_type x) { m_ /= x; return *this; }
  self& operator%=(value_type x) { m_ %= x; return *this; }
  self operator-( ) { return -m_; }
  self operator+( ) { return +m_; }
  self operator!( ) { return !m_; }
  self operator~( ) { return ~m_; }

  // friend operators
  friend self operator+(const self& x, const self& y) { return self(x) += y; }
  friend self operator-(const self& x, const self& y) { return self(x) -= y; }
  friend self operator+(const self& x, value_type y) { return self(x) += y; }
  friend self operator-(const self& x, value_type y) { return self(x) -= y; }
  friend self operator*(const self& x, value_type y) { return self(x) *= y; }
  friend self operator/(const self& x, value_type y) { return self(x) /= y; }
  friend self operator%(const self& x, value_type y) { return self(x) %= y; }

 private:
  int num_rows_;
  int num_cols_;
  mutable valarray<T> m_;
};

// Resets all entries in matrix to the provided value or default if  none given.
template<class T>
inline void Reset(matrix<T>* m, T value = T()) {
  const int num_rows = m->num_rows();
  const int num_cols = m->num_cols();

  for (int i = 0; i < num_rows; ++i)
    for (int j = 0; j < num_cols; ++j)
      (*m)[i][j] = value;
}


template <class T>
class Matrix {
 public:
  typedef T value_type;

  Matrix();
  Matrix(size_t n, size_t m);
  Matrix(size_t n, size_t m, const T &a);
  Matrix(size_t n, size_t m, const T *a);
  Matrix(const Matrix &rhs);
  ~Matrix();

  Matrix & operator=(const Matrix &rhs);
  T* operator[](const size_t i);
  const T* operator[](const size_t i) const;
  size_t nrows() const;
  size_t ncols() const;
  void Resize(size_t newn, size_t newm);
  void Assign(size_t newn, size_t newm, const T &a);
  T* begin() { return *v; }
  const T* begin() const { return *v; }

 private:
  size_t nn;
  size_t mm;
  T **v;
};

template <class T>
Matrix<T>::Matrix() : nn(0), mm(0), v(NULL) {}

template <class T>
Matrix<T>::Matrix(size_t n, size_t m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL) {
  size_t i,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
Matrix<T>::Matrix(size_t n, size_t m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL) {
  size_t i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
Matrix<T>::Matrix(size_t n, size_t m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL) {
  size_t i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
Matrix<T>::Matrix(const Matrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL) {
  size_t i,j,nel=mm*nn;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs) {
  if (this != &rhs) {
    size_t i,j,nel;
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
        delete[] (v[0]);
        delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v = nn>0 ? new T*[nn] : NULL;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : NULL;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}

template <class T>
inline T* Matrix<T>::operator[](const size_t i) {
#ifdef CHECKBOUNDS
  if (i<0 || i>=nn) {
    throw Exception("Matrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline const T* Matrix<T>::operator[](const size_t i) const {
#ifdef CHECKBOUNDS
  if (i<0 || i>=nn) {
    throw Excpetion("Matrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline size_t Matrix<T>::nrows() const {
  return nn;
}

template <class T>
inline size_t Matrix<T>::ncols() const {
  return mm;
}

template <class T>
void Matrix<T>::Resize(size_t newn, size_t newm) {
  size_t i,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
}

template <class T>
void Matrix<T>::Assign(size_t newn, size_t newm, const T& a) {
  size_t i,j,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
Matrix<T>::~Matrix() {
  if (v != NULL) {
    delete[] (v[0]);
    delete[] (v);
  }
}

// Assigns given constant value or default to all entries in matrix
template<class T>
inline void Assign(Matrix<T>& m, T val = T()) {
  for (size_t i = 0; i < m.nrows(); ++i)
    for (size_t j = 0; j < m.ncols(); ++j)
      m[i][j] = val;
}

// Prints profile in human-readable format for debugging.
template<class T>
std::ostream& operator<< (std::ostream& out, const Matrix<T>& m) {
    out << "Matrix" << std::endl;
    for (size_t j = 0; j < m.ncols(); ++j)
        out << "\t" << j;
    out << std::endl;
    for (size_t i = 0; i < m.nrows(); ++i) {
        out << i;
        for (size_t j = 0; j < m.ncols(); ++j)
            out << "\t" << std::setprecision(4) << m[i][j];
        out << std::endl;
    }
    return out;
}

#endif  // CS_MATRIX_H_
