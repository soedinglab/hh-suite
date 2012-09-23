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

#ifndef CS_VECTOR_H_
#define CS_VECTOR_H_

template <class T>
class Vector {
 public:
  typedef T value_type;

  Vector();
  explicit Vector(size_t n);
  Vector(size_t n, const T &a);
  Vector(size_t n, const T *a);
  Vector(const Vector &rhs);
  ~Vector();

  Vector & operator=(const Vector &rhs);

  inline T & operator[](const size_t i);
  inline const T & operator[](const size_t i) const;
  inline size_t size() const;
  void Resize(size_t newn);
  void Assign(size_t newn, const T &a);
  inline void Assign(const T &a);

 private:
  size_t nn;
  T *v;
};

// Vector definitions
template <class T>
Vector<T>::Vector() : nn(0), v(NULL) {}

template <class T>
Vector<T>::Vector(size_t n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
Vector<T>::Vector(size_t n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL) {
  for(size_t i=0; i<n; i++) v[i] = a;
}

template <class T>
Vector<T>::Vector(size_t n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL) {
  for(size_t i=0; i<n; i++) v[i] = *a++;
}

template <class T>
Vector<T>::Vector(const Vector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL) {
  for(size_t i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
Vector<T> & Vector<T>::operator=(const Vector<T> &rhs) {
  if (this != &rhs) {
    if (nn != rhs.nn) {
      if (v != NULL) delete [] (v);
      nn=rhs.nn;
      v= nn>0 ? new T[nn] : NULL;
    }
    for (size_t i=0; i<nn; i++)
      v[i]=rhs[i];
  }
  return *this;
}

template <class T>
inline T & Vector<T>::operator[](const size_t i) {
#ifdef CHECKBOUNDS
  if (i<0 || i>=nn) {
    throw Exception("Vector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline const T & Vector<T>::operator[](const size_t i) const {
#ifdef CHECKBOUNDS
  if (i<0 || i>=nn) {
    throw Exception("Vector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline size_t Vector<T>::size() const {
  return nn;
}

template <class T>
void Vector<T>::Resize(size_t newn) {
  if (newn != nn) {
    if (v != NULL) delete[] (v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
}

template <class T>
void Vector<T>::Assign(size_t newn, const T& a) {
  if (newn != nn) {
    if (v != NULL) delete[] (v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
  for (size_t i=0;i<nn;i++) v[i] = a;
}

template <class T>
inline void Vector<T>::Assign(const T& a) {
  Assign(nn, a);
}

template <class T>
Vector<T>::~Vector() {
  if (v != NULL) delete[] (v);
}

// Assigns given constant value or default to all entries in vector
template<class T>
inline void Assign(Vector<T>& v, T val = T()) {
  for (size_t i = 0; i < v.size(); ++i) v[i] = val;
}

#endif  // CS_VECTOR_H_
