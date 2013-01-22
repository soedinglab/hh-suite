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

#ifndef CS_SCOPED_PTR_H_
#define CS_SCOPED_PTR_H_

#include "globals.h"

// This implementation of scoped_ptr is PARTIAL - it only contains
// enough stuff to satisfy our needs.
template <typename T>
class scoped_ptr {
 public:
  explicit scoped_ptr(T* p = NULL) : ptr_(p) {}
  ~scoped_ptr() { reset(); }

  T& operator*() const { return *ptr_; }
  T* operator->() const { return ptr_; }
  T* get() const { return ptr_; }
  operator bool () const { return ptr_ != NULL; }

  T* release() {
    T* const ptr = ptr_;
    ptr_ = NULL;
    return ptr;
  }

  void reset(T* p = NULL) {
    if (p != ptr_) {
      if (sizeof(T) > 0) {  // Makes sure T is a complete type.
        delete ptr_;
      }
      ptr_ = p;
    }
  }

 private:
  T* ptr_;

  DISALLOW_COPY_AND_ASSIGN(scoped_ptr);
};

#endif  // CS_SCOPED_PTR_H_
