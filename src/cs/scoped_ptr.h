// Copyright 2009, Andreas Biegert

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
