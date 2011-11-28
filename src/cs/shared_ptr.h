// Copyright 2009, Andreas Biegert

#ifndef CS_SHARED_PTR_H_
#define CS_SHARED_PTR_H_

#include <stddef.h>

// The shared_ptr class template stores a pointer to a dynamically allocated
// object, typically with a C++ new-expression. The object pointed to is
// guaranteed to be deleted when the last shared_ptr pointing to it is
// destroyed or reset.
template <typename T>
class shared_ptr {
 public:
  shared_ptr() : p_(NULL), c_(new long(0)) {}
  explicit shared_ptr(T* p) : p_(p), c_(new long(1)) {}
  ~shared_ptr() { if(!--*c_) { delete c_; delete p_; } }

  // Copy constructor
  shared_ptr(const shared_ptr& other) : p_(other.p_), c_(other.c_) { ++*c_; }

  // Generalized copy constructor
  template <typename U>
  shared_ptr(const shared_ptr<U>& other)
      : p_(other.get()), c_(other.use_count()) { ++*c_; }

  // Copy assignment
  shared_ptr& operator= (const shared_ptr &rhs) {
    if (this != &rhs) {
      if (!--*c_) { delete c_; delete p_; }
      p_ = rhs.p_;
      ++*(c_ = rhs.c_);
    }
    return *this;
  }

  // Generalized copy assignment
  template <typename U>
  shared_ptr& operator= (const shared_ptr<U> &rhs) {
    if (this != &rhs) {
      if (!--*c_) { delete c_; delete p_; }
      p_ = rhs.get();
      ++*(c_ = rhs.use_count());
    }
    return *this;
  }

  T& operator*() const { return *p_; }
  T* operator->() const { return p_; }
  T* get() const { return p_; }
  long use_count() const { return *c_; }
  bool unique() const { return *c_ == 1; }
  operator bool () const { return p_ != NULL; }

  // Template function for implicit conversion (see More Effectove C++
  // page 176 for details).
  // template <class U>
  // operator shared_ptr<U>() { return shared_ptr<U>(p_); }

 private:
  T* p_;
  long* c_;
};

#endif  // CS_SHARED_PTR_H_
