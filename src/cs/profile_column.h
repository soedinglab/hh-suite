// Copyright 2009, Andreas Biegert

#ifndef CS_PROFILE_COLUMN_H_
#define CS_PROFILE_COLUMN_H_

namespace cs {

template<class Abc>
class ProfileColumn {
 public:
  ProfileColumn();
  explicit ProfileColumn(const double &a);
  explicit ProfileColumn(const double *a);
  ProfileColumn(const ProfileColumn &rhs);
  ~ProfileColumn();

  ProfileColumn & operator=(const ProfileColumn &rhs);

  inline double & operator[](const size_t i);
  inline const double & operator[](const size_t i) const;
  inline size_t size() const;

 private:
  double *v;
};

template<class Abc>
ProfileColumn<Abc>::ProfileColumn() : v(new double[Abc::kSizeAny]) {}

template<class Abc>
ProfileColumn<Abc>::ProfileColumn(const double& a) : v(new double[Abc::kSizeAny]) {
  for(size_t i=0; i<Abc::kSizeAny; i++) v[i] = a;
}

template<class Abc>
ProfileColumn<Abc>::ProfileColumn(const double *a) : v(new double[Abc::kSizeAny]) {
  for(size_t i=0; i<Abc::kSizeAny; i++) v[i] = *a++;
}

template<class Abc>
ProfileColumn<Abc>::ProfileColumn(const ProfileColumn<Abc> &rhs)
  : v(new double[Abc::kSizeAny]) {
  for(size_t i=0; i<Abc::kSizeAny; i++) v[i] = rhs[i];
}

template<class Abc>
ProfileColumn<Abc> & ProfileColumn<Abc>::operator=(const ProfileColumn<Abc> &rhs) {
  if (this != &rhs)
    for (size_t i=0; i<Abc::kSizeAny; i++) v[i]=rhs[i];
  return *this;
}

template<class Abc>
inline double & ProfileColumn<Abc>::operator[](const size_t i) {
#ifdef CHECKBOUNDS
  if (i<0 || i>=nn) {
    throw Exception("Profile column subscript out of bounds");
  }
#endif
  return v[i];
}

template<class Abc>
inline const double & ProfileColumn<Abc>::operator[](const size_t i) const {
#ifdef CHECKBOUNDS
  if (i<0 || i>=nn) {
    throw Exception("Profile column subscript out of bounds");
  }
#endif
  return v[i];
}

template<class Abc>
inline size_t ProfileColumn<Abc>::size() const {
  return Abc::kSizeAny;
}

template<class Abc>
ProfileColumn<Abc>::~ProfileColumn() {
  if (v != NULL) delete[] (v);
}

// Assigns given constant value or default to all entries in vector
template<class Abc>
inline void Assign(ProfileColumn<Abc>& v, double val) {
  for (size_t i = 0; i < Abc::kSizeAny; ++i) v[i] = val;
}

// Prints profile column in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const ProfileColumn<Abc>& col) {
  for (size_t a = 0; a < Abc::kSizeAny; ++a)
    out << strprintf("   %c  \t", Abc::kIntToChar[a]);
  out << std::endl;
  for (size_t a = 0; a < Abc::kSizeAny; ++a)
    out << strprintf("%6.4f\t", col[a]);
  out << std::endl;
  return out;
}

// Normalizes all profile column to fixed value. Iff 'incl_any' is true,
// normalization also includes value of ANY letter
template <class Abc>
inline void Normalize(ProfileColumn<Abc>& col, double val, bool incl_any = false) {
  const size_t abc_size = incl_any ? Abc::kSizeAny : Abc::kSize;
  double sum = 0;
  for (size_t a = 0; a < abc_size; ++a) sum += col[a];
  double fac = val / sum;
  for (size_t a = 0; a < abc_size; ++a) col[a] *= fac;
}

// Calculates entropy of given profile column using logarithm base 2.
template <class Abc>
inline double Entropy(const ProfileColumn<Abc>& col) {
  double rv = 0.0;
  for (size_t a = 0; a < Abc::kSize; ++a)
    if (col[a] > FLT_MIN) rv -= col[a] * log2(col[a]);
  return rv;
}

}  // namespace cs

#endif  // CS_PROFILE_COLUMN_H_
