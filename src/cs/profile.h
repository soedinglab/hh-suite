// Copyright 2009, Andreas Biegert

#ifndef CS_PROFILE_H_
#define CS_PROFILE_H_

namespace cs {

// A simple container class acting pretty much like a matrix but with column size
// fixed to Abc::kSizeAny
template <class Abc>
class Profile {
  public:
    Profile();
    explicit Profile(size_t n);
    Profile(size_t n, const double &a);
    Profile(size_t n, const double *a);
    Profile(const Profile &rhs);
    ~Profile();

    Profile & operator=(const Profile &rhs);
    double* operator[](const size_t i);
    const double* operator[](const size_t i) const;
    size_t length() const { return nn; }
    void Resize(size_t newn);
    void Assign(size_t newn, const double &a);

  private:
    size_t nn;
    double **v;
};

// Assigns given constant value or default to all entries in matrix
template <class Abc>
void Assign(Profile<Abc>& p, double val = 0.0);

// Normalizes all profile columns to fixed value. Iff 'incl_any' is true,
// normalization also includes values for ANY letter
template <class Abc>
void Normalize(Profile<Abc>& p, double val, bool incl_any = false);

// Normalizes all profile columns to corresponding value in vector 'norm'.
// Iff 'incl_any' is true, normalization also includes values for ANY letter.
template <class Abc>
void Normalize(Profile<Abc>& p, const Vector<double>& norm, bool incl_any = false);

// Prints profile in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const Profile<Abc>& p);

}  // namespace cs

#endif  // CS_PROFILE_H_
