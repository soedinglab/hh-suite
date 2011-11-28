// Copyright 2009, Andreas Biegert

#ifndef CS_PROFILE_INL_H_
#define CS_PROFILE_INL_H_

#include "profile.h"

namespace cs {

template <class Abc>
Profile<Abc>::Profile() : nn(0), v(NULL) {}

template <class Abc>
Profile<Abc>::Profile(size_t n) : nn(n), v(n>0 ? new double*[n] : NULL) {
    size_t i,nel=n*Abc::kSizeAny;
    if (v) v[0] = nel>0 ? new double[nel] : NULL;
    for (i=1;i<n;i++) v[i] = v[i-1] + Abc::kSizeAny;
}

template <class Abc>
Profile<Abc>::Profile(size_t n, const double &a)
        : nn(n), v(n>0 ? new double*[n] : NULL) {
    size_t i,j,nel=n*Abc::kSizeAny;
    if (v) v[0] = nel>0 ? new double[nel] : NULL;
    for (i=1; i< n; i++) v[i] = v[i-1] + Abc::kSizeAny;
    for (i=0; i< n; i++) for (j=0; j<Abc::kSizeAny; j++) v[i][j] = a;
}

template <class Abc>
Profile<Abc>::Profile(size_t n, const double *a)
        : nn(n), v(n>0 ? new double*[n] : NULL) {
    size_t i,j,nel=n*Abc::kSizeAny;
    if (v) v[0] = nel>0 ? new double[nel] : NULL;
    for (i=1; i< n; i++) v[i] = v[i-1] + Abc::kSizeAny;
    for (i=0; i< n; i++) for (j=0; j<Abc::kSizeAny; j++) v[i][j] = *a++;
}

template <class Abc>
Profile<Abc>::Profile(const Profile &rhs)
        : nn(rhs.nn), v(nn>0 ? new double*[nn] : NULL) {
    size_t i,j,nel=nn*Abc::kSizeAny;
    if (v) v[0] = nel>0 ? new double[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + Abc::kSizeAny;
    for (i=0; i< nn; i++) for (j=0; j<Abc::kSizeAny; j++) v[i][j] = rhs[i][j];
}

template <class Abc>
Profile<Abc> & Profile<Abc>::operator=(const Profile<Abc> &rhs) {
    if (this != &rhs) {
        size_t i,j,nel;
        if (nn != rhs.nn) {
            if (v != NULL) {
                delete[] v[0];
                delete[] v;
            }
            nn=rhs.nn;
            v = nn>0 ? new double*[nn] : NULL;
            nel = nn*Abc::kSizeAny;
            if (v) v[0] = nel>0 ? new double[nel] : NULL;
            for (i=1; i< nn; i++) v[i] = v[i-1] + Abc::kSizeAny;
        }
        for (i=0; i< nn; i++) for (j=0; j<Abc::kSizeAny; j++) v[i][j] = rhs[i][j];
    }
    return *this;
}

template <class Abc>
inline double* Profile<Abc>::operator[](const size_t i) {
#ifdef CHECKBOUNDS
    if (i<0 || i>=nn) {
        throw Exception("Profile subscript out of bounds");
    }
#endif
    return v[i];
}

template <class Abc>
inline const double* Profile<Abc>::operator[](const size_t i) const {
#ifdef CHECKBOUNDS
    if (i<0 || i>=nn) {
        throw Exception("Profile subscript out of bounds");
    }
#endif
    return v[i];
}

template <class Abc>
void Profile<Abc>::Resize(size_t newn) {
    size_t i,nel;
    if (newn != nn) {
        if (v != NULL) {
            delete[] v[0];
            delete[] v;
        }
        nn = newn;
        v = nn>0 ? new double*[nn] : NULL;
        nel = nn*Abc::kSizeAny;
        if (v) v[0] = nel>0 ? new double[nel] : NULL;
        for (i=1; i< nn; i++) v[i] = v[i-1] + Abc::kSizeAny;
    }
}

template <class Abc>
void Profile<Abc>::Assign(size_t newn, const double& a) {
    size_t i,j,nel;
    if (newn != nn) {
        if (v != NULL) {
            delete[] v[0];
            delete[] v;
        }
        nn = newn;
        v = nn>0 ? new double*[nn] : NULL;
        nel = nn*Abc::kSizeAny;
        if (v) v[0] = nel>0 ? new double[nel] : NULL;
        for (i=1; i< nn; i++) v[i] = v[i-1] + *Abc::kSizeAny;
    }
    for (i=0; i< nn; i++) for (j=0; j<*Abc::kSizeAny; j++) v[i][j] = a;
}

template <class Abc>
Profile<Abc>::~Profile() {
    if (v != NULL) {
        delete[] v[0];
        delete[] v;
    }
}

// Assigns given constant value or default to all entries in matrix
template <class Abc>
inline void Assign(Profile<Abc>& p, double val) {
    for (size_t i = 0; i < p.length(); ++i)
        for (size_t j = 0; j < Abc::kSizeAny; ++j)
            p[i][j] = val;
}

// Normalizes all profile columns to fixed value. Iff 'incl_any' is true,
// normalization also includes values for ANY letter
template <class Abc>
inline void Normalize(Profile<Abc>& p, double val, bool incl_any) {
    const size_t abc_size = incl_any ? Abc::kSizeAny : Abc::kSize;
    for (size_t i = 0; i < p.length(); ++i) {
        double sum = 0.0;
        for (size_t a = 0; a < abc_size; ++a) sum += p[i][a];
        if (sum != 0.0) {
            double fac = val / sum;
            for (size_t a = 0; a < abc_size; ++a) p[i][a] *= fac;
        }
    }
}

// Normalizes all profile columns to corresponding value in vector 'norm'.
// Iff 'incl_any' is true, normalization also includes values for ANY letter.
template <class Abc>
inline void Normalize(Profile<Abc>& p, const Vector<double>& norm, bool incl_any) {
    const size_t abc_size = incl_any ? Abc::kSizeAny : Abc::kSize;
    for (size_t i = 0; i < p.length(); ++i) {
        double sum = 0.0;
        for (size_t a = 0; a < abc_size; ++a) sum += p[i][a];
        if (sum != 0.0) {
            double fac = norm[i] / sum;
            for (size_t a = 0; a < abc_size; ++a) p[i][a] *= fac;
        }
    }
}

// Prints profile in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const Profile<Abc>& p) {
    out << "Profile" << std::endl;
    for (size_t a = 0; a < Abc::kSizeAny; ++a)
        out << "\t" << Abc::kIntToChar[a];
    out << std::endl;
    for (size_t i = 0; i < p.length(); ++i) {
        out << i+1;
        for (size_t a = 0; a < Abc::kSizeAny; ++a)
            out << strprintf("\t%6.4f", p[i][a]);
        out << std::endl;
    }
    return out;
}

}  // namespace cs

#endif  // CS_PROFILE_INL_H_
