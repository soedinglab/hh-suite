/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

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

#ifndef CS_PROFILE_H_
#define CS_PROFILE_H_


namespace cs {

// Forward declaration
template <class Abc>
struct CountProfile;

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
    Profile(const CountProfile<Abc>& cp);
    ~Profile();

    Profile & operator=(const Profile &rhs);
    double* operator[](const size_t i);
    const double* operator[](const size_t i) const;
    size_t length() const { return nn; }
    void Resize(size_t newn);
    void Assign(size_t newn, const double &a);
    void Insert(size_t idx, const Profile<Abc>& other);

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

// Calculates the entropy of the given profile using logarithm base 2.
template <class Abc>
inline double Entropy(const Profile<Abc>& p);

// Calculates the Neff in the given profile.
template <class Abc>
inline double Neff(const Profile<Abc>& p);

}  // namespace cs

#endif  // CS_PROFILE_H_
