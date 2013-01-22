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

#ifndef CS_CONTEXT_PROFILE_H_
#define CS_CONTEXT_PROFILE_H_

#include "profile-inl.h"
#include "count_profile-inl.h"
#include "profile_column.h"

namespace cs {

// Simple struct for working with RGB color values.
struct Color {
    Color(double r = 1.0, double g = 1.0, double b = 1.0)
	    : red(r), green(g), blue(b) {}

    Color(std::string coldef) {
        std::vector<std::string> tokens;
        Tokenize(coldef, ',', &tokens);
        red   = atof(tokens[0].c_str());
        green = atof(tokens[1].c_str());
        blue  = atof(tokens[2].c_str());
    }

    std::string ToString() const {
        return strprintf("%4.2f,%4.2f,%4.2f", red, green, blue);
    }

    bool operator< (const Color& rhs) const {
        if (red != rhs.red)     return (red < rhs.red);
        if (green != rhs.green) return (green < rhs.green);
        return (blue < rhs.blue);
    }

    friend std::ostream& operator<< (std::ostream& out, const Color& c) {
        out << c.ToString();
        return out;
    }

    double red, green, blue;
};

template<class Abc>
struct ContextProfile {
    // Default construction
    ContextProfile();

    // Constructs a context profile with 'len' columns
    explicit ContextProfile(size_t len);

    // Construction from serialized profile read from input stream.
    explicit ContextProfile(FILE* fin);

    // Constructs a context profile from normalized values in given profile 'p'.
    explicit ContextProfile(const Profile<Abc>& p);

    // Construct a context profile by copying subprofile starting at index 'start'
    // for 'len' columns and normalizing afterwards.
    ContextProfile(const Profile<Abc>& p, size_t start, size_t len);

    // Compares context profiles
    bool operator < (const ContextProfile<Abc>& other) const {
      return prior < other.prior;
    }

    // Initializes count profile with a serialized profile read from stream.
    void Read(FILE* fin);

    // Initializes count profile with a serialized profile read from stream.
    void Write(FILE* fin) const;

    // Returns number of columns.
    size_t length() const;

    std::string name;               // optional name of this context profile
    double prior;                   // prior probability
    bool is_log;                    // flag indicating if profile is in log-space
    Profile<Abc> probs;             // profile probabilities (includes ANY)
    ProfileColumn<Abc> pc;          // predicted pseudocounts (always in lin-space)
    Color color;                    // RGB color for profile logos
};


// Prints context profile probabilities in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const ContextProfile<Abc>& cp);

// Transforms probabilites in context profile to log-space and sets 'is_log' flag.
template<class Abc>
void TransformToLog(ContextProfile<Abc>& cp);

// Transforms probabilites in context profile to lin-space and sets 'is_log' flag.
template<class Abc>
void TransformToLin(ContextProfile<Abc>& cp);

}  // namespace cs

#endif  // CS_CONTEXT_PROFILE_H_
