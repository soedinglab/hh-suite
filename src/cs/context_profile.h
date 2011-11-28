// Copyright 2009, Andreas Biegert

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
    ContextProfile() : prior(0.0), probs(), pc() {};

    // Constructs a context profile with 'len' columns
    explicit ContextProfile(size_t len)
	    : prior(0.0), is_log(false), probs(len), pc() {
	assert(len & 1);
    }

    // Construction from serialized profile read from input stream.
    explicit ContextProfile(FILE* fin) { Read(fin); }

    // Constructs a context profile from normalized values in given profile 'p'.
    explicit ContextProfile(const Profile<Abc>& p)
	    : prior(0.0),
	      is_log(false),
	      probs(p),
	      pc(&p[(p.length() - 1) / 2][0]) {
	assert(p.length() & 1);
	for (size_t i = 0; i < probs.length(); ++i)
	    probs[i][Abc::kAny] = 1.0;
	pc[Abc::kAny] = 1.0;

	Normalize(probs, 1.0);
	Normalize(pc, 1.0);
    }

    // Construct a context profile by copying subprofile starting at index 'start'
    // for 'len' columns and normalizing afterwards.
    ContextProfile(const Profile<Abc>& p, size_t start, size_t len)
	    : prior(0.0),
	      is_log(false),
	      probs(len),
	      pc(&p[start + (len - 1) / 2][0]) {
	assert(len & 1);
	for (size_t i = 0; i < len; ++i) {
	    for (size_t a = 0; a < Abc::kSize; ++a)
		probs[i][a] = p[start + i][a];
	    probs[i][Abc::kAny] = 1.0;
	}
	pc[Abc::kAny] = 1.0;

	Normalize(probs, 1.0);
	Normalize(pc, 1.0);
    }

    // Initializes count profile with a serialized profile read from stream.
    void Read(FILE* fin);

    // Initializes count profile with a serialized profile read from stream.
    void Write(FILE* fin) const;

    // Returns number of columns.
    size_t length() const { return probs.length(); }

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
