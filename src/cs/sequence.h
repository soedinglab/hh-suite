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

#include <stdint.h>
#include <sstream>

#ifndef CS_SEQUENCE_H_
#define CS_SEQUENCE_H_

namespace cs {

// A container class representing a sequence consisting of letters over a
// sequence alphabet.
template<class Abc>
class Sequence {
  public:
    typedef uint8_t value_type;
    typedef value_type* Iter;
    typedef const value_type* ConstIter;

    // Constructs sequence with specified length.
    Sequence(size_t length = 0);

    // Copy constructor
    Sequence(const Sequence& other);

    // Subsequence constructor
    Sequence(const Sequence& other, size_t idx, size_t len);

    // Constructs sequence from serialized sequence in FASTA format
    explicit Sequence(FILE* fin);

    // Constructs sequence with given header and sequence string of characters.
    Sequence(const std::string& sequence, const std::string& header = "");

    // Deallocates sequence array
    ~Sequence() { delete[] seq_; }

    // Assignment operator
    Sequence& operator= (const Sequence& rhs);

    // Accessors for integer at position i of the sequence.
    value_type& operator[](size_t i) { return seq_[i]; }
    const value_type& operator[](size_t i) const { return seq_[i]; }
    value_type& at(size_t i) { return seq_[i]; }
    const value_type& at(size_t i) const { return seq_[i]; }

    // Returns the character at position i of the sequence.
    int chr(size_t i) const { return Abc::kIntToChar[seq_[i]]; }

    // Returns the sequence length.
    size_t length() const { return length_; }

    // Sets the header to given string.
    void set_header(const std::string& header) { header_ = header; }

    // Returns the header information of this sequence.
    const std::string& header() const { return header_; }

    // Returns a const iterator to the first integer element of the sequence.
    ConstIter begin() const { return seq_; }

    // Returns a const iterator just past the end of the sequence.
    ConstIter end() const { return begin() + length(); }

    // Returns an iterator to the first integer element of the sequence.
    Iter begin() { return seq_; }

    // Returns an iterator just past the end of the sequence.
    Iter end() { return begin() + length(); }

    // Initializes the sequence object with a sequence in FASTA format
    void Read(FILE* in);

    // Prints the sequence in FASTA format to output stream.
    void Write(FILE* fout, size_t width = 100) const;

    void Write(std::stringstream& ss, size_t width = 100) const;


    // Resizes current sequence to new length (Note: old sequence is NOT preserved!)
    void Resize(size_t newlen);

    // Returns sequence as character string.
    std::string ToString() const;

  protected:
    // Convert the sequence in character representation to integer representation.
    void Init(std::string sequence, std::string header);

    // Length of sequence
    size_t length_;
    // The sequence itself in integer representation
    value_type* seq_;
    // The header without leading '>'
    std::string header_;
};  // Sequence


// Prints the Alignment in A2M format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const Sequence<Abc>& s) {
    const size_t kWidth = 100;
    out << '>' << s.header() << std::endl;
    for (size_t i = 0; i < s.length(); ++i) {
        out << s.chr(i);
        if ((i+1) % kWidth == 0) out << std::endl;
    }
    if (s.length() % kWidth != 0) out << std::endl;
    return out;
}

// Predicate indicating if character is a gap character, that is '-' or '.'
inline bool isgap(char c) {
    return (c == '-' || c == '.');
}

}  // namespace cs

#endif  // CS_SEQUENCE_H_
