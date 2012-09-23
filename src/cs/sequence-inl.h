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

#ifndef CS_SEQUENCE_INL_H_
#define CS_SEQUENCE_INL_H_

#include "sequence.h"

namespace cs {

template<class Abc>
Sequence<Abc>::Sequence(size_t length)
        : length_(length),
          seq_(length_ > 0 ? new value_type[length_] : NULL) {}

template<class Abc>
Sequence<Abc>::Sequence(FILE* in)
        : length_(0),
          seq_(NULL) {
    Read(in);
}

template<class Abc>
Sequence<Abc>::Sequence(const std::string& sequence, const std::string& header)
        : length_(0),
          seq_(NULL) {
    Init(sequence, header);
}

template<class Abc>
Sequence<Abc>::Sequence(const Sequence& other)
        : length_(other.length_),
          seq_(length_ > 0 ? new value_type[length_] : NULL),
          header_(other.header_) {
    for (size_t i = 0; i < length_; ++i) seq_[i] = other[i];
}

template<class Abc>
Sequence<Abc>::Sequence(const Sequence& other, size_t idx, size_t len)
        : length_(len),
          seq_(len > 0 ? new value_type[len] : NULL),
          header_(other.header_) {
    for (size_t i = 0; i < length_; ++i) seq_[i] = other[idx + i];
}

template<class Abc>
Sequence<Abc>& Sequence<Abc>::operator= (const Sequence& rhs) {
    if (this == &rhs) return *this;  // handle self assignment
    delete[] seq_;

    length_ = rhs.length_;
    if (length_ > 0) {
        seq_ = new value_type[length_];
        for (size_t i = 0; i < length_; ++i) seq_[i] = rhs[i];
    } else {
        seq_ = NULL;
    }
    header_ = rhs.header_;

    return *this;
}

template<class Abc>
void Sequence<Abc>::Init(std::string sequence, std::string header) {
    assert(seq_ == NULL);

    // Init header using swap trick to trim excess capacity
    std::string(header).swap(header_);

    // Strip whitespace and newlines from sequence.
    sequence.erase(remove_if(sequence.begin(), sequence.end(), isspace), sequence.end());

    // Strip gap characters from sequence.
    sequence.erase(remove_if(sequence.begin(), sequence.end(), isgap), sequence.end());

    // First validate each character before copying
    const size_t seqlen = sequence.length();
    for (size_t i = 0; i < seqlen; ++i) {
        if (!Abc::kValidChar[static_cast<int>(sequence[i])])
            throw Exception("Invalid character with ASCII number %i at position %zu of sequence '%s'",
                            static_cast<int>(sequence[i]), i, header_.c_str());
    }
    // Copy character sequence in packed format into sequence array
    length_ = seqlen;
    seq_ = new value_type[length_];
    for (size_t i = 0; i < seqlen; ++i) seq_[i] = Abc::kCharToInt[static_cast<int>(sequence[i])];
}

template<class Abc>
void Sequence<Abc>::Read(FILE* fin) {
    delete [] seq_;
    const size_t kBuffSize = MB;
    char buffer[kBuffSize];
    int c = '\0';
    std::string header;
    std::string sequence;

    // Read header
    while (fgetline(buffer, kBuffSize, fin)) {
        if (!strscn(buffer)) continue;
        if (buffer[0] == '>') {
            header.append(buffer + 1);
            break;
        } else {
            throw Exception("Sequence header does not start with '>'!");
        }
    }
    // Read sequence and stop if either a new header or delimiter is found
    while (fgetline(buffer, kBuffSize, fin) && buffer[0] != '/' && buffer[1] != '/') {
        if (strscn(buffer))
            sequence.append(buffer);

        c = getc(fin);
        if (c == EOF) break;
        ungetc(c, fin);
        if (static_cast<char>(c) == '>') break;
    }
    Init(sequence, header);
}

template<class Abc>
void Sequence<Abc>::Write(FILE* fout, size_t width) const {
    fprintf(fout, ">%s\n", header_.c_str());
    for (size_t i = 0; i < length(); ++i) {
        // fprintf(stdout, "%i\n", chr(i));
        fputc(chr(i), fout);
        if ((i+1) % width == 0) fputc('\n', fout);
    }
    if (length() % width != 0) fputc('\n', fout);
}

template<class Abc>
std::string Sequence<Abc>::ToString() const {
    std::string s(length_, '\0');
    for (size_t i = 0; i < length_; ++i)
        s[i] = Abc::kIntToChar[seq_[i]];
    return s;
}

template<class Abc>
void Sequence<Abc>::Resize(size_t newlen) {
    if (newlen != length_) {
        if (seq_ != NULL) delete[] seq_;
        length_ = newlen;
        seq_ = length_ > 0 ? new value_type[length_] : NULL;
    }
}

}  // namespace cs

#endif  // CS_SEQUENCE_INL_H_
