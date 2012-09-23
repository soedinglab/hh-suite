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

#ifndef CS_BITSET_H_
#define CS_BITSET_H_

namespace cs {

class Bitset {
 public:
  Bitset(uint32_t sz) {
    uint32_t nwords = (sz >> 5)+1;
    words_ = new uint32_t[nwords];

    assert(words_ != NULL);
    memset(words_, 0, nwords * 4);
    sz_ = nwords << 5;
  }

  ~Bitset() {
    delete[] words_;
  }

  bool test(uint32_t i) const {
    bool ret = false;
    if(i < sz_) {
      ret = ((words_[i >> 5] >> (i & 0x1f)) & 1) != 0;
    }
    return ret;
  }

  void set(uint32_t i) {
    while(i >= sz_) { expand(); }
    // Fast path
    words_[i >> 5] |= (1 << (i & 0x1f));
    assert(((words_[i >> 5] >> (i & 0x1f)) & 1) == 1);
  }

  void flip(uint32_t i) {
    while(i >= sz_) { expand(); }
    // Fast path
    words_[i >> 5] ^= (1 << (i & 0x1f));
  }

  size_t size() const { return sz_; }

private:
  void expand() {
    uint32_t *newwords = realloc(sz_, words_);
    delete[] words_;   // delete old array
    words_ = newwords; // install new array
  }

  uint32_t* realloc(uint32_t& sz, uint32_t* words) {
    uint32_t oldsz = sz;
    if(sz > 0) {
      sz += (sz >> 1) + 31; // Add 50% more elements, plus a bit
                sz &= ~31;            // Make sure it's 32-aligned
    } else {
      sz = 1024; // Start off at 1024 bits to avoid many expansions
    }
    uint32_t *newwords;
    newwords = new uint32_t[sz >> 5 /* convert to words */];
    if(oldsz > 0) {
      // Move old values into new array
      memcpy(newwords, words, oldsz >> 3 /* convert to bytes */);
    }
    // Initialize all new words to 0
    memset(newwords + (oldsz >> 5 /*convert to words*/), 0,
           (sz - oldsz) >> 3 /* convert to bytes */);
    return newwords; // return new array
  }

  uint32_t sz_;        // size as # of bits
  uint32_t *words_;    // storage

  DISALLOW_COPY_AND_ASSIGN(Bitset);
};

// Equality operator for bitsets: two bitsets are equal iff they have same size
// and same bitpattern.
inline bool operator== (const Bitset& lhs, const Bitset& rhs) {
  if (lhs.size() != rhs.size()) return false;
  bool rv = true;
  for (size_t b = 0; b < lhs.size(); ++b)
    if (lhs.test(b) != rhs.test(b)) { rv = false; break; }
  return rv;
}

}  // namespace cs

#endif  // CS_BITSET_H_
