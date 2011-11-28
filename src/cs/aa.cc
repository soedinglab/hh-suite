// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "aa.h"

namespace cs {

// Size of alphabet excluding wildcard character ANY
const size_t AA::kSize = 20;

// Size of alphabet includding wildcard character ANY
const size_t AA::kSizeAny = 21;

// Integer code of ANY character
const uint8_t AA::kAny = 20;

// Integer code of GAP
const uint8_t AA::kGap = 21;

// Integer code of ENDGAP
const uint8_t AA::kEndGap = 22;

// For converting from ASCII to the amino acid code
const uint8_t AA::kCharToInt[] = {
  /*   0 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*  16 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*  32 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 21, 21,  0,
  /*                                                             -   . */
  /*  48 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /*  64 */  0,  0,  3,  4,  3,  6, 13,  7,  8,  9, 20, 11, 10, 12,  2, 20,
  /*             A   B   C   D   E   F   G   H   I   J   K   L   M   N   O */
  /*  80 */ 14,  5,  1, 15, 16,  4, 19, 17, 20, 18,  6,  0,  0,  0,  0,  0,
  /*         P   Q   R   S   T   U   V   W   X   Y   Z */
  /*  96 */  0,  0,  3,  4,  3,  6, 13,  7,  8,  9, 20, 11, 10, 12,  2, 20,
  /*             a   b   c   d   e   f   g   h   i   j   k   l   m   n   o */
  /* 112 */ 14,  5,  1, 15, 16,  4, 19, 17, 20, 18,  6,  0,  0,  0,  0,  0,
  /*         p   q   r   s   t   u   v   w   x   y   z */
  /* 128 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 144 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 160 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 176 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 192 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 208 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 224 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* 240 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

// For converting from integer code back to ASCII character
const char AA::kIntToChar[] = {
  'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-', '-'
};

// For testing if ASCII character is from amino acid code
const bool AA::kValidChar[] = {
  /*   0 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  16 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  32 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,   true,   true,  false,
  /*  48 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /*  64 */  false,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,
  /*  80 */   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,  false,  false,  false,  false,  false,
  /*  96 */  false,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,
  /* 112 */   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,   true,  false,  false,  false,  false,  false,
  /* 128 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 144 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 160 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 176 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 192 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 208 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 224 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,
  /* 240 */  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false,  false
};

// Functional groups of amino acid alphabet needed for coloring of profile logos
const int AA::kFuncGroup[] = {
  4, 8, 6, 5, 3, 6, 5, 4, 7, 1, 1, 8, 1, 2, 9, 4, 4, 2, 2, 1, 0
};

// Shorthand name for this amino acid alphabet
const char AA::kName[] = "aa";

}  // namespace cs
