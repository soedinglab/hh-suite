// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "as.h"

namespace cs {

// Size of AS62 alphabet
const size_t AS62::kSize = 62;

// Size of AS62 alphabet including wildcard character ANY
const size_t AS62::kSizeAny = 63;

// Integer code of ANY character
const uint8_t AS62::kAny = 62;

// Integer code of GAP
const uint8_t AS62::kGap = 63;

// Integer code of ENDGAP
const uint8_t AS62::kEndGap = 64;

// For converting from ASCII to integer-based AS62-code
const uint8_t AS62::kCharToInt[] = {
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  62,   0,   0,  63,  63,   0,
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   0,   0,   0,   0,   0,   0,
   0,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,
  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,   0,   0,   0,   0,   0,
   0,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,
  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
};

// For converting from integer code back to ASCII character
const char AS62::kIntToChar[] = {
 '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '*', '-',
 '-'
};

// For testing if ASCII character is from AS62 alphabet
const bool AS62::kValidChar[] = {
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false,  true, false, false,  true,  true, false,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false, false,
 false,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false,
 false,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
};



// Size of AS90 alphabet
const size_t AS90::kSize = 90;

// Size of AS90 alphabet including wildcard character ANY
const size_t AS90::kSizeAny = 91;

// Integer code of ANY character
const uint8_t AS90::kAny = 90;

// Integer code of GAP
const uint8_t AS90::kGap = 91;

// Integer code of ENDGAP
const uint8_t AS90::kEndGap = 02;

// For converting from ASCII to integer-based AS90-code
const uint8_t AS90::kCharToInt[] = {
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   1,   2,   3,   4,   5,   6,   7,   8,  90,   9,  10,  91,  91,  11,
  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,   0,  26,
  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,
  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,
  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
};

// For converting from integer code back to ASCII character
const char AS90::kIntToChar[] = {
 '!', '"', '#', '$', '%', '&', '\'', '(', ')', '+', ',', '/', '0', '1', '2', '3',
 '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '?', '@', 'A', 'B', 'C', 'D',
 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd',
 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
 'u', 'v', 'w', 'x', 'y', 'z', '{', '|', '}', '~', '*', '-', '-'
};

// For testing if ASCII character is from AS90 alphabet
const bool AS90::kValidChar[] = {
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, false,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
};




// Size of AS219 alphabet
const size_t AS219::kSize = 219;

// Size of AS219 alphabet including wildcard character ANY
const size_t AS219::kSizeAny = 220;

// Integer code of ANY character
const uint8_t AS219::kAny = 219;

// Integer code of GAP
const uint8_t AS219::kGap = 220;

// Integer code of ENDGAP
const uint8_t AS219::kEndGap = 221;

// For converting from ASCII to integer-based AS219-code
const uint8_t AS219::kCharToInt[] = {
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   1,   2,   3,   4,   5,   6,   7,   8, 219,   9,  10, 220, 220,  11,
  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,   0,  26,
  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,
  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,
  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106,
 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154,
 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186,
 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202,
 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218
};

// For converting from integer code back to ASCII character
const int AS219::kIntToChar[] = {
  33,  34,  35,  36,  37,  38,  39,  40,  41,  43,  44,  47,  48,  49,  50,  51,
  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  63,  64,  65,  66,  67,  68,
  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100,
 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132,
 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228,
 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244,
 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,  42,  45,  45
 // '\x21', '\x22', '\x23', '\x24', '\x25', '\x26', '\x27', '\x28', '\x29', '\x2B', '\x2C', '\x2F', '\x30', '\x31', '\x32', '\x33',
 // '\x34', '\x35', '\x36', '\x37', '\x38', '\x39', '\x3A', '\x3B', '\x3C', '\x3D', '\x3F', '\x40', '\x41', '\x42', '\x43', '\x44',
 // '\x45', '\x46', '\x47', '\x48', '\x49', '\x4A', '\x4B', '\x4C', '\x4D', '\x4E', '\x4F', '\x50', '\x51', '\x52', '\x53', '\x54',
 // '\x55', '\x56', '\x57', '\x58', '\x59', '\x5A', '\x5B', '\x5C', '\x5D', '\x5E', '\x5F', '\x60', '\x61', '\x62', '\x63', '\x64',
 // '\x65', '\x66', '\x67', '\x68', '\x69', '\x6A', '\x6B', '\x6C', '\x6D', '\x6E', '\x6F', '\x70', '\x71', '\x72', '\x73', '\x74',
 // '\x75', '\x76', '\x77', '\x78', '\x79', '\x7A', '\x7B', '\x7C', '\x7D', '\x7E', '\x7F', '\x80', '\x81', '\x82', '\x83', '\x84',
 // '\x85', '\x86', '\x87', '\x88', '\x89', '\x8A', '\x8B', '\x8C', '\x8D', '\x8E', '\x8F', '\x90', '\x91', '\x92', '\x93', '\x94',
 // '\x95', '\x96', '\x97', '\x98', '\x99', '\x9A', '\x9B', '\x9C', '\x9D', '\x9E', '\x9F', '\xA0', '\xA1', '\xA2', '\xA3', '\xA4',
 // '\xA5', '\xA6', '\xA7', '\xA8', '\xA9', '\xAA', '\xAB', '\xAC', '\xAD', '\xAE', '\xAF', '\xB0', '\xB1', '\xB2', '\xB3', '\xB4',
 // '\xB5', '\xB6', '\xB7', '\xB8', '\xB9', '\xBA', '\xBB', '\xBC', '\xBD', '\xBE', '\xBF', '\xC0', '\xC1', '\xC2', '\xC3', '\xC4',
 // '\xC5', '\xC6', '\xC7', '\xC8', '\xC9', '\xCA', '\xCB', '\xCC', '\xCD', '\xCE', '\xCF', '\xD0', '\xD1', '\xD2', '\xD3', '\xD4',
 // '\xD5', '\xD6', '\xD7', '\xD8', '\xD9', '\xDA', '\xDB', '\xDC', '\xDD', '\xDE', '\xDF', '\xE0', '\xE1', '\xE2', '\xE3', '\xE4',
 // '\xE5', '\xE6', '\xE7', '\xE8', '\xE9', '\xEA', '\xEB', '\xEC', '\xED', '\xEE', '\xEF', '\xF0', '\xF1', '\xF2', '\xF3', '\xF4',
 // '\xF5', '\xF6', '\xF7', '\xF8', '\xF9', '\xFA', '\xFB', '\xFC', '\xFD', '\xFE', '\xFF', '*', '-', '-'
};

// For testing if ASCII character is from AS219 alphabet
const bool AS219::kValidChar[] = {
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
 false,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true, false,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true
};

}  // namespace cs
