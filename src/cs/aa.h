// Copyright 2009, Andreas Biegert

#ifndef CS_AA_H_
#define CS_AA_H_

namespace cs {

class AA {
 public:
  // Size of alphabet excluding wildcard character ANY
  static const size_t kSize;

  // Size of alphabet includding wildcard character ANY
  static const size_t kSizeAny;

  // Integer code of ANY character
  static const uint8_t kAny;

  // Integer code of GAP
  static const uint8_t kGap;

  // Integer code of ENDGAP
  static const uint8_t kEndGap;

  // For converting from ASCII to the amino acid code
  static const uint8_t kCharToInt[];

  // For converting from integer code back to ASCII character
  static const char kIntToChar[];

  // For testing if ASCII character is from amino acid code
  static const bool kValidChar[];

  // Functional groups of amino acid alphabet needed for coloring of profile logos
  static const int kFuncGroup[];

  // Name of this alphabet
  static const char kName[];

 private:
  DISALLOW_COPY_AND_ASSIGN(AA);
};

}  // namespace cs

#endif  // CS_AA_H_
