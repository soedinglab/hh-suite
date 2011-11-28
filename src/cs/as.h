// Copyright 2009, Andreas Biegert

#ifndef CS_AS_H_
#define CS_AS_H_

namespace cs {

class AS62 {
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

  // For converting from ASCII to the integer code
  static const uint8_t kCharToInt[];

  // For converting from integer code back to ASCII character
  static const char kIntToChar[];

  // For testing if ASCII character is valid
  static const bool kValidChar[];

private:
  DISALLOW_COPY_AND_ASSIGN(AS62);
};


class AS90 {
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

  // For converting from ASCII to the integer code
  static const uint8_t kCharToInt[];

  // For converting from integer code back to ASCII character
  static const char kIntToChar[];

  // For testing if ASCII character is valid
  static const bool kValidChar[];

 private:
  DISALLOW_COPY_AND_ASSIGN(AS90);
};


class AS219 {
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

  // For converting from ASCII to the integer code
  static const uint8_t kCharToInt[];

  // For converting from integer code back to ASCII character
  static const int kIntToChar[];

  // For testing if ASCII character is valid
  static const bool kValidChar[];

 private:
  DISALLOW_COPY_AND_ASSIGN(AS219);
};

}  // namespace cs

#endif  // CS_AS_H_
