// Utility subroutines

#include "util.h"

/////////////////////////////////////////////////////////////////////////////////////
// Arithmetics
/////////////////////////////////////////////////////////////////////////////////////


// copies substring str[a,b] into substr and returns substr
char *substr(char* substr, char* str, int a, int b) {
  if (b < a) {
    int i = b;
    b = a;
    a = i;
  }
  if (b - a > 1000) {
    printf("Function substr: >1000 chars to copy. Exiting.\n");
    exit(6);
  }
  char* dest = substr;
  char* source = str + a;
  char* send = str + b;
  while (*source != '\0' && source <= send)
    *(dest++) = *(source++);
  *dest = '\0';
  return substr;
}


// Similar to Perl's tr/abc/ABC/: Replaces all chars in str found in one list with characters from the second list
// Returns the number of replaced characters
int strtr(char* str, const char oldchars[], const char newchars[]) {
  char* ptr;
  const char *plist;
  int ntr = 0;
  for (ptr = str; *ptr != '\0'; ptr++)
    for (plist = oldchars; *plist != '\0'; plist++)
      if (*ptr == *plist) {
        *ptr = newchars[plist - oldchars];
        ntr++;
        break;
      }
  return ntr;
}

// Similar to Perl's tr/abc//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, const char chars[]) {
  char* ptr0 = str;
  char* ptr1 = str;
  const char *plist;
  while (*ptr1 != '\0') {
    for (plist = chars; *plist != '\0'; plist++)
      if (*ptr1 == *plist)
        break;
    if (*plist == '\0') {
      *ptr0 = *ptr1;
      ptr0++;
    }
    ptr1++;
  }
  *ptr0 = *ptr1;
  return ptr1 - ptr0;
}

// Similar to Perl's tr/a-z//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, char char1, char char2) {
  char* ptr0 = str;
  char* ptr1 = str;
  while (*ptr1 != '\0') {
    if (*ptr1 >= char1 && *ptr1 <= char2) {
      *ptr0 = *ptr1;
      ptr0++;
    }
    ptr1++;
  }
  *ptr0 = *ptr1;
  return ptr1 - ptr0;
}

// Counts the number of characters in str that are in range between char1 and char2
int strcount(char* str, char char1, char char2) {
  char* ptr = str;
  int count = 0;
  while (*ptr != '\0') {
    if (*ptr >= char1 && *ptr <= char2)
      count++;
    ptr++;
  }
  return count;
}

// transforms str into an all uppercase string
char* uprstr(char* str) {
  char* s = str;
  while (*s != '\0') {
    if (*s >= 'a' && *s <= 'z')
      *s += 'A' - 'a';
    s++;
  }
  return (str);
}

// transforms str into an all uppercase string
char* lwrstr(char* str) {
  char* s = str;
  while (*s != '\0') {
    if (*s >= 'A' && *s <= 'Z')
      *s += 'a' - 'A';
    s++;
  }
  return (str);
}

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets pt to NULL
int strint(char*& ptr) {
  int i;
  char* ptr0 = ptr;
  if (!ptr)
    return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9'))
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return INT_MIN;
  }
  if (ptr > ptr0 && *(ptr - 1) == '-')
    i = -atoi(ptr);
  else
    i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9')
    ptr++;
  return i;
}

// Same as strint, but interpretes '*' as default
int strinta(char*& ptr, int deflt) {
  int i;
  if (!ptr)
    return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9') && *ptr != '*')
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return INT_MIN;
  }
  if (*ptr == '*') {
    ptr++;
    return deflt;
  }
  if (*(ptr - 1) == '-')
    i = atoi(ptr - 1);
  else
    i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9')
    ptr++;
  return i;
}

// Returns leftmost float in ptr and sets the pointer to first char after
// the float. If no float is found, returns FLT_MIN and sets pt to NULL
float strflt(char*& ptr) {
  float i;
  char* ptr0 = ptr;
  if (!ptr)
    return FLT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9'))
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return FLT_MIN;
  }
  if (ptr > ptr0 && *(ptr - 1) == '-')
    i = -atof(ptr);
  else
    i = atof(ptr);
  while ((*ptr >= '0' && *ptr <= '9') || *ptr == '.')
    ptr++;
  return i;
}

// Same as strint, but interpretes '*' as default
float strflta(char*& ptr, float deflt) {
  float i;
  if (!ptr)
    return FLT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9') && *ptr != '*')
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return FLT_MIN;
  }
  if (*ptr == '*') {
    ptr++;
    return deflt;
  }
  if (*(ptr - 1) == '-')
    i = -atof(ptr);
  else
    i = atof(ptr);
  while ((*ptr >= '0' && *ptr <= '9') || *ptr == '.')
    ptr++;
  return i;
}


// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortInt(int v[], int k[], int left, int right, int up) {
  int i;
  int last;   // last element to have been swapped

  if (left >= right)
    return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k, left, (left + right) / 2);
  last = left;
  // ... and swap all elements i SMALLER than the pivot
  // with an element that is LARGER than the pivot (element last+1):
  if (up == 1) {
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] < v[k[left]])
        swapi(k, ++last, i);
  }
  else
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] > v[k[left]])
        swapi(k, ++last, i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k, left, last);

  // Sort the elements left from the pivot and right from the pivot
  QSortInt(v, k, left, last - 1, up);
  QSortInt(v, k, last + 1, right, up);
}


// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortFloat(float v[], int k[], int left, int right, int up) {
  int i;
  int last;   // last element to have been swapped
  void swapi(int k[], int i, int j);

  if (left >= right)
    return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k, left, (left + right) / 2);
  last = left;
  // ... and swap all elements i SMALLER than the pivot
  // with an element that is LARGER than the pivot (element last+1):
  if (up == 1) {
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] < v[k[left]])
        swapi(k, ++last, i);
  }
  else
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] > v[k[left]])
        swapi(k, ++last, i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k, left, last);

  // Sort the elements left from the pivot and right from the pivot
  QSortFloat(v, k, left, last - 1, up);
  QSortFloat(v, k, last + 1, right, up);
}

// Fast SSE2 log2 for four floats
// Calculate integer of log2 for four floats in parallel with SSE2
// Maximum deviation: +/- 2.1E-5
// Run time: ~5.6ns on Intel core2 2.13GHz.
// For a negative argument, nonsense is returned. Otherwise, when <1E-38, a value
// close to -126 is returned and when >1.7E38, +128 is returned.
// The function makes use of the representation of 4-byte floating point numbers:
// seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign, eee eee e gives the exponent + 127 (in hex: 0x7f).
// The following 23 bits give the mantisse, the binary digits after the decimal
// point:  x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeeee-127)
// Therefore,  log2(x) = eeeeeeee-127 + log2(1.mmmmmm...)
//                     = eeeeeeee-127 + log2(1+y),  where y = 0.mmmmmm...
//                     ~ eeeeeeee-127 + ((a*y+b)*y+c)*y
// The coefficients a, b  were determined by a least squares fit, and c=1-a-b to get 1 at y=1.
// Lower/higher order polynomials may be used for faster or more precise calculation:
// Order 1: log2(1+y) ~ y
// Order 2: log2(1+y) = (a*y + 1-a)*y, a=-0.3427
//  => max dev = +/- 8E-3, run time ~ 3.8ns
// Order 3: log2(1+y) = ((a*y+b)*y + 1-a-b)*y, a=0.1564, b=-0.5773
//  => max dev = +/- 1E-3, run time ~ 4.4ns
// Order 4: log2(1+y) = (((a*y+b)*y+c)*y + 1-a-b-c)*y, a=-0.0803 b=0.3170 c=-0.6748
//  => max dev = +/- 1.4E-4, run time ~ 5.0ns?
// Order 5: log2(1+y) = ((((a*y+b)*y+c)*y+d)*y + 1-a-b-c-d)*y, a=0.0440047 b=-0.1903190 c=0.4123442 d=-0.7077702
//  => max dev = +/- 2.1E-5, run time ~ 5.6ns?

#ifdef HH_SSE2
__m128 _mm_flog2_ps(__m128 X)
{
  const __m128i CONST32_0x7f = _mm_set_epi32(0x7f,0x7f,0x7f,0x7f);
  const __m128i CONST32_0x7fffff = _mm_set_epi32(0x7fffff,0x7fffff,0x7fffff,0x7fffff);
  const __m128i CONST32_0x3f800000 = _mm_set_epi32(0x3f800000,0x3f800000,0x3f800000,0x3f800000);
  const __m128 CONST32_1f = _mm_set_ps(1.0,1.0,1.0,1.0);
  // const float a=0.1564, b=-0.5773, c=1.0-a-b;  // third order
  const float a=0.0440047, b=-0.1903190, c=0.4123442, d=-0.7077702, e=1.0-a-b-c-d;// fifth order
  const __m128 CONST32_A = _mm_set_ps(a,a,a,a);
  const __m128 CONST32_B = _mm_set_ps(b,b,b,b);
  const __m128 CONST32_C = _mm_set_ps(c,c,c,c);
  const __m128 CONST32_D = _mm_set_ps(d,d,d,d);
  const __m128 CONST32_E = _mm_set_ps(e,e,e,e);
  __m128i E;// exponents of X
  __m128 R;//  result

  E = _mm_srli_epi32((__m128i) X, 23);// shift right by 23 bits to obtain exponent+127
  E = _mm_sub_epi32(E, CONST32_0x7f);// subtract 127 = 0x7f
  X = (__m128) _mm_and_si128((__m128i) X, CONST32_0x7fffff);// mask out exponent => mantisse
  X = (__m128) _mm_or_si128((__m128i) X, CONST32_0x3f800000);// set exponent to 127 (i.e., 0)
  X = _mm_sub_ps(X, CONST32_1f);// subtract one from mantisse
  R = _mm_mul_ps(X, CONST32_A);// R = a*X
  R = _mm_add_ps(R, CONST32_B);// R = a*X+b
  R = _mm_mul_ps(R, X);// R = (a*X+b)*X
  R = _mm_add_ps(R, CONST32_C);// R = (a*X+b)*X+c
  R = _mm_mul_ps(R, X);// R = ((a*X+b)*X+c)*X
  R = _mm_add_ps(R, CONST32_D);// R = ((a*X+b)*X+c)*X+d
  R = _mm_mul_ps(R, X);// R = (((a*X+b)*X+c)*X+d)*X
  R = _mm_add_ps(R, CONST32_E);// R = (((a*X+b)*X+c)*X+d)*X+e
  R = _mm_mul_ps(R, X);// R = ((((a*X+b)*X+c)*X+d)*X+e)*X ~ log2(1+X) !!
  R = _mm_add_ps(R, _mm_cvtepi32_ps(E));// convert integer exponent to float and add to mantisse
  return R;
}
#endif


void readU16(char** ptr, uint16_t &result) {
  unsigned char array[2];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8);
}

void readU32(char** ptr, uint32_t &result) {
  unsigned char array[4];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;
  array[2] = (unsigned char) (**ptr);
  (*ptr)++;
  array[3] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);
}
