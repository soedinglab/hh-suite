/*
 * util.h
 *
 *  Created on: Mar 28, 2014
 *      Author: meiermark
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <cassert>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <cstring>
#include <cmath>
#include <climits>
#include <float.h>

#include "util-inl.h"

#ifdef HH_SSE41
#include <tmmintrin.h>   // SSSE3
#include <smmintrin.h>   // SSE4.1
#define HH_SSE3
#endif

#ifdef HH_SSE3
#include <pmmintrin.h>   // SSE3
#define HH_SSE2
#endif

#ifdef HH_SSE2
#ifndef __SUNPRO_C
#include <emmintrin.h>   // SSE2
#else
#include <sunmedia_intrin.h>
#endif
#endif

char *substr(char* substr, char* str, int a, int b);

// Similar to Perl's tr/abc/ABC/: Replaces all chars in str found in one list with characters from the second list
// Returns the number of replaced characters
int strtr(char* str, const char oldchars[], const char newchars[]);

// Similar to Perl's tr/abc//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, const char chars[]);

// Similar to Perl's tr/a-z//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, char char1, char char2);

// Counts the number of characters in str that are in range between char1 and char2
int strcount(char* str, char char1, char char2);

// transforms str into an all uppercase string
char* uprstr(char* str);

// transforms str into an all uppercase string
char* lwrstr(char* str);

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets pt to NULL
int strint(const char* ptr);

// Same as strint, but interpretes '*' as default
int strinta(char*& ptr, int deflt = 99999);

// Returns leftmost float in ptr and sets the pointer to first char after
// the float. If no float is found, returns FLT_MIN and sets pt to NULL
float strflt(char*& ptr);

// Same as strint, but interpretes '*' as default
float strflta(char*& ptr, float deflt = 99999);

void QSortInt(int v[], int k[], int left, int right, int up = +1);

// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortFloat(float v[], int k[], int left, int right, int up = +1);

void readU16(char** ptr, uint16_t &result);

void readU32(char** ptr, uint32_t &result);

#ifdef HH_SSE2
__m128 _mm_flog2_ps(__m128 X);
#endif


#endif /* UTIL_H_ */
