/*
 * hhutil-inl.h
 *
 *  Created on: Mar 28, 2014
 *      Author: meiermark
 */

#ifndef HHUTIL_INL_H_
#define HHUTIL_INL_H_

#include "simd.h"

/////////////////////////////////////////////////////////////////////////////////////
// Transform a character to lower case and '.' to '-' and vice versa
/////////////////////////////////////////////////////////////////////////////////////
inline char MatchChr(char c)  {
	return ((c>='a' && c<='z')? c-'a'+'A' : (c=='.'? '-':c) );
}

inline char InsertChr(char c) {
	return ((c>='A' && c<='Z')? c+'a'-'A' : ((c>='0' && c<='9') || c=='-')? '.':c );
}

inline int  WordChr(char c) {
	return (int)((c>='A' && c<='Z') || (c>='a' && c<='z'));
}

// Compute the sum of bits of one or two integers
inline int NumberOfSetBits(int i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//TODO: check
//inline int NumberOfSetBits(int i)
//{
//  return _mm_popcnt_u32(i);
//}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms the one-letter amino acid code into an integer between 0 and 22
/////////////////////////////////////////////////////////////////////////////////////
inline char aa2i(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case 'A': return 0;
    case 'R': return 1;
    case 'N': return 2;
    case 'D': return 3;
    case 'C': return 4;
    case 'Q': return 5;
    case 'E': return 6;
    case 'G': return 7;
    case 'H': return 8;
    case 'I': return 9;
    case 'L': return 10;
    case 'K': return 11;
    case 'M': return 12;
    case 'F': return 13;
    case 'P': return 14;
    case 'S': return 15;
    case 'T': return 16;
    case 'W': return 17;
    case 'Y': return 18;
    case 'V': return 19;
    case 'X': return ANY;
    case 'J': return ANY;
    case 'O': return ANY;
    case 'U': return 4;  //Selenocystein -> Cystein
    case 'B': return 3;  //D (or N)
    case 'Z': return 6;  //E (or Q)
    case '-': return GAP;
    case '.': return GAP;
    case '_': return GAP;
    }
  if (c>=0 && c<=32) return -1; // white space and control characters
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 22 into the one-letter amino acid code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2aa(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  switch (c)
    {
    case 0: return 'A';
    case 1: return 'R';
    case 2: return 'N';
    case 3: return 'D';
    case 4: return 'C';
    case 5: return 'Q';
    case 6: return 'E';
    case 7: return 'G';
    case 8: return 'H';
    case 9: return 'I';
    case 10: return 'L';
    case 11: return 'K';
    case 12: return 'M';
    case 13: return 'F';
    case 14: return 'P';
    case 15: return 'S';
    case 16: return 'T';
    case 17: return 'W';
    case 18: return 'Y';
    case 19: return 'V';
    case ANY: return 'X';
    case GAP: return '-';
    case ENDGAP: return '-';
    }
  return '?';
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms the dssp/psipred secondary structure code into an integer number
/////////////////////////////////////////////////////////////////////////////////////
inline char ss2i(char c)
{
  //- H E C S T G B
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'X': return 0;
    case 'H': return 1;
    case 'E': return 2;
    case 'C': return 3;
    case '~': return 3;
    case 'S': return 4;
    case 'T': return 5;
    case 'G': return 6;
    case 'B': return 7;
    case 'I': return 3;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 8 into the dssp/psipred secondary structure code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2ss(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'H';
    case 2: return 'E';
    case 3: return 'C';
    case 4: return 'S';
    case 5: return 'T';
    case 6: return 'G';
    case 7: return 'B';
    case 8: return 'I';
    }
  return '?';
}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms the solvend accessiblity code into an integer number
/////////////////////////////////////////////////////////////////////////////////////
inline char sa2i(char c)
{
  //- A B C D E
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'A': return 1;
    case 'B': return 2;
    case 'C': return 3;
    case 'D': return 4;
    case 'E': return 5;
    case 'F': return 6;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 5 into the solvent accessibility code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2sa(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'A';
    case 2: return 'B';
    case 3: return 'C';
    case 4: return 'D';
    case 5: return 'E';
    case 6: return 'F';
    }
  return '?';
}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms alternative secondary structure symbols into symbols
/////////////////////////////////////////////////////////////////////////////////////
inline char ss2ss(char c)
{
  //- H E C S T G B
  switch (c)
    {
    case '~': return 'C';
    case 'I': return 'C';
    case 'i': return 'c';
    case 'H':
    case 'E':
    case 'C':
    case 'S':
    case 'T':
    case 'G':
    case 'B':
    case 'h':
    case 'e':
    case 'c':
    case 's':
    case 't':
    case 'g':
    case 'b':
    case '.':
      return c;
    }
  return '-';
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms confidence values of psipred into internal code
/////////////////////////////////////////////////////////////////////////////////////
inline char cf2i(char c)
{
  switch (c)
    {
    case '-': return 0;
    case '.': return 0;
    case '0': return 1;
    case '1': return 2;
    case '2': return 3;
    case '3': return 4;
    case '4': return 5;
    case '5': return 6;
    case '6': return 7;
    case '7': return 8;
    case '8': return 9;
    case '9': return 10;
    }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms internal representation of psipred confidence values into printable chars
/////////////////////////////////////////////////////////////////////////////////////
inline char i2cf(char c)
{
  switch (c)
    {
    case 0: return '-';
    case 1: return '0';
    case 2: return '1';
    case 3: return '2';
    case 4: return '3';
    case 5: return '4';
    case 6: return '5';
    case 7: return '6';
    case 8: return '7';
    case 9: return '8';
    case 10: return '9';
    }
  return '-';
}


/////////////////////////////////////////////////////////////////////////////////////
// Fast lookup of log2(1+2^(-x)) for x>=0 (precision < 0.35%)
/////////////////////////////////////////////////////////////////////////////////////
inline float fast_addscore(float x)
{
  static float val[2001];         // val[i]=log2(1+2^(-x))
  static char initialized=0;
  if (x>20) return 0.0;
  if (x<0)
    {
      fprintf(stderr,"Error in function fast_addscore: argument %g is negative\n",x);
      exit(7);
    }
  if (!initialized)   //First fill in the log2-vector
    {
      for (int i=0; i<=2000; i++) val[i]=log2(1.0+pow(2,-0.01*(i+0.5)));
      initialized=1;
    }
  return val[(int)(100.0*x)];
}



/////////////////////////////////////////////////////////////////////////////////////
// Little utilities for output
/////////////////////////////////////////////////////////////////////////////////////
inline void fout(FILE* outf, int d)
{
  if (d>=99999) fprintf(outf,"*\t"); else fprintf(outf,"%i\t",d);
  return;
}


inline void sout(std::stringstream& out, int d) {
  if (d>=99999)
    out << "*\t";
  else
    out << d << "\t";
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Takes family code (eg. a.1.2.3) and returns strings 'a', 'a.1', and 'a.1.2'
/////////////////////////////////////////////////////////////////////////////////////
inline void  ScopID(char cl[], char fold[], char sfam[], const char fam[])
{
  char* ptr;

  //get scop class ID
  strcpy(cl,fam);
  ptr = strchr(cl,'.');               //return adress of next '.' in name
  if(ptr) ptr[0]='\0';

  //get scop fold ID
  strcpy(fold,fam);
  ptr = strchr(fold,'.');             //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');    //return adress of next '.' in name
  if(ptr) ptr[0]='\0';

  //get scop superfamily ID
  strcpy(sfam,fam);
  ptr = strchr(sfam,'.');            //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr[0]='\0';
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// SIMD 2^x for four floats
// Calculate float of 2pow(x) for four floats in parallel with SSE2
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 4.6E-6  (< 2.3E-7 with 5'th order polynomial)
//
// Internal representation of float number according to IEEE 754 (__m128 --> 4x):
//   1bit sign, 8 bits exponent, 23 bits mantissa: seee eeee emmm mmmm mmmm mmmm mmmm mmmm
//                                    0x4b400000 = 0100 1011 0100 0000 0000 0000 0000 0000
//   In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
/////////////////////////////////////////////////////////////////////////////////////
inline simd_float simdf32_fpow2(simd_float X) {

    simd_int* xPtr = (simd_int*) &X;    // store address of float as pointer to int

    const simd_float CONST32_05f       = simdf32_set(0.5f); // Initialize a vector (4x32) with 0.5f
    // (3 << 22) --> Initialize a large integer vector (shift left)
    const simd_int CONST32_3i          = simdi32_set(3);
    const simd_int CONST32_3shift22    = simdi32_slli(CONST32_3i, 22);
    const simd_float CONST32_1f        = simdf32_set(1.0f);
    const simd_float CONST32_FLTMAXEXP = simdf32_set(FLT_MAX_EXP);
    const simd_float CONST32_FLTMAX    = simdf32_set(FLT_MAX);
    const simd_float CONST32_FLTMINEXP = simdf32_set(FLT_MIN_EXP);
    // fifth order
    const simd_float CONST32_A = simdf32_set(0.00187682f);
    const simd_float CONST32_B = simdf32_set(0.00898898f);
    const simd_float CONST32_C = simdf32_set(0.0558282f);
    const simd_float CONST32_D = simdf32_set(0.240153f);
    const simd_float CONST32_E = simdf32_set(0.693153f);

    simd_float tx;
    simd_int lx;
    simd_float dx;
    simd_float result    = simdf32_set(0.0f);
    simd_float maskedMax = simdf32_set(0.0f);
    simd_float maskedMin = simdf32_set(0.0f);

    // Check wheter one of the values is bigger or smaller than FLT_MIN_EXP or FLT_MAX_EXP
    // The correct FLT_MAX_EXP value is written to the right place
    maskedMax = simdf32_gt(X, CONST32_FLTMAXEXP);
    maskedMin = simdf32_gt(X, CONST32_FLTMINEXP);
    maskedMin = simdf32_xor(maskedMin, maskedMax);
    // If a value is bigger than FLT_MAX_EXP --> replace the later result with FLTMAX
    maskedMax = simdf32_and(CONST32_FLTMAX, simdf32_gt(X, CONST32_FLTMAXEXP));

    tx = simdf32_add((simd_float ) CONST32_3shift22, simdf32_sub(X, CONST32_05f)); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
                                                                             // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
                                                                             // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)

    lx = simdf32_f2i(tx);                                       // integer value of x

    dx = simdf32_sub(X, simdi32_i2f(lx));                       // float remainder of x

    //   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
    //            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
    //            + dx*(0.0558282f            // Speed: 2.3E-8s
    //            + dx*(0.00898898f
    //            + dx* 0.00187682f ))));
    X = simdf32_mul(dx, CONST32_A);
    X = simdf32_add(CONST32_B, X);  // add constant B
    X = simdf32_mul(dx, X);
    X = simdf32_add(CONST32_C, X);  // add constant C
    X = simdf32_mul(dx, X);
    X = simdf32_add(CONST32_D, X);  // add constant D
    X = simdf32_mul(dx, X);
    X = simdf32_add(CONST32_E, X);  // add constant E
    X = simdf32_mul(dx, X);
    X = simdf32_add(X, CONST32_1f); // add 1.0f

    simd_int lxExp = simdi32_slli(lx, 23); // add integer power of 2 to exponent

    *xPtr = simdi32_add(*xPtr, lxExp); // add integer power of 2 to exponent

    // Add all Values that are greater than min and less than max
    result = simdf32_and(maskedMin, X);
    // Add MAX_FLT values where entry values were > FLT_MAX_EXP
    result = simdf32_or(result, maskedMax);

    return result;
}

// Fast SIMD log2 for four floats
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

inline simd_float simdf32_flog2(simd_float X) {

  const simd_int CONST32_0x7f = simdi32_set(0x7f);
  const simd_int CONST32_0x7fffff = simdi32_set(0x7fffff);
  const simd_int CONST32_0x3f800000 = simdi32_set(0x3f800000);
  const simd_float  CONST32_1f = simdf32_set(1.0);
  // const float a=0.1564, b=-0.5773, c=1.0-a-b;  // third order
  const float a=0.0440047f, b=-0.1903190f, c=0.4123442f, d=-0.7077702f, e=1.0-a-b-c-d; // fifth order
  const simd_float  CONST32_A = simdf32_set(a);
  const simd_float  CONST32_B = simdf32_set(b);
  const simd_float  CONST32_C = simdf32_set(c);
  const simd_float  CONST32_D = simdf32_set(d);
  const simd_float  CONST32_E = simdf32_set(e);
  simd_int E; // exponents of X
  simd_float R; //  result
  E = simdi32_srli((simd_int) X, 23);    // shift right by 23 bits to obtain exponent+127
  E = simdi32_sub(E, CONST32_0x7f);     // subtract 127 = 0x7f
  X = (simd_float) simdi_and((simd_int) X, CONST32_0x7fffff);  // mask out exponent => mantisse
  X = (simd_float) simdi_or ((simd_int) X, CONST32_0x3f800000); // set exponent to 127 (i.e., 0)
  X = simdf32_sub(X, CONST32_1f);          // subtract one from mantisse
  R = simdf32_mul(X, CONST32_A);           // R = a*X
  R = simdf32_add(R, CONST32_B);           // R = a*X+b
  R = simdf32_mul(R, X);                   // R = (a*X+b)*X
  R = simdf32_add(R, CONST32_C);           // R = (a*X+b)*X+c
  R = simdf32_mul(R, X);                   // R = ((a*X+b)*X+c)*X
  R = simdf32_add(R, CONST32_D);           // R = ((a*X+b)*X+c)*X+d
  R = simdf32_mul(R, X);                   // R = (((a*X+b)*X+c)*X+d)*X
  R = simdf32_add(R, CONST32_E);           // R = (((a*X+b)*X+c)*X+d)*X+e
  R = simdf32_mul(R, X);                   // R = ((((a*X+b)*X+c)*X+d)*X+e)*X ~ log2(1+X) !!
  R = simdf32_add(R, simdi32_i2f(E));  // convert integer exponent to float and add to mantisse
  return R;

}

#define LOG_POLY_DEGREE 4
#define POLY0(x, c0) simdf32_set(c0)
#define POLY1(x, c0, c1) simdf32_add(simdf32_mul(POLY0(x, c1), x), simdf32_set(c0))
#define POLY2(x, c0, c1, c2) simdf32_add(simdf32_mul(POLY1(x, c1, c2), x), simdf32_set(c0))
#define POLY3(x, c0, c1, c2, c3) simdf32_add(simdf32_mul(POLY2(x, c1, c2, c3), x), simdf32_set(c0))
#define POLY4(x, c0, c1, c2, c3, c4) simdf32_add(simdf32_mul(POLY3(x, c1, c2, c3, c4), x), simdf32_set(c0))
#define POLY5(x, c0, c1, c2, c3, c4, c5) simdf32_add(simdf32_mul(POLY4(x, c1, c2, c3, c4, c5), x), simdf32_set(c0))

inline simd_float log2f4(simd_float x)
{
    simd_int exp  = simdi32_set(0x7F800000);
    simd_int mant = simdi32_set(0x007FFFFF);

    simd_float one = simdf32_set( 1.0f);

    simd_int i = simdf_f2icast(x);

    simd_float e = simdi32_i2f(simdi32_sub(simdi32_srli(simdi_and(i, exp), 23), simdi32_set(127)));

    simd_float m = simdf32_or(simdi_i2fcast(simdi_and(i, mant)), one);

    simd_float p;

    /* Minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ */
#if LOG_POLY_DEGREE == 6
    p = POLY5( m, 3.1157899f, -3.3241990f, 2.5988452f, -1.2315303f,  3.1821337e-1f, -3.4436006e-2f);
#elif LOG_POLY_DEGREE == 5
    p = POLY4(m, 2.8882704548164776201f, -2.52074962577807006663f, 1.48116647521213171641f, -0.465725644288844778798f, 0.0596515482674574969533f);
#elif LOG_POLY_DEGREE == 4
    p = POLY3(m, 2.61761038894603480148f, -1.75647175389045657003f, 0.688243882994381274313f, -0.107254423828329604454f);
#elif LOG_POLY_DEGREE == 3
    p = POLY2(m, 2.28330284476918490682f, -1.04913055217340124191f, 0.204446009836232697516f);
#else
#error
#endif

    /* This effectively increases the polynomial degree by one, but ensures that log2(1) == 0*/
    p = simdf32_mul(p, simdf32_sub(m, one));

    return simdf32_add(p, e);
}

// Perform log-sum-exp calculation with six SIMD variables
//              result = x1 + x2 + x3 + x4 + x5 + x6
//  -->     result = log2( 2^(x1) + 2^(x2) + 2^(x3) + 2^(x4) + 2^(x5) + 2^(x6))
//      // to prevent overflows apply log sum of exp trick
//      --> xMax = max(x1, x2, x3, x4, x5, x6)
//  -->     result = log2(2^(xMax) * (2^(x1 - xMax) + 2^(x2 - xMax) + 2^(x3 - xMax) + 2^(x4 - xMax) + 2^(x5 - xMax) + 2^(x6 - xMax)))
//  -->     result = log2(2^(xMax)) + log2(2^(x1 - xMax) + 2^(x2 - xMax) + 2^(x3 - xMax) + 2^(x4 - xMax) + 2^(x5 - xMax) + 2^(x6 - xMax))
//  -->     result = xMax + log2(2^(x1 - xMax) + 2^(x2 - xMax) + 2^(x3 - xMax) + 2^(x4 - xMax) + 2^(x5 - xMax) + 2^(x6 - xMax))
// WHERE x1, x2, x3, x4, x5 and x6 contain the log-values respectively!
inline simd_float simd_flog2_sum_fpow2(simd_float x1, simd_float x2,
        simd_float x3, simd_float x4, simd_float x5, simd_float x6) {

    // Calculate the maximum out of the six variables
    simd_float x_max0 = simdf32_max(x1,x2);
    simd_float x_max1 = simdf32_max(x3,x4);
    simd_float x_max2 = simdf32_max(x5,x6);
    x_max0 = simdf32_max(x_max0, x_max1);
    simd_float x_max  = simdf32_max(x_max2, x_max0);

    simd_float max_comp_vec = simdf32_set(-FLT_MAX);
    simd_float term0 = simdf32_fpow2(simdf32_sub(x1, x_max));
    simd_float term1 = simdf32_fpow2(simdf32_sub(x2, x_max));
    simd_float term2 = simdf32_fpow2(simdf32_sub(x3, x_max));
    simd_float term3 = simdf32_fpow2(simdf32_sub(x4, x_max));
    simd_float term4 = simdf32_fpow2(simdf32_sub(x5, x_max));
    simd_float term5 = simdf32_fpow2(simdf32_sub(x6, x_max));
    
    term0 = simdf32_add(term0, term1);          //      2^(x1 - x_max)
    term2 = simdf32_add(term2, term3);          // +    2^(x2 - x_max)
    term4 = simdf32_add(term4, term5);          // +    2^(x3 - x_max)
    term0 = simdf32_add(term0, term2);          // +    2^(x4 - x_max)
    term4 = simdf32_add(term4, term0);          // +    2^(x5 - x_max)
     //      max(-FLT_MAX, x_max + log2(sum(terms)))
    return simdf32_max(max_comp_vec, simdf32_add(x_max, simdf32_flog2(term4)));
}

// Perform log-sum-exp calculation with three SIMD variables
//              result = x1 + x2 + x3
//  -->     result = log2( 2^(x1) + 2^(x2) + 2^(x3) )
//      // to prevent overflows apply log sum of exp trick
//      --> xMax = max(x1, x2, x3)
//  -->     result = log2(2^(xMax) * (2^(x1 - xMax) + 2^(x2 - xMax) + 2^(x3 - xMax)))
//  -->     result = log2(2^(xMax)) + log2(2^(x1 - xMax) + 2^(x2 - xMax) + 2^(x3 - xMax))
//  -->     result = xMax + log2(2^(x1 - xMax) + 2^(x2 - xMax) + 2^(x3 - xMax))
// WHERE x1, x2 and x3 contain the log-values respectively!
inline simd_float simd_flog2_sum_fpow2(simd_float x1, simd_float x2, simd_float x3) {

    // Calculate the maximum out of the six variables
    simd_float x_max = simdf32_max(x1, simdf32_max(x2, x3));

    simd_float max_comp_vec = simdf32_set(-FLT_MAX);


    return simdf32_max(max_comp_vec, simdf32_add(x_max, simdf32_flog2(  //      x_max + log2(
                    simdf32_add(simdf32_fpow2(simdf32_sub(x1, x_max)),                  //      2^(x1 - x_max)
                            simdf32_add(simdf32_fpow2(simdf32_sub(x2, x_max)),          // +    2^(x2 - x_max)
                                    simdf32_fpow2(simdf32_sub(x3, x_max)))))));                 // +    2^(x3 - x_max))

}

// Perform log-sum-exp calculation with two SIMD variables
//              result = x1 + x2
//  -->     result = log2(2^(x1) + 2^(x2))
//      // to prevent overflows apply log sum of exp trick
//      --> xMax = max(x1, x2)
//  -->     result = log2(2^(xMax) * (2^(x1 - xMax) + 2^(x2 - xMax)))
//  -->     result = log2(2^(xMax)) + log2(2^(x1 - xMax) + 2^(x2 - xMax))
//  -->     result = xMax + log2(2^(x1 - xMax) + 2^(x2 - xMax))
// WHERE x1, x2 and x3 contain the log-values respectively!
inline simd_float simd_flog2_sum_fpow2(simd_float x1, simd_float x2) {

    // Calculate the maximum out of the six variables
    simd_float x_max = simdf32_max(x1, x2);

    simd_float max_comp_vec = simdf32_set(-FLT_MAX);

    return simdf32_max(max_comp_vec, simdf32_add(x_max, simdf32_flog2(  //      fmax(-FLT_MAX, x_max + log2(
                    simdf32_add(simdf32_fpow2(simdf32_sub(x1, x_max)),                  //      2^(x1 - x_max)
                                        simdf32_fpow2(simdf32_sub(x2, x_max))))));              // +    2^(x2 - x_max))

}






#endif /* HHUTIL_INL_H_ */
