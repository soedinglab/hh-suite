/*
 * hhhit-inl.h
 *
 *  Created on: Mar 28, 2014
 *      Author: meiermark
 */

#ifndef HHHIT_INL_H_
#define HHHIT_INL_H_

#include "util.h"

// /////////////////////////////////////////////////////////////////////////////////////
// //// Function for Viterbi()
// /////////////////////////////////////////////////////////////////////////////////////
// inline float max2(const float& xMM, const float& xX, char& b)
// {
//   if (xMM>xX) { b=MM; return xMM;} else { b=SAME;  return xX;}
// }
inline float max2(const float& xMM, const float& xSAME, char& b,
    const unsigned char bit) {
  if (xMM > xSAME) {
    b |= bit;
    return xMM;
  }
  else { /* b |= 0x00!*/
    return xSAME;
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Functions for StochasticBacktrace()
/////////////////////////////////////////////////////////////////////////////////////

inline int pickprob2(const double& xMM, const double& xX, const int& state) {
  if ((xMM + xX) * frand() < xMM)
    return MM;
  else
    return state;
}
inline int pickprob3_GD(const double& xMM, const double& xDG,
    const double& xGD) {
  double x = (xMM + xDG + xGD) * frand();
  if (x < xMM)
    return MM;
  else if (x < xMM + xDG)
    return DG;
  else
    return GD;
}
inline int pickprob3_IM(const double& xMM, const double& xMI,
    const double& xIM) {
  double x = (xMM + xMI + xIM) * frand();
  if (x < xMM)
    return MM;
  else if (x < xMM + xMI)
    return MI;
  else
    return IM;
}
inline int pickprob6(const double& x0, const double& xMM, const double& xGD,
    const double& xIM, const double& xDG, const double& xMI) {
  double x = (x0 + xMM + xGD + xIM + xDG + xMI) * frand();
  x -= xMM;
  if (x < 0)
    return MM;
  x -= x0;
  if (x < 0)
    return STOP;
  x -= xGD;
  if (x < 0)
    return GD;
  x -= xIM;
  if (x < 0)
    return IM;
  if (x < xDG)
    return DG;
  else
    return MI;
}

inline int pickmax2(const double& xMM, const double& xX, const int& state) {
  if (xMM > xX)
    return MM;
  else
    return state;
}
inline int pickmax3_GD(const double& xMM, const double& xDG,
    const double& xGD) {
  char state;
  double x;
  if (xMM > xDG) {
    state = MM;
    x = xMM;
  }
  else {
    state = DG;
    x = xDG;
  }
  if (xGD > x) {
    state = GD;
    x = xGD;
  }
  return state;
}
inline int pickmax3_IM(const double& xMM, const double& xMI,
    const double& xIM) {
  char state;
  double x;
  if (xMM > xMI) {
    state = MM;
    x = xMM;
  }
  else {
    state = MI;
    x = xMI;
  }
  if (xIM > x) {
    state = IM;
    x = xIM;
  }
  return state;
}
inline int pickmax6(const double& x0, const double& xMM, const double& xGD,
    const double& xIM, const double& xDG, const double& xMI) {
  char state;
  double x;
  if (x0 > xMM) {
    state = STOP;
    x = x0;
  }
  else {
    state = MM;
    x = xMM;
  }
  if (xGD > x) {
    state = GD;
    x = xGD;
  }
  if (xIM > x) {
    state = IM;
    x = xIM;
  }
  if (xDG > x) {
    state = DG;
    x = xDG;
  }
  if (xMI > x) {
    state = MI;
    x = xMI;
  }
  return state;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Functions that calculate P-values and probabilities
/////////////////////////////////////////////////////////////////////////////////////

//// Evaluate the CUMULATIVE extreme value distribution at point x
//// p(s)ds = lamda * exp{ -exp[-lamda*(s-mu)] - lamda*(s-mu) } ds = exp( -exp(-x) - x) dx = p(x) dx
//// => P(s>S) = integral_-inf^inf {p(x) dx}  = 1 - exp{ -exp[-lamda*(S-mu)] }
inline double Pvalue(double x, double a[]) {
  //a[0]=lamda, a[1]=mu
  double h = a[0] * (x - a[1]);
  return (h > 10) ? exp(-h) : double(1.0) - exp(-exp(-h));
}

inline double Pvalue(float x, float lamda, float mu) {
  double h = lamda * (x - mu);
  return (h > 10) ? exp(-h) : (double(1.0) - exp(-exp(-h)));
}

inline double logPvalue(float x, float lamda, float mu) {
  double h = lamda * (x - mu);
  return (h > 10) ? -h :
         (h < -2.5) ? -exp(-exp(-h)) : log((double(1.0) - exp(-exp(-h))));
}

inline double logPvalue(float x, double a[]) {
  double h = a[0] * (x - a[1]);
  return (h > 10) ? -h :
         (h < -2.5) ? -exp(-exp(-h)) : log((double(1.0) - exp(-exp(-h))));
}

inline float ScalarProd20(const float* qi, const float* tj) {
    
#ifdef AVX
  float __attribute__((aligned(ALIGN_FLOAT))) res;
  __m256 P; // query 128bit SSE2 register holding 4 floats
  __m256 S; // aux register
  __m256 R; // result
  __m256* Qi = (__m256*) qi;
  __m256* Tj = (__m256*) tj;

  R = _mm256_mul_ps(*(Qi++),*(Tj++));
  P = _mm256_mul_ps(*(Qi++),*(Tj++));
  S = _mm256_mul_ps(*Qi,*Tj); // floats A, B, C, D, ?, ?, ? ,?
  R = _mm256_add_ps(R,P);     // floats 0, 1, 2, 3, 4, 5, 6, 7
  P = _mm256_permute2x128_si256(R, R, 0x01); // swap hi and lo 128 bits: 4, 5, 6, 7, 0, 1, 2, 3
  R = _mm256_add_ps(R,P);     // 0+4, 1+5, 2+6, 3+7, 0+4, 1+5, 2+6, 3+7
  R = _mm256_add_ps(R,S);     // 0+4+A, 1+5+B, 2+6+C, 3+7+D, ?, ?, ? ,?
  R = _mm256_hadd_ps(R,R);    // 04A15B, 26C37D, ?, ?, 04A15B, 26C37D, ?, ?
  R = _mm256_hadd_ps(R,R);    // 01234567ABCD, ?, 01234567ABCD, ?, 01234567ABCD, ?, 01234567ABCD, ?
  _mm256_store_ps(&res, R);
  return res;
#elif SSE
//#ifdef SSE
    float __attribute__((aligned(16))) res;
    __m128 P; // query 128bit SSE2 register holding 4 floats
    __m128 R;// result
    __m128* Qi = (__m128*) qi;
    __m128* Tj = (__m128*) tj;
    
    R = _mm_mul_ps(*(Qi++),*(Tj++));
    P = _mm_mul_ps(*(Qi++),*(Tj++));
    R = _mm_add_ps(R,P);
    P = _mm_mul_ps(*(Qi++),*(Tj++));
    R = _mm_add_ps(R,P);
    P = _mm_mul_ps(*(Qi++),*(Tj++));
    R = _mm_add_ps(R,P);
    P = _mm_mul_ps(*Qi,*Tj);
    R = _mm_add_ps(R,P);

    R = _mm_hadd_ps(R,R);
    R = _mm_hadd_ps(R,R);
//    P = _mm_shuffle_ps(R,R, _MM_SHUFFLE(2,0,2,0));
//    R = _mm_shuffle_ps(R,R, _MM_SHUFFLE(3,1,3,1));
//    R = _mm_add_ps(R,P);
//    P = _mm_shuffle_ps(R,R, _MM_SHUFFLE(2,0,2,0));
//    R = _mm_shuffle_ps(R,R, _MM_SHUFFLE(3,1,3,1));
//    R = _mm_add_ps(R,P);
    _mm_store_ss(&res, R);
    return res;
#endif
    return tj[0] * qi[0] + tj[1] * qi[1] + tj[2] * qi[2] + tj[3] * qi[3]
         + tj[4] * qi[4] + tj[5] * qi[5] + tj[6] * qi[6] + tj[7] * qi[7]
         + tj[8] * qi[8] + tj[9] * qi[9] + tj[10] * qi[10] + tj[11] * qi[11]
         + tj[12] * qi[12] + tj[13] * qi[13] + tj[14] * qi[14]
         + tj[15] * qi[15] + tj[16] * qi[16] + tj[17] * qi[17]
         + tj[18] * qi[18] + tj[19] * qi[19];
}

// Calculate score between columns i and j of two HMMs (query and template)
inline float ProbFwd(float* qi, float* tj) {
  return ScalarProd20(qi, tj); //
}

//Calculate score between columns i and j of two HMMs (query and template)
inline float Score(float* qi, float* tj) {
  return fast_log2(ScalarProd20(qi, tj));
}

// Calculate secondary structure score between columns i and j of two HMMs (query and template)
static inline float ScoreSS(const HMM* q, const HMM* t, const int i,
    const int j, const int ssm) {
  //TODO
  return 0.0;
}


//// Calculate score between columns i and j of two HMMs (query and template)
//inline float ProbFwd(float* qi, float* tj) {
//  return ScalarProd20(qi, tj); //
//}
//
////Calculate score between columns i and j of two HMMs (query and template)
//inline float Score(float* qi, float* tj) {
//  return fast_log2(ProbFwd(qi, tj));
//}

#endif /* HHHIT_INL_H_ */
