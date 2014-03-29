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

#endif /* HHHIT_INL_H_ */
