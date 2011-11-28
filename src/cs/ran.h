// Copyright 2009, Andreas Biegert

#include "utils.h"

namespace cs {

// Implementation of the highest quality recommended random number generator
// from Numerical Recipes. The period of the generator is 3.138e57.
struct Ran {
  typedef unsigned long long int Ullong;

  Ran(Ullong j) : v(4101842887655102017LL), w(1) {
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
  }

  inline Ullong int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }

  inline double doub() { return 5.42101086242752217E-20 * int64(); }

  inline unsigned int int32() { return (unsigned int)int64(); }

  inline unsigned int operator() (unsigned int n) { return int32() % n; }

  Ullong u,v,w;
};


// Normal diistribution generator also from Numerical Recipes.
struct Gaussian : public Ran {
  typedef unsigned long long int Ullong;

  Gaussian(double mmu, double ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig) {}

  double operator() () {
    double u,v,x,y,q;
    do {
      u = doub();
      v = 1.7156*(doub()-0.5);
      x = u - 0.449871;
      y = std::abs(v) + 0.386595;
      q = SQR(x) + y*(0.19600*y-0.25472*x);
    } while (q > 0.27597
             && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
    return mu + sig*v/u;
  }

  double mu,sig;
};

}  // namespace cs
