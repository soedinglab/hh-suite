/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_PSEUDOCOUNTS_INL_H_
#define CS_PSEUDOCOUNTS_INL_H_
#include "pseudocounts.h"

namespace cs {


// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const Sequence<Abc>& seq, Admix& admix) const {
    Profile<Abc> p(seq.length());
    AddToSequence(seq, p);
    if (target_neff_ >= 1.0) {
      AdmixToTargetNeff(seq, p, admix);
    } else {
      AdmixTo(seq, p, admix);
    }
    for(size_t i = 0; i < seq.length(); ++i) p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const CountProfile<Abc>& cp, Admix& admix) const {
    Profile<Abc> p(cp.counts.length());
    AddToProfile(cp, p);
    if (target_neff_ >= 1.0) {
      AdmixToTargetNeff(cp, p, admix);
    } else {
      AdmixTo(cp, p, admix);
    }
    for(size_t i = 0; i < cp.counts.length(); ++i)
        p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

template<class Abc>
void Pseudocounts<Abc>::AdmixTo(const Sequence<Abc>& q, Profile<Abc>& p, const Admix& admix) const {
    double tau = admix(1.0);
    double t = 1 - tau;
    for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] *= tau;
        p[i][q[i]] += t;
    }
}

template<class Abc>
void Pseudocounts<Abc>::AdmixTo(const CountProfile<Abc>& q, Profile<Abc>& p, const Admix& admix) const {
    for (size_t i = 0; i < p.length(); ++i) {
        double tau = admix(q.neff[i]);
        double t = 1 - tau;
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = tau * p[i][a] + t * q.counts[i][a] / q.neff[i];
    }
}

// Adjusts the Neff in 'p' to 'neff' by admixing q and returns tau.
template<class Abc>
template<class T>
double Pseudocounts<Abc>::AdmixToTargetNeff(const T& q, Profile<Abc>& p, Admix& admix) const {

    double l = kTargetNeffParamMin;
    double r = kTargetNeffParamMax;
    admix.SetTargetNeffParam(kTargetNeffParamInit);
    Profile<Abc> pp;
    while (l < kTargetNeffParamMax - kTargetNeffEps && r > kTargetNeffParamMin + kTargetNeffEps) {
        pp = p;
        AdmixTo(q, pp, admix);
        double ne = Neff(pp);
        if (fabs(ne - target_neff_) <= target_neff_delta_) {
            break;
        } else {
            if (ne < target_neff_) l = admix.GetTargetNeffParam();
            else r = admix.GetTargetNeffParam();
        }
        admix.SetTargetNeffParam(0.5 * (l + r));
    }
    if (l > kTargetNeffParamMax - kTargetNeffEps) {
        admix.SetTargetNeffParam(kTargetNeffParamMax);
        AdmixTo(q, p, admix);
    } else if (r < kTargetNeffParamMin + kTargetNeffEps) {
        admix.SetTargetNeffParam(kTargetNeffParamMin);
        AdmixTo(q, p, admix);
    } else {
        p = pp;
    }
    return admix.GetTargetNeffParam();
}


}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_INL_H_
