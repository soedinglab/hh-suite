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

#ifndef CS_PSEUDOCOUNTS_H_
#define CS_PSEUDOCOUNTS_H_

namespace cs {

// Forward declarations
template<class Abc>
class Sequence;
template<class Abc>
class Profile;
template<class Abc>
class CountProfile;

// Calculates pseudocount admixture for profile column.
struct Admix {
  public:
    Admix() {}
    virtual ~Admix() {}

    // Returns the admixture fraction tau given the Neff.
    virtual double operator() (double neff) const = 0;
    // Gets the parameter value adjusted to obtain a certain target Neff.
    virtual double GetTargetNeffParam() const = 0; 
    // Sets the parameter value adjusted to obtain a certain target Neff.
    virtual void SetTargetNeffParam(double p) = 0;
};

// Calculates constant pseudocount admixture independent of number of effective
// sequences.
struct ConstantAdmix : public Admix {
    ConstantAdmix(double a) : pca(a) {}
    virtual ~ConstantAdmix() {}

    virtual double operator() (double) const { 
      return pca; 
    }
    
    virtual double GetTargetNeffParam() const {
      return pca;
    }

    virtual void SetTargetNeffParam(double p) {
      pca = p;
    }

    double pca;
};

// Calculates divergence-dependent pseudocount admixture as in CS-BLAST.
// tau = A * (B + 1) / (B + Neff)
struct CSBlastAdmix : public Admix {
    CSBlastAdmix(double a, double b) : pca(a), pcb(b) {}
    virtual ~CSBlastAdmix() {}

    virtual double operator() (double neff) const {
        return MIN(1.0, pca * (pcb + 1.0) / (pcb + neff));
    }
    
    virtual double GetTargetNeffParam() const {
      return pca;
    }

    virtual void SetTargetNeffParam(double p) {
      pca = p;
    }

    double pca, pcb;
};

// Calculates divergence-dependent pseudocount admixture as in HHsearch.
struct HHsearchAdmix : public Admix {
    HHsearchAdmix(double a, double b, double c = 1.0) : pca(a), pcb(b), pcc(c) {}
    virtual ~HHsearchAdmix() {}

    virtual double operator() (double neff) const {
        double rv = 0.0;
        if (pcc == 1.0)
            rv = MIN(1.0, pca / (1.0 + neff / pcb));
        else
            rv = MIN(1.0, pca / (1.0 + pow(neff / pcb, pcc)));
        return rv;
    }
    
    virtual inline double GetTargetNeffParam() const {
      return pca;
    }

    virtual inline void SetTargetNeffParam(double p) {
      pca = p;
    }

    double pca, pcb, pcc;
};



// An abstract base class for pseudocount factories.
template<class Abc>
class Pseudocounts {
  public:
    Pseudocounts() : target_neff_(0.0), target_neff_delta_(0.025) {}

    virtual ~Pseudocounts() {}

    // Adds pseudocounts to sequence using admixture and returns normalized profile.
    Profile<Abc> AddTo(const Sequence<Abc>& seq, Admix& admix) const;

    // Adds pseudocounts to sequence using admixture and returns normalized profile.
    Profile<Abc> AddTo(const CountProfile<Abc>& cp, Admix& admix) const;

    // Gets the target Neff in the resulting profile after admixing pseudocounts.
    double GetTargetNeff() const {
      return target_neff_;
    }

    // Sets the target Neff in the resulting profile after admixing pseudocounts.
    void SetTargetNeff(double target_neff) {
      target_neff_ = target_neff;
    }

    // Gets maximal deviation from the target Neff.
    double GetTargetNeffDelta() const {
      return target_neff_delta_;
    }

    // Sets maximal deviation from the target Neff.
    void SetTargetNeffDelta(double target_neff_delta) {
      target_neff_delta_ = target_neff_delta;
    }

  private:
    // Adds pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const = 0;

    // Adds pseudocounts to alignment derived profile.
    virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const = 0;

    // Admixes Sequence q to Profile p.
    void AdmixTo(const Sequence<Abc>& q, Profile<Abc>& p, const Admix& admix) const;

    // Admixes CountProfile q to Profile q.
    void AdmixTo(const CountProfile<Abc>& q, Profile<Abc>& p, const Admix& admix) const;

    // Admixes q to Profile p such that the Neff in p converges to the target Neff.
    template<class T>
    double AdmixToTargetNeff(const T& q, Profile<Abc>& p, Admix& admix) const;

  private:
    double target_neff_;       // Target Neff in the resulting profile.
    double target_neff_delta_; // Maximal deviation from the target Neff.

  private:
    static const double kNormalize           = 1e-5; // Normalization threshold.
    static const double kTargetNeffParamMin  = 0.0;  // Minimal paramater value for adjusting to the target Neff.
    static const double kTargetNeffParamMax  = 1.0;  // Maximal parameter value for adjusting to the target Neff.
    static const double kTargetNeffParamInit = 0.5;  // Initial parameter value for adjusting to the target Neff.
    static const double kTargetNeffEps       = 0.01; // Convergence threshold for adjusting to the target Neff.

    DISALLOW_COPY_AND_ASSIGN(Pseudocounts);
};  // Pseudocounts

}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_H_
