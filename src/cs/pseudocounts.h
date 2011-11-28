// Copyright 2009, Andreas Biegert

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

    virtual double operator() (double neff) const = 0;
};

// An abstract base class for pseudocount factories.
template<class Abc>
class Pseudocounts {
  public:
    Pseudocounts() {}
    virtual ~Pseudocounts() {}

    // Adds pseudocounts to sequence and returns normalized profile.
    Profile<Abc> AddTo(const Sequence<Abc>& seq, const Admix& pca) const {
        Profile<Abc> rv(seq.length());
        AddToSequence(seq, pca, rv);
        for(size_t i = 0; i < seq.length(); ++i) rv[i][Abc::kAny] = 0.0;
        Normalize(rv, 1.0);
        return rv;
    }

    // Adds pseudocounts to sequence and returns normalized profile.
    Profile<Abc> AddTo(const CountProfile<Abc>& cp, const Admix& pca) const {
        Profile<Abc> rv(cp.counts.length());
        AddToProfile(cp, pca, rv);
        for(size_t i = 0; i < cp.counts.length(); ++i) rv[i][Abc::kAny] = 0.0;
        Normalize(rv, 1.0);
        return rv;
    }

  private:
    // Adds pseudocounts to sequence and stores resulting frequencies in given
    // profile.
    virtual void AddToSequence(const Sequence<Abc>& seq, const Admix& pca, Profile<Abc>& p) const = 0;

    // Adds pseudocounts to alignment derived profile.
    virtual void AddToProfile(const CountProfile<Abc>& cp, const Admix& pca, Profile<Abc>& p) const = 0;

    DISALLOW_COPY_AND_ASSIGN(Pseudocounts);
};  // Pseudocounts


// Calculates constant pseudocount admixture independent of number of effective
// sequences.
struct ConstantAdmix : public Admix {
    ConstantAdmix(double a) : pca(a) {}
    virtual ~ConstantAdmix() {}

    virtual double operator() (double) const { return pca; }

    const double pca;
};

// Calculates divergence-dependent pseudocount admixture as in CS-BLAST
// tau = A * (B + 1) / (B + Neff)
struct CSBlastAdmix : public Admix {
    CSBlastAdmix(double a, double b) : pca(a), pcb(b) {}
    virtual ~CSBlastAdmix() {}

    virtual double operator() (double neff) const {
        return MIN(1.0, pca * (pcb + 1.0) / (pcb + neff));
    }

    double pca, pcb;
};

// Calculates divergence-dependent pseudocount admixture as in HHsearch
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

    double pca, pcb, pcc;
};

}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_H_
