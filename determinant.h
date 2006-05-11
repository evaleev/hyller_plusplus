
#ifndef _hyller_determinant_h_
#define _hyller_determinant_h_

#include <vector>
#include "orbital.h"

namespace hyller {

  typedef enum {SpinAlpha = 1, SpinBeta = -1} SpinCase1;
  typedef enum {SpinAlphaAlpha = 1, SpinAlphaBeta = 0} SpinCase2;
  inline SpinCase1 spin1(SpinCase2 sc2) { return SpinAlpha; }
  inline SpinCase1 spin2(SpinCase2 sc2) { return sc2==SpinAlphaBeta ? SpinBeta : SpinAlpha; }
  inline bool samespin(SpinCase2 sc2) { return sc2 != SpinAlphaBeta; }

  /// 2-e Slater determinant
  struct SD {
    SD(const Orbital& i, const Orbital& j, SpinCase2 s);
    const Orbital* o1;
    const Orbital* o2;
    SpinCase2 spin;
  };

  bool operator==(const SD& A, const SD& B);

  /// Basis of 2-e Slater determinants supported by an 1-e basis
  class SDBasisSet {
  public:
    SDBasisSet(SpinCase2 spin,
	       const OrbitalBasisSet& obs);
    ~SDBasisSet() {}
    
    SpinCase2 spin() const { return spin_; }
    const OrbitalBasisSet& obs() const { return obs_; }
    int num_bf() const { return bfs_.size(); }
    const SD& bf(int i) const { return bfs_.at(i); }
    /// find this basis function and return its index. Throw, if not found
    int find(const SD& bf) const;

  private:
    SpinCase2 spin_;
    const OrbitalBasisSet& obs_;
    std::vector<SD> bfs_;
  };

};

#endif
