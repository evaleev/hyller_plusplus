
#ifndef _hyller_determinant_h_
#define _hyller_determinant_h_

#include <vector>
#include "spin.h"
#include "orbital.h"

namespace hyller {

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
