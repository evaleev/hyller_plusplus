
#ifndef _hyller_csf_h_
#define _hyller_csf_h_

#include <vector>
#include "determinant.h"

namespace hyller {

  typedef enum {SpinSinglet = 1, SpinTriplet = 3} Spin2;

  /// Configuration State Function for 2 electrons
  struct CSF {
    CSF(Spin2 s, const Orbital& o1, const Orbital& o2, const SDBasisSet& sdbasis);
    Spin2 spin;
    int nd;
    const SD* d[2];
    // coefficients defined for determinants ordered by the spins, not orbital indices (i.e. singlet = ab + ab, not ab - ba)
    double c[2];
  };

  bool operator==(const CSF& I, const CSF& J);

  /// Basis of CSFs
  class CSFBasisSet {
  public:
    CSFBasisSet(Spin2 spin, const SDBasisSet& sdbasis);

    Spin2 spin() const { return spin_; }
    const SDBasisSet& sdbasis() const { return sdbasis_; }
    int num_bf() const { return bfs_.size(); }
    const CSF& bf(int i) const { return bfs_.at(i); }
    /// find this basis function and return its index. Can throw BasisFunctionNotFound
    int find(const CSF& bf) const;

  private:
    std::vector<CSF> bfs_;
    Spin2 spin_;
    const SDBasisSet& sdbasis_;
  };

  class HylleraasBasisSet;
  class HylleraasBasisFunction;
  /// Given a CSF, return coefficients for its expansion in terms of hbs. Throw if hbs does not span I exactly
  std::vector<double>
    CSF_to_hylleraas(const CSF& csf, const HylleraasBasisSet& hbs);
  /// Given a Hylleraas function, return coefficients for its expansion in terms of cbs. Throw if bf has m=0 (cannot be expanded in terms of s-type orbitals)
  std::vector<double>
    hylleraas_to_CSF(const HylleraasBasisFunction& bf, const CSFBasisSet& cbs);

}

#endif
