
#ifndef _hyller_fock_h_
#define _hyller_fock_h_

#include <smartptr.h>
#include <orbital.h>
#include <basis.h>

namespace hyller {

  /// Density-fit Fock operator
  class DFFockOperator {
  public:
    /// Defined by the HF wave function and the fitting basis
    DFFockOperator(const Ptr<OrbitalWfn>& hfwfn,
		   const Ptr<OrbitalBasisSet>& fbs);
    ~DFFockOperator();

    /// Compute density fit of the Coulomb operator
    const Ptr<OrbitalWfn>& J_df() const;
    /// Compute density fit of the exchange operator
    const Ptr<OrbitalWfn>& K_df() const;
    /// Compute the two-electron Coulomb operator between the two basis sets. The result is allocated.
    template <typename F> double** J12(const BasisSet<F>& brabs,
				       const BasisSet<F>& ketbs) const;

  private:
    Ptr<OrbitalWfn> hfwfn_;
    Ptr<OrbitalBasisSet> fbs_;
    mutable Ptr<OrbitalWfn> J_df_;
    mutable Ptr<OrbitalWfn> K_df_;

  };

};

#endif
