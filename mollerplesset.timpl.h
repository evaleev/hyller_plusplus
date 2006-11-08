
#ifndef _hyller_mollerplessettimpl_h_
#define _hyller_mollerplessettimpl_h_

#include <fock.h>
#include <fock.timpl.h>
#include <mollerplesset.h>

namespace hyller {

  template <class B>
  void
  MollerPlessetSeries<B>::compute() {

    // the one-particle exchange will be fit using this basis
    fbs_ = Ptr<OrbitalBasisSet>(new OrbitalBasisSet(10,0,1.8));
    DFFockOperator Fock(hfwfn_,fbs_);

    // compute the exchange
    double** J12 = Fock.J12(*bs_,*bs_);

  }

};

#endif
