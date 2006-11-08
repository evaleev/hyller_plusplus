
#ifndef _hyller_focktimpl_h_
#define _hyller_focktimpl_h_

#include <fock.h>
#include <matrix.h>

namespace hyller {

template <typename F>
double**
DFFockOperator::J12(const BasisSet<F>& brabs,
		    const BasisSet<F>& ketbs) const
{
  const unsigned int nbra = brabs.nbf();
  const unsigned int nket = ketbs.nbf();
  const unsigned int nf = fbs_->nbf();
  
  // compute the density-fit one-particle J
  J_df();
  
  // compute the two-particle matrix
  double** J = block_matrix(nbra,nket);
  for(unsigned int b=0; b<nbra; ++b) {
    for(unsigned int k=0; k<nket; ++k) {
      double J_bk = 0.0;
      
      for(unsigned int i=0; i<nf; ++i) {
	const Orbital& f = fbs_->bf(i);
	const double J_bik_1 = hyller::gen_mult_oper(f.n,0,0,f.zeta,0.0,0.0,brabs.bf(b),ketbs.bf(k));
	const double J_bik_2 = hyller::gen_mult_oper(0,f.n,0,0.0,f.zeta,0.0,brabs.bf(b),ketbs.bf(k));
	J_bk += J_bik_1 + J_bik_2;
      }

      J[b][k] = J_bk;
    }
  }

  return J;
}

};

#endif
