
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
  const unsigned int npbra = brabs.nprim();
  const unsigned int npket = ketbs.nprim();
  const unsigned int nf = fbs_->nbf();
  
  // compute the density-fit one-particle J
  J_df();
  const std::vector<double>& Jcoefs = J_df_->coefs();
  
  // compute the two-particle matrix in the primitive basis
  double** Jpp = block_matrix(npbra,npket);
  for(unsigned int b=0; b<npbra; ++b) {
    for(unsigned int k=0; k<npket; ++k) {
      double J_bk = 0.0;
      
      for(unsigned int i=0; i<nf; ++i) {
	const Orbital& f = fbs_->bf(i);
	const double oofourpi = 1.0 / (4.0 * M_PI);
	const double norm_pfac = normConst(f) * oofourpi;
	const double J_bik_1 = hyller::gen_mult_oper(f.n,0,0,f.zeta,0.0,0.0,brabs.prim(b),ketbs.prim(k));
	const double J_bik_2 = hyller::gen_mult_oper(0,f.n,0,0.0,f.zeta,0.0,brabs.prim(b),ketbs.prim(k));
	J_bk += (J_bik_1 + J_bik_2) * norm_pfac * Jcoefs[i] * (8.0 * M_PI * M_PI);
      }

      Jpp[b][k] = J_bk;
    }
  }

  // Convert to a contracted basis
  double** p2c_bra = brabs.coefs();
  double** p2c_ket = ketbs.coefs();
  double** J = UtVX(p2c_bra,npbra,nbra,Jpp,p2c_ket,npket,nket);
  free_block(p2c_bra);
  free_block(p2c_ket);
  free_block(Jpp);

  return J;
}

};

#endif
