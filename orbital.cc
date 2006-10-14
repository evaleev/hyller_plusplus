
#include <cmath>
#include "misc.h"
#include "orbital.h"
#include "hylleraas.h"
#include "matrix.h"
#include "except.h"
#include <basisfn.h>

using namespace hyller;

Orbital::Orbital() : n(0), zeta(0.0) {}

Orbital::Orbital(int nn, double zz) :
  n(nn), zeta(zz)
{
  if (zeta <= 0.0) {
    throw std::runtime_error("Orbital::Orbital -- zeta must be positive");
  }
}

bool
hyller::operator==(const Orbital& I, const Orbital& J)
{
  return I.n==J.n && I.zeta==J.zeta;
}

////

OrbitalBasisSet::OrbitalBasisSet(int n_max, int n_min, double zeta) :
  n_max_(n_max), n_min_(n_min), zeta_(zeta)
{
  for(int n=n_min_; n<=n_max_; n++) {
    Orbital bf(n,zeta_);
    bfs_.push_back(bf);
  }
}

void
OrbitalBasisSet::set_zeta(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("OrbitalBasisSet::set_zeta -- zeta must be positive");
  zeta_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].zeta = zeta_;
}

////

OrbitalWfn::OrbitalWfn(const OrbitalBasisSet& bs, const std::vector<double>& coefs) :
  bs_(bs), coefs_(coefs)
{
  if (bs_.num_bf() != coefs_.size())
    throw std::runtime_error("OrbitalWfn::OrbitalWfn -- size of basis set and coefficient vector do not match");
}

////

std::vector<double>
hyller::orbital_to_hylleraas(const Orbital& I, const Orbital& J, const HylleraasBasisSet& hbs)
{
  const bool singlet = hbs.singlet();
  if (I.zeta != J.zeta)
    throw std::runtime_error("orbital_to_hylleraas -- I and J have different orbital exponents");
  const double zeta = I.zeta;
  const double oofourpi = 1.0 / (4.0 * M_PI);

  std::vector<double> X(hbs.num_bf(),0.0);

  for(int m=0; m<=I.n; m++) {
    for(int l=0; l<=J.n; l++) {
      const int y = I.n+J.n-m-l;
      if ( (y%2 == 0) == singlet) {
	// Note the factor of 2 -- my implementation of Hylleraas function has an extra factor of 2
	HylleraasBasisFunction bf(m+l,y,0,2.0*zeta);
	try {
	  const int ii = hbs.find(bf);
	  X[ii] += oofourpi * minus1_to_l(J.n-l) * binomial(I.n,m) * binomial(J.n,l) * normConst(I) * normConst(J) / ( normConst(bf) * pow(2.0*zeta,bf.nlm()+3) * pow(2.0,I.n+J.n));
	}
	catch (BasisFunctionNotFound&) {}
      }
    }
  }

  return X;
}


ContractedBasisFunction<Orbital>
hyller::contract(const OrbitalWfn& wfn)
{
  const OrbitalBasisSet& basis = wfn.basis();
  const std::vector<double>& coefs = wfn.coefs();
  const unsigned int nbf = basis.nbf();
  ContractedBasisFunction<Orbital> result;
  for(int f=0; f<nbf; ++f) {
    const Orbital& O = basis.bf(f);
    // The product is in terms of non-normalized functions
    result.add(O,coefs[f]*normConst(O));
  }
  return result;
}
