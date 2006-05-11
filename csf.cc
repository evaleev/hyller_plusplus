
#include <cmath>
#include "misc.h"
#include "csf.h"
#include "hylleraas.h"
#include "matrix.h"
#include "except.h"

using namespace hyller;

CSF::CSF(Spin2 s, const Orbital& o1, const Orbital& o2, const SDBasisSet& sdbasis) :
  spin(s)
{
  if (s == SpinSinglet && samespin(sdbasis.spin()))
    throw std::runtime_error("CSF::CSF -- cannot construct CSF of requested spin with the given set of determinants");

  if (o1 == o2) {
    if (s != SpinSinglet)
      throw std::runtime_error("CSF::CSF -- cannot put 2 electrons into same spin-orbital");
    nd = 1;
    c[0] = 1.0;
    SD det(o1,o2,SpinAlphaBeta);
    const int i = sdbasis.find(det);
    d[0] = &(sdbasis.bf(i));
  }
  else {
    // only triplet possible
    if (samespin(sdbasis.spin())) {
      nd = 1;
      c[0] = 1.0;
      SD det(o1,o2,sdbasis.spin());
      const int i = sdbasis.find(det);
      d[0] = &(sdbasis.bf(i));
    }
    // can form either
    else {
      nd = 2;
      c[0] = M_SQRT1_2;
      c[1] = (spin == SpinSinglet) ? +M_SQRT1_2 : -M_SQRT1_2;
      SD det1(o1,o2,sdbasis.spin());
      const int i1 = sdbasis.find(det1);
      d[0] = &(sdbasis.bf(i1));
      SD det2(o2,o1,sdbasis.spin());
      const int i2 = sdbasis.find(det2);
      d[1] = &(sdbasis.bf(i2));
    }
  }
}

bool
hyller::operator==(const CSF& I, const CSF& J)
{
  bool c0 = I.spin==J.spin && I.nd==J.nd;
  if (c0) {
    bool c1 = I.d[0][0]==J.d[0][0] && I.c[0]==J.c[0];
    if (c1 && I.nd > 1)
      c1 = c1 && I.d[1][0]==J.d[1][0] && I.c[1]==J.c[1];
    bool c2 = false;
    if (!c1 && I.nd > 1) {
      c2 = I.d[0][0]==J.d[1][0] && I.c[0]==J.c[1];
      if (c2)
	c2 = c2 && I.d[1][0]==J.d[0][0] && I.c[1]==J.c[0];
    }
    if (c1 || c2)
      return true;
  }
  return false;
}

////

CSFBasisSet::CSFBasisSet(Spin2 spin, const SDBasisSet& sdbasis) :
  spin_(spin), sdbasis_(sdbasis)
{
  if (spin == SpinSinglet && samespin(sdbasis.spin()))
    throw std::runtime_error("CSFBasisSet::CSFBasisSet -- cannot construct CSFs of requested spin with the given set of determinants");
  
  const OrbitalBasisSet& obs = sdbasis.obs();
  const int norbs = obs.num_bf();

  for(int i=0; i<norbs; i++) {
    for(int j=0; j<=i; j++) {
      if (spin == SpinTriplet && i == j) continue;

      CSF csf(spin,obs.bf(i),obs.bf(j),sdbasis);
      bfs_.push_back(csf);
    }
  }
}

int
CSFBasisSet::find(const CSF& bf) const
{
  std::vector<CSF>::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
  if (result == bfs_.end())
    throw BasisFunctionNotFound();
  return result - bfs_.begin();
}

////

std::vector<double>
hyller::CSF_to_hylleraas(const CSF& C, const HylleraasBasisSet& hbs)
{
  const bool singlet = hbs.singlet();
  const Spin2 spin2 = singlet ? SpinSinglet : SpinTriplet;
  if (C.spin != spin2)
    throw std::runtime_error("CSF_to_hylleraas() -- spin cases do not match");

  const int ndets = C.nd;
  const Orbital& I = C.d[0]->o1[0];
  const Orbital& J = C.d[0]->o2[0];
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
	  double sign = (ndets == 1) ?  minus1_to_l(J.n-l) :  minus1_to_l(J.n-l)*C.c[0] +  minus1_to_l(I.n-m)*C.c[1];
	  X[ii] += oofourpi * sign * binomial(I.n,m) * binomial(J.n,l) * normConst(I) * normConst(J) / ( normConst(bf) * pow(2.0*zeta,bf.nlm()+3) * pow(2.0,I.n+J.n));
	}
	catch (BasisFunctionNotFound&) {}
      }
    }
  }

  return X;
}


std::vector<double>
hyller::hylleraas_to_CSF(const HylleraasBasisFunction& bf, const CSFBasisSet& cbs)
{
  if (bf.m != 0)
    throw std::runtime_error("hylleraas_to_CSF() -- Hylleraas functions with nonzero m cannot be expanded yet");

  const bool singlet = (bf.l%2 == 0);
  const Spin2 spin2 = singlet ? SpinSinglet : SpinTriplet;
  if (cbs.spin() != spin2)
    throw std::runtime_error("hylleraas_to_CSF() -- spin cases do not match");

  const SDBasisSet& sbs = cbs.sdbasis();
  const OrbitalBasisSet& obs = sbs.obs();
  // Note the factor of 2 -- my implementation of Hylleraas function has an extra factor of 2
  const double zeta = bf.zeta/2.0;
  const double oofourpi = 1.0 / (4.0 * M_PI);

  std::vector<double> X(cbs.num_bf(),0.0);

  for(int n=0; n<=bf.n; n++) {
    for(int l=0; l<=bf.l; l++) {
      const int x = n + l;
      const int y = bf.n + bf.l - x;

      if (x == y && !singlet) continue;

      Orbital o1(x,zeta);
      Orbital o2(y,zeta);
      try {
	CSF csf(spin2,o1,o2,sbs);
	const int ind = cbs.find(csf);
	const CSF& C = cbs.bf(ind);
	const Orbital& I = C.d[0]->o1[0];
	const Orbital& J = C.d[0]->o2[0];
	bool swapped_orbs = !(o1 == I);
	
	const int ndets = C.nd;
	double perm_pfac = (ndets == 1) ? 1.0 : M_SQRT1_2;
	double norm_pfac = pow(2.0*zeta,bf.nlm()+3) * normConst(bf) / (normConst(I) * normConst(J) * oofourpi);
	double sign = !swapped_orbs ? minus1_to_l(bf.l-l) : minus1_to_l(l);
	X[ind] +=  (!singlet ? 2.0 : 1.0) * perm_pfac * norm_pfac * sign * binomial(bf.n,n) * binomial(bf.l,l);
      }
      catch (BasisFunctionNotFound&) {}
    }
  }

  return X;
}
