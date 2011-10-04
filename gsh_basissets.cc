
#include <gsh_basissets.h>
#include <orbital.h>
#include <except.h>

using namespace hyller;

GSHBasisSet::GSHBasisSet(double alpha, double beta, double gamma,
			 bool mutable_alpha, bool mutable_beta, bool mutable_gamma) :
  GSHBasisSet::BaseBasisSet()
{
  Ptr<PSet> params(new PSet(3));
  params->param(0,Parameter(alpha,mutable_alpha));
  params->param(1,Parameter(beta,mutable_beta));
  params->param(2,Parameter(gamma,mutable_gamma));
  BaseBasisSet::params(params);
}

GSHBasisSet::~GSHBasisSet()
{
}

namespace gsh {
namespace {
  typedef GSHBasisSet::PSet PSet;
  typedef GSHBasisSet::PrimBF PrimBF;
  inline bool primbf_params_ok(const Ptr<PSet>& params, const PrimBF& bf) {
    const double alpha = params->param(0).value();
    const double beta = params->param(1).value();
    const double gamma = params->param(2).value();

    bool result = ( (bf.gamma == gamma || bf.gamma == 0.0) && ( (bf.alpha == alpha && bf.beta == beta) ||
								(bf.alpha == beta && bf.beta == alpha) ));
    return result;
  }
  inline void set_primbf_params(PrimBF& bf, double alpha, double beta, double gamma) {
    bf.alpha = alpha;
    bf.beta  = beta;
    // do not change gamma is it's zero
    if (bf.gamma != 0.0)
      bf.gamma = gamma;
  }
}
}

void
GSHBasisSet::add(const ContrBF& bf)
{
  // loop over all primitives and make sure alpha, beta, and gamma are correct
  const unsigned int np = bf.n();
  for(unsigned int p=0; p<np; ++p) {
    const ContrTerm& t = bf.term(p);
#define ALLOW_DIFFERENT_EXPONENTS 1
#if !ALLOW_DIFFERENT_EXPONENTS
    if (!gsh::primbf_params_ok(params(),t.first))
      throw InvalidBasisFunctionParams();
#endif
  }
  TopBaseBasisSet::add(bf);
}

void
GSHBasisSet::param(unsigned int i, const Parameter& p)
{
  // update the parameter
  PSet::param(i,p);
  const double alpha = params()->param(0).value();
  const double beta = params()->param(1).value();
  const double gamma = params()->param(2).value();

  // loop over all primitives and adjust their alpha, beta, and gamma
  const unsigned int np = nprim();
  for(unsigned int p=0; p<np; ++p) {
    gsh::set_primbf_params(prim(p),alpha,beta,gamma);
  }

  // update contractions with the new primitives
  update_contr_with_prims();
}


////

SymmGSHBasisSet::SymmGSHBasisSet(double alpha, double gamma,
			 bool mutable_alpha, bool mutable_gamma) :
  SymmGSHBasisSet::BaseBasisSet(alpha,alpha,gamma,mutable_alpha,mutable_alpha,mutable_gamma)
{
  Ptr<PSet> params(new PSet(2));
  params->param(0,Parameter(alpha,mutable_alpha));
  params->param(1,Parameter(gamma,mutable_gamma));
  BaseBasisSet::params(params);
}

SymmGSHBasisSet::~SymmGSHBasisSet()
{
}

namespace {
namespace symmgsh{
  typedef SymmGSHBasisSet::PSet PSet;
  typedef SymmGSHBasisSet::PrimBF PrimBF;
  inline bool primbf_params_ok(const Ptr<PSet>& params, const PrimBF& bf) {
    const double alpha = params->param(0).value();
    const double gamma = params->param(1).value();

    bool result = ((bf.gamma == gamma || bf.gamma == 0.0) && (bf.alpha == alpha && bf.beta == alpha));
    return result;
  }
}
}

void
SymmGSHBasisSet::add(const ContrBF& bf)
{
  // loop over all primitives and make sure alpha, beta, and gamma are correct
  const unsigned int np = bf.n();
  for(unsigned int p=0; p<np; ++p) {
    const ContrTerm& t = bf.term(p);
#define ALLOW_DIFFERENT_EXPONENTS 1
#if !ALLOW_DIFFERENT_EXPONENTS
    if (!symmgsh::primbf_params_ok(params(),t.first))
      throw InvalidBasisFunctionParams();
#endif
  }
  TopBaseBasisSet::add(bf);
}

void
SymmGSHBasisSet::add(int ijk_max, int ijk_min, int ij_max, int k_max)
{
  const double alpha = params()->param(0).value();
  const double gamma = params()->param(1).value();

  /// iterate over basis functions, given the constraints
  for(int i=0;i<=ijk_max;i++)
    for(int j=0;j<=ijk_max-i;j++)
      for(int k=0;k<=ijk_max-i-j;k++) {
	if (i > ij_max || j > ij_max || k > k_max || (i+j+k < ijk_min))
	  continue;
        
	PrimBF pbf(i,j,k,alpha,alpha,gamma);
	ContrBF cbf(pbf);
	add(cbf);
      }

}

void
SymmGSHBasisSet::param(unsigned int i, const Parameter& p)
{
  // update the parameter
  PSet::param(i,p);
  const double alpha = params()->param(0).value();
  const double gamma = params()->param(1).value();

  // loop over all primitives and adjust their alpha and gamma
  const unsigned int np = nprim();
  for(unsigned int p=0; p<np; ++p) {
    gsh::set_primbf_params(prim(p),alpha,alpha,gamma);
  }

  // update contractions with the new primitives
  update_contr_with_prims();
}


////

Ptr<SymmGSHBasisSet>
hyller::operator^(const SymmGSHBasisSet& bs1,
		  const SymmGSHBasisSet& bs2)
{
  typedef SymmGSHBasisSet Basis;
  const double alpha = bs1.params()->param(0).value() + bs2.params()->param(0).value();
  const double gamma = bs1.params()->param(1).value() + bs2.params()->param(1).value();

  Ptr<SymmGSHBasisSet> bs(new SymmGSHBasisSet(alpha,gamma,false,false));
  const unsigned int nbf1 = bs1.nbf();
  const unsigned int nbf2 = bs2.nbf();
  for(unsigned int bf1=0; bf1<nbf1; ++bf1) {
    const Basis::ContrBF& f1 = bs1.bf(bf1);
    for(unsigned int bf2=0; bf2<nbf2; ++bf2) {
      const Basis::ContrBF f2 = bs2.bf(bf2);
      bs->add(f1 * f2);
    }
  }

  return bs;
}

Ptr<GSHBasisSet>
hyller::operator^(const OrbitalBasisSet& bs1,
		  const OrbitalBasisSet& bs2)
{
  const double alpha = bs1.zeta();
  const double beta = bs2.zeta();
  const double gamma = 0.0;
  Ptr<GSHBasisSet> bs(new GSHBasisSet(alpha,beta,gamma,false,false,false));

  const unsigned int norbs1 = bs1.nbf();
  const unsigned int norbs2 = bs2.nbf();
  for(unsigned int o1=0; o1<norbs1; ++o1) {
    const Orbital& orb1 = bs1.bf(o1);
    for(unsigned int o2=0; o2<norbs2; ++o2) {
      const Orbital& orb2 = bs2.bf(o2);

      typedef GSHBasisSet::ContrBF ContrBF;
      bs->add(ContrBF(orb1^orb2));

    }
  }

  return bs;
}
