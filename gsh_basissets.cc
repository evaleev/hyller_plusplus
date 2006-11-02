
#include <gsh_basissets.h>
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
    if (!gsh::primbf_params_ok(params(),t.first))
      throw InvalidBasisFunctionParams();
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
    if (!symmgsh::primbf_params_ok(params(),t.first))
      throw InvalidBasisFunctionParams();
  }
  TopBaseBasisSet::add(bf);
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

SymmGSHBasisSet
hyller::operator^(const SymmGSHBasisSet& bs1,
		  const SymmGSHBasisSet& bs2)
{
  const double alpha = bs1.params().param(0).value() + bs2.params().param(0).value();
  const double gamma = bs1.params().param(1).value() + bs2.params().param(1).value();
  SymmGSHBasisSet bs(alpha,gamma);
}
