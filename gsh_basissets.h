
#ifndef _hyller_gshbasissets_h_
#define _hyller_gshbasissets_h_

#include <basis.h>
#include <slaterhylleraas.h>

namespace hyller {

  /** A simple BasisSet of (contracted) GenSlaterHylleraasBasisFunctions which allows functions to be added.
      This basis has 3 parameters: alpha, beta, gamma. Any basis function Phi(n,l,m; a,b,c) added to this set must
      have a=alpha, b=beta, c=gamma or a=beta, b=alpha, c=gamma or a=alpha, b=beta, c=0.0 or a=beta, b=alpha, c=0.0.
      The basis functions with c=0.0 are just Hylleraas basis functions and are an important special case.
      When mutable_gamma is true, only nonzero parameters c are be varied.

      Meant to serve as a base to more restrictive basis sets. */
  class GSHBasisSet : public BasisSet<GenSlaterHylleraasBasisFunction,double> {
  public:
    typedef BasisSet<GenSlaterHylleraasBasisFunction,double> BaseBasisSet;
#include <basebasisset_typedefs.h>

    GSHBasisSet(double alpha, double beta, double gamma, bool mutable_alpha = true, bool mutable_beta = true, bool mutable_gamma = true);
    virtual ~GSHBasisSet();

    // overloads BasisSet::add. Applies parameter constraints described above
    void add(const ContrBF& bf);

  protected:
    // overloads PSet::param()
    void param(unsigned int i, const Parameter& p);

  };

  /** Specialization of GSHBasisSet which assumes that alpha==beta.
      This basis has 2 parameters: alpha, gamma. Any basis function Phi(n,l,m; a,b,c) added to this set must
      have a=b=alpha, c=gamma. */
  class SymmGSHBasisSet : public GSHBasisSet {
  public:
    typedef GSHBasisSet BaseBasisSet;
#include <basebasisset_typedefs.h>

    SymmGSHBasisSet(double alpha, double gamma, bool mutable_alpha = true, bool mutable_gamma = true);
    ~SymmGSHBasisSet();

    // overloads GSHBasisSet::add. Applies parameter constraints described above
    void add(const ContrBF& bf);

  protected:
    // overloads PSet::param()
    void param(unsigned int i, const Parameter& p);

  };

  /// Outer product of 2 basis sets is a set of every pairwise product
  SymmGSHBasisSet operator^(const SymmGSHBasisSet& bs1,
			    const SymmGSHBasisSet& bs2);

};

#endif
