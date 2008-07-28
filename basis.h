
#ifndef _hyller_basis_h_
#define _hyller_basis_h_

#include <utility>
extern "C" {
#include <libciomr/libciomr.h>
}
#include <except.h>
#include <basisfn.h>
#include <params.h>
#include <smartptr.h>

namespace hyller {

  /** Generic basis set of basis functions which are linear combinations ("contractions") of BF.
      By default all basis functions are independent of each other, i.e. there are no parameters
      which determine the composition of the basis. If a parameterized basis is needed, a derived class
      should be created.

      Note the public inheritance here. It helps derived BasisSets to define what happens when
      parameters are changed.
  */
  template <class PBF, class P = double> class BasisSet : public ParameterSet< Variable<P> >, public EnablePtrFromThis< BasisSet<PBF,P> > {
  public:
    typedef PBF PrimBF;
    typedef P ParamValue;
    typedef ContractedBasisFunction<PrimBF> ContrBF;
    typedef typename ContrBF::ContrTerm ContrTerm;
    typedef std::vector<PrimBF> PrimBFSet;
    typedef std::vector<ContrBF> ContrBFSet;
    typedef Variable<ParamValue> Parameter;
    typedef ParameterSet<Parameter> PSet;
    typedef RefParameterSet<Parameter> RefPSet;
    typedef PSet PSet_parent;
    typedef EnablePtrFromThis< BasisSet<PBF,P> > EnablePtr_parent;
    typedef BasisSet<PrimBF,ParamValue> TopBaseBasisSet;
  
    /// Constructs an empty, parameter-free basis set
    BasisSet() : PSet_parent(0), bfs_(0) {
      bfs_.reserve(expected_num_bf);
    }
    /** If possible, derived basis set should use this constructor to provide parameters upfront.
	If constructing parameters is too involved, use the default constructor and the protected params() to reset. */
    BasisSet(const Ptr<PSet>& pars) : bfs_(0) { bfs_.reserve(expected_num_bf); PSet_parent::params(pars->params()); }
    virtual ~BasisSet() {}

    Ptr<PSet> params() const {
      Ptr<BasisSet const> thisptr = EnablePtr_parent::Ptr_from_this();
      Ptr<PSet const> pset = static_pointer_cast<PSet const,BasisSet const>(thisptr);
      return const_pointer_cast<PSet,PSet const>(pset);
    }

    /// Adds a contracted basis function
    virtual void add(const ContrBF& bf) {
      // Adding contracted basis function is easy
      typename ContrBFSet::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
      typename ContrBFSet::const_iterator bend = bfs_.end();
      if (result != bend)
	return;

      bfs_.push_back(bf);

      // Now do the same for primitives, and also keep track of coefficients

      // loop over all terms in the contraction
      const unsigned int np = bf.n();
      ContrData cdata;
      for(unsigned int p=0; p<np; ++p) {
	const ContrTerm& t = bf.term(p);
	PrimBF p = t.first;
	typename PrimBFSet::const_iterator pend = prims_.end();
	typename PrimBFSet::const_iterator piter = std::find(prims_.begin(),prims_.end(),p);
	unsigned int pabs;
	if (piter == pend) {
	  pabs = prims_.size();
	  prims_.push_back(p);
	}
	else {
	  pabs = piter - prims_.begin();
	}
	cdata.push_back(std::make_pair(pabs,t.second));
      }
      coefs_.push_back(cdata);
    }

    /// Adds a unit-length contraction
    void add(const PrimBF& bf) {
      add(ContrBF(bf));
    }

    /// The number of contracted basis functions
    int nbf() const { return bfs_.size(); }
    /// The number of primitive basis functions
    int nprim() const { return prims_.size(); }
    /// ith contracted basis function
    const ContrBF& bf(int i) const { return bfs_.at(i); }
    /// ith primitive basis function
    const PrimBF& prim(int i) const { return prims_.at(i); }
    /// ith contracted basis function
    ContrBF& bf(int i) { return bfs_.at(i); }
    /// ith primitive basis function
    PrimBF& prim(int i) { return prims_.at(i); }
    /// find this contracted basis function and return its index. Throw, if not found
    int find(const ContrBF& bf) const
    {
      typename ContrBFSet::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
      if (result == bfs_.end())
	throw BasisFunctionNotFound();
      return result - bfs_.begin();
    }
    /// find this primitive basis function and return its index. Throw, if not found
    int find(const PrimBF& bf) const
    {
      typename PrimBFSet::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
      if (result == bfs_.end())
	throw BasisFunctionNotFound();
      return result - bfs_.begin();
    }

    /// Compute the primitive basis -> contracted basis coefficient matrix. The dimension of the matrix is nprims by nbf. Dense block matrix is returned.
    double** coefs() const {
      const int np = nprim();
      const int nb = nbf();

      // allocate the matrix
      double** C = block_matrix(np,nb);
#if 0
      new double*[np];
      C[0] = new double[np*nb];
      for(unsigned int i=1; i<np; i++)
	C[i] = C[i-1] + nb;
      memset((void*)C[0],0,np*nb*sizeof(double));
#endif

      // convert coefs_ to C
      std::vector<ContrData>::const_iterator end = coefs_.end();
      unsigned int bfindex = 0;
      for(std::vector<ContrData>::const_iterator bf=coefs_.begin(); bf!=end; ++bf, ++bfindex) {
	typename ContrData::const_iterator bend = bf->end();
	for(typename ContrData::const_iterator c=(*bf).begin(); c!=bend; ++c) {
	  C[c->first][bfindex] = c->second;
	}
      }

      return C;
    }

    std::string to_string() const {
      std::ostringstream oss;
      oss << "  -Basis Set:" << std::endl;
      const unsigned int np = nprim();
      for(unsigned int p=0; p<np; ++p) {
	const unsigned int np = nprim();
	oss << "    Prim #" << p << std::endl << prim(p).to_string();
      }
      return oss.str();
    }

  protected:
    /// Derived class may need to use a nontrivial parameter set. Need to make a deep copy here because the parameter set is owned by this (as reflected by this class being derived from PSet)
    void params(const Ptr<PSet>& p) {
      PSet_parent::params(p->params());
    }
    /// update contracted basis functions with the (presumably updated) parameters of the primitives
    void update_contr_with_prims() {
      typedef typename ContrBFSet::const_iterator citer;
      typedef typename ContrBFSet::iterator iter;
      citer cend = bfs_.end();
      unsigned int i=0;
      for(iter c=bfs_.begin(); c!=cend; ++c, ++i) {
	const unsigned int np = c->n();
	for(unsigned int p=0; p<np; ++p) {
	  const unsigned int pp = coefs_[i][p].first;
	  ContrTerm ct = c->term(p);
	  ct.first = prims_[pp];
	  c->term(p,ct);
	}
      }
    }

  private:

    /** maximum number of basis functions to expect. More can be handled without a problem,
	but memory will be reallocated */
    const static unsigned int expected_num_bf = 10000;
    // Primitive basis functions
    PrimBFSet prims_;
    // Contracted basis functions
    ContrBFSet bfs_;

    //
    // It is also convenient to store metadata to make construction of the transformation matrix easier 
    // 1) basis functions expansion coefficients in terms of primitives
    //

    // Contracted basis function is a vector of pairs of numbers: 1) index of the primitive 2) coefficient 
    typedef std::vector<std::pair<unsigned int, double> > ContrData;
    // Contraction coefficients
    std::vector<ContrData> coefs_;
  };

  template <typename BS> class Wavefunction;
  class GenSlaterHylleraasBasisSet;

  /// Describes traits of various BasisSet<T> instantiations
  template <typename BS>
    struct BasisSetTraits {
      typedef hyller::Wavefunction<BS> Wavefunction;
      /// No parameters by default
      static const int nparam = 0;
    };
  template <>
    struct BasisSetTraits< BasisSet<GenSlaterHylleraasBasisSet> > {
      typedef hyller::Wavefunction< BasisSet<GenSlaterHylleraasBasisSet> > Wavefunction;
      /// 2 parameters: alpha and gamma
      static const int nparam = 2;
    };

};

#endif
