
#ifndef _hyller_basis_h_
#define _hyller_basis_h_

#include <utility>
#include <except.h>
#include <basisfn.h>

namespace hyller {

  /// Generic basis set of basis functions which are linear combinations ("contractions") of BF
  template <class T> class BasisSet {
  public:
    typedef T BF;
    typedef BF PrimBF;
    typedef std::pair<PrimBF,double> ContrTerm;
    typedef std::vector<ContrTerm> ContrBF;
    typedef std::vector<PrimBF> PrimBFSet;
    typedef std::vector<ContrBF> ContrBFSet;

    /// Constructs an empty basis set, given the parameters and restrictions
    BasisSet() : bfs_(0) { bfs_.reserve(expected_num_bf); }
    ~BasisSet() {}

    /// Adds a contracted basis function
    void add(const ContrBF& bf) {
      // Adding contracted basis function is easy
      typename ContrBFSet::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
      typename ContrBFSet::const_iterator bend = bfs_.end();
      if (result != bend)
	return;

      bfs_.push_back(bf);

      // Now do the same for primitives, and also keep track of coefficients
      typename ContrBF::const_iterator end = bf.end();
      ContrData cdata(bf.size());
      unsigned int loc_pindex = 0;
      for(typename ContrBF::const_iterator t=bf.begin(); t!=end; ++t, ++loc_pindex) {
	PrimBF p = t->first;
	typename PrimBFSet::const_iterator pend = prims_.end();
	typename PrimBFSet::const_iterator piter = std::find(prims_.begin(),prims_.end(),p);
	unsigned int abs_pindex;
	if (piter == pend) {
	  abs_pindex = prims_.size();
	  prims_.push_back(p);
	}
	else {
	  abs_pindex = piter - prims_.begin();
	}
	cdata[loc_pindex].first = abs_pindex;
	cdata[loc_pindex].second = t->second;
      }
      coefs_.push_back(cdata);
    }

    /// Adds a unit-length contraction
    void add(const PrimBF& bf) {
      add(ContrBF(1,std::make_pair(bf,1.0)));
    }

    /// The number of contracted basis functions
    int nbf() const { return bfs_.size(); }
    /// The number of primitive basis functions
    int nprim() const { return prims_.size(); }
    /// ith contracted basis function
    const ContrBF& bf(int i) const { return bfs_.at(i); }
    /// ith primitive basis function
    const PrimBF& prim(int i) const { return prims_.at(i); }
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
      double** C = new double*[np];
      C[0] = new double[np*nb];
      for(unsigned int i=1; i<np; i++)
	C[i] = C[i-1] + nb;
      memset((void*)C[0],0,np*nb*sizeof(double));

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
    
  private:


    /** maximum number of basis functions to expect. More can be handled without a problem,
	but memory will be reallocated */
    const static unsigned int expected_num_bf = 10000;
    // Primitive basis fucntions
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

};

#endif
