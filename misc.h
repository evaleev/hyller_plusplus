
#ifndef _hyller_misc_h_
#define _hyller_misc_h_

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "defines.h"

namespace hyller {
  
  template <bool Safe> class Precomputed {
  public:
    Precomputed(unsigned long maxindex) :
      maxindex_(maxindex), fac_(maxindex_+1), facdfac_(maxindex_+1), binomial_(maxindex_+1)
	{
	  fac_[0] = 1.0;
	  for(unsigned long i=1; i<=maxindex_; i++)
	    fac_[i] = i*fac_[i-1];
	  
	  for(unsigned long i=0; i<=maxindex_; i++) {
	    facdfac_[i].resize(maxindex_+1);
	    for(unsigned long j=0; j<=maxindex_; j++) {
	      facdfac_[i][j] = fac_[i] / fac_[j];
	    }
	  }

	  for(unsigned long i=0; i<=maxindex_; i++) {
	    binomial_[i].resize(i+1);
	    for(unsigned long j=0; j<=i; j++) {
	      binomial_[i][j] = facdfac_[i][j] / fac_[i-j];
	    }
	  }

	}

    ~Precomputed() {}
    
    /// Return i!
    double fac(unsigned int long i) const {
      return Safe ? fac_.at(i) : fac_[i];
    }
    /// Return i!/j!
    double facdfac(unsigned int long i, unsigned long j) const {
      return Safe ? facdfac_.at(i).at(j) : facdfac_[i][j];
    }
    /// Return i!/ ( (i-j)! j!)
    double binomial(unsigned long i, unsigned long j) const {
      return Safe ? binomial_.at(i).at(j) : binomial_[i][j];
    }
  private:
    unsigned long maxindex_;
    std::vector<double> fac_;         // factorial
    std::vector< std::vector<double> > facdfac_;     // factorial divided by factorial
    std::vector< std::vector<double> > binomial_;     // binomial coefficient
  };

#if CORRECT_RATHER_THAN_FAST
  extern Precomputed<true> cache;
#else
  extern Precomputed<false> cache;
#endif

  /// Computes i!
  inline double fac(unsigned long i) { return cache.fac(i); }
  /// Computes i!/j!
  inline double facdfac(unsigned long i, unsigned long j) { return cache.facdfac(i,j); }
  /// Computes i!/( j! (i-j)!) (i >= j must hold)
  inline double binomial(unsigned long i, unsigned long j) { return cache.binomial(i,j); }

  /// Computes (-1)^l
  inline double minus1_to_l(int l) { return l%2 ? -1.0 : 1.0; }

  /** This function computes Matrix^(-1/2),
      result stored in *Result (dim x dim-num_dep) (must not be allocated yet).
      dim is dimentionality, returns number of linearly dependent vectors */
  int sq_sqrt(double **Matrix, double ***Result, int dim);
  /// Invert Matrix
  double** sq_inverse(double **Matrix, int dim);

  /** (allocates and) computes U^t.V.U -- U is i by o, V is i by i */
  double** UtVU(double** U, double**V, int i, int o);
  /** (allocates and) computes U.V.U^t -- U is o by i, V is i by i */
  double** UVUt(double** U, double**V, int i, int o);
  /** (allocates and) computes U^t.v -- U is i by o, v is i long */
  double* Utv(double** U, double* v, int i, int o);

  /// Converts a vector<T> to T* (memory is allocated using new)
  template <typename T> T* to_C_array(const std::vector<T>& v) {
    const unsigned int rank = v.size();
    if (rank) {
      T* d = new T[rank];
      std::copy(v.begin(),v.end(),d);
      return d;
    }
    return (T*)0;
  }
  /// Converts a vector< vector<T> > to T** (memory is allocated as a block (see block_matrix()) using new)
  template <typename T> T** to_C_array(const std::vector< std::vector<T> >& v) {
    const unsigned int nrow = v.size();
    if (nrow) {
      const unsigned int ncol = v[0].size();
      if (ncol) {
	// check that all rows have the same number of columns
	for(unsigned int r=0; r<nrow; ++r) {
	  if (v[r].size() != ncol)
	    throw std::runtime_error("to_C_array() -- all rows must have the same number of columns");
	}
	// allocate a block matrix
	T** d = new T*[nrow];
	d[0] = new T[nrow*ncol];
	for(unsigned int r=1; r<nrow; ++r)
	  d[r] = d[r-1] + ncol;
	for(unsigned int r=0; r<nrow; ++r) {
	  std::copy(v[r].begin(),v[r].end(),d[r]);
	}
	return d;
      }
      return (T**)0;
    }
    return (T**)0;
  }

  /// Norm of a vector<T>
  template <typename T>
  T norm(const std::vector<T>& v) {
    T norm2 = 0.0;
    typedef typename std::vector<T>::const_iterator citer;
    citer vend = v.end();
    for(citer i=v.begin(); i!=vend; ++i) {
      norm2 += (*i)*(*i);
    }
    T norm = sqrt(norm2);
    return norm;
  }

};

#endif
