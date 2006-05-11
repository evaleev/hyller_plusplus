
#ifndef _hyller_misc_h_
#define _hyller_misc_h_

#include <vector>
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
};

#endif
