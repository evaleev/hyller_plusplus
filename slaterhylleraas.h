
#ifndef _hyller_slaterhylleraas_h_
#define _hyller_slaterhylleraas_h_

#include "coords.h"
#include <vector>

namespace hyller {

  /** A Slater-Hylleraas-type basis function (exponentially correlated Hylleraas):
      \Phi_{nlm}(\zeta,\gamma) = e^{-\zeta s} e^{-\gamma u} s^n t^l u^m
  */
  struct SlaterHylleraasBasisFunction {
    /// Construct a dummy function
    SlaterHylleraasBasisFunction();
    /// Construct a Hylleraas basis function
    SlaterHylleraasBasisFunction(int n, int l, int m, double zeta, double gamma);
    ~SlaterHylleraasBasisFunction() {}
    int n;
    int l;
    int m;
    // zeta is not necessary at the moment, but may be useful in the future
    double zeta;
    double gamma;

    int nlm() const { return n+l+m; }
  };

  /// Comparison operator
  bool operator==(const SlaterHylleraasBasisFunction& A, const SlaterHylleraasBasisFunction& B);

  /// A Slater-Hylleraas-type basis set
  class SlaterHylleraasBasisSet {
  public:
    /// Constructs a Slater-Hylleraas-type basis set, given the parameters and restrictions
    SlaterHylleraasBasisSet(bool singlet, int nlm_max, int nlm_min, int n_max, int l_max, int m_max, double zeta, double gamma);
    ~SlaterHylleraasBasisSet() {}

    /// Return true is singlet
    bool singlet() const { return singlet_; }
    /// maximum n+l+m value
    int nlm_max() const { return nlm_max_; }
    int n_max() const { return n_max_; }
    int l_max() const { return l_max_; }
    int m_max() const { return m_max_; }
    /// minimum n+l+m value
    int nlm_min() const { return nlm_min_; }
    double zeta() const { return zeta_; }
    double gamma() const { return gamma_; }
    void set_zeta(double zeta);
    void set_gamma(double gamma);

    /// The number of basis functions
    int num_bf() const { return bfs_.size(); }
    /// ith basis function
    const SlaterHylleraasBasisFunction& bf(int i) const { return bfs_.at(i); }
    /// find this basis function and return its index. Throw, if not found
    int find(const SlaterHylleraasBasisFunction& bf) const;
    
  private:
    /** maximum number of basis functions to expect. More can be handled without a problem,
	but memory will be reallocated */
    const static unsigned int expected_num_bf = 10000;
    std::vector<SlaterHylleraasBasisFunction> bfs_;
    bool singlet_;
    int nlm_max_;
    int n_max_;
    int l_max_;
    int m_max_;
    int nlm_min_;
    double zeta_;
    double gamma_;
  };

  /**
     Is a wave vector in subspace of Hilbert space supported by the basis
   */
  class SlaterHylleraasWfn {
  public:
    SlaterHylleraasWfn(const SlaterHylleraasBasisSet& bs, const std::vector<double>& coefs);
    ~SlaterHylleraasWfn() {}

    const SlaterHylleraasBasisSet& basis() const { return bs_; }
    const std::vector<double> coefs() const { return coefs_; }
  private:
    const SlaterHylleraasBasisSet& bs_;
    std::vector<double> coefs_;
  };


  //////////

  /// 
  /** A generalized Slater-Hylleraas-type basis function (exponentially correlated Hylleraas):
      \Phi_{ijk}(\alpha,\beta,\gamma) = e^{-\alpha r_1 -\beta r_2 - \gamma r_{12}} r_1^i r_2^j r_{12}^k

      This function is not symmetry-adapted.
  */
  struct GenSlaterHylleraasBasisFunction {
    /// Regular coordinates are used
    static const int Coords = Coordinates::R1R2R12;
    /// The identity function
    static GenSlaterHylleraasBasisFunction Identity;

    /// Construct a dummy function
    GenSlaterHylleraasBasisFunction();
    /// Construct a Hylleraas basis function
    GenSlaterHylleraasBasisFunction(int i, int j, int k, double alpha, double beta, double gamma);
    ~GenSlaterHylleraasBasisFunction() {}
    int i;
    int j;
    int k;
    double alpha;
    double beta;
    double gamma;

    std::string to_string() const;
  };

  /// Comparison operator
  bool operator==(const GenSlaterHylleraasBasisFunction& A, const GenSlaterHylleraasBasisFunction& B);
  /// Return A*B
  GenSlaterHylleraasBasisFunction
    operator*(const GenSlaterHylleraasBasisFunction& A, const GenSlaterHylleraasBasisFunction& B);
  /// Returns the Jacobian * the function
  GenSlaterHylleraasBasisFunction
    apply_jacobian(const GenSlaterHylleraasBasisFunction& A);
  /// Expresses r_1^i * r_2^j * r_{12}^k operator as a GenSlaterHylleraasBasisFunction
  GenSlaterHylleraasBasisFunction
    gen_r1r2r12_oper(int i, int j, int k);

  class Orbital;
  /// a product of 2 1-particle functions can be represented as a GenSlaterHylleraasBasisFunction
  GenSlaterHylleraasBasisFunction
    operator^(const Orbital& f1,
	      const Orbital& f2);

  template <typename T> class Wavefunction;

  /// Basis set of generalized Slater-Hylleraas-type functions
  class GenSlaterHylleraasBasisSet {
  public:
    /// Type of the basis functions this set contains
    typedef GenSlaterHylleraasBasisFunction BasisFunction;
    /// Wave function type based on this basis
    typedef Wavefunction<GenSlaterHylleraasBasisSet> Wfn;

    /// Constructs a generalized Slater-Hylleraas-type basis set, given the parameters and restrictions
    GenSlaterHylleraasBasisSet(int ijk_max, int ijk_min, int ij_max, int k_max, double alpha, double beta, double gamma);
    ~GenSlaterHylleraasBasisSet() {}

    /// maximum i+j+k value
    int ijk_max() const { return ijk_max_; }
    /// maximum value of i or j
    int ij_max() const { return ij_max_; }
    int k_max() const { return k_max_; }
    /// minimum i+j+k value
    int ijk_min() const { return ijk_min_; }
    double alpha() const { return alpha_; }
    double beta() const { return beta_; }
    double gamma() const { return gamma_; }
    void set_alpha(double alpha);
    void set_beta(double beta);
    void set_gamma(double gamma);

    /// The number of basis functions
    int num_bf() const { return bfs_.size(); }
    /// The number of basis functions
    int nbf() const { return bfs_.size(); }
    /// ith basis function
    const GenSlaterHylleraasBasisFunction& bf(int i) const { return bfs_.at(i); }
    /// find this basis function and return its index. Throw, if not found
    int find(const GenSlaterHylleraasBasisFunction& bf) const;
    
  private:
    /** maximum number of basis functions to expect. More can be handled without a problem,
	but memory will be reallocated */
    const static unsigned int expected_num_bf = 10000;
    std::vector<GenSlaterHylleraasBasisFunction> bfs_;
    int ijk_max_;
    int ij_max_;
    int k_max_;
    int ijk_min_;
    double alpha_;
    double beta_;
    double gamma_;
  };

};

#endif
