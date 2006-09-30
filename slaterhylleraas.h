
#ifndef _hyller_slaterhylleraas_h_
#define _hyller_slaterhylleraas_h_

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

};

#endif
