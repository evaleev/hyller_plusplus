
#ifndef _hyller_hylleraas_h_
#define _hyller_hylleraas_h_

#include <vector>

namespace hyller {

  /// A Hylleraas-type basis function
  struct HylleraasBasisFunction {
    /// Construct a dummy function
    HylleraasBasisFunction();
    /// Construct a Hylleraas basis function
    HylleraasBasisFunction(int n, int l, int m, double zeta);
    ~HylleraasBasisFunction() {}
    int n;
    int l;
    int m;
    // zeta is not necessary at the moment, but may be useful in the future
    double zeta;

    int nlm() const { return n+l+m; }
  };

  /// Comparison operator
  bool operator==(const HylleraasBasisFunction& A, const HylleraasBasisFunction& B);

  /// A Hylleraas-type basis set
  class HylleraasBasisSet {
  public:
    /// Constructs a Hylleraas-type basis set, given the parameters and restrictions
    HylleraasBasisSet(bool singlet, int nlm_max, int nlm_min, int n_max, int l_max, int m_max, double zeta);
    ~HylleraasBasisSet() {}

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
    void set_zeta(double zeta);

    /// The number of basis functions
    int num_bf() const { return bfs_.size(); }
    /// ith basis function
    const HylleraasBasisFunction& bf(int i) const { return bfs_.at(i); }
    /// find this basis function and return its index. Throw, if not found
    int find(const HylleraasBasisFunction& bf) const;
    
  private:
    /** maximum number of basis functions to expect. More can be handled without a problem,
	but memory will be reallocated */
    const static unsigned int expected_num_bf = 10000;
    std::vector<HylleraasBasisFunction> bfs_;
    bool singlet_;
    int nlm_max_;
    int n_max_;
    int l_max_;
    int m_max_;
    int nlm_min_;
    double zeta_;
  };

  /**
     Is a wave vector in subspace of Hilbert space supported by the basis
   */
  class HylleraasWfn {
  public:
    HylleraasWfn(const HylleraasBasisSet& bs, const std::vector<double>& coefs);
    ~HylleraasWfn() {}

    const HylleraasBasisSet& basis() const { return bs_; }
    const std::vector<double> coefs() const { return coefs_; }
  private:
    const HylleraasBasisSet& bs_;
    std::vector<double> coefs_;
  };

};

#endif
