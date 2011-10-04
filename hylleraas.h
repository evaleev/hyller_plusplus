
#ifndef _hyller_hylleraas_h_
#define _hyller_hylleraas_h_

#include <vector>
#include <string>
#include "spin.h"

namespace hyller {

  template <typename TwoElectronBasisFunction>
  double value_at(const TwoElectronBasisFunction& bf, double r1, double r2, double r12);

  /// A single-exponent Hylleraas-type basis function:
  /// \f$ (r_1+r_2)^n (r_1-r_2)^l r_{12}^m e^{- \zeta (r_1 + r_2)} \f$
  struct HylleraasBasisFunction {
    /// Construct a dummy function
    HylleraasBasisFunction();
    /// Construct a Hylleraas basis function
    HylleraasBasisFunction(int n, int l, int m, double zeta);
    ~HylleraasBasisFunction() {}
    int n;
    int l;
    int m;
    double zeta;

    std::string to_string() const;
    std::string to_C_string() const;

    int nlm() const { return n+l+m; }

  };

  /// Comparison operator
  bool operator==(const HylleraasBasisFunction& A, const HylleraasBasisFunction& B);

  /// A Hylleraas-type basis set
  class HylleraasBasisSet {
  public:
    typedef HylleraasBasisFunction BasisFunction;
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
    typedef HylleraasBasisSet Basis;
    HylleraasWfn(const HylleraasBasisSet& bs, const std::vector<double>& coefs);
    ~HylleraasWfn() {}

    const HylleraasBasisSet& basis() const { return bs_; }
    const std::vector<double> coefs() const { return coefs_; }

    std::string to_string() const;
    std::string to_C_string() const;

  private:
    const HylleraasBasisSet& bs_;
    std::vector<double> coefs_;
  };

  ////

  /// A double-exponent Hylleraas-type basis function: (r1^i r2^j exp(-zeta1 r1) exp(-zeta2 r2) +- r1^j r2^i exp(-zeta2 r1) exp(-zeta1 r2)) r12^k, where '+' is for singlets, '-' is for triplet
  struct GenHylleraasBasisFunction {
    /// Construct a dummy function
    GenHylleraasBasisFunction();
    /// Construct a GenHylleraas basis function
    GenHylleraasBasisFunction(Spin2 s, int i, int j, int k, double zeta1, double zeta2);
    ~GenHylleraasBasisFunction() {}
    Spin2 spin;
    int i;
    int j;
    int k;
    // zeta is not necessary at the moment, but may be useful in the future
    double zeta1;
    double zeta2;

  };

  /// Comparison operator
  bool operator==(const GenHylleraasBasisFunction& A, const GenHylleraasBasisFunction& B);

  /// A GenHylleraas basis set
  class GenHylleraasBasisSet {
  public:
    /// Constructs a GenHylleraas-type basis set, given the parameters and restrictions
    GenHylleraasBasisSet(Spin2 spin, int ijk_max, int i_max, int j_max, int k_max, double zeta1, double zeta2);
    ~GenHylleraasBasisSet() {}

    /// Return spin
    Spin2 spin() const { return spin_; }
    /// maximum i+j+k value
    int ijk_max() const { return ijk_max_; }
    int i_max() const { return i_max_; }
    int j_max() const { return j_max_; }
    int k_max() const { return k_max_; }
    double zeta1() const { return zeta1_; }
    double zeta2() const { return zeta2_; }
    void set_zeta1(double zeta);
    void set_zeta2(double zeta);

    /// The number of basis functions
    int num_bf() const { return bfs_.size(); }
    /// ith basis function
    const GenHylleraasBasisFunction& bf(int i) const { return bfs_.at(i); }
    /// find this basis function and return its index. Throw, if not found
    int find(const GenHylleraasBasisFunction& bf) const;
    
  private:
    /** maximum number of basis functions to expect. More can be handled without a problem,
	but memory will be reallocated */
    const static unsigned int expected_num_bf = 10000;
    std::vector<GenHylleraasBasisFunction> bfs_;
    Spin2 spin_;
    int ijk_max_;
    int i_max_;
    int j_max_;
    int k_max_;
    double zeta1_;
    double zeta2_;
  };
  

};

#endif
