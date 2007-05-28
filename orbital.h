
#ifndef _hyller_orbital_h_
#define _hyller_orbital_h_

#include <vector>
#include <stdexcept>
#include <cmath>

namespace hyller {

  /// An orbital function of type r^n e^(-zeta*r). Spin not included (see class SD)
  struct Orbital {
  public:
    static Orbital Identity;
    /// Construct an orbital
    Orbital(int n, double zeta);
    ~Orbital() {}

    int n;
    double zeta;

    /// value at (r,theta,phi)
    double value_at(double r) const { return std::pow(r,(double)n) * std::exp(-zeta*r); }

  private:
    /// an identity function
    Orbital();
  };

  bool operator==(const Orbital& I, const Orbital& J);

  /** return I(1) * J(1) */
  Orbital operator*(const Orbital& I, const Orbital& J);

  /// Orbital basis set -- includes orbitals with same zeta and n in (n_min, n_max)
  class OrbitalBasisSet {
  public:
    OrbitalBasisSet(int n_max, int n_min, double zeta);
    ~OrbitalBasisSet() {}

    int n_max() const { return n_max_; }
    int n_min() const { return n_min_; }
    /// zeta
    double zeta() const { return zeta_; }
    void set_zeta(double zeta);

    /// number of basis functions
    int num_bf() const { return bfs_.size(); }
    /// number of basis functions
    int nbf() const { return bfs_.size(); }
    /// ith basis function
    const Orbital& bf(int i) const { return bfs_.at(i); }

  private:
    std::vector<Orbital> bfs_;
    int n_max_;
    int n_min_;
    double zeta_;
  };

  /**
     Is a wave vector in subspace of Hilbert space supported by the basis
   */
  class OrbitalWfn {
  public:
    OrbitalWfn(const OrbitalBasisSet& bs, const std::vector<double>& coefs, double E);
    ~OrbitalWfn() {}

    const OrbitalBasisSet& basis() const { return bs_; }
    const std::vector<double>& coefs() const { return coefs_; }
    double E() const { return E_; }
  private:
    const OrbitalBasisSet& bs_;
    std::vector<double> coefs_;
    double E_;
  };

  class HylleraasBasisSet;
  /// Given a product of Orbitals I and J, return coefficients for its expansion in terms of hbs. Throw if hbs does not span I exactly
  std::vector<double>
    orbital_to_hylleraas(const Orbital& I, const Orbital& J, const HylleraasBasisSet& hbs);

  template <typename BF> class ContractedBasisFunction;
  /// Convert an OrbitalWfn to a contracted basis function
  ContractedBasisFunction<Orbital> contract(const OrbitalWfn&);

};

#endif
