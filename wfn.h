
#ifndef _hyller_wfn_h_
#define _hyller_wfn_h_

namespace hyller {

  /**
     Is a wave vector in subspace of Hilbert space supported by the basis
   */
  template <typename B>
  class Wavefunction {
  public:
    /// Basis set type
    typedef B Basis;

    Wavefunction(const Basis& bs, const std::vector<double>& coefs) : bs_(bs), coefs_(coefs)
      {
	if (bs_.num_bf() != coefs_.size())
	  throw std::runtime_error("Wavefunction::Wavefunction -- number of functions in basis does not match the size of the coefficient vector.");
      }
    ~Wavefunction() {}

    const Basis& basis() const { return bs_; }
    const std::vector<double> coefs() const { return coefs_; }
  private:
    const Basis& bs_;
    std::vector<double> coefs_;
  };

};

#endif
