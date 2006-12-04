
#ifndef _hyller_mollerplesset_h_
#define _hyller_mollerplesset_h_

namespace hyller {

  /// Evaluated MP wavefunctions and energies
  template <typename B>
  class MollerPlessetSeries {
  public:
    typedef B Basis;

    /// Z is the nuclear charge, lambda scales the perturbation operator
    MollerPlessetSeries(const Ptr<OrbitalWfn>& hfwfn,
			const Ptr<Basis>& bs,
			double Z,
			double lambda) :
      hfwfn_(hfwfn), bs_(bs), Z_(Z), lambda_(lambda)
      {
      }
    ~MollerPlessetSeries() {
    }

    /// Sets the maximum order for the energy
    void nmax(unsigned int n) { nmax_ = n; }
    /// what do you think this does?
    void compute();
      
  private:
    Ptr<OrbitalWfn> hfwfn_;
    Ptr<Basis> bs_;
    Ptr<OrbitalBasisSet> fbs_;
    double Z_;
    double lambda_;
    unsigned int nmax_;

  };

};

#endif
