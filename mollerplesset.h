
#ifndef _hyller_mollerplesset_h_
#define _hyller_mollerplesset_h_

namespace hyller {

  /// Evaluated MP wavefunctions and energies
  template <typename B>
  class MollerPlessetSeries {
  public:
    typedef B Basis;

    MollerPlessetSeries(const Ptr<OrbitalWfn>& hfwfn,
			const Ptr<Basis>& bs) :
      hfwfn_(hfwfn), bs_(bs)
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
    unsigned int nmax_;

  };

};

#endif
