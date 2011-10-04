
#ifndef _hyller_wfn_h_
#define _hyller_wfn_h_

#include <iomanip>

// These must be global -- required by Psi3
extern "C" FILE *outfile;

namespace hyller {

  /**
     Is a wave vector in subspace of Hilbert space supported by the basis
   */
  template <typename B>
  class Wavefunction {
  public:
    /// Basis set type
    typedef B Basis;

    Wavefunction(const Basis& bs, const std::vector<double>& coefs,
                 double energy = 0.0) : bs_(bs), coefs_(coefs), energy_(energy)
      {
      if (bs_.nbf() != coefs_.size())
        throw std::runtime_error("Wavefunction::Wavefunction -- number of functions in basis does not match the size of the coefficient vector.");
      }
    ~Wavefunction() {}

    const Basis& basis() const { return bs_; }
    const std::vector<double>& coefs() const { return coefs_; }
    double energy() const { return energy_; }

    std::string to_string() const {
      const size_t nbf = bs_.nbf();
      std::ostringstream oss;
      if (nbf) {
        oss << std::setprecision(15) << coefs_[0] << " * "
            << bs_.bf(0).to_string();
        for(int f=1; f<nbf; ++f) {
          oss << " ";
          if (coefs_[f] > 0.0)
            oss << " + ";
          oss << std::setprecision(15) << coefs_[f] << " * "
              << bs_.bf(f).to_string() << " " << std::endl;
        }
      }
      return oss.str();
    }

    std::string to_C_string(bool skip_exponential = false) const {
      const size_t nbf = bs_.nbf();
      std::ostringstream oss;
      if (nbf) {
        oss << std::setprecision(15) << coefs_[0] << " * "
            << bs_.bf(0).to_C_string(skip_exponential);
        for(int f=1; f<nbf; ++f) {
          oss << " ";
          if (coefs_[f] > 0.0)
            oss << " + ";
          oss << std::setprecision(15) << coefs_[f] << " * "
              << bs_.bf(f).to_C_string(skip_exponential) << " " << std::endl;
        }
      }
      return oss.str();
    }

  private:
    const Basis& bs_;
    std::vector<double> coefs_;
    double energy_;
  };

  template <typename TwoElectronWavefunction>
  void plot_angular_profile(const TwoElectronWavefunction& wfn, int npoints, double radius, bool psinorm)
  {
    FILE *fpcusp;

    /*---------------------------------------------------------------------------
      Plotting wave function vs. r12 on a circle
     ---------------------------------------------------------------------------*/
    fprintf(outfile,"\tThe wavefunction near the cusp is plotted\n");
    fprintf(outfile,"\tRadius = %lf\n",radius);

    typedef typename TwoElectronWavefunction::Basis Basis;
    typedef typename Basis::BasisFunction BasisFunction;

    const Basis& basis = wfn.basis();
    const std::vector<double>& coeff = wfn.coefs();
    const int nbf = basis.num_bf();
    double psi00 = 0.0;
    for(int i=0;i<nbf;i++) {
      const BasisFunction& bf = basis.bf(i);
      psi00 += value_at(bf, radius, radius, 0.0) * coeff[i];
    }
    fprintf(outfile,"\tPsi(r12=0) = %11.10e\n",psi00);
    fpcusp = fopen("cusp_plot.out","w");
    //fprintf(fpcusp,"#%d %d %d %d %d %d %lf %25.15lf %s\n",npoints,basis.nlm_max(),basis.nlm_min(),
    // basis.n_max(),basis.l_max(),basis.m_max(),
   // radius,psi00,(psinorm ? "true" : "false"));
    if (psinorm)
      fprintf(fpcusp,"0.000000\t%25.15lf\n",1.0);
    else
      fprintf(fpcusp,"0.000000\t%25.15lf\n",psi00);
    for(int i=1;i<npoints;i++) {
      const double phi = M_PI * i / (npoints-1);
      const double r12 = 2.0 * sin(phi/2.0) * radius;
      double psi = 0.0;
      for(int j=0;j<nbf;j++) {
        const BasisFunction& bf = basis.bf(j);
        psi += value_at(bf, radius, radius, r12) * coeff[j];
      }
      if (psinorm)
        psi /= psi00;
      fprintf(fpcusp,"%lf\t%25.15lf\n",phi,psi);
    }
    fclose(fpcusp);
  }

};

#endif
