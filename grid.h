
#ifndef _hyller_grid_h_
#define _hyller_grid_h_

#include <utility>
#include <chemistry/qc/dft/integrator.h>

// These must be global -- required by Psi3
extern "C" FILE *outfile;

namespace hyller {

  /**
     Integrates a 2-electron wavefunction. Assumed that this is an l=0 (S) state.
   */
  template <typename TwoElectronWavefunction>
  class NIntegrate {
  public:
    /// Wavefunction type
    typedef TwoElectronWavefunction Wfn;
    /// Basis set type
    typedef typename TwoElectronWavefunction::Basis Basis;
    typedef typename Basis::BasisFunction BasisFunction;

    NIntegrate(const Wfn& wfn,
               int npts_rad = 50,
               int npts_ang = 194) : wfn_(wfn),
               npts_rad_(npts_rad), npts_ang_(npts_ang),
               integrator_rad_(new sc::EulerMaclaurinRadialIntegrator(npts_rad)),
               integrator_ang_(new sc::LebedevLaikovIntegrator(npts_ang))
      {
      }
    ~NIntegrate() {}

    const Wfn& wfn() const { return wfn_; }
    std::pair<int,int> npts() const { return std::make_pair(npts_rad_,npts_ang_); }

    // compute the integral
    void compute(double& norm) {

      // prepare
      const Basis& basis = wfn_.basis();
      const std::vector<double>& coeff = wfn_.coefs();
      const int nbf = basis.nbf();

      // integrate in a 10 bohr sphere
      const double radius = 10.0;
      norm = 0.0;

      // electron 1 (this is an s-state so I can put it on the z axis and multiply the result by 4*pi)
      for(int ir1=0; ir1<npts_rad_; ++ir1) {
        double weight1_rad;
        const double r1 = integrator_rad_->radial_value(ir1, npts_rad_, radius,
                                                        weight1_rad);
        const int npts_ang = integrator_ang_->num_angular_points(r1, ir1);
        for(int ia1=0; ia1<npts_ang; ++ia1) {
          sc::SCVector3 xyz1;
          const double weight1_ang = integrator_ang_->angular_point_cartesian(ia1, r1,
                                                                              xyz1);

          const double weight1 = weight1_rad * weight1_ang;

          // electron 2
          for(int ir2=0; ir2<npts_rad_; ++ir2) {
            double weight2_rad;
            const double r2 = integrator_rad_->radial_value(ir2, npts_rad_, radius,
                                                            weight2_rad);
            const int npts_ang = integrator_ang_->num_angular_points(r2, ir2);
            for(int ia2=0; ia2<npts_ang; ++ia2) {
              sc::SCVector3 xyz2;
              const double weight2_ang = integrator_ang_->angular_point_cartesian(ia2, r2,
                                                                                  xyz2);

              const double weight2 = weight2_rad * weight2_ang;

              sc::SCVector3 r12_vec(xyz1);
              r12_vec -= xyz2;
              const double r12 = r12_vec.norm();

              double psi = 0.0;
              for(int j=0;j<nbf;j++) {
                const BasisFunction& bf = basis.bf(j);
                psi += value_at(bf, r1, r2, r12) * coeff[j];
              }

              norm += weight1 * weight2 * psi * psi;

            }
          } // electron 2
        }
      } // electron 3

    }

  private:
    const Wfn& wfn_;
    int npts_rad_, npts_ang_;
    sc::Ref<sc::EulerMaclaurinRadialIntegrator> integrator_rad_;
    sc::Ref<sc::LebedevLaikovIntegrator> integrator_ang_;
  };

};

#endif
