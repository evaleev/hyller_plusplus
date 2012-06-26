/*
 * r12geminal.h
 *
 *  Created on: Oct 4, 2011
 *      Author: Edward Valeev (valeyev.net)
 *      Edward Valeev (C) 2011
 */

#ifndef _hyller_r12geminal_h_
#define _hyller_r12geminal_h_

#include <vector>
#include <wfn.h>

namespace hyller {

  /** Computes \f$ \hat{Q}_{12} f(r_1,r_2,r_{12}) |\chi \chi\rangle \f$
   *  where \f$ \hat{Q}_{12} = (1 - \hat{Q}_1) (1 - \hat{Q}_2) = 1 - \hat{Q}_1 - \hat{Q}_2 + \hat{Q}_1 \hat{Q}_2 \f$
   *  and \f$ f(r_1,r_2,r_{12}) \f$ is a correlation factor
  */
  template <typename B>
  class R12Geminal {
  public:
    typedef B Basis;
    typedef typename Basis::PrimBF PrimBF;
    typedef typename Basis::ContrBF ContrBF;
    typedef PrimBF CorrelationFactor;
    typedef Wavefunction<Basis> WFN;

    /// @param corrfac the correlation factor \f$ \equiv f(r_1,r_2,r_{12}) \f$
    /// @param corrfac_prefactor scaling coefficient for the correlation factor
    R12Geminal(const Ptr<OrbitalWfn>& hfwfn,
               const CorrelationFactor& corrfac,
               double corrfac_prefactor) :
      hfwfn_(hfwfn), corrfac_(corrfac), corrfac_prefactor_(corrfac_prefactor),
      fbs_(new OrbitalBasisSet(std::min(hfwfn->basis().n_max() * 2, 8), 0, hfwfn->basis().zeta()))
      {
        // convert HF wfn into GSH basis
        typedef ContractedBasisFunction<Orbital> COrbital;
        const COrbital uorb(contract(*hfwfn, false)); // in terms of non-normalized primitives
        const COrbital orb(contract(*hfwfn, true)); // in terms of non-normalized primitives
        ContrBF phi0(operator^<PrimBF,Orbital,Orbital>(uorb,uorb));
        phi0 *= 1.0 / (4.0 * M_PI);
        std::cout << "Orbital product wave function (phi0):" << std::endl
                  << phi0.to_string() << std::endl;

        // wrap |0> into a basis
        Ptr<Basis> hfbs(new Basis);
        hfbs->add(phi0);

        // compute R12 geminal x |0> and wrap it into a basis
        Ptr<Basis> r12bs (new Basis);
        ContrBF r12phi0 = phi0 * corrfac;
        r12phi0 *= corrfac_prefactor_;
        r12bs->add(r12phi0);

        // wrap 1/r12 |0> into a basis
        Ptr<Basis> v12bs(new Basis);
        v12bs->add(PrimBF(0,0,-1,0.0,0.0,0.0) * phi0);

        // verify that |0> is normalized
        Ptr<Overlap<Basis> > S00(new Overlap<Basis>(hfbs));
        double** s00_mat = (*S00)();
        assert(fabs(s00_mat[0][0] - 1.0) < 1e-12);

        // compute <R12|0>
        Ptr<NSOverlap<Basis,Basis> > SR0(new NSOverlap<Basis,Basis>(r12bs,hfbs));
        double** sR0_mat = (*SR0)();
        double sR0 = sR0_mat[0][0];
        std::cout << "<0|R12> = " << sR0 << std::endl;

        // compute \hat{O}_i|R12> (i=1,2) using RI/density fitting
        // \hat{O}_1 |R12> \approx \sum_i  |o_1 f_2> <o_1 f_2 | R12>
        // f_i is the fitting one-electron basis
        {
          Ptr<Basis> ofbs(new Basis); // {|of>}
          Ptr<Basis> fobs(new Basis); // {|fo>}
          const unsigned int nfbs = fbs_->nbf();
          for(unsigned int f=0; f<nfbs; ++f) {
            ContrBF of_bf(operator^<PrimBF,Orbital,Orbital>(uorb,
                                                            COrbital(fbs_->bf(f),
                                                                     normConst(fbs_->bf(f)))
                                                           ) );
            of_bf *= 1.0 / (4.0 * M_PI);
            ofbs->add( of_bf );
            ContrBF fo_bf(operator^<PrimBF,Orbital,Orbital>(COrbital(fbs_->bf(f),
                                                                     normConst(fbs_->bf(f))
                                                                    ),
                                                            uorb) );
            fo_bf *= 1.0 / (4.0 * M_PI);
            fobs->add( fo_bf );
          }

          Ptr<NSOverlap<Basis, Basis> > S_of_R(
              new NSOverlap<Basis, B>(ofbs, r12bs));
          Ptr<NSOverlap<Basis, Basis> > S_fo_R(
              new NSOverlap<Basis, Basis>(fobs, r12bs));
          double** ofR1 = (*S_of_R)();
          double** ofR2 = (*S_fo_R)();

          typedef BasisSet<Orbital> OrbitalBasis;
          Ptr<OrbitalBasis> orbbs(new OrbitalBasis);
          orbbs->add(orb);
          Ptr<OrbitalBasis> forbbs(new OrbitalBasis);
          for(unsigned int f=0; f<nfbs; ++f)
            forbbs->add(COrbital(fbs_->bf(f)));

          // get the inverse metric for |f>
          double** Sinvf;
          {
            // get overlap in the fitting basis
            Ptr<Overlap<OrbitalBasis> > Sff(new Overlap<OrbitalBasis>(forbbs));
            double** Sff_mat = (*Sff)();
#if 1
            fprintf(stdout,"  <f|f> matrix\n");
            print_mat(Sff_mat,nfbs,nfbs,stdout);
#endif
            Sinvf = sq_inverse(Sff_mat, nfbs);
          }

          double** f1 = block_matrix(nfbs,1); // = |f> <fo|R12>
          double** f2 = block_matrix(nfbs,1); // = |f> <of|R12>
          mmult(Sinvf,0,ofR1,0,f2,0,nfbs,nfbs,1,0);
          mmult(Sinvf,0,ofR2,0,f1,0,nfbs,nfbs,1,0);
          std::vector<double> f1coefs(nfbs); std::copy(f1[0], f1[0]+nfbs, f1coefs.begin());
          std::vector<double> f2coefs(nfbs); std::copy(f2[0], f2[0]+nfbs, f2coefs.begin());

          OrbitalWfn f1wfn(*fbs_, f1coefs, 0.0);
          OrbitalWfn f2wfn(*fbs_, f2coefs, 0.0);
          const COrbital f1orb(contract(f1wfn));
          const COrbital f2orb(contract(f2wfn));
          ContrBF O2r12wfn(operator^<PrimBF,Orbital,Orbital>(f1orb,uorb));
          ContrBF O1r12wfn(operator^<PrimBF,Orbital,Orbital>(uorb,f2orb));
          O1r12wfn *= 1.0 / (4.0 * M_PI);
          O2r12wfn *= 1.0 / (4.0 * M_PI);

          std::cout << "\\hat{O}_1 |R12> function:" << std::endl
                    << O1r12wfn.to_string() << std::endl;
          std::cout << "\\hat{O}_2 |R12> function:" << std::endl
                    << O2r12wfn.to_string() << std::endl;

#if 1
          { // test these functions <R12|\hat{O}_1 \hat{O}_2|R12> = <0|R12>^2
            Ptr<Basis> O1r12wfn_bs(new Basis);
            O1r12wfn_bs->add(O1r12wfn);
            Ptr<Basis> O2r12wfn_bs(new Basis);
            O2r12wfn_bs->add(O2r12wfn);

            {
              Ptr<NSOverlap<Basis,Basis> > Stmp(new NSOverlap<Basis,Basis>(O1r12wfn_bs, O2r12wfn_bs));
              double** Stmp_mat = (*Stmp)();
              fprintf(stdout,"  <R12|\\hat{O}_1 \\hat{O}_2|R12> = <0|R12>^2 ?\n");
              print_mat(Stmp_mat,1,1,stdout);
            }
            {
              Ptr<NSOverlap<Basis,Basis> > Stmp(new NSOverlap<Basis,Basis>(O1r12wfn_bs, hfbs));
              double** Stmp_mat = (*Stmp)();
              fprintf(stdout,"  <R12|\\hat{O}_1 |0> = <R12|0>?\n");
              print_mat(Stmp_mat,1,1,stdout);
            }
          }

#endif

          // compute <0|1/r12 \\hat{Q}_{12} |R12>
          {
            // compute <0|V12|0>
            double s0V0 = 0.0;
            {
              Ptr<NSOverlap<Basis, Basis> > s(
                  new NSOverlap<Basis, Basis>(hfbs, v12bs));
              double** s_mat = (*s)();
              s0V0 = s_mat[0][0];
              std::cout << "<0|1/r12|0> = " << s0V0 << std::endl;
            }

            // compute <R12|V12|0>
            double sRV0 = 0.0;
            {
              Ptr<NSOverlap<Basis, Basis> > s(
                  new NSOverlap<Basis, Basis>(r12bs, v12bs));
              double** s_mat = (*s)();
              sRV0 = s_mat[0][0];
              std::cout << "<0|1/r12|R12> = " << sRV0 << std::endl;
            }

            // compute <R12|\hat{O}_1 V12|0>
            double sRVQ10 = 0.0;
            {
              Ptr<Basis> O1r12wfn_bs(new Basis);
              O1r12wfn_bs->add(O1r12wfn);
              Ptr<NSOverlap<Basis, Basis> > s(
                  new NSOverlap<Basis, Basis>(O1r12wfn_bs, v12bs));
              double** s_mat = (*s)();
              sRVQ10 = s_mat[0][0];
              std::cout << "<0|1/r12 \\hat{Q}_1 |R12> = " << sRVQ10 << std::endl;
            }

            std::cout << "<0|1/r12 \\hat{Q}_{12} |R12> = " << (sRV0 + sR0 * s0V0 - 2.0 * sRVQ10) << std::endl;
          }

        }

      }

    ~R12Geminal() {
    }

    /// computes the result
    Ptr<WFN> operator()() {
    }

  private:
    Ptr<OrbitalWfn> hfwfn_;
    CorrelationFactor corrfac_;
    double corrfac_prefactor_;
    Ptr<WFN> result_;

    Ptr<OrbitalBasisSet> fbs_;

  };

};

#endif
