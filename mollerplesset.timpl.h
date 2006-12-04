
#ifndef _hyller_mollerplessettimpl_h_
#define _hyller_mollerplessettimpl_h_

#define STRONGLY_ORTHOGONAL 0
#define TEST_JFIT 0

#include <fock.h>
#include <fock.timpl.h>
#include <hamiltonian.h>
#include <mollerplesset.h>

namespace hyller {

  template <class B>
  void
  MollerPlessetSeries<B>::compute() {

    typedef typename B::ContrBF CBF;
    typedef typename B::PrimBF PBF;

    // the one-particle exchange will be fit using this basis
    //fbs_ = Ptr<OrbitalBasisSet>(new OrbitalBasisSet(12,0,1.6875));
    fbs_ = Ptr<OrbitalBasisSet>(new OrbitalBasisSet(10,0,3.0));
    const unsigned int nfbf = fbs_->nbf();
    DFFockOperator Fock(hfwfn_,fbs_);
    const int nbf = bs_->nbf();

    // compute the core Hamiltonian
    Ptr< TwoBodyHamiltonian<B> > hcore(new TwoBodyHamiltonian<B>(bs_,Z_,1.0,0.0));
    double** Hc = (*hcore)();
    double** S = (*(hcore->overlap()))();
    fprintf(outfile,"  MollerPlesset::compute -- S:\n");
    print_mat(S,nbf,nbf,outfile);

#if TEST_JFIT
    {
      const unsigned int np = bs_->nprim();
      // compute the two-particle matrix in the primitive basis
      double** Jpp = block_matrix(np,np);
      for(unsigned int b=0; b<np; ++b) {
	for(unsigned int k=0; k<np; ++k) {
	  double J_bk = 0.0;
	  
	  const Orbital& f = fbs_->bf(0);
	  const double J_bik_1 = hyller::gen_mult_oper(0,0,0,f.zeta,0.0,0.0,bs_->prim(b),bs_->prim(k));
	  const double J_bik_2 = hyller::gen_mult_oper(0,0,0,0.0,f.zeta,0.0,bs_->prim(b),bs_->prim(k));
	  J_bk += (J_bik_1 + J_bik_2);

	  Jpp[b][k] = J_bk;
	}
	  
      }

      fprintf(outfile,"  MollerPlesset::compute -- Jtest:\n");
      print_mat(Jpp,np,np,outfile);
    }
#endif


    // Orthogonalizer (X = S^-1/2) and inverse (Xi = S^1/2 = SX)
    double **X;
    int rank;  // rank of the two-electron basis
    {
      const int numDep = sq_sqrt(S,&X,nbf);
      rank = nbf - numDep;
    }
    double** Xi = block_matrix(nbf,rank);
    mmult(S,0,X,0,Xi,0,nbf,nbf,rank,0);

    double** Hch = UtVU(X,Hc,nbf,rank);

    // compute the fock matrix
    double** J12 = Fock.J12(*bs_,*bs_);
#if TEST_JFIT
    fprintf(outfile,"  MollerPlesset::compute -- J12:\n");
    print_mat(J12,nbf,nbf,outfile);
#endif
    double** J12h = UtVU(X,J12,nbf,rank);
    // H0 = Hc + J12
    double** H0h = add(rank,rank,1.0,Hch,1.0,J12h);

    // compute the perturbation operator H1 = 1/r12 - F(1) - F(2)
    Ptr< TwoBodyHamiltonian<B> > hee(new TwoBodyHamiltonian<B>(bs_,Z_,0.0,1.0));
    double** Hee = (*hee)();
    double** H1 = add(nbf,nbf,1.0,Hee,-1.0,J12);
    double** H1h = UtVU(X,H1,nbf,rank);

    // Construct the HF vector in the two-particle basis. bs_ must span the product space of OrbitalBasisSet
    const OrbitalBasisSet& obs = hfwfn_->basis();
    const std::vector<double>& hfcoefs = hfwfn_->coefs();
    std::vector<double> hfcoefs12(nbf,0.0);
    const unsigned int norbs = obs.nbf();
    const double oofourpi = 1.0/(4.0 * M_PI);
    for(unsigned int o1=0; o1<norbs; ++o1) {
      const Orbital& orb1 = obs.bf(o1);
      const double c1 = hfcoefs[o1];
      const double n1 = normConst(orb1);
      for(unsigned int o2=0; o2<norbs; ++o2) {
	const Orbital& orb2 = obs.bf(o2);
	const double c2 = hfcoefs[o2];
	const double n2 = normConst(orb2);

	const PBF pbf12(orb1^orb2);
	const CBF cbf12(pbf12);
	unsigned int i = bs_->find(cbf12);

	hfcoefs12[i] += c1*c2 * n1*n2 * oofourpi;
      }
    }
    double* V0 = to_C_array(hfcoefs12);
    double* V0h = Utv(Xi,V0,nbf,rank);

    double** HE0h;
    {
      HE0h = block_matrix(rank,rank);
      const double epsilon = hfwfn_->E();
      for(unsigned int r=0; r<rank; ++r) {
	for(unsigned int c=0; c<rank; ++c) {
	  HE0h[r][c] = H0h[r][c];
	  if (r == c)
	    HE0h[r][c] -= 2*epsilon;
	}
      }
    }
    double** HE0_inv_h = sq_inverse(HE0h,rank);

    double* H1V0h = Utv(H1h,V0h,rank,rank);
    double E1; dot_arr(V0h,H1V0h,rank,&E1);
    double** HE1h;
    {
      HE1h = block_matrix(rank,rank);
      for(unsigned int r=0; r<rank; ++r) {
	for(unsigned int c=0; c<rank; ++c) {
	  HE1h[r][c] = H1h[r][c];
	  if (r == c)
	    HE1h[r][c] -= E1;
	}
      }
    }
    double* HE1V0h = Utv(HE1h,V0h,rank,rank);


    //
    // Strongly orthogonal projector P = (1 - O1) (1 - O2) = 1 - O1 - O2 + O1O2
    //
    // 1 + O1O2
    double** P0h = block_matrix(rank,rank);
    for(unsigned int r=0; r<rank; ++r) {
      for(unsigned int c=0; c<rank; ++c) {
#if STRONGLY_ORTHOGONAL
	P0h[r][c] =  + V0h[r] * V0h[c];
#else
	P0h[r][c] =  - V0h[r] * V0h[c];
#endif
	if (r == c)
	  P0h[r][c] += 1.0;
      }
    }
#if STRONGLY_ORTHOGONAL
    // - O1 - O2
    {
      typedef typename B::BaseBasisSet BaseBS;
      Ptr<GSHBasisSet> ofbs = obs^(*fbs_); 
      Ptr<GSHBasisSet> fobs = (*fbs_)^obs;
      Ptr< NSOverlap<GSHBasisSet,B> > of_S_g(new NSOverlap<GSHBasisSet,B>(ofbs,bs_));
      Ptr< NSOverlap<GSHBasisSet,B> > fo_S_g(new NSOverlap<GSHBasisSet,B>(fobs,bs_));
      double** ofg1 = (*of_S_g)();
      double** ofg2 = (*fo_S_g)();

      // contract orbital o in <o f | g > to the MO coefficients to get <mo f | g >
      double** fg1 = block_matrix(nfbf,nbf);
      double** fg2 = block_matrix(nfbf,nbf);
      for(unsigned int o=0, of=0; o<norbs; ++o) {
	const double oofourpi = 1.0/(4.0 * M_PI);
	const double hfcoef = hfcoefs[o] * normConst(obs.bf(o)) * oofourpi;
	for(unsigned int f=0; f<nfbf; ++f, ++of) {
	  const unsigned int fo = f*norbs + o;
	  const double ncf = normConst(fbs_->bf(f));
	  const double normcoef = ncf * hfcoef;
	  for(unsigned int g=0; g<nbf; ++g) {
	    fg1[f][g] += ofg1[of][g] * normcoef;
	    fg2[f][g] += ofg2[fo][g] * normcoef;
	  }
	}
      }
      fprintf(outfile,"  <x(1) f|g>\n");
      print_mat(fg1,nfbf,nbf,outfile);
      fprintf(outfile,"  <x(2) f|g>\n");
      print_mat(fg2,nfbf,nbf,outfile);

      // transform |g > to the orthonormal |g'>
      double** fg1h = block_matrix(nfbf,rank);
      double** fg2h = block_matrix(nfbf,rank);
      mmult(fg1,0,X,0,fg1h,0,nfbf,nbf,rank,0);
      mmult(fg2,0,X,0,fg2h,0,nfbf,nbf,rank,0);

      // orthogonalizer for |f>
      double** Xf;
      int frank;
      {
	// get overlap in the fitting basis
	double** S = block_matrix(nfbf,nfbf);
	for(unsigned int i=0, ij=0; i<nfbf; ++i) {
	  for(unsigned int j=0; j<=i; ++j, ++ij) {
	    double S_ij = hyller::S(fbs_->bf(i), fbs_->bf(j));
	    S[i][j] = S[j][i] = S_ij;
	  }
	}
	const int nflindep = sq_sqrt(S,&Xf,nfbf);
	frank = nfbf - nflindep;
      }

      // transform <f| to the orthonormal <f'|
      double** hfg1h = block_matrix(frank,rank);
      double** hfg2h = block_matrix(frank,rank);
      mmult(Xf,1,fg1h,0,hfg1h,0,frank,nfbf,rank,0);
      mmult(Xf,1,fg2h,0,hfg2h,0,frank,nfbf,rank,0);

      fprintf(outfile,"  <x(1) f'|g'>\n");
      print_mat(hfg1h,frank,rank,outfile);
      fprintf(outfile,"  <x(2) f'|g'>\n");
      print_mat(hfg2h,frank,rank,outfile);

#if 1
      for(unsigned int r=0; r<rank; ++r) {
	for(unsigned int c=0; c<rank; ++c) {
	  for(unsigned int f=0; f<frank; ++f) {
	    P0h[r][c] +=  - hfg1h[f][r] * hfg1h[f][c];
	    P0h[r][c] +=  - hfg2h[f][r] * hfg2h[f][c];
	  }
	}
      }
#endif
      
    }
#endif

    double* H0V0h = Utv(H0h,V0h,rank,rank);
    double* P0H0V0h = Utv(P0h,H0V0h,rank,rank);
    double* P0H1V0h = Utv(P0h,H1V0h,rank,rank);
    double* P0HE1V0h = Utv(P0h,HE1V0h,rank,rank);

    // test 0: <0|P0|0> = 0
    {
      double* P0V0h = Utv(P0h,V0h,rank,rank);
      double S; dot_arr(V0h,P0V0h,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 0):  <0|P0|0> = %12.8lf\n", S);
    }
    // test 1: <0|H0|0> = E0
    {
      double Ehf; dot_arr(V0h,H0V0h,rank,&Ehf);
      fprintf(outfile,"Moller-Plesset (test 1):  <0|H0|0> = %12.8lf\n", Ehf);
    }
    // test 2: <0|H1|0> = E1
    {
      double E1; dot_arr(V0h,H1V0h,rank,&E1);
      fprintf(outfile,"Moller-Plesset (test 2):  <0|H1|0> = %12.8lf\n", E1);
    }
    // test 3: <0|H0-E0|0> = 0
    double* HE0V0h = Utv(HE0h,V0h,rank,rank);
    {
      double E; dot_arr(V0h,HE0V0h,rank,&E);
      fprintf(outfile,"Moller-Plesset (test 3):  <0|H0-E0|0> = %12.8lf\n", E);
    }
    // test 4: <0|H1-E1|0> = 0
    {
      double E; dot_arr(V0h,HE1V0h,rank,&E);
      fprintf(outfile,"Moller-Plesset (test 4):  <0|H1-E1|0> = %12.8lf\n", E);
    }
    // test 5: <0|P0 H1|0> = 1
    {
      double S; dot_arr(V0h,P0H1V0h,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 5):  <0|P0 H1|0> = %12.8lf\n", S);
    }
    // test 6: <0|H1 P0 (H0-E0)-1 P0 H1|0> = 0
    double* tmp = Utv(HE0_inv_h,P0H1V0h,rank,rank);
    {
      double* tmp = Utv(HE0_inv_h,P0H1V0h,rank,rank);
      double E2 = 0.0; dot_arr(P0H1V0h,tmp,rank,&E2);
      fprintf(outfile,"Moller-Plesset (test 6):  <0| H1 P0 (H0-E0)-1 P0 H1 |0> = %12.8lf\n", E2);
    }
    // test 7: <0|HE1 P0 HE1|0> = 1
    {
      double S; dot_arr(HE1V0h,P0HE1V0h,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 7):  <0|HE1 P0 HE1|0> = %12.8lf\n", S);
    }
    // test 8: <0|HE1 HE1|0> = 1
    {
      double S; dot_arr(HE1V0h,HE1V0h,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 8):  <0|HE1 HE1|0> = %12.8lf\n", S);
    }
    // test 9: <0|HE0 HE0-1 HE1|0> = 1
    {
      double* tmp = Utv(HE0_inv_h,HE0V0h,rank,rank);
      double S; dot_arr(HE0V0h,tmp,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 9):  <0|HE0 HE0^-1 HE0|0> = %12.8lf\n", S);
    }
    // test 10: <0|HE0 P0 HE0-1 P0 HE1|0> = 1
    {
      double* P0HE0V0h = Utv(P0h,HE0V0h,rank,rank);
      double* tmp = Utv(HE0_inv_h,P0HE0V0h,rank,rank);
      double S; dot_arr(P0HE0V0h,tmp,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 10):  <0|HE0 P0 HE0^-1 P0 HE0|0> = %12.8lf\n", S);
    }
    // test 11: <0|H0 P0 H0|0> = E0
    {
      double S; dot_arr(H0V0h,P0H0V0h,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 1):  <0|H0 P0 H0|0> = %12.8lf\n", S);
    }
    // test 12: spectrum of H0
    {
      double** evecs = block_matrix(rank,rank);
      double* evals = new double[rank];
      sq_rsp(rank,rank,H0h,evals,1,evecs,1.0E-20);
      fprintf(outfile,"Moller-Plesset (test 12):  eigenvalues of H0\n");
      for(int i=0;i<rank;i++)
	fprintf(outfile,"  %d  E = %12.8lf\n",i,evals[i]);
      free_block(evecs);
      delete[] evals;
    }
    // test 13: spectrum of H0 + H1
    {
      double** evecs = block_matrix(rank,rank);
      double* evals = new double[rank];
      double** Hh = add(rank,rank,1.0,H0h,1.0,H1h);
      sq_rsp(rank,rank,Hh,evals,1,evecs,1.0E-20);
      fprintf(outfile,"Moller-Plesset (test 13):  eigenvalues of H0 + H1\n");
      for(int i=0;i<rank;i++)
	fprintf(outfile,"  %d  E = %12.8lf\n",i,evals[i]);
      free_block(evecs);
      delete[] evals;
    }
    // test 14: spectrum of H0 + P0H1P0
    {
      double** evecs = block_matrix(rank,rank);
      double* evals = new double[rank];
      double** P0H1P0h = UtVU(P0h,H1h,rank,rank);
      double** Hh = add(rank,rank,1.0,H0h,1.0,P0H1P0h);
      sq_rsp(rank,rank,Hh,evals,1,evecs,1.0E-20);
      fprintf(outfile,"Moller-Plesset (test 14):  eigenvalues of H0 + P0H1P0\n");
      for(int i=0;i<rank;i++)
	fprintf(outfile,"  %d  E = %12.8lf\n",i,evals[i]);
      free_block(evecs);
      delete[] evals;
    }
    // test 15: <0|0> = 1
    {
      double S; dot_arr(V0h,V0h,rank,&S);
      fprintf(outfile,"Moller-Plesset (test 15):  <0|0> = %12.8lf\n", S);
    }
    // test 16: print P0 H1|0>
    {
      double* P0H1V0 = Uv(Xi,tmp,rank,nbf);
      fprintf(outfile,"Moller-Plesset (test 16):  P0 H1 V0\n");
      for(int i=0;i<rank;i++)
	fprintf(outfile,"  %40s  E = %12.8lf\n",bs_->bf(i).to_string().c_str(),P0H1V0[i]);
    }

  }

};

#endif
