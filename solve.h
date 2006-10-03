
#ifndef _hyller_solve_h_
#define _hyller_solve_h_

// These must be global -- required by Psi3
extern "C" {
  extern FILE *infile,*outfile;
  extern char *psi_file_prefix;
}

namespace hyller {

  template <typename Basis>
    typename Basis::Wfn&
    solveCi(Basis& basis, double Z, int root_num, bool opt, double threshold, int maxiter)
    {
      typedef typename Basis::Wfn Wfn;
      typedef typename Basis::BasisFunction BF;

      double **H,**S,**X,**temp,*evals;
      double E = 1.0;
      double E_old = 0.0;
      double Ed = 0.0;
      double Ed2 = 0.0;
      double tmp,nuc_dens,alpha,gamma,u,x;
      int iter = 1;
      int indDim = 0;
      
      double alpha_new = basis.alpha();
      double gamma_new = basis.gamma();
      const int nbf = basis.num_bf();
      std::vector<double> coeff(basis.num_bf());
      
      while(iter == 1 || fabs(E-E_old) > threshold) {
	E_old = E;
	alpha = alpha_new;
	gamma = gamma_new;
	if (opt) {
	  fprintf(outfile,"\t Test message : iteration #%d\n",iter);
	  fprintf(outfile,"\t Test message : alpha = %3.12lf\n",alpha);
	  fprintf(outfile,"\t Test message : gamma = %3.12lf\n",gamma);
	  fflush(outfile);
	}
	if (iter == 1) {
	  S = block_matrix(nbf,nbf);
	  for(int i=0;i<nbf;i++) {
	    const BF& bfi = basis.bf(i);
	    for(int j=0;j<=i;j++) {
	      const BF& bfj = basis.bf(j);
	      S[i][j] = S[j][i] = hyller::S(bfi,bfj) * NormConst(bfi) * NormConst(bfj);
	    }
	  }
	  fprintf(outfile,"\t Overlap matrix\n");
	  print_mat(S,nbf,nbf,outfile);
	  indDim = nbf - sq_sqrt(S,&X,nbf);
	  fprintf(outfile,"\t%d basis functions are linearly independent\n",indDim);
	  free_block(S);
	  fflush(outfile);
	}
	H = block_matrix(nbf,nbf);
	for(int i=0;i<nbf;i++) {
	  const BF& bfi = basis.bf(i);
	  for(int j=0;j<=i;j++) {
	    const BF& bfj = basis.bf(j);
	    
	    const double norm_pfac = NormConst(bfi) * NormConst(bfj);

	    double Hij = T(bfi,bfj);
	    Hij += Z*V_ne(bfi,bfj);
	    //Hij += Z*(-DeltaR1(bfi,bfj)-DeltaR2(bfi,bfj));
	    Hij += V_ee(bfi,bfj);
	    //Hij += DeltaR12(bfi,bfj);

	    Hij *= norm_pfac;
	    H[j][i] = H[i][j] = Hij;
	  }
	}
	fprintf(outfile,"\tHamiltonian matrix\n");
	print_mat(H,nbf,nbf,outfile);
	
	temp = block_matrix(nbf,indDim);
	
	mmult(H,0,X,0,temp,0,nbf,nbf,indDim,0);
	free_block(H);
	H = block_matrix(indDim,indDim);
	mmult(X,1,temp,0,H,0,indDim,nbf,indDim,0);
	free_block(temp);
	temp = block_matrix(indDim,indDim);
	
	evals = init_array(indDim);
	sq_rsp(indDim,indDim,H,evals,1,temp,1.0E-20);
	free_block(H);
	/* Test - prints out all roots */
	for(int i=0;i<=root_num;i++)
	  fprintf(outfile,"\tState #%d  E = %3.12lf\n",i+1,evals[i]);
	E = evals[root_num];
	
	for(int i=0;i<indDim;i++) /* temporarily putting eigenvector to evals array */
	  evals[i] = temp[i][root_num];
	free_block(temp);
	
	for(int i=0;i<nbf;i++) {
	  double c = 0.0;
	  for(int j=0;j<indDim;j++)
	    c += X[i][j]*evals[j];
	  coeff[i] = c;
	}
	
	if (!opt)
	  free_block(X);
	
	tmp = 0;
	for(int i=0;i<nbf;i++) {
	  const BF& bfi = basis.bf(i);
	  for(int j=0;j<nbf;j++) {
	    const BF& bfj = basis.bf(j);
	    tmp =  T(bfi,bfj);
	    Ed += tmp*coeff[i]*coeff[j]; 
	  }
	}
#if 0
	Ed2 = Ed*2/(zeta*zeta);
	Ed = (Ed+E)/zeta;
	//zeta_new = zeta-100000*Ed/Ed2;
	zeta_new = zeta-Ed/Ed2;
	if (opt) {
	  fprintf(outfile,"\tdE/d(zeta) = %3.12lf\n",Ed);
	  fprintf(outfile,"\t Test message : new zeta1S/2 = %3.12lf\n",zeta_new/2);
	  fprintf(outfile,"\t Test message : delta(E) = %3.12lf\n",E-E_old);
	}
#endif
	iter++;
	
	if (!opt)
	  break;     /* if optimization flag is not set - quit 'while' loop */
	else
	  fflush(outfile);
      }
      
#if TEST_NORMALIZATION
      /* verify c^T . S . c = 1 */
      {
	double norm = 0.0;
	for(int i=0;i<nbf;i++) {
	  const BF& bfi = basis.bf(i);
	  for(int j=0;j<nbf;j++) {
	    const BF& bfj = basis.bf(j);
	    norm += coeff[i] * coeff[j] *
	      NormConst(bfi) * NormConst(bfj) * hyller::S(bfi,bfj);
	  }
	}
	fprintf(outfile,"\tTEST: c^T * S * c = %10.9lf (should be 1.0)\n",norm);
      }
#endif

      Wfn* wfn = new Wfn(basis,coeff);
      return *wfn;
    }

};

#endif
