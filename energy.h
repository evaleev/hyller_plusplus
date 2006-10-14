
#ifndef _hyller_energy_h_
#define _hyller_energy_h_

#include <smartptr.h>

namespace hyller {

  /// This class computes eigenvalue (and eigenfunction) of a given Hamiltonian. Hamiltonian must be a PFunction<double**, double>.
  template <typename H>
  class EigenEnergy : public PFunction<double,double> {
  public:
    typedef H Hamiltonian;
    typedef typename Hamiltonian::BasisSet BasisSet;
    typedef typename BasisSet::Wavefunction Wavefunction;
    typedef typename BasisSet::BasisFunction BF;
    
    /// solve for r-th root
    EigenEnergy(unsigned int r, const Ptr<Hamiltonian>& h) : PFunction(h->nparam()), root_(r), h_(h) {}
    ~EigenEnergy();

    const Ptr<Wavefunction>& wfn() {
      double E = (*this)();
      return wfn_;
    }
    
  private:
    unsigned int root_;
    Ptr<Hamiltonian> h_;
    Ptr<Wavefunction> wfn_;

    /// Solves the generalized eigenvalue problem
    void compute() {

      // Pass on all parameters to the Hamiltonian
      const unsigned int np = nparam();
      for(unsigned int p=0; p<np; ++p)
	h_->param(p,param(p));

      BasisSet basis = h_->basis();
      const int nbf = basis.nbf();
      const int nprim = basis.nprim();

      // contraction coefficient matrix
      double **C = basis.coefs();
      // eigenvector in the contracted basis
      std::vector<double> evec(basis.nbf());

      //
      // compute the orthogonalizer X
      //

      typedef typename Hamiltonian::Overlap Overlap;
      Ptr<Overlap> ov = H_->overlap();
      double** S = ov();
      fprintf(outfile,"  -Overlap matrix in contracted basis\n");
      print_mat(S,nbf,nbf,outfile);

      // get the orthogonalizer X ...
      double** X;
      int indDim = nbf - sq_sqrt(S,&X,nbf);
      fprintf(outfile,"\t%d basis functions are linearly independent\n",indDim);
      // ... and its inverse Xinv = X^t * S
      double** Xinv = block_matrix(indDim,nbf);
      mmult(X,1,S,0,Xinv,0,indDim,nbf,nbf,0);
      fflush(outfile);

      //
      // Compute and diagonalize the Hamiltonian
      //

      double** Hc = h_();
      fprintf(outfile,"  -Hamiltonian matrix in contracted basis\n");
      print_mat(Hc,nbf,nbf,outfile);

      // convert to orthogonal basis
      double** H = UtVU(X,Hc,nbf,indDim);
      fprintf(outfile,"  -Hamiltonian matrix in orthogonal basis\n");
      print_mat(H,indDim,indDim,outfile);

      // diagonalize
      double* evals = init_array(indDim);
      double** evecs = block_matrix(indDim,indDim);
      sq_rsp(indDim,indDim,H,evals,1,evecs,1.0E-20);
      free_block(H);
      /* Test - prints out all roots */
      for(int i=0;i<=root_;i++)
	fprintf(outfile,"\tState #%d  E = %3.12lf\n",i+1,evals[i]);
      const double E = evals[root_];
      
      for(int i=0;i<indDim;i++) /* temporarily putting eigenvector to evals array */
	evals[i] = evecs[i][root_];
      free_block(evecs);
	
      for(int i=0;i<nbf;i++) {
	double c = 0.0;
	for(int j=0;j<indDim;j++)
	  c += Xinv[j][i]*evals[j];
	evec[i] = c;
      }
	
      free_block(X);
      free_block(Xinv);
	
      Ptr<Wavefunction> wfn(new Wfn(basis,evec));

      wfn_ = wfn;
      value_ = E;
    }

  };

};

#endif
