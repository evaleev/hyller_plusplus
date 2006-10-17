
#ifndef _hyller_hamiltonian_h_
#define _hyller_hamiltonian_h_

#include <overlap.h>
#include <pfunction.h>
#include <io.h>

namespace hyller {

  template <typename B>
  class TwoBodyHamiltonian : public PFunction<double**, typename B::RefPSet > {
  public:
    typedef B BasisSet;
    typedef typename BasisSet::PrimBF BF;
    typedef Overlap<BasisSet> Overlap;
    typedef typename BasisSet::RefPSet ParamSet;
    typedef PFunction<double**, ParamSet> parent;

    TwoBodyHamiltonian(const Ptr<BasisSet>& bs, double Z) :
      parent(Ptr<ParamSet>(new ParamSet(bs->params()))),
      basis_(bs), overlap_(new Overlap(bs)), Z_(Z), matrix_(0)
      {
      }
    ~TwoBodyHamiltonian() {
      if (matrix_) free_block(matrix_);
    }

    const Ptr<BasisSet>& basis() const { return basis_; }
    double Z() const { return Z_; }
    const Ptr<Overlap>& overlap() { return overlap_; }

  private:
    Ptr<BasisSet> basis_;
    Ptr<Overlap> overlap_;
    double Z_;
    double** matrix_;

    void compute() {

      const int nbf = basis_->nbf();
      const int nprim = basis_->nprim();
      // contraction coefficient matrix
      double **C = basis_->coefs();

      // compute Hamiltonian in the primitive basis
      double** Hp = block_matrix(nprim,nprim);
      for(int i=0;i<nprim;i++) {
	const BF& bfi = basis_->prim(i);
	for(int j=0;j<=i;j++) {
	  const BF& bfj = basis_->prim(j);
	  
	  double Hij = T(bfi,bfj);
	  Hij += Z_*V_ne(bfi,bfj);
	  Hij += V_ee(bfi,bfj);
	  
	  Hp[j][i] = Hp[i][j] = Hij;
	}
      }
      fprintf(outfile,"  -Hamiltonian matrix in primitive basis\n");
      print_mat(Hp,nprim,nprim,outfile);

      // convert to contracted basis
      if (matrix_)
	free_block(matrix_);
      matrix_ = UtVU(C,Hp,nprim,nbf);
      free_block(Hp);
      fprintf(outfile,"  -Hamiltonian matrix in contracted basis\n");
      print_mat(matrix_,nbf,nbf,outfile);

      parent::set_value(matrix_);
    }

  };

};

#endif

