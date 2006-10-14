
#ifndef _hyller_hamiltonian_h_
#define _hyller_hamiltonian_h_

#include <io.h>
#include <pfunction.h>
#include <smartptr.h>

namespace hyller {

  template <typename B>
  class Overlap : PFunction<double**, double> {
  public:
    typedef B BasisSet;
    typedef typename BasisSet::PrimBF BF;
    typedef PFunction<double**, double> parent;

    Overlap(const Ptr<BasisSet>& bs) : parent(BasisSet::nparam),
      basis_(bs), matrix_(0)
      {
      }
    ~Overlap() {
      if (matrix_) free_block(matrix_);
    }

    const Ptr<BasisSet>& basis() const { return basis_; }

  private:
    Ptr<BasisSet> basis_;
    double** matrix_;

    void compute() {

      const int nbf = basis_.nbf();
      const int nprim = basis_.nprim();
      // contraction coefficient matrix
      double **C = basis_.coefs();

      // Compute S matrix in primitive basis
      double** Sp = block_matrix(nprim,nprim);
      for(int i=0;i<nprim;i++) {
	const BF& bfi = basis.prim(i);
	for(int j=0;j<=i;j++) {
	  const BF& bfj = basis.prim(j);
	  Sp[i][j] = Sp[j][i] = hyller::S(bfi,bfj);
	}
      }
      fprintf(outfile,"  -Overlap matrix in primitive basis\n");
      print_mat(Sp,nprim,nprim,outfile);

      if (matrix_) free_block(matrix_);
      // convert to contracted basis
      matrix_ = UtVU(C,Sp,nprim,nbf);
      fprintf(outfile,"  -Overlap matrix in contracted basis\n");
      print_mat(matrix_,nbf,nbf,outfile);

      value_ = matrix_;
    }

  };

};

#endif

