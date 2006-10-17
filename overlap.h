
#ifndef _hyller_overlap_h_
#define _hyller_overlap_h_

#include <io.h>
#include <pfunction.h>
#include <smartptr.h>

namespace hyller {

  template <typename B>
    class Overlap : public PFunction<double**, typename B::RefPSet> {
  public:
    typedef B BasisSet;
    typedef typename BasisSet::PrimBF BF;
    typedef typename BasisSet::RefPSet ParamSet;
    typedef PFunction<double**, ParamSet > parent;

    Overlap(const Ptr<BasisSet>& bs) :
      parent(Ptr<ParamSet>(new ParamSet(bs->params()))), basis_(bs), matrix_(0)
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

#if 1
      const std::string basisstring = basis_->to_string();
      fprintf(outfile,"  Basis:\n");
      fprintf(outfile,"%s\n",basisstring.c_str());
#endif

      const int nbf = basis_->nbf();
      const int nprim = basis_->nprim();
      // contraction coefficient matrix
      double **C = basis_->coefs();

      // Compute S matrix in primitive basis
      double** Sp = block_matrix(nprim,nprim);
      for(int i=0;i<nprim;i++) {
	const BF& bfi = basis_->prim(i);
	for(int j=0;j<=i;j++) {
	  const BF& bfj = basis_->prim(j);
	  Sp[i][j] = Sp[j][i] = hyller::S(bfi,bfj);
	}
      }
#if 0
      fprintf(outfile,"  -Overlap matrix in primitive basis\n");
      print_mat(Sp,nprim,nprim,outfile);
#endif

      if (matrix_) free_block(matrix_);
      // convert to contracted basis
      matrix_ = UtVU(C,Sp,nprim,nbf);
#if 0
      fprintf(outfile,"  -Overlap matrix in contracted basis\n");
      print_mat(matrix_,nbf,nbf,outfile);
#endif

      parent::set_value(matrix_);
    }

  };

};

#endif

