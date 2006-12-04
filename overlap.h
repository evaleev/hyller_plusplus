
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
      free_block(Sp);

      parent::set_value(matrix_);
    }

  };

  template <typename BraB, typename KetB>
    class NSOverlap : public PFunction<double**, typename BraB::RefPSet> {
  public:
    typedef BraB BraBasis;
    typedef KetB KetBasis;
    typedef BraBasis BasisSet;
    typedef typename BasisSet::PrimBF BF;
    typedef typename BasisSet::RefPSet ParamSet;
    typedef PFunction<double**, ParamSet > parent;

    NSOverlap(const Ptr<BraBasis>& brabs, const Ptr<KetBasis>& ketbs) :
      parent(Ptr<ParamSet>(new ParamSet(brabs->params()))), brabs_(brabs), ketbs_(ketbs), matrix_(0)
      {
      }
    ~NSOverlap() {
      if (matrix_) free_block(matrix_);
    }

    const Ptr<BraBasis>& bra_basis() const { return brabs_; }
    const Ptr<KetBasis>& ket_basis() const { return ketbs_; }

  private:
    Ptr<BraBasis> brabs_;
    Ptr<KetBasis> ketbs_;
    double** matrix_;

    void compute() {

      const int nbra = brabs_->nbf();
      const int nket = ketbs_->nbf();
      const int npbra = brabs_->nprim();
      const int npket = ketbs_->nprim();
      // contraction coefficient matrix
      double **Cbra = brabs_->coefs();
      double **Cket = ketbs_->coefs();

      // Compute S matrix in primitive basis
      double** Sp = block_matrix(npbra,npket);
      for(int i=0;i<npbra;i++) {
	const BF& bfi = brabs_->prim(i);
	for(int j=0;j<npket;j++) {
	  const BF& bfj = ketbs_->prim(j);
	  Sp[i][j] = hyller::S(bfi,bfj);
	}
      }

      if (matrix_) free_block(matrix_);
      // convert to contracted basis
      matrix_ = UtVX(Cbra,npbra,nbra,Sp,Cket,npket,nket);
      free_block(Sp);

      parent::set_value(matrix_);
    }

  };

  /// Computes overlap of an OrbitalBasisSet with basis B
  template <typename B>
    class OrbitalOverlap : public PFunction<std::pair<double**,double**>, typename B::RefPSet> {
  public:
    typedef B BasisSet;
    typedef typename BasisSet::PrimBF BF;
    typedef typename BasisSet::RefPSet ParamSet;
    typedef std::pair<double**,double**> Value;
    typedef PFunction<Value,ParamSet> parent;

    OrbitalOverlap(const OrbitalBasisSet& obs, const Ptr<BasisSet>& bs) :
      parent(Ptr<ParamSet>(new ParamSet(bs->params()))), obs_ (obs), basis_(bs), overlaps_()
      {
      }
    ~OrbitalOverlap() {
      if (overlaps_.first) free_block(overlaps_.first);
      if (overlaps_.second) free_block(overlaps_.second);
    }

    const OrbitalBasisSet& obs() const { return obs_; }
    const Ptr<BasisSet>& basis() const { return basis_; }

  private:
    OrbitalBasisSet obs_;
    Ptr<BasisSet> basis_;
    Value overlaps_;

    void compute() {

      const int norbs = obs_.nbf();
      const int nbf = basis_->nbf();
      const int nprim = basis_->nprim();
      // contraction coefficient matrix
      double **C = basis_->coefs();

      // Compute S matrix in primitive basis
      double** S1p = block_matrix(norbs,nprim);
      double** S2p = block_matrix(norbs,nprim);
      for(int i=0;i<norbs;i++) {
	const Orbital& bfi = obs_.bf(i);
	const double nci = normConst(bfi);
	for(int j=0;j<nprim;j++) {
	  const BF& bfj = basis_->prim(j);
	  S1p[i][j] = nci * hyller::S1(bfi,bfj);
	  S2p[i][j] = nci * hyller::S2(bfi,bfj);
	}
      }
#if 1
      fprintf(outfile,"  -Overlap matrix 1 in primitive basis\n");
      print_mat(S1p,norbs,nprim,outfile);
      fprintf(outfile,"  -Overlap matrix 2 in primitive basis\n");
      print_mat(S2p,norbs,nprim,outfile);
#endif

      if (overlaps_.first)  free_block(overlaps_.first);
      if (overlaps_.second) free_block(overlaps_.second);
      // convert to contracted basis
      overlaps_.first  = block_matrix(norbs,nbf);
      overlaps_.second = block_matrix(norbs,nbf);
      mmult(S1p,0,C,0,overlaps_.first ,0,norbs,nprim,nbf,0);
      mmult(S2p,0,C,0,overlaps_.second,0,norbs,nprim,nbf,0);
#if 1
      fprintf(outfile,"  -Overlap matrix 1 in contracted basis\n");
      print_mat(overlaps_.first, norbs,nbf,outfile);
      fprintf(outfile,"  -Overlap matrix 2 in contracted basis\n");
      print_mat(overlaps_.second,norbs,nbf,outfile);
#endif
      free_block(S1p);
      free_block(S2p);

      parent::set_value(overlaps_);
    }

  };

};

#endif

