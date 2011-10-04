
#include "includes.h"
#include "projector.h"
#include "matrix.h"
#include "except.h"
#include "misc.h"

extern "C" FILE* outfile;

using namespace hyller;

SD_2_Hylleraas::SD_2_Hylleraas(const SDBasisSet& sbs, const HylleraasBasisSet& hbs) :
  sbs_(sbs), hbs_(hbs)
{
  const Spin2 spin2 = hbs_.singlet() ? SpinSinglet : SpinTriplet;
  if (samespin(sbs_.spin()) && hbs_.singlet())
    throw std::runtime_error("SD_2_Hylleraas::SD_2_Hylleraas -- spin case of the determinant basis incompatible with the spin of Hylleraas basis");

  const int nsd = sbs_.num_bf();
  const int nhyll = hbs_.num_bf();
  fprintf(outfile,"\n\tNumber of slater determinants = %d\n", nsd);

  // expansion coefficient for ij product in terms of Hylleraas functions
  X_ = block_matrix(nsd,nhyll);
  for(int d=0; d<nsd; ++d) {
    // expansion coefficient for ij product in terms of Hylleraas functions
    std::vector<double> xij = orbital_to_hylleraas(sbs.bf(d).o1[0],sbs.bf(d).o2[0],hbs);
    for(int h=0; h<nhyll; h++)
      X_[d][h] = xij[h];
  }
  fprintf(outfile,"\t Slater determinants (in rows) expressed in terms of Hylleraas functions\n");
  print_mat(X_,nsd,nhyll,outfile);

  // compute overlap between hylleraas functions
  double** Sh = block_matrix(nhyll,nhyll);
  for(int i=0; i<nhyll; i++) {
    const HylleraasBasisFunction& bfi = hbs_.bf(i);
    for(int j=0; j<=i; j++) {
      const HylleraasBasisFunction& bfj = hbs_.bf(j);
      Sh[i][j] = Sh[j][i] = S(bfi,bfj);
    }
  }

  // compute overlap between slater determinants via Hylleraas basis
  double** tmp1 = block_matrix(nsd,nhyll);
  mmult(X_,0,Sh,0,tmp1,0,nsd,nhyll,nhyll,0);
  double** Sd = block_matrix(nsd,nsd);
  mmult(tmp1,0,X_,1,Sd,0,nsd,nhyll,nsd,0);
  fprintf(outfile,"\t Overlap matrix in slater determinant basis computed via Hylleraas basis\n");
  print_mat(Sd,nsd,nsd,outfile);

  // compute overlap between slater determinants directly
  for(int d1=0; d1<nsd; ++d1) {
    for(int d2=0; d2<nsd; ++d2) {
      Sd[d1][d2] = S(sbs.bf(d1),sbs.bf(d2));
    }
  }
  fprintf(outfile,"\t Overlap matrix in slater determinant basis computed directly\n");
  print_mat(Sd,nsd,nsd,outfile);

  free_block(tmp1);
  free_block(Sd);
  free_block(Sh);
}

SD_2_Hylleraas::~SD_2_Hylleraas()
{
  free_block(X_);
}

////

CSF_2_Hylleraas::CSF_2_Hylleraas(const CSFBasisSet& cbs, const HylleraasBasisSet& hbs) :
  cbs_(cbs), hbs_(hbs)
{
  const Spin2 spin2 = hbs_.singlet() ? SpinSinglet : SpinTriplet;
  if (cbs_.spin() != spin2)
    throw std::runtime_error("CSF_2_Hylleraas::CSF_2_Hylleraas -- spin cases incompatible");

  const int ncsf = cbs_.num_bf();
  const int nhyll = hbs_.num_bf();
  fprintf(outfile,"\n\tNumber of CSF = %d\n", ncsf);

  // expansion coefficient for Hylleraas fucntions in terms of CSF
  X_ = block_matrix(nhyll,ncsf);
  for(int h=0; h<nhyll; ++h) {
    // WARNING only functions with m=0 can be expanded in terms of s-type CSFs
    const HylleraasBasisFunction& bf = hbs_.bf(h);
    if (bf.m !=0 ) continue;
    // expansion coefficient for Hylleraas function in terms of CSFs
    std::vector<double> x = hylleraas_to_CSF(bf,cbs_);
    for(int c=0; c<ncsf; c++)
      X_[h][c] = x[c];
  }
  fprintf(outfile,"\t Hylleraas functions (in rows) expressed in terms of CSFs\n");
  print_mat(X_,nhyll,ncsf,outfile);

  // expansion coefficient for CSF in terms of Hylleraas functions
  Xinv_ = block_matrix(ncsf,nhyll);
  for(int c=0; c<ncsf; ++c) {
    // expansion coefficient for CSF in terms of Hylleraas functions
    std::vector<double> x = CSF_to_hylleraas(cbs.bf(c),hbs);
    for(int h=0; h<nhyll; h++)
      Xinv_[c][h] = x[h];
  }
  fprintf(outfile,"\t CSFs (in rows) expressed in terms of Hylleraas functions\n");
  print_mat(Xinv_,ncsf,nhyll,outfile);

  // compute overlap between hylleraas functions
  double** Sh = block_matrix(nhyll,nhyll);
  for(int i=0; i<nhyll; i++) {
    const HylleraasBasisFunction& bfi = hbs_.bf(i);
    for(int j=0; j<=i; j++) {
      const HylleraasBasisFunction& bfj = hbs_.bf(j);
      Sh[i][j] = Sh[j][i] = S(bfi,bfj);
    }
  }

  // compute overlap between CSFs via Hylleraas basis
  double** tmp1 = block_matrix(ncsf,nhyll);
  mmult(Xinv_,0,Sh,0,tmp1,0,ncsf,nhyll,nhyll,0);
  double** Sd = block_matrix(ncsf,ncsf);
  mmult(tmp1,0,Xinv_,1,Sd,0,ncsf,nhyll,ncsf,0);
  fprintf(outfile,"\t Overlap matrix in CSF basis computed via Hylleraas basis\n");
  print_mat(Sd,ncsf,ncsf,outfile);

  // compute overlap between CSFs directly
  for(int c1=0; c1<ncsf; ++c1) {
    for(int c2=0; c2<ncsf; ++c2) {
      Sd[c1][c2] = S(cbs.bf(c1),cbs.bf(c2));
    }
  }
  fprintf(outfile,"\t Overlap matrix in CSF basis computed directly\n");
  print_mat(Sd,ncsf,ncsf,outfile);


  fprintf(outfile,"\t Overlap matrix in Hylleraas basis basis computed directly\n");
  print_mat(Sh,nhyll,nhyll,outfile);

  // compute overlap between Hylleraas functions via CSFs
  double** tmp2 = block_matrix(nhyll,ncsf);
  mmult(X_,0,Sd,0,tmp2,0,nhyll,ncsf,ncsf,0);
  mmult(tmp2,0,X_,1,Sh,0,nhyll,ncsf,nhyll,0);
  fprintf(outfile,"\t Overlap matrix in Hylleraas basis computed via CSFs\n");
  print_mat(Sh,nhyll,nhyll,outfile);

  free_block(tmp1);
  free_block(tmp2);
  free_block(Sd);
  free_block(Sh);

  // test Xinv.X
  double** XinvX = block_matrix(ncsf,ncsf);
  mmult(Xinv_,0,X_,0,XinvX,0,ncsf,nhyll,ncsf,0);
  fprintf(outfile,"\tShould be 1\n");
  print_mat(XinvX,ncsf,ncsf,outfile);
  free_block(XinvX);
  // test X.Xinv
  double** XXinv = block_matrix(nhyll,nhyll);
  mmult(X_,0,Xinv_,0,XXinv,0,nhyll,ncsf,nhyll,0);
  fprintf(outfile,"\tShould be 1\n");
  print_mat(XXinv,nhyll,nhyll,outfile);
  free_block(XXinv);
}

CSF_2_Hylleraas::~CSF_2_Hylleraas()
{
  free_block(X_);
  free_block(Xinv_);
}

////

SD_2_CSF::SD_2_CSF(const SDBasisSet& sbs, const CSFBasisSet& cbs) :
  sbs_(sbs), cbs_(cbs)
{
  if (samespin(sbs_.spin()) && cbs_.spin() == SpinSinglet)
    throw std::runtime_error("SD_2_CSF::SD_2_CSF -- spin cases incompatible");

  const int nsd = sbs_.num_bf();
  const int ncsf = cbs_.num_bf();
  fprintf(outfile,"\n\tNumber of SD  = %d\n", nsd);
  fprintf(outfile,"\n\tNumber of CSF = %d\n", ncsf);

  // expansion coefficient for CSF in terms of SD
  X_ = block_matrix(ncsf,nsd);
  for(int c=0; c<ncsf; ++c) {
    const CSF& csf = cbs_.bf(c);
    for(int i=0; i<csf.nd; ++i) {
      const SD& sd = csf.d[i][0];
      try {
	const int ii = sbs_.find(sd);
	X_[c][ii] = csf.c[i];
      }
      catch (BasisFunctionNotFound&) {}
    }
  }
  fprintf(outfile,"\t CSFs (in rows) expressed in terms of Slater determinants\n");
  print_mat(X_,ncsf,nsd,outfile);
}

SD_2_CSF::~SD_2_CSF()
{
  free_block(X_);
}

////

GenHylleraas_2_Hylleraas::GenHylleraas_2_Hylleraas(const GenHylleraasBasisSet& gbs, const HylleraasBasisSet& hbs) :
  gbs_(gbs), hbs_(hbs)
{
  if ((gbs_.spin() == SpinSinglet)^ hbs_.singlet())
    throw std::runtime_error("GenHylleraas_2_Hylleraas::GenHylleraas_2_Hylleraas -- different spin cases");

  const int nghyll = gbs_.num_bf();
  const int nhyll = hbs_.num_bf();
  fprintf(outfile,"\n\tNumber of GenHylleraas  = %d\n", nghyll);
  fprintf(outfile,"\n\tNumber of Hylleraas = %d\n", nhyll);

  // expansion coefficient for Hylleraas in terms of GenHylleraas
  X_ = block_matrix(nhyll,nghyll);

  // contributions from i,j functions
  double** ijcoeff = new double*[hbs_.nlm_max()];
  for(int i=0; i<hbs_.nlm_max(); i++) {
    ijcoeff[i] = new double[hbs_.nlm_max()];
  }
  
  for(int k=0; k<nhyll; ++k) {
    const HylleraasBasisFunction& h = hbs_.bf(k);
    for(int i=0; i<h.n+h.l; i++) {
      for(int j=0; j<h.n+h.l; j++) {
	ijcoeff[i][j] = 0.0;
      }
    }
    const double hnorm = normConst(h);
    for(int i1=0; i1<h.n; i1++) {
      for(int i2=0; i2<h.l; i2++) {
	GenHylleraasBasisFunction g(gbs_.spin(),i1+i2,h.n+h.l-i1-i2,h.m,h.zeta,h.zeta);
	ijcoeff[i1+i2][h.n+h.l-i1-i2] += binomial(h.n,i1) * binomial(h.l,i2) * pow(-1.0,h.l-i2) * normConst(g) / hnorm;
      }
    }

    for(int i=0; i<h.n+h.l; i++) {
      for(int j=0; j<=i; j++) {
	if (fabs(fabs(ijcoeff[i][j]) - fabs(ijcoeff[i][j])) > 1.0e-12)
	  throw std::runtime_error("GenHylleraas_2_Hylleraas::GenHylleraas_2_Hylleraas -- coefficients of spin-components are nor equal");
	GenHylleraasBasisFunction g(gbs_.spin(),i,j,h.m,h.zeta,h.zeta);
	const int ii = gbs_.find(g);
	X_[k][ii] = ijcoeff[i][j];
      }
    }

  }
  fprintf(outfile,"\t Hylleraas functions (in rows) expressed in terms of GenHylleraas functions\n");
  print_mat(X_,nhyll,nghyll,outfile);
}

GenHylleraas_2_Hylleraas::~GenHylleraas_2_Hylleraas() {
}

////

GSH_2_Hylleraas::GSH_2_Hylleraas(const HylleraasBasisSet& hbs) :
  gshbs_(make_gshbs(hbs)), hbs_(hbs)
{
  const int ngsh = gshbs_.num_bf();
  const int nhyll = hbs_.num_bf();
  fprintf(outfile, "\n\tNumber of GenSlaterHylleraas  = %d\n", ngsh);
  fprintf(outfile, "\n\tNumber of Hylleraas = %d\n", nhyll);

  // expansion coefficient for Hylleraas in terms of GSH
  X_ = block_matrix(nhyll, ngsh);

  // contributions from i,j functions
  std::vector< std::vector<double> > ijcoeff(hbs_.nlm_max()+1);
  for (int i = 0; i <= hbs_.nlm_max(); i++) {
    ijcoeff.at(i).resize(hbs_.nlm_max()+1);
  }

  for (int k = 0; k < nhyll; ++k) {
    const HylleraasBasisFunction& h = hbs_.bf(k);
    for (int i = 0; i <= h.n + h.l; i++) {
      for (int j = 0; j <= h.n + h.l; j++) {
        ijcoeff.at(i).at(j) = 0.0;
      }
    }
    const double hnorm = normConst(h) * pow(h.zeta, 2.5 + 0.5 * h.nlm()) * pow(h.zeta/2.0, -1.5 *  (h.nlm())) * pow(2.0 * M_PI, h.nlm());
    abort();
    for (int i1 = 0; i1 <= h.n; i1++) {
      for (int i2 = 0; i2 <= h.l; i2++) {
        GenSlaterHylleraasBasisFunction g(i1 + i2, h.n + h.l - i1 - i2,
                                          h.m, h.zeta/2, h.zeta/2, 0.0);
        ijcoeff.at(i1 + i2).at(h.n + h.l - i1 - i2) += binomial(h.n, i1)
            * binomial(h.l, i2) * pow(-1.0, h.l - i2) * NormConst(g) / hnorm;
      }
    }

    for (int i = 0; i <= h.n + h.l; i++) {
      for (int j = 0; j <= h.n + h.l; j++) {
        GenSlaterHylleraasBasisFunction g(i, j, h.m, h.zeta/2.0, h.zeta/2.0, 0.0);
        const int ii = gshbs_.find(g);
        X_[k][ii] = ijcoeff.at(i).at(j);
      }
    }

  }
  fprintf(
      outfile,
      "\t Hylleraas functions (in rows) expressed in terms of GenHylleraas functions\n");
  print_mat(X_, nhyll, ngsh, outfile);
}

GSH_2_Hylleraas::~GSH_2_Hylleraas() {
}

GenSlaterHylleraasBasisSet
GSH_2_Hylleraas::make_gshbs(const HylleraasBasisSet& hbs) {
  const int nhyll = hbs.num_bf();
  int ijk_max = 0;
  int ijk_min = 1000000000;
  int ij_max  = 0;
  int k_max  = 0;
  for(int bf=0; bf<nhyll; ++bf) {
    const int ij = hbs.bf(bf).n + hbs.bf(bf).l;
    const int k = hbs.bf(bf).m;
    ij_max = std::max(ij_max, ij);
    k_max = std::max(k_max, k);
    ijk_min = std::min(ijk_min, ij+k);
  }
  ijk_max = ij_max + k_max;
  GenSlaterHylleraasBasisSet result(ijk_max, ijk_min, ij_max, k_max, hbs.zeta()/2.0, hbs.zeta()/2.0, 0.0);
  return result;
}


