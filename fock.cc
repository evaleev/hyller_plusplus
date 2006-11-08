
#include <matrix.h>
#include <misc.h>
#include <fock.h>

#define F77_DSPSVX dspsvx_
extern "C" {
extern void F77_DSPSVX(const char* fact, const char* uplo, const int* n, const int* nrhs,
                       const double* AP, double* AFP, int* ipiv, const double* BB, const int* nb,
                       double* XX, const int* nx, double* rcond, double* ferr, double* berr,
                       double* work, int* iwork, int* info);
};

using namespace hyller;

DFFockOperator::DFFockOperator(const Ptr<OrbitalWfn>& hfwfn,
			       const Ptr<OrbitalBasisSet>& fbs) :
  hfwfn_(hfwfn), fbs_(fbs)
{
}

DFFockOperator::~DFFockOperator()
{
}

const Ptr<OrbitalWfn>&
DFFockOperator::J_df() const
{
  const OrbitalBasisSet& obs = hfwfn_->basis();
  const unsigned int nbf = obs.nbf();
  const unsigned int nbf2 = nbf*nbf;
  const unsigned int nfbf = fbs_->nbf();

  const std::vector<double>& coefs = hfwfn_->coefs();
  std::vector<double> J_f(nfbf,0.0);

  // Compute the one-center J matrix using three-center electron repulsion integral and the density.
  for(unsigned int i=0; i<nbf; ++i) {
    const double coef_i = coefs[i];
    for(unsigned int j=0; j<=i; ++j) {
      const double coef_j = coefs[j];

      // assuming 1s^2 state!
      const double ij_perm_pfac = (i == j) ? 1.0 : 2.0;
      const double D_ij = 2.0 * coef_i * coef_j * ij_perm_pfac;

      for(unsigned int k=0; k<nfbf; ++k) {

	const double J_kij = V_ee_3ct(fbs_->bf(k), obs.bf(i), obs.bf(j));
	J_f[k] += J_kij * D_ij;

      }
    }
  }

  // get overlap in the fitting basis
  const unsigned int nftri = nfbf * (nfbf+1)/2;
  double* S = new double[nftri];
  double* SF = new double[nftri];
  for(unsigned int i=0, ij=0; i<nfbf; ++i) {
    for(unsigned int j=0; j<=i; ++j, ++ij) {
      double S_ij = hyller::S(fbs_->bf(i), fbs_->bf(j));
      S[ij] = S_ij;
    }
  }

  const double* J_f_v = to_C_array(J_f);
  double* X_v = new double[nfbf];

  // solve the linear system J_f = S * X, where X are the desired fitting coefficient vector 
  int *ipiv = new int[nfbf];
  char fact = 'N';
  char uplo = 'U';
  int nrhs = 1;
  double rcond = 0.0;
  double* ferr = new double[1];
  double* berr = new double[1];
  double* work = new double[3*nfbf];
  int* iwork = new int[nfbf];
  int info = 0;
  const int n = nfbf;
  F77_DSPSVX(&fact, &uplo, &n, &nrhs, S, SF, ipiv, J_f_v, &n, X_v, &n, &rcond, ferr, berr, work, iwork, &info);

  if (info) {
    if (info > 0 && info <= nfbf)
      throw std::runtime_error("DFFockOperator::J_df() -- matrix S has factors which are exactly zero");
    if (info == nfbf + 1)
      throw std::runtime_error("DFFockOperator::J_df() -- matrix S is singular");
    if (info < 0)
      throw std::runtime_error("DFFockOperator::J_df() -- one of the arguments to DSPSVX is invalid");
  }
  delete[] ferr;
  delete[] berr;
  delete[] work;
  delete[] iwork;

  // convert the solution to OrbitalWfn and return
  std::vector<double> fit(nfbf);
  std::copy(X_v,X_v+nfbf+1,fit.begin());
  Ptr<OrbitalWfn> result(new OrbitalWfn(*fbs_,fit));
  J_df_ = result;

  delete[] X_v;
  delete[] J_f_v;
  delete[] S;
  delete[] SF;

  return J_df_;
}


