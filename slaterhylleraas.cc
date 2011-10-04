
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "slaterhylleraas.h"
#include "except.h"
#include "matrix.h"
#include <orbital.h>
#include <matrix.timpl.h>

using namespace hyller;

SlaterHylleraasBasisFunction::SlaterHylleraasBasisFunction() :
  n(0), l(0), m(0), zeta(0.0), gamma(0.0) {}

SlaterHylleraasBasisFunction::SlaterHylleraasBasisFunction(int nn, int ll, int mm, double zz, double gg) :
  n(nn), l(ll), m(mm), zeta(zz), gamma(gg)
{
  if (zeta < 0.0)
    throw std::runtime_error("SlaterHylleraasBasisFunction::SlaterHylleraasBasisFunction -- nonpositive zeta");
#define ALLOW_NEGATIVE_GAMMA 1
#if ALLOW_NEGATIVE_GAMMA
  if (gamma < 0.0)
    throw std::runtime_error("SlaterHylleraasBasisFunction::SlaterHylleraasBasisFunction -- nonpositive gamma");
#endif
}

bool
hyller::operator==(const SlaterHylleraasBasisFunction& A, const SlaterHylleraasBasisFunction& B)
{
  return (A.n==B.n && A.l==B.l && A.m==B.m && A.zeta==B.zeta && A.gamma==B.gamma);
}

////

SlaterHylleraasBasisSet::SlaterHylleraasBasisSet(bool singlet, int nlm_max, int nlm_min, int n_max, int l_max, int m_max, double zeta, double gamma) :
  bfs_(0), singlet_(singlet), nlm_max_(nlm_max), nlm_min_(nlm_min), n_max_(n_max), l_max_(l_max), m_max_(m_max), zeta_(zeta), gamma_(gamma)
{
  if (nlm_min < -1)
    throw std::runtime_error("SlaterHylleraasBasisSet::SlaterHylleraasBasisSet -- nlm_min must be greater or equal -1");
  bfs_.reserve(expected_num_bf);

  /// iterate over basis functions, given the constraints
  if (nlm_min >= 0)
    // Standard SlaterHylleraas
    for(int n=nlm_min;n<=nlm_max;n++)
      for(int l=nlm_min;l<=nlm_max-n;l++)
        for(int m=nlm_min;m<=nlm_max-n-l;m++) {
          if (n > n_max || l > l_max || m > m_max)
            continue;
          
	  const bool singlet = (l%2 == 0);
	  if (singlet == singlet_) {
	    SlaterHylleraasBasisFunction bf(n,l,m,zeta_, gamma_);
	    bfs_.push_back(bf);
	  }
        }
  else {
    // Kinoshita case
    if (singlet_) {

      for(int n=-1;n<=nlm_max;n++)
	for(int l=0;l<=nlm_max+2;l+=2)
	  for(int m=-1;m<=nlm_max+1;m++) {
	    
	    if (n > n_max || l > l_max || m > m_max)
	      continue;
	    
	    if ((n+l+m <= nlm_max) && (n+l+m >= -1) && (l+m >= 0)) {
	      SlaterHylleraasBasisFunction bf(n,l,m,zeta_,gamma_);
	      bfs_.push_back(bf);
	    }
	  }
    }
    else {

      for(int n=-1;n<=nlm_max;n++)
	for(int l=1;l<=nlm_max+2;l+=2)
	  for(int m=-1;m<=nlm_max;m++) {
	    if (n > n_max || l > l_max || m > m_max)
	      continue;
          
	    if (n+l+m <= nlm_max && n+l+m >= -1) {
	      SlaterHylleraasBasisFunction bf(n,l,m,zeta_,gamma_);
	      bfs_.push_back(bf);
	    }
	  }
    }
  }
}

void
SlaterHylleraasBasisSet::set_zeta(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("SlaterHylleraasBasisSet::set_zeta -- zeta must be positive");
  zeta_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].zeta = zeta_;
}

void
SlaterHylleraasBasisSet::set_gamma(double gamma)
{
  if (gamma <= 0.0)
    throw std::runtime_error("SlaterHylleraasBasisSet::set_gamma -- gamma must be positive");
  gamma_ = gamma;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].gamma = gamma_;
}

int
SlaterHylleraasBasisSet::find(const SlaterHylleraasBasisFunction& bf) const
{
  std::vector<SlaterHylleraasBasisFunction>::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
  if (result == bfs_.end())
    throw BasisFunctionNotFound();
  return result - bfs_.begin();
}

////

SlaterHylleraasWfn::SlaterHylleraasWfn(const SlaterHylleraasBasisSet& bs, const std::vector<double>& coefs) :
  bs_(bs), coefs_(coefs)
{
  if (bs_.num_bf() != coefs_.size())
    throw std::runtime_error("SlaterHylleraasWfn::SlaterHylleraasWfn -- size of basis set and coefficient vector do not match");
}


////////////////////////

GenSlaterHylleraasBasisFunction::GenSlaterHylleraasBasisFunction() :
  i(-1), j(-1), k(-1), alpha(-1.0), beta(-1.0), gamma(-1.0)
{
}

GenSlaterHylleraasBasisFunction::GenSlaterHylleraasBasisFunction(int ii, int jj, int kk, double a, double b, double g) :
  i(ii), j(jj), k(kk), alpha(a), beta(b), gamma(g)
{
  if (alpha < 0.0)
    throw std::runtime_error("GenSlaterHylleraasBasisFunction::GenSlaterHylleraasBasisFunction -- negative alpha");
  if (beta < 0.0)
    throw std::runtime_error("GenSlaterHylleraasBasisFunction::GenSlaterHylleraasBasisFunction -- negative beta");
  if (gamma < 0.0)
    throw std::runtime_error("GenSlaterHylleraasBasisFunction::GenSlaterHylleraasBasisFunction -- negative gamma");
}

GenSlaterHylleraasBasisFunction
GenSlaterHylleraasBasisFunction::Identity(0,0,0,0.0,0.0,0.0);

std::string
GenSlaterHylleraasBasisFunction::to_string() const {
  std::ostringstream oss;
  if (i != 0) {
    oss << "r_1^" << i << " * ";
  }
  if (j != 0) {
    oss << "r_2^" << j << " * ";
  }
  if (k != 0) {
    oss << "r_{12}^" << k << " * ";
  }
  
  if (alpha == 0.0 && beta == 0.0 && gamma == 0.0)
    return oss.str();
  else {
    oss << "e^{";
  }
  
  if (alpha != 0.0)
    oss << "-" << alpha << "*r_1";
  if (beta != 0.0)
    oss << "-" << beta << "*r_2";
  if (gamma != 0.0)
    oss << "-" << gamma << "*r_{12}";
  oss << "}";
  return oss.str();
}

std::string
GenSlaterHylleraasBasisFunction::to_C_string(bool skip_exponential) const {
  std::ostringstream oss;
  bool have_preexp = false;
  if (i != 0) {
    oss << "r1_"  << i;
    have_preexp = true;
  }
  if (j != 0) {
    if (have_preexp) oss << " * ";
    oss << "r2_"  << j;
    have_preexp = true;
  }
  if (k != 0) {
    if (have_preexp) oss << " * ";
    oss << "r12_" << k;
    have_preexp = true;
  }

  if (alpha == 0.0 && beta == 0.0 && gamma == 0.0 || skip_exponential) {
    if (have_preexp == false) oss << "1";
    return oss.str();
  }
  else {
    if (have_preexp) oss << " * ";
    oss << "exp(";
  }

  if (alpha != 0.0)
    oss << "-" << std::setprecision(15) << alpha << "*r1_1";
  if (beta != 0.0)
    oss << "-" << std::setprecision(15) << beta << "*r2_1";
  if (gamma != 0.0)
    oss << "-" << std::setprecision(15) << gamma << "*r12_1";
  oss << ")";
  return oss.str();
}

bool
hyller::operator==(const GenSlaterHylleraasBasisFunction& A, const GenSlaterHylleraasBasisFunction& B)
{
  return (A.i==B.i && A.j==B.j && A.k==B.k && A.alpha==B.alpha && A.beta==B.beta && A.gamma==B.gamma);
}

GenSlaterHylleraasBasisFunction
hyller::operator*(const GenSlaterHylleraasBasisFunction& A, const GenSlaterHylleraasBasisFunction& B)
{
  return GenSlaterHylleraasBasisFunction(A.i+B.i,A.j+B.j,A.k+B.k,A.alpha+B.alpha,A.beta+B.beta,A.gamma+B.gamma);
}

GenSlaterHylleraasBasisFunction
hyller::apply_jacobian(const GenSlaterHylleraasBasisFunction& A)
{
  return GenSlaterHylleraasBasisFunction(A.i+1,A.j+1,A.k+1,A.alpha,A.beta,A.gamma);
}

GenSlaterHylleraasBasisFunction
hyller::gen_r1r2r12_oper(int i, int j, int k)
{
  return GenSlaterHylleraasBasisFunction(i,j,k,0.0,0.0,0.0);
}

GenSlaterHylleraasBasisFunction
hyller::operator^(const Orbital& f1,
		  const Orbital& f2)
{
  return GenSlaterHylleraasBasisFunction(f1.n,f2.n,0,f1.zeta,f2.zeta,0.0);
}

namespace hyller {
  template <typename TwoElectronBasisFunction>
  double value_at(const TwoElectronBasisFunction& bf, double r1, double r2, double r12);

  template<>
  double value_at<GenSlaterHylleraasBasisFunction>(const GenSlaterHylleraasBasisFunction& bf,
                                                   double r1, double r2, double r12) {
    const double value = NormConst(bf) *
        pow(r1,bf.i) * pow(r2,bf.j) * pow(r12,bf.k) *
        exp(- bf.alpha * r1
            - bf.beta * r2
            - bf.gamma * r12);
    return value;
  }
}

GenSlaterHylleraasBasisSet::GenSlaterHylleraasBasisSet(int ijk_max, int ijk_min, int ij_max, int k_max, double alpha, double beta, double gamma,
                                                       bool odd_k) :
  bfs_(0), ijk_max_(ijk_max), ij_max_(ij_max), k_max_(k_max), alpha_(alpha), beta_(beta), gamma_(gamma)
{
  bfs_.reserve(expected_num_bf);

  const bool alpha_eq_beta = (alpha_ == beta_);
  /// iterate over basis functions, given the constraints
  for (int i = 0; i <= ijk_max; i++)
    for (int j = 0; j <= ijk_max - i; j++)
      for (int k = 0; k <= ijk_max - i - j; k++) {
        if (i > ij_max ||
            j > ij_max ||
            k > k_max ||
            (i + j + k < ijk_min) ||
            (odd_k==false && k%2==1)
           )
          continue;

        GenSlaterHylleraasBasisFunction bf(i, j, k, alpha_, beta_, gamma_);
        bfs_.push_back(bf);
        // Include the function where alpha and beta are permuted also
        if (i == j && !alpha_eq_beta) {
          GenSlaterHylleraasBasisFunction bf(i, j, k, beta_, alpha_, gamma_);
          bfs_.push_back(bf);
        }
      }
}

GenSlaterHylleraasBasisSet::GenSlaterHylleraasBasisSet(const GenSlaterHylleraasBasisSet & gbs) :
    bfs_(gbs.bfs_), ijk_max_(gbs.ijk_max_), ij_max_(gbs.ij_max_), k_max_(gbs.k_max_),
    alpha_(gbs.alpha_), beta_(gbs.beta_), gamma_(gbs.gamma_)
{
  std::cout << "GenSlaterHylleraasBasisSet cctor -- " << bfs_.size() << " " << gbs.bfs_.size() << std::endl;
}

void
GenSlaterHylleraasBasisSet::set_alpha(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("GenSlaterHylleraasBasisSet::set_alpha -- alpha must be positive");
  alpha_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].alpha = alpha_;
}

void
GenSlaterHylleraasBasisSet::set_beta(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("GenSlaterHylleraasBasisSet::set_beta -- beta must be positive");
  beta_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].beta = beta_;
}

void
GenSlaterHylleraasBasisSet::set_gamma(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("GenSlaterHylleraasBasisSet::set_gamma -- gamma must be positive");
  gamma_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].gamma = gamma_;
}

int
GenSlaterHylleraasBasisSet::find(const GenSlaterHylleraasBasisFunction& bf) const
{
  std::vector<GenSlaterHylleraasBasisFunction>::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
  if (result == bfs_.end())
    throw BasisFunctionNotFound();
  return result - bfs_.begin();
}
