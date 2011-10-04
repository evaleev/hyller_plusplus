
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "hylleraas.h"
#include "except.h"
#include "matrix.h"

using namespace hyller;

HylleraasBasisFunction::HylleraasBasisFunction() :
  n(-1), l(-1), m(-1), zeta(-1.0) {}

HylleraasBasisFunction::HylleraasBasisFunction(int nn, int ll, int mm, double zz) :
  n(nn), l(ll), m(mm), zeta(zz) { if (zeta <= 0.0) throw std::runtime_error("HylleraasBasisFunction::HylleraasBasisFunction -- nonpositive zeta"); }

std::string
HylleraasBasisFunction::to_string() const {
  std::ostringstream oss;
  if (n != 0) {
    oss << "(r1+r2)^" << n << " * ";
  }
  if (l != 0) {
    oss << "(r1-r2)^" << l << " * ";
  }
  if (m != 0) {
    oss << "r12^" << m;
    if (zeta != 0.0)
      oss << " * ";
  }

  if (zeta == 0.0)
    return oss.str();
  else {
    oss << "Exp[ -" << std::setprecision(15) << zeta/2.0 << "*(r1+r2)]";
  }

  return oss.str();
}

std::string
HylleraasBasisFunction::to_C_string() const {
  std::ostringstream oss;
  if (n != 0) {
    oss << "pow(r1+r2," << n << ")*";
  }
  if (l != 0) {
    oss << "pow(r1-r2," << l << ")*";
  }
  if (m != 0) {
    oss << "pow(r12," << m << ")";
    if (zeta != 0.0)
      oss << " * ";
  }

  if (zeta == 0.0)
    return oss.str();
  else {
    oss << "exp( -" << std::setprecision(15) << zeta/2.0 << "*(r1+r2))";
  }

  return oss.str();
}


bool
hyller::operator==(const HylleraasBasisFunction& A, const HylleraasBasisFunction& B)
{
  return (A.n==B.n && A.l==B.l && A.m==B.m && A.zeta==B.zeta);
}

namespace hyller {
template<>
double value_at<HylleraasBasisFunction>(const HylleraasBasisFunction& bf,
             double r1, double r2, double r12) {
  const double s = r1+r2;
  const double t = r1-r2;
  const double u = r12;
  double value;
  if (bf.m == -1 && u == 0.0) // special case: m=-1 && r12=0 cancel out 0 (l is guaranteed nonzero in this case)
    value = normConst(bf) * pow(s,bf.n) * pow(t,bf.l-1) * exp(- bf.zeta * s / 2.0);
  else // default case
    value = normConst(bf) * pow(s,bf.n) * pow(t,bf.l) * pow(u,bf.m) * exp(- bf.zeta * s / 2.0);

  return value;
}
}

////

HylleraasBasisSet::HylleraasBasisSet(bool singlet, int nlm_max, int nlm_min, int n_max, int l_max, int m_max, double zeta) :
  bfs_(0), singlet_(singlet), nlm_max_(nlm_max), nlm_min_(nlm_min), n_max_(n_max), l_max_(l_max), m_max_(m_max), zeta_(zeta)
{
  if (nlm_min < -1)
    throw std::runtime_error("HylleraasBasisSet::HylleraasBasisSet -- nlm_min must be greater or equal -1");
  bfs_.reserve(expected_num_bf);

  /// iterate over basis functions, given the constraints
  if (nlm_min >= 0)
    // Standard Hylleraas
    for(int n=nlm_min;n<=nlm_max;n++)
      for(int l=nlm_min;l<=nlm_max-n;l++)
        for(int m=nlm_min;m<=nlm_max-n-l;m++) {
          if (n > n_max || l > l_max || m > m_max)
            continue;
          
	  const bool singlet = (l%2 == 0);
	  if (singlet == singlet_) {
	    HylleraasBasisFunction bf(n,l,m,zeta_);
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
	      HylleraasBasisFunction bf(n,l,m,zeta_);
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
	      HylleraasBasisFunction bf(n,l,m,zeta_);
	      bfs_.push_back(bf);
	    }
	  }
    }
  }
}

void
HylleraasBasisSet::set_zeta(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("HylleraasBasisSet::set_zeta -- zeta must be positive");
  zeta_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].zeta = zeta_;
}

int
HylleraasBasisSet::find(const HylleraasBasisFunction& bf) const
{
  std::vector<HylleraasBasisFunction>::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
  if (result == bfs_.end())
    throw BasisFunctionNotFound();
  return result - bfs_.begin();
}

////

HylleraasWfn::HylleraasWfn(const HylleraasBasisSet& bs, const std::vector<double>& coefs) :
  bs_(bs), coefs_(coefs)
{
  if (bs_.num_bf() != coefs_.size())
    throw std::runtime_error("HylleraasWfn::HylleraasWfn -- size of basis set and coefficient vector do not match");
}

std::string
HylleraasWfn::to_string() const {
  const size_t nbf = bs_.num_bf();
  std::ostringstream oss;
  if (nbf) {
    oss << std::setprecision(15) << coefs_[0] * normConst(bs_.bf(0)) * pow(bs_.bf(0).zeta, bs_.bf(0).nlm()+3) << " * "
        << bs_.bf(0).to_string();
    for(int f=1; f<nbf; ++f) {
      oss << " ";
      if (coefs_[f] > 0.0)
        oss << " + ";
      oss << std::setprecision(15) << coefs_[f] * normConst(bs_.bf(f)) * pow(bs_.bf(f).zeta, bs_.bf(f).nlm()+3)   << " * "
          << bs_.bf(f).to_string() << " " << std::endl;
    }
  }
  return oss.str();
}

std::string
HylleraasWfn::to_C_string() const {
  const size_t nbf = bs_.num_bf();
  std::ostringstream oss;
  if (nbf) {
    oss << std::setprecision(15) << coefs_[0] << " * "
        << bs_.bf(0).to_C_string();
    for(int f=1; f<nbf; ++f) {
      oss << " ";
      if (coefs_[f] > 0.0)
        oss << " + ";
      oss << std::setprecision(15) << coefs_[f] << " * "
          << bs_.bf(f).to_C_string() << " " << std::endl;
    }
  }
  return oss.str();
}

////

GenHylleraasBasisFunction::GenHylleraasBasisFunction() :
  i(-1), j(-1), k(-1), zeta1(-1.0), zeta2(-1.0), spin(SpinSinglet)
{
}

GenHylleraasBasisFunction::GenHylleraasBasisFunction(Spin2 s, int ii, int jj, int kk, double zz1, double zz2) :
  i(ii), j(jj), k(kk), zeta1(zz1), zeta2(zz2), spin(s)
{
  if (zeta1 <= 0.0)
    throw std::runtime_error("GenHylleraasBasisFunction::GenHylleraasBasisFunction -- nonpositive zeta1");
  if (zeta2 <= 0.0)
    throw std::runtime_error("GenHylleraasBasisFunction::GenHylleraasBasisFunction -- nonpositive zeta2");
}

bool
hyller::operator==(const GenHylleraasBasisFunction& A, const GenHylleraasBasisFunction& B)
{
  return (A.spin==B.spin && A.i==B.i && A.j==B.j && A.k==B.k && A.zeta1==B.zeta1 && A.zeta2==B.zeta2);
}

GenHylleraasBasisSet::GenHylleraasBasisSet(Spin2 spin, int ijk_max, int i_max, int j_max, int k_max, double zeta1, double zeta2) :
  bfs_(0), spin_(spin), ijk_max_(ijk_max), i_max_(i_max), j_max_(j_max), k_max_(k_max), zeta1_(zeta1), zeta2_(zeta2)
{
  bfs_.reserve(expected_num_bf);

  /// iterate over basis functions, given the constraints
  for(int i=0;i<=ijk_max;i++)
    for(int j=0;j<=ijk_max-i;j++)
      for(int k=0;k<=ijk_max-i-j;k++) {
	if (i > i_max || j > j_max || k > k_max || i < j)
            continue;
        
	GenHylleraasBasisFunction bf(spin_,i,j,k,zeta1_,zeta2_);
	bfs_.push_back(bf);
      }
}

void
GenHylleraasBasisSet::set_zeta1(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("GenHylleraasBasisSet::set_zeta1 -- zeta must be positive");
  zeta1_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].zeta1 = zeta1_;
}

void
GenHylleraasBasisSet::set_zeta2(double zeta)
{
  if (zeta <= 0.0)
    throw std::runtime_error("GenHylleraasBasisSet::set_zeta2 -- zeta must be positive");
  zeta2_ = zeta;
  const int nbf = bfs_.size();
  for(int i=0; i<nbf; i++)
    bfs_[i].zeta2 = zeta2_;
}

int
GenHylleraasBasisSet::find(const GenHylleraasBasisFunction& bf) const
{
  std::vector<GenHylleraasBasisFunction>::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
  if (result == bfs_.end())
    throw BasisFunctionNotFound();
  return result - bfs_.begin();
}
