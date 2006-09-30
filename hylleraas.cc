
#include <stdexcept>
#include "hylleraas.h"
#include "except.h"

using namespace hyller;

HylleraasBasisFunction::HylleraasBasisFunction() :
  n(-1), l(-1), m(-1), zeta(-1.0) {}

HylleraasBasisFunction::HylleraasBasisFunction(int nn, int ll, int mm, double zz) :
  n(nn), l(ll), m(mm), zeta(zz) { if (zeta <= 0.0) throw std::runtime_error("HylleraasBasisFunction::HylleraasBasisFunction -- nonpositive zeta"); }

bool
hyller::operator==(const HylleraasBasisFunction& A, const HylleraasBasisFunction& B)
{
  return (A.n==B.n && A.l==B.l && A.m==B.m && A.zeta==B.zeta);
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
