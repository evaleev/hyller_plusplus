
#ifndef _hyller_basisfn_h_
#define _hyller_basisfn_h_

#include <vector>
#include <sstream>
#include <algorithm>
#include <utility>
#include <iomanip>

namespace hyller {

  namespace {
    template <typename BF, typename X> struct first_equal {
      typedef std::pair<BF,X> first_argument_type;
      typedef BF second_argument_type;
      typedef bool result_type;
      bool operator()(const std::pair<BF,X>& p, const BF& b) const { return p.first == b; }
    };
  };

  /// A contracted basis function
  template <typename BF>
  class ContractedBasisFunction {
    static const int safe = 0;

  public:
    typedef BF PrimBF;
    typedef std::pair<PrimBF,double> ContrTerm;
    
    /// Empty contraction
    ContractedBasisFunction() {}
    /// Contraction contains a single function
    ContractedBasisFunction(const PrimBF& a) : contr_(ContrTermSet(1,std::make_pair(a,1.0))) {}
    /// Contraction contains a single function, with a coefficient
    ContractedBasisFunction(const PrimBF& a, double coef) : contr_(ContrTermSet(1,std::make_pair(a,coef))) {}
    ~ContractedBasisFunction() {}

    /// Number of terms in the contraction
    unsigned int n() const { return contr_.size(); }
    /// Return the n-th term
    const ContrTerm& term(unsigned int n) const { return (safe ? contr_.at(n) : contr_[n]); }
    /// Sets the n-th term
    void term(unsigned int n, const ContrTerm& t) { safe ? (contr_.at(n) = t): (contr_[n] = t); }

    /// Add a term. If this primitive exists -- simply add the coefficient (if SafeAddPolicy is false) or throw (if SafeAddPolicy is true)
    void add(const ContrTerm& c) {
      typename ContrTermSet::iterator result = std::find_if(contr_.begin(),contr_.end(),std::bind2nd(first_equal<PrimBF,double>(),c.first));
      if (result == contr_.end()) {
	contr_.push_back(c);
      }
      else {
	if (safe)
	  throw std::runtime_error("ContractedBasisFunction::add() -- this primitive function already exists");
	else {
	  result->second += c.second;
	}
      }
    }

    /// Add a term. See above.
    void add(const PrimBF& p, double coef) {
      add(std::make_pair(p,coef));
    }

    void operator*=(double c) {
      const unsigned int nt = n();
      for(unsigned int t=0; t<nt; ++t)
	contr_[t].second *= c;
    }

    /// Will return false if the order of primitives is not the same
    bool operator==(const ContractedBasisFunction& A) const
    {
      const unsigned int nt = n();
      if (nt != A.n()) return false;
      for(int t=0; t<nt; ++t) {
        if (contr_[t] != A.contr_[t]) return false;
      }
      return true;
    }

    std::string to_string() const {
      std::ostringstream oss;
      const unsigned int nt = n();
      if (nt) {
        oss << "  " << std::fixed << std::setprecision(20) << contr_[0].second << " * " << contr_[0].first.to_string();
        for(unsigned int t=1; t<nt; ++t) {
          oss << " + " << std::fixed << std::setprecision(20) << contr_[t].second << " * " << contr_[t].first.to_string();
        }
        oss << std::endl;
      }
      return oss.str();
    }

    std::string to_C_string(bool skip_exponential) const {
      std::ostringstream oss;
      const unsigned int nt = n();
      if (nt) {
        oss << " ( ";
        if (contr_[0].second != 1.0)
          oss << std::setprecision(15) << contr_[0].second << " * ";
        oss << contr_[0].first.to_C_string(skip_exponential) << std::endl;
        for(unsigned int t=1; t<nt; ++t) {
          oss << "+ ";
          if (contr_[t].second != 1.0)
            oss << std::setprecision(15) << contr_[t].second << " * ";
          oss << contr_[t].first.to_C_string(skip_exponential) << std::endl;
        }
        oss << " ) " << std::endl;
      }
      return oss.str();
    }
    
    private:
    typedef std::vector<ContrTerm> ContrTermSet;
    ContrTermSet contr_;
    
  };

  /// Returns a sum of 2 contracted functions. A difference can be computed using this operator combined with "*="
  template <typename BF>
    ContractedBasisFunction<BF> operator+(const ContractedBasisFunction<BF>& A,
					  const ContractedBasisFunction<BF>& B)
    {
      typedef ContractedBasisFunction<BF> CBF;
      CBF result;
      const unsigned int na = A.n();
      const unsigned int nb = B.n();
      for(unsigned int i=0; i<na; ++i)
	result.add(A.term(i));
      for(unsigned int i=0; i<nb; ++i)
	result.add(B.term(i));
      return result;
    }
  /// Returns a product of 2 contracted functions
  template <typename BF>
    ContractedBasisFunction<BF> operator*(const ContractedBasisFunction<BF>& A,
					  const ContractedBasisFunction<BF>& B)
    {
      typedef ContractedBasisFunction<BF> CBF;
      typedef typename CBF::ContrTerm Term;
      CBF result;
      const unsigned int na = A.n();
      const unsigned int nb = B.n();
      for(unsigned int i=0; i<na; ++i) {
	const Term& ta = A.term(i);
	for(unsigned int j=0; j<nb; ++j) {
	  const Term& tb = B.term(j);
	  result.add(ta.first*tb.first,ta.second*tb.second);
	}
      }
      return result;
    }
  /// Returns a product of 1 contracted function and 1 primitive function
  template <typename BF>
    ContractedBasisFunction<BF> operator*(const BF& A,
					  const ContractedBasisFunction<BF>& B)
    {
      typedef ContractedBasisFunction<BF> CBF;
      typedef typename CBF::ContrTerm Term;
      CBF result;
      const unsigned int nb = B.n();
      for(unsigned int j=0; j<nb; ++j) {
	const Term& tb = B.term(j);
	result.add(A*tb.first,tb.second);
      }
      return result;
    }
  /// Returns a product of 1 contracted function and 1 primitive function
  template <typename BF>
    ContractedBasisFunction<BF> operator*(const ContractedBasisFunction<BF>& B,
					  const BF& A) {
    return A*B;
  }
  /// Returns an outer product of 2 contracted functions
  template <typename BF, typename bf1, typename bf2>
    ContractedBasisFunction<BF> operator^(const ContractedBasisFunction<bf1>& A,
					  const ContractedBasisFunction<bf2>& B)
    {
      typedef ContractedBasisFunction<BF> CBF; 
      typedef typename ContractedBasisFunction<bf1>::ContrTerm TermA;
      typedef typename ContractedBasisFunction<bf2>::ContrTerm TermB;
      typedef typename CBF::ContrTerm Term;
      CBF result;
      const unsigned int na = A.n();
      const unsigned int nb = B.n();
      for(unsigned int i=0; i<na; ++i) {
	const TermA& ta = A.term(i);
	const bf1& fa = ta.first;
	for(unsigned int j=0; j<nb; ++j) {
	  const TermB& tb = B.term(j);
	  const bf2& fb = tb.first;
	  result.add(fa^fb,ta.second*tb.second);
	}
      }
      return result;
    }

    template<typename BF>
    double value_at(const ContractedBasisFunction<BF>& bf,
                    double r1, double r2, double r12) {
      typedef ContractedBasisFunction<BF> CBF;
      typedef typename CBF::ContrTerm Term;
      double value = 0.0;
      for(unsigned int i=0; i<bf.n(); ++i) {
        const Term& term = bf.term(i);
        value += term.second * value_at(term.first, r1, r2, r12) / NormConst(term.first);
      }
      return value;
    }

};

#endif
