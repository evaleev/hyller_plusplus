
#ifndef _hyller_basisfn_h_
#define _hyller_basisfn_h_

#include <sstream>
#include <algorithm>
#include <utility>

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
    ~ContractedBasisFunction() {}

    /// Number of terms in the contraction
    unsigned int n() const { return contr_.size(); }
    /// Return the n-th term
    const ContrTerm& term(unsigned int n) const { return (safe ? contr_.at(n) : contr_[n]); }

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

    std::string to_string() const {
      std::ostringstream oss;
      const unsigned int nt = n();
      if (nt) {
	oss << " " << contr_[0].second << " * " << contr_[0].first.to_string() << std::endl;	
	for(unsigned int t=1; t<nt; ++t) {
	  oss << "+" << contr_[t].second << " * " << contr_[t].first.to_string() << std::endl;
	}
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

};

#endif
