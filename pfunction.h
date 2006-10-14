
#ifndef _hyller_pfunction_h_
#define _hyller_pfunction_h_

namespace hyller {

  /// Traits for value types
  template <typename T>
    struct Traits {};
  /// Traits for double
  template <>
    struct Traits<double> {
      double invalid_value() { return -6.66e66; }
    };

  /// PFunction class is an abstract interface to functions which depend on 0 or more parameters
  template <typename V, typename P>
  class PFunction {
  public:
    /// data type for function values
    typedef V Value;
    /// data type for parameter values
    typedef P Parameter;

    /// Function depends on n parameters
    PFunction(unsigned int n) : params_(n), obsolete_(true), value_(Traits<Value>::invalid_value()) {
    }
    ~PFunction() {
    }

    /// Number of parameters
    unsigned int nparam() const { return params_.size(); }
    /// Returns i-th parameter
    const Parameter& param(unsigned int i) const { return params_.at(i); }
    /// Sets i-th parameter
    void param(unsigned int i, const Parameter& p) {
      if (p != params_.at(i)) {
	params_[i] = p;
	obsolete_ = true;
      }
    }
    /// Obsolete
    void obsolete() { obsolete_ = true; }
    /// Compute the value
    const Value& operator()() {
      if (obsolete_)
	compute();
      obsolete_ = false;
      return value_;
    }

  protected:
    std::vector<Parameter> params_;
    bool obsolete_;
    Value value_;

  private:
    /// compute() must update value_!
    virtual void compute() =0;       // compute() is private so that children cannot force a compute without changing obsolete_ to false

  };

};

#endif
