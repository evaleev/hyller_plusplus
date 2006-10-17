
#ifndef _hyller_pfunction_h_
#define _hyller_pfunction_h_

#include <smartptr.h>
#include <params.h>
#include <traits.h>
#include <clock.h>

namespace hyller {

  /// PFunction class is an abstract interface to functions which depend on a set of parameters
  template <typename V, typename PS>
  class PFunction {
  public:
    /// data type for function values
    typedef V Value;
    /// data type for parameter set
    typedef PS PSet;
    /// data type for parameter values
    typedef typename PSet::Parameter Parameter;
    /// timestamps
    typedef TheClockType::Timestamp Timestamp;

    /// Default
    PFunction() : params_(), obsolete_(true), timestamp_(TheClock->nulltimestamp()),
      value_(Traits<Value>::invalid_value()) {
    }
    /// Parameter set passed explicitly
      PFunction(const Ptr<PSet>& params) : params_(params), obsolete_(true), timestamp_(TheClock->nulltimestamp()),
	value_(Traits<Value>::invalid_value()) {
    }
    virtual ~PFunction() {
    }

    /// Parameter set
    const Ptr<PSet>& params() const { return params_; }
    /// Number of parameters
    unsigned int nparam() const { return params_->n(); }
    /// Returns i-th parameter
    const Parameter& param(unsigned int i) const { return params_->param(i); }
    /// Sets i-th parameter
    void param(unsigned int i, const Parameter& p) {
      if (!(p == params_->param(i))) {
	params_->param(i,p);
      }
    }
    /// To force recomputation
    void obsolete() { obsolete_ = true; }
    /// Timestamp
    const Timestamp& timestamp() const { return timestamp_; }
    /// Compute the value
    const Value& operator()() {
      if ( (timestamp() < params_->timestamp()) || obsolete_)
	compute();
      obsolete_ = false;
      timestamp_ = TheClock->timestamp();
      return value_;
    }

  protected:
    /// Needed if at the time of construction parameter set does not exist (e.g. construction of an energy object starts before geometry has been read)
    void set_params(const Ptr<PSet>& p) { params_ = p; }
    void set_value(const Value& v) { value_ = v; }

  private:
    /// compute() must call value()!
    virtual void compute() =0;       // compute() is private so that children cannot force a compute without changing obsolete_ to false

    Ptr< PSet > params_;
    bool obsolete_;
    Timestamp timestamp_;
    Value value_;
  };

};

#endif
