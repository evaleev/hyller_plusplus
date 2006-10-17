
#ifndef _hyller_params_h_
#define _hyller_params_h_

#include <vector>

#include <traits.h>
#include <smartptr.h>
#include <clock.h>

namespace hyller {

  /// Variable represents a quantity which may only be changed if mutable is set
  template <typename T>
  class Variable {
  public:
    typedef T Type;

    Variable() :
      value_(Traits<Type>::invalid_value()), mutable_(false) {
    }
    Variable(const Type& value,
	     bool mut) :
      value_(value), mutable_(mut) {
    }
    ~Variable() {
    }

    class ChangingNonMutableVariable : public std::runtime_error {
    public:
      ChangingNonMutableVariable() : std::runtime_error("Attempted to change a non-mutable Variable") {
      }
    };

    const Type& value() const { return value_; }
    void value(const Type& v) {
      if (mutable_)
	value_ = v;
      else
	throw ChangingNonMutableVariable();
    }
    bool mut_able() const { return mutable_; }

    Variable& operator+=(const Type& d) {
      if (mutable_)
	value_ += d;
      else
	throw ChangingNonMutableVariable();
    }

  private:
    Type value_;
    bool mutable_;
  };

  template <typename T>
    bool operator==(const Variable<T>& A, const Variable<T>& B) {
    return (A.value() == B.value() && A.mut_able() == B.mut_able());
  }

  /// ParameterSet owns the parameters
  template <typename P>
  class ParameterSet {
  public:
    /// data type for parameter values
    typedef P Parameter;
    typedef std::vector<Parameter> PSet;
    typedef TheClockType::Timestamp Timestamp;

    ParameterSet(unsigned int n) : params_(new PSet(n)), timestamp_(TheClock->nulltimestamp()) {
    }
    virtual ~ParameterSet() {
    }

    unsigned int n() const { return params_->size(); }
    const Ptr<PSet>& params() const { return params_; }
    /// Parameter access 
    const Parameter& param(unsigned int i) const { return params_->at(i); }
    /** Set parameter.
	This member function is virtual so that the owning class can define the behavior appropriately. For example,
	when a parameter is changed, often the onwing object must do something to update its state.
	In order to do this, the ownership must be expressed as public (or private, or protected) inheritance, not as containment.
	the default behavior is to change the parameters only.

	Update timestamp so that users of this parameter set can check their timestamp to determine whether
	need to recompute themselves.
    */
    virtual void param(unsigned int i, const Parameter& p) {
      params_->at(i) = p;
      timestamp_ = TheClock->timestamp();
    }
    const Timestamp& timestamp() const { return timestamp_; }

  protected:
    /** Default constructor needed if number of parameters is not available until later, e.g. if reading from file.
	Only makes sense if calling from derived class' constructor. */
    ParameterSet() : params_() {
    }
    /** Derived class the also needs to be able to "create" the parameter set. Deep copy is needed here, hence the following is provided */
    void params(const Ptr<PSet>& pset) {
      params_ = pset;
    }

  private:
    Ptr<PSet> params_;
    Timestamp timestamp_;

  };

  /** RefParameterSet refers to parameters owned by someone else. managing parameters via param() member functions is
      someone else's. */
  template <typename P>
  class RefParameterSet {
  public:
    /// data type for parameter values
    typedef P Parameter;
    typedef TheClockType::Timestamp Timestamp;

    RefParameterSet(const Ptr< ParameterSet<Parameter> >& paramset) : paramset_(paramset) {
    }
    virtual ~RefParameterSet() {
    }

    unsigned int n() const { return paramset_->n(); }
    /// Parameter access 
    const Parameter& param(unsigned int i) const {
      return paramset_->param(i);
    }
    /// Set parameter
    void param(unsigned int i, const Parameter& p) {
      return paramset_->param(i,p);
    }
    const Timestamp& timestamp() const { return paramset_->timestamp(); }

  private:
    Ptr< ParameterSet<Parameter> > paramset_;
  };

};

#endif
