
#ifndef _hyller_clock_h_
#define _hyller_clock_h_

#include <smartptr.h>

namespace hyller {

  /// Clock issues timestamps in the increasing order. It is a Singleton.
  template <typename T>
  class Clock {
  public:
    typedef T Time;
    typedef Time Timestamp;

    static const Ptr<Clock>& Instance() {
      if (!default_clock_) {
	Ptr<Clock> ptr(new Clock<T>);
	default_clock_ = ptr;
      }
      return default_clock_;
    }
    virtual ~Clock() {
    }

    /// Must be an atomic operation
    virtual Timestamp timestamp() const {
      // no locking yet
      Time result = time_;
      ++time_;
      return result;
    }

    static Timestamp nulltimestamp() { return -1; }

  private:
    static Ptr<Clock> default_clock_;
    Clock() : time_((T)0) {
    }

    mutable Time time_;
  };

  template <typename T> Ptr<Clock<T> > Clock<T>::default_clock_(new Clock<T>);

  /// By default use a long int
  typedef Clock<long int> TheClockType;
  extern Ptr<TheClockType> TheClock;
};

#endif
