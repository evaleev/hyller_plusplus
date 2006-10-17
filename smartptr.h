
#ifndef _hyller_smartptr_h_
#define _hyller_smartptr_h_

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace hyller {

  // For now I'll do a cheat since templated typedefs are not standard
  // Should probably at least derive SafePtr from shared_ptr
#define Ptr boost::shared_ptr
#define EnablePtrFromThis boost::enable_shared_from_this
#define Ptr_from_this shared_from_this

  using boost::static_pointer_cast;
  using boost::const_pointer_cast;

  /** Can be used to determine whether a type is a Ptr */
  template <typename T>
    struct IsPtr {
      enum { result = false };
    };
  
  template <typename T>
    struct IsPtr< Ptr<T> > {
      enum { result = true };
    };
  template <typename T>
    struct IsPtr< const Ptr<T> > {
      enum { result = true };
    };
  template <typename T>
    struct IsPtr< Ptr<T>& > {
      enum { result = true };
    };
  template <typename T>
    struct IsPtr< const Ptr<T>& > {
      enum { result = true };
    };

};

#endif
