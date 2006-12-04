
#ifndef _hyller_traits_h_
#define _hyller_traits_h_

namespace hyller {

  /// Traits for value types
  template <typename T>
    struct Traits {};
  /// Traits for double
  template <>
    struct Traits<double> {
      // norm of a double is a double 
      typedef double Norm;
      static double invalid_value() { return -6.66e66; }
    };
  /// Traits for any pointer
  template <typename T>
    struct Traits<T*> {
      // norm of a vector or a matrix is a double 
      typedef double Norm;
      static T* invalid_value() { return (T*)-1; }
    };
  /// Traits for std::pair<T,T>
  template <typename T>
    struct Traits< std::pair<T,T> > {
      typedef typename Traits<T>::Norm Norm;
      static std::pair<T,T> invalid_value() { return std::make_pair(Traits<T>::invalid_value(),Traits<T>::invalid_value()); }
    };

};

#endif
