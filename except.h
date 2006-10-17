
#ifndef _hyller_except_h_
#define _hyller_except_h_

#include <stdexcept>

namespace hyller {

  /// basis function not found -- this may be OK
  class BasisFunctionNotFound : public std::exception {
  public:
    BasisFunctionNotFound();
    ~BasisFunctionNotFound() throw();
    const char* what() const throw();
  };

  /// basis function parameters do not match the basis set parameters
  class InvalidBasisFunctionParams : public std::exception {
  public:
    InvalidBasisFunctionParams();
    ~InvalidBasisFunctionParams() throw();
    const char* what() const throw();
  };

};

#endif
