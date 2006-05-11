
#include "except.h"

using namespace hyller;

BasisFunctionNotFound::BasisFunctionNotFound() :
  std::exception()
{
}

BasisFunctionNotFound::~BasisFunctionNotFound() throw()
{
}

const char*
BasisFunctionNotFound::what() const throw()
{
  return "BasisFunctionNotFound::what() -- basis function not found";
}
