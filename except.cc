
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

////

InvalidBasisFunctionParams::InvalidBasisFunctionParams() :
  std::exception()
{
}

InvalidBasisFunctionParams::~InvalidBasisFunctionParams() throw()
{
}

const char*
InvalidBasisFunctionParams::what() const throw()
{
  return "InvalidBasisFunctionParams::what() -- basis function parameters do not match those of the basis set";
}
