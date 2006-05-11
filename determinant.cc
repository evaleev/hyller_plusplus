
#include "determinant.h"
#include "except.h"

using namespace hyller;

SD::SD(const Orbital& i, const Orbital& j, SpinCase2 s) :
  o1(&i), o2(&j), spin(s)
{
}

bool
hyller::operator==(const SD& A, const SD& B)
{
  return (A.o1[0]==B.o1[0] && A.o2[0]==B.o2[0] && A.spin==B.spin);
}

////

SDBasisSet::SDBasisSet(SpinCase2 spin,
		       const OrbitalBasisSet& obs) :
  spin_(spin), obs_(obs)
{
  const int norbs = obs_.num_bf();

  for(int i=0; i<norbs; ++i) {
    for(int j=0; j<norbs; ++j) {

      if (spin_ != SpinAlphaBeta && j >= i)
	continue;

      SD det(obs.bf(i),obs.bf(j),spin_);
      bfs_.push_back(det);
    }
  }

}

int
SDBasisSet::find(const SD& bf) const
{
  std::vector<SD>::const_iterator result = std::find(bfs_.begin(),bfs_.end(),bf);
  if (result == bfs_.end())
    throw BasisFunctionNotFound();
  return result - bfs_.begin();
}
