
#ifndef _hyller_matrix_h_
#define _hyller_matrix_h_

#include "hylleraas.h"
#include "orbital.h"
#include "determinant.h"
#include "csf.h"

namespace hyller {
  /**
     matrix elements in Hylleraas basis
   */
  double normConst(const HylleraasBasisFunction& bf);
  double Overlap(const HylleraasBasisFunction& bfi,
		 const HylleraasBasisFunction& bfj);
  double V_en(const HylleraasBasisFunction& bfi,
	      const HylleraasBasisFunction& bfj);
  double V_ee(const HylleraasBasisFunction& bfi,
	      const HylleraasBasisFunction& bfj);
  double T(const HylleraasBasisFunction& bfi,
	   const HylleraasBasisFunction& bfj);
  double density_Nuc(const HylleraasBasisFunction& bfi,
		     const HylleraasBasisFunction& bfj);

  /**
     matrix elements in Slater-Hylleraas basis
   */
  /// The elementary integral is the integral over whole space of the given function without the Jacobian!
  double I(const SlaterHylleraasBasisFunction& bf);
  double normConst(const SlaterHylleraasBasisFunction& bf);
  double Overlap(const SlaterHylleraasBasisFunction& bfi,
		 const SlaterHylleraasBasisFunction& bfj);
  double V_en(const SlaterHylleraasBasisFunction& bfi,
	      const SlaterHylleraasBasisFunction& bfj);
  double V_ee(const SlaterHylleraasBasisFunction& bfi,
	      const SlaterHylleraasBasisFunction& bfj);
  double T(const SlaterHylleraasBasisFunction& bfi,
	   const SlaterHylleraasBasisFunction& bfj);

  /**
     matrix elements in orbital basis phi(r,i,z) = r^i e^(-zr)
   */
  /// Normalization constant for phi(i,z)
  double normConst(const Orbital& bf);
  /// Overlap between phi(i,z) and phi(j,z)
  double Overlap(const Orbital& bfi, const Orbital& bfj);
  double V_en(const Orbital& bfi, const Orbital& bfj);
  double T(const Orbital& bfi, const Orbital& bfj);

  /**
     matrix elements in determinant basis
  */
  // normConst is not provided for determinants because they are normalized automatically
  /// overlap
  double Overlap(const SD& bfi, const SD& bfj);

  /**
     matrix elements in CSF basis
  */
  // normConst is not provided for CSFs because they cannot be normalized easily
  /// overlap
  double Overlap(const CSF& bfi, const CSF& bfj);
}

#endif
