
#ifndef _hyller_matrix_h_
#define _hyller_matrix_h_

#include "hylleraas.h"
#include "slaterhylleraas.h"
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
     matrix elements in nonsymmetric Hylleraas basis (i.e. in terms of r1, r2, and r12)
  */
  double normConst(const GenHylleraasBasisFunction& bf);
  double Overlap(const GenHylleraasBasisFunction& bfi,
		 const GenHylleraasBasisFunction& bfj);
  double V_en(const GenHylleraasBasisFunction& bfi,
	      const GenHylleraasBasisFunction& bfj);
  double V_ee(const GenHylleraasBasisFunction& bfi,
	      const GenHylleraasBasisFunction& bfj);
  double T(const GenHylleraasBasisFunction& bfi,
	   const GenHylleraasBasisFunction& bfj);

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



  /////////////////////////////////////////
  // Generic functions for matrix elements
  /////////////////////////////////////////

  // The fundamental operations:

  /** The "radial" integral:
      I = \int f(r_1,r_2,r_{12}) dr_1 dr_2 dr_{12}.
      
      The integral is elementary: it does not include the Jacobian and angles are not integrated.
      All other integrals depend on this function.
   */
  template <typename F>
    double I(const F& f);

  /** Overlap between 2 unnormalized functions. Implementation of this function depends on the details of F and the Jacobian.

      Hence operator*(class F&, class F&) and apply_jacobian(class F&) must be defined for class F.
  */
  template <typename F>
    double S(const F& bra, const F& ket);
  /** Computes the normalization constant. Uses S(). */
  template <typename F>
    double NormConst(const F& bf);
  /** Matrix element of the r_1^i r_2^j r_{12}^k operator over unnormalized functions. Uses S(). */
  template <typename F>
    double gen_mult_oper(int i, int j, int k, const F& bra, const F& ket);
  /** Matrix element of the electron repulsion operator (r_{12}^{-1}) over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double V_ee(const F& bra, const F& ket);
  /** Matrix element of the nuclear attraction operator (- r_1^{-1} - r_2^{-1}) over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double V_ne(const F& bra, const F& ket);
  /** Matrix element of the kinetic energy operator over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double T(const F& bra, const F& ket);

}

#endif
