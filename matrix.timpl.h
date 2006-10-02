
#ifndef _hyller_matrixtimpl_h_
#define _hyller_matrixtimpl_h_

namespace hyller {

  /** Overlap between 2 unnormalized functions. Implementation of this function depends on the Jacobian, hence class F
      must define F::coord which specifies the coordinate system. Currently only the standard (r1,r2,r12) and Hylleraas (r1+r2,r1-r2,r12)
      coordinates are known. */
  template <typename F>
    double S(const F& bra, const F& ket) {
    const double result = I(apply_jacobian(bra*ket));
    return result;
  }

  /** Computes the normalization constant. Uses S(). */
  template <typename F>
    double NormConst(const F& bf) {
    const double selfoverlap = I(apply_jacobian(bf*bf));
    const double result = 1.0/std::sqrt(selfoverlap);
    return result;
  }

  /** Matrix element of the r_1^i r_2^j r_{12}^k operator over unnormalized functions. Can only work for functions expressed in terms of r1, r2, and r12.
      Must supply "F gen_r1r2r12_oper(i,j,k)". Uses S(). */
  template <typename F>
    double gen_mult_oper(int i, int j, int k, const F& bra, const F& ket) {
    const F oper(gen_r1r2r12_oper(i,j,k));
    const double result = S(bra,oper*ket);
    return result;
  }
  /** Matrix element of the electron repulsion operator (r_{12}^{-1}) over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double V_ee(const F& bra, const F& ket) {
    const double result = gen_mult_oper(0,0,-1,bra,ket);
    return result;
  }
  /** Matrix element of the nuclear attraction operator (- r_1^{-1} - r_2^{-1}) over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double V_ne(const F& bra, const F& ket) {
    double result = gen_mult_oper(-1,0,0,bra,ket);
    result += gen_mult_oper(0,-1,0,bra,ket);
    result *= -1.0;
    return result;
  }
  /** Matrix element of the kinetic energy operator over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double T(const F& bra, const F& ket) {
    return 0.0;
  }

};

#endif
