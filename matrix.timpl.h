
#ifndef _hyller_matrixtimpl_h_
#define _hyller_matrixtimpl_h_

namespace hyller {

  /** Overlap between 2 unnormalized functions. Implementation of this function depends on the Jacobian, hence class F
      must define F::coord which specifies the coordinate system. Currently only the standard (r1,r2,r12) and Hylleraas (r1+r2,r1-r2,r12)
      coordinates are known. */
  template <typename F>
    double S(const F& bra, const F& ket) {
    const double result = I(apply_jacobian(bra*ket)) * (8.0 * M_PI * M_PI);
    return result;
  }

  /** Computes the normalization constant. Uses S(). */
  template <typename F>
    double NormConst(const F& bf) {
    const double selfoverlap = I(apply_jacobian(bf*bf)) * (8.0 * M_PI * M_PI);
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

  /** Matrix element of the kinetic energy operator over unnormalized functions expressed in (r1,r2,r12) coordinates */
  template <typename F>
    double T_R1R2R12(const F& bra, const F& ket) {

    const F _m1_0_0(gen_r1r2r12_oper(-1,0,0));
    const F _m2_0_0(gen_r1r2r12_oper(-2,0,0));
    const F _0_m1_0(gen_r1r2r12_oper(0,-1,0));
    const F _0_m2_0(gen_r1r2r12_oper(0,-2,0));
    const F _0_0_m1(gen_r1r2r12_oper(0,0,-1));
    const F _0_0_m2(gen_r1r2r12_oper(0,0,-2));
    
    double T1 = 0.0;
    double T2 = 0.0;

    T1 -= ket.alpha * S(bra,_m1_0_0*ket);
    if (ket.i) {
      T1 += ket.i*S(bra,_m2_0_0*ket);
    }
    T1 -= ket.beta * S(bra,_0_m1_0*ket);
    if (ket.j) {
      T1 += ket.j*S(bra,_0_m2_0*ket);
    }
    T1 -= ket.gamma * S(bra,_0_0_m1*ket);
    if (ket.k) {
      T1 += ket.k*S(bra,_0_0_m2*ket);
    }

    T1 *= -2.0;

    
    T2 -= (2.0*ket.alpha*ket.alpha + 2.0*ket.beta*ket.beta + ket.gamma*ket.gamma) * S(bra,ket);
    if (ket.i) {
      T2 += 4.0*ket.i*ket.alpha * S(bra,_m1_0_0*ket);
      if (ket.i > 1) {
	T2 -= 2.0*ket.i*(ket.i-1) * S(bra,_m2_0_0*ket);
      }
    }
    if (ket.j) {
      T2 += 4.0*ket.j*ket.beta * S(bra,_0_m1_0*ket);
      if (ket.j > 1) {
	T2 -= 2.0*ket.j*(ket.j-1) * S(bra,_0_m2_0*ket);
      }
    }
    if (ket.k) {
      T2 += 2.0*ket.k*ket.gamma * S(bra,_0_0_m1*ket);
      if (ket.k > 1) {
	T2 -= ket.k*(ket.k-1) * S(bra,_0_0_m2*ket);
      }
    }

    T2 -= ket.alpha*ket.gamma * (gen_mult_oper(-1,0,1,bra,ket) - gen_mult_oper(-1,2,-1,bra,ket));
    T2 -= ket.beta*ket.gamma * (gen_mult_oper(0,-1,1,bra,ket) - gen_mult_oper(2,-1,-1,bra,ket));

    if (ket.i) {
      T2 += ket.i * ket.gamma * (gen_mult_oper(-2,0,1,bra,ket) - gen_mult_oper(-3,2,1,bra,ket));
    }
    if (ket.k) {
      T2 += ket.k * ket.alpha * (gen_mult_oper(-1,0,0,bra,ket) - gen_mult_oper(-1,2,-2,bra,ket));
    }
    if (ket.i && ket.k) {
      T2 -= ket.i * ket.k * (gen_mult_oper(-2,0,0,bra,ket) - gen_mult_oper(-2,2,-2,bra,ket));
    }

    if (ket.j) {
      T2 += ket.j * ket.gamma * (gen_mult_oper(0,-2,1,bra,ket) - gen_mult_oper(2,-3,1,bra,ket));
    }
    if (ket.k) {
      T2 += ket.k * ket.beta * (gen_mult_oper(0,-1,0,bra,ket) - gen_mult_oper(2,-1,-2,bra,ket));
    }
    if (ket.j && ket.k) {
      T2 -= ket.j * ket.k * (gen_mult_oper(0,-2,0,bra,ket) - gen_mult_oper(2,-2,-2,bra,ket));
    }

    const double result = T1 + T2;
    return result;
  }

  /** Matrix elements of the kinetic energy in Hylleraas coordinates */
  template <typename F>
    double T_Hylleraas(const F& bra, const F& ket);

  /** Matrix element of the kinetic energy operator over unnormalized functions. Uses S() and other functions. */
  template <typename F>
    double T(const F& bra, const F& ket) {
    if (F::Coords == Coordinates::R1R2R12) {
      return T_R1R2R12(bra,ket);
    }
    if (F::Coords == Coordinates::Hylleraas) {
      return T_Hylleraas(bra,ket);
    }
  }

};

#endif
