
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
  /** Matrix element of the r_1^i r_2^j r_{12}^k e^{- alpha * r_1 - beta * r_2 - gamma * r_{12} }
      operator over unnormalized functions. Can only work for functions expressed in terms of r1, r2, and r12.
      Must supply operator*(const GenSlaterHylleraasBasisFunction&, const F&).
      Uses S(). */
  template <typename F>
    double gen_mult_oper(int i, int j, int k,
			 double alpha, double beta, double gamma,
			 const F& bra, const F& ket) {
    const GenSlaterHylleraasBasisFunction oper(i,j,k,alpha,beta,gamma);
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

  /** Matrix element of the kinetic energy operator over unnormalized functions expressed in (r1,r2,r12) coordinates. See IJQC 101, 246 (2005). */
  template <typename F>
    double T_R1R2R12(const F& bra, const F& ket) {

    double T1 = 0.0;
    double T2 = 0.0;
    const double cij = 2.0;

#if 1
    T1 -= ket.alpha * gen_mult_oper(-1,0,0,bra,ket);
    if (ket.i) {
      T1 += ket.i*gen_mult_oper(-2,0,0,bra,ket);
    }
    T1 -= ket.beta * gen_mult_oper(0,-1,0,bra,ket);
    if (ket.j) {
      T1 += ket.j*gen_mult_oper(0,-2,0,bra,ket);
    }
    T1 -= cij * ket.gamma * gen_mult_oper(0,0,-1,bra,ket);
    if (ket.k) {
      T1 += cij * ket.k*gen_mult_oper(0,0,-2,bra,ket);
    }
    T1 *= -1.0;
#endif

#if 1
    T2 += -0.5*(ket.alpha*ket.alpha + ket.beta*ket.beta + cij*ket.gamma*ket.gamma) * S(bra,ket);
    if (ket.i) {
      T2 += ket.i*ket.alpha * gen_mult_oper(-1,0,0,bra,ket);
      if (ket.i > 1) {
	T2 += -0.5*ket.i*(ket.i-1) * gen_mult_oper(-2,0,0,bra,ket);
      }
    }
    if (ket.j) {
      T2 += ket.j*ket.beta * gen_mult_oper(0,-1,0,bra,ket);
      if (ket.j > 1) {
	T2 += -0.5*ket.j*(ket.j-1) * gen_mult_oper(0,-2,0,bra,ket);
      }
    }
    if (ket.k) {
      T2 += cij * ket.k*ket.gamma * gen_mult_oper(0,0,-1,bra,ket);
      if (ket.k > 1) {
	T2 += -0.5* cij * ket.k*(ket.k-1) * gen_mult_oper(0,0,-2,bra,ket);
      }
    }
#endif

#define SYMMETRIC_MIXED_T2 1
#if SYMMETRIC_MIXED_T2
    const double c1 = -0.5;
#else
    const double c1 = -1.0;
#endif

#if 1
    T2 += c1* ket.alpha*ket.gamma * (gen_mult_oper(-1,0,1,bra,ket) + gen_mult_oper(1,0,-1,bra,ket) - gen_mult_oper(-1,2,-1,bra,ket));
#if SYMMETRIC_MIXED_T2
    T2 += c1* ket.beta*ket.gamma  * (gen_mult_oper(0,-1,1,bra,ket) + gen_mult_oper(0,1,-1,bra,ket) - gen_mult_oper(2,-1,-1,bra,ket));
#endif
#endif

#if 1
    if (ket.i) {
      T2 -= c1* ket.i*ket.gamma * (gen_mult_oper(-2,0,1,bra,ket) + gen_mult_oper(0,0,-1,bra,ket) - gen_mult_oper(-2,2,-1,bra,ket));
    }
    if (ket.k) {
      T2 -= c1* ket.k*ket.alpha * (gen_mult_oper(-1,0,0,bra,ket) + gen_mult_oper(1,0,-2,bra,ket) - gen_mult_oper(-1,2,-2,bra,ket));
    }
    if (ket.i && ket.k) {
      T2 += c1* ket.i*ket.k * (gen_mult_oper(-2,0,0,bra,ket) + gen_mult_oper(0,0,-2,bra,ket) - gen_mult_oper(-2,2,-2,bra,ket));
    }

#if SYMMETRIC_MIXED_T2
    if (ket.j) {
      T2 -= c1* ket.j*ket.gamma * (gen_mult_oper(0,-2,1,bra,ket) + gen_mult_oper(0,0,-1,bra,ket) - gen_mult_oper(2,-2,-1,bra,ket));
    }
    if (ket.k) {
      T2 -= c1* ket.k*ket.beta * (gen_mult_oper(0,-1,0,bra,ket) + gen_mult_oper(0,1,-2,bra,ket) - gen_mult_oper(2,-1,-2,bra,ket));
    }
    if (ket.j && ket.k) {
      T2 += c1* ket.j*ket.k * (gen_mult_oper(0,-2,0,bra,ket) + gen_mult_oper(0,0,-2,bra,ket) - gen_mult_oper(2,-2,-2,bra,ket));
    }
#endif
#endif

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

  template <typename O, typename G>
    double S1(const O& bra, const G& ket) {

    G braG(bra^(O::Identity));
    return S(braG,ket);
  }

  template <typename O, typename G>
    double S2(const O& bra, const G& ket) {

    G braG((O::Identity)^bra);
    return S(braG,ket);
  }

};

#endif
