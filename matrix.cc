/* This file contains definitions for matrix elements and other 
   miscellaneous functions */
#include "includes.h"
#include "matrix.h"
#include "misc.h"

extern "C" FILE *outfile;

namespace hyller {

/* Note: Normalization constant does include factor 1/sqrt(8*Pi^2)
   arising from integration over angles. */
double normConst(const HylleraasBasisFunction& bfi)
{
  const double oosqrt8pi2 = M_SQRT2/(4.0*M_PI);
  double result = (2*bfi.l+1)*(2*bfi.l+3)*(2*bfi.l+2*bfi.m+3)*(2*bfi.l+2*bfi.m+5);
  result = sqrt(result/((2*bfi.l+bfi.m+3)*fac(2*bfi.nlm()+5)));
  result *= oosqrt8pi2;
  return result;
}

/*
   Matrix elements are defined as in Wesley Allen's homework, except for
   the angular part of the Jacobian is also included (see 8pi^2). The
   diagonal elements of S therefore equal 1.
 */
double Overlap(const HylleraasBasisFunction& bfi,
	       const HylleraasBasisFunction& bfj)
{
  const double eightpi2 = 8.0*M_PI*M_PI;
  double result = (0.25*(1+pow(-1,bfi.l+bfj.l))*(bfi.m+bfj.m+2*(bfi.l+bfj.l)+6)*fac(bfi.nlm()+bfj.nlm()+5)*
		   normConst(bfi)*normConst(bfj));
  result = result/((bfi.l+bfj.l+1)*(bfi.l+bfj.l+3)*(bfi.l+bfj.l+bfi.m+bfj.m+3)*(bfi.l+bfj.l+bfi.m+bfj.m+5));
  result *= eightpi2;
  return result;
}

double V_en(const HylleraasBasisFunction& bfi,
	    const HylleraasBasisFunction& bfj)
{
  const double zeta = bfi.zeta;
  const double eightpi2 = 8.0*M_PI*M_PI;
  double result = -zeta*(1+pow(-1,bfi.l+bfj.l))*normConst(bfi)*normConst(bfj)*
    fac(bfi.nlm()+bfj.nlm()+4);
  result = result/((bfi.l+bfj.l+1)*(bfi.l+bfj.l+bfi.m+bfj.m+3));
  result *= eightpi2;
  return result;
}


double V_ee(const HylleraasBasisFunction& bfi,
	    const HylleraasBasisFunction& bfj)
{
  const double zeta = bfi.zeta;
  const double eightpi2 = 8.0*M_PI*M_PI;
  double result = 0.25*zeta*(1+pow(-1,bfi.l+bfj.l))*(bfi.m+bfj.m+2*bfi.l+2*bfj.l+5)*fac(bfi.nlm()+bfj.nlm()+4)*
    normConst(bfi)*normConst(bfj);
  result = result/((bfi.l+bfj.l+1)*(bfi.l+bfj.l+3)*(bfi.l+bfj.l+bfi.m+bfj.m+2)*(bfi.l+bfj.l+bfi.m+bfj.m+4));
  result *= eightpi2;
  return result;
}


double T(const HylleraasBasisFunction& bfi,
	 const HylleraasBasisFunction& bfj)
{
  const double zeta = bfi.zeta;
  const double eightpi2 = 8.0*M_PI*M_PI;
  double nom1,den1,nom2,den2,nom3,den3,nom4,den4,nom5,den5;
  double result,tmp;
  double tmp1,tmp2,tmp3,tmp4,tmp5;
  result = 0.125*zeta*zeta*(1 + pow(-1,bfi.l + bfj.l))*fac(bfi.nlm()+bfj.nlm()+3);
  result *= normConst(bfi)*normConst(bfj);
  /* first term in the big brackets : nominator at first, then denominator */
  nom1 = (4*(bfi.n+bfi.l+2*bfi.m+1)*(bfi.l-bfi.n+2)+(bfi.nlm()+bfj.nlm())*(4*bfi.n+4*bfi.m-(bfi.nlm()+bfj.nlm())-1)+4);
  den1 = (4*(bfi.l+bfj.l+1)*(bfi.l+bfj.l+bfi.m+bfj.m+3));
  tmp1 = nom1 / den1;

  /* .. second and third terms .. */
  nom2 = (bfi.m*(bfi.m+2*bfi.l+1));
  den2 = ((bfi.l+bfj.l+1)*(bfi.l+bfj.l+bfi.m+bfj.m+1));
  tmp2 = nom2 / den2;
  nom3 = (bfi.l*(bfi.l-1));
  den3 = ((bfi.l+bfj.l-0.5*(1+pow(-1,bfi.l+bfj.l)))*(bfi.l+bfj.l+bfi.m+bfj.m+1));
  tmp3 = nom3 / den3;

  /* .. nominator of the fourth term, then its denominator .. */
  nom4 = (4*bfi.n*(bfi.n-(bfi.nlm()+bfj.nlm())-5)+(bfi.nlm()+bfj.nlm()+4)*(bfi.nlm()+bfj.nlm()+5));
  den4 = (4*(bfi.l+bfj.l+3)*(bfi.l+bfj.l+bfi.m+bfj.m+5));
  tmp4 = nom4 / den4;

  /* .. the fifth term .. */
  nom5 = (bfi.m*(bfi.m+2*bfi.n-(bfi.nlm()+bfj.nlm())-3));
  den5 = ((bfi.l+bfj.l+3)*(bfi.l+bfj.l+bfi.m+bfj.m+3));
  tmp5 = nom5 / den5;
  tmp = tmp1 - tmp2 - tmp3 + tmp4 + tmp5;
  /* printf("Pre-bracket term = %lf\n",result);
     printf("First nominator = %lf\n",nom1);
     printf("First denominator = %lf\n",den1);
     printf("Second nominator = %lf\n",nom2);
     printf("Second denominator = %lf\n",den2);
     printf("Third nominator = %lf\n",nom3);
     printf("Third denominator = %lf\n",den3);
     printf("Fourth nominator = %lf\n",nom4);
     printf("Fourth denominator = %lf\n",den4);
     printf("Fifth nominator = %lf\n",nom5);
     printf("Fifth denominator = %lf\n",den5); */

  result = result*tmp;
  result *= eightpi2;
  return result;
}

/* This function computes Integral of phi(n,l,m)*phi(n1,l1,m1)
   from 0 to Infinity over R2 */
double density_Nuc(const HylleraasBasisFunction& bfi,
		   const HylleraasBasisFunction& bfj)
{
  const double zeta = bfi.zeta;
  return 8*M_PI*normConst(bfi)*normConst(bfj)*pow(zeta,3)*pow(-1,bfi.l+bfj.l)*
    fac(bfi.nlm()+bfj.nlm()+2);
}

////////////////////

// evaluates elementary integral \int_0^\infty dr r^i e^(-zr)
namespace {
  double S(double z, int i) {
    double s = fac(i)/pow(z,i+1);
    return s;
  }
};

// normalization constant
double normConst(const Orbital& bf)
{
  double N = 1.0/sqrt( S(2.0*bf.zeta,2*bf.n+2) );
  return N;
}

// overlap of two one-electron basis functions
double Overlap(const Orbital& bfi, const Orbital& bfj)
{
  const double zeta = bfi.zeta;
  double overlap = normConst(bfi) * normConst(bfj) * S(2.0*zeta,bfi.n+bfj.n+2);
  return overlap;
}

// nuclear attraction over two one-electron basis functions
double V_en(const Orbital& bfi, const Orbital& bfj)
{
  const double zeta = bfi.zeta;
  double v_ne = normConst(bfi) * normConst(bfj) * S(2.0*zeta,bfi.n+bfj.n+1);
  return v_ne;
}

// kinetic energy over two one-electron basis functions
double T(const Orbital& bfi, const Orbital& bfj)
{
  const double zeta = bfi.zeta;
  const double twoz = 2.0 * zeta;
  const int ipj = bfi.n + bfj.n;
  double t = zeta*zeta*S(twoz,ipj+2) - 2.0*(bfj.n+1)*zeta*S(twoz,ipj+1) + (bfj.n+1)*bfj.n*S(twoz,ipj);
  t *= -0.5 * normConst(bfi) * normConst(bfj);
  return t;
}

////

double Overlap(const SD& bfi, const SD& bfj)
{
  if (bfi.spin != bfj.spin) {
    return 0.0;
  }
  else {
    double S11 = Overlap(*bfi.o1,*bfj.o1);
    double S22 = Overlap(*bfi.o2,*bfj.o2);
    double S = S11*S22;
    if (samespin(bfi.spin)) {
      double S12 = Overlap(*bfi.o1,*bfj.o2);
      double S21 = Overlap(*bfi.o2,*bfj.o1);
      S -= S12*S21;
    }
    return S;
  }
}

////

double Overlap(const CSF& bfi, const CSF& bfj)
{
  if (bfi.spin != bfj.spin) {
    return 0.0;
  }
  else {
    double S = 0.0;
    for(int di=0; di<bfi.nd; ++di) {
      for(int dj=0; dj<bfj.nd; ++dj) {
	S += bfi.c[di] * bfj.c[dj] * Overlap(*bfi.d[di],*bfj.d[dj]);
      }
    }

    return S;
  }
}

};
