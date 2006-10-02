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
  double I(double z, int i) {
    double integral = fac(i)/pow(z,i+1);
    return integral;
  }
};

// normalization constant
double normConst(const Orbital& bf)
{
  double N = 1.0/sqrt( I(2.0*bf.zeta,2*bf.n+2) );
  return N;
}

// overlap of two one-electron basis functions
double Overlap(const Orbital& bfi, const Orbital& bfj)
{
  const double zeta = bfi.zeta;
  double overlap = normConst(bfi) * normConst(bfj) * I(2.0*zeta,bfi.n+bfj.n+2);
  return overlap;
}

// nuclear attraction over two one-electron basis functions
double V_en(const Orbital& bfi, const Orbital& bfj)
{
  const double zeta = bfi.zeta;
  double v_ne = normConst(bfi) * normConst(bfj) * I(2.0*zeta,bfi.n+bfj.n+1);
  return v_ne;
}

// kinetic energy over two one-electron basis functions
double T(const Orbital& bfi, const Orbital& bfj)
{
  const double zeta = bfi.zeta;
  const double twoz = 2.0 * zeta;
  const int ipj = bfi.n + bfj.n;
  double t = zeta*zeta*I(twoz,ipj+2) - 2.0*(bfj.n+1)*zeta*I(twoz,ipj+1) + (bfj.n+1)*bfj.n*I(twoz,ipj);
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

/////////////

namespace {

  // computes integral of function r1^i r2^j r12^k exp(-zeta1*r1) exp(-zeta2*r2)
  double I(double zeta1, double zeta2, int i, int j, int k)
  {
    const double sixteenpi2 = 16.0 * M_PI * M_PI;
    const int ki = k+i;
    const int kij = k+i+j;
    const double zeta12 = zeta1+zeta2;

    double result = 0.0;
    for(int p=1; p<=k+2; p+=2) {
      const double pfac = binomial(k+2,p);
      const int jp = j+p+1;
      const int kip = ki-p+3;
      const double term1 = fac(jp+1)*fac(ki+3-p)/(pow(zeta1,ki+4-p)*pow(zeta2,jp+2));

      double term2 = 0.0;
      for(int l=0; l<=jp; l++) {
	term2 += fac(kij+4-l)/(pow(zeta12,kij+5-l) * fac(jp-l) * pow(zeta2,l+1));
      }
      term2 *= fac(jp);

      double term3 = 0.0;
      for(int l=0; l<=kip; l++) {
	term3 += fac(kij+4-l) / (pow(zeta12,kij+5-l) * fac(kip-l) * pow(zeta2,l+1));
      }
      term3 *= fac(kip);

      result += pfac * (term1 - term2 + term3);
    }

    result *= sixteenpi2 / (k+2);
    return result;
  }
};

double normConst(const GenHylleraasBasisFunction& bfi)
{
  const double Sijij = I(2.0*bfi.zeta1,
			 2.0*bfi.zeta2,
			 2*bfi.i,
			 2*bfi.j,
			 2*bfi.k);
  const double Sijji = I(bfi.zeta1+bfi.zeta2,
			 bfi.zeta1+bfi.zeta2,
			 bfi.i+bfi.j,
			 bfi.i+bfi.j,
			 2*bfi.k);
  const double S = (bfi.spin == SpinSinglet) ? 2.0*(Sijij+Sijji) : 2.0*(Sijij-Sijji);
  return 1.0/std::sqrt(S);
}

double
Overlap(const GenHylleraasBasisFunction& bfi,
	const GenHylleraasBasisFunction& bfj)
{
  if (bfi.spin != bfj.spin) 
    return 0.0;
  const double Sijij = I(bfi.zeta1+bfj.zeta1,
			 bfi.zeta2+bfj.zeta2,
			 bfi.i+bfj.i,
			 bfi.j+bfj.j,
			 bfi.k+bfj.k);
  const double Sijji = I(bfi.zeta1+bfj.zeta2,
			 bfi.zeta2+bfj.zeta1,
			 bfi.i+bfj.j,
			 bfi.j+bfj.i,
			 bfi.k+bfj.k);
  const double S = (bfi.spin == SpinSinglet) ? 2.0*(Sijij+Sijji) : 2.0*(Sijij-Sijji);
  const double norm_pfac = normConst(bfi)*normConst(bfj);
  return S * norm_pfac;
}

double
V_en(const GenHylleraasBasisFunction& bfi,
     const GenHylleraasBasisFunction& bfj)
{
  if (bfi.spin != bfj.spin)
    return 0.0;
  double V = 0.0;
  const GenHylleraasBasisFunction* bf1;
  const GenHylleraasBasisFunction* bf2;
  if (bfi.i < bfj.i) {
    bf1 = &bfi; bf2 = &bfj;
  }
  else {
    bf1 = &bfj; bf2 = &bfi;
  }
  {
    GenHylleraasBasisFunction bf22(bf2->spin,bf2->i-1,bf2->j,bf2->k,bf2->zeta1,bf2->zeta2);
    V += normConst(*bf2) * Overlap(*bf1,bf22) / normConst(bf22);
  }
  if (bfi.j < bfj.j) {
    bf1 = &bfi; bf2 = &bfj;
  }
  else {
    bf1 = &bfj; bf2 = &bfi;
  }
  {
    GenHylleraasBasisFunction bf22(bf2->spin,bf2->i,bf2->j-1,bf2->k,bf2->zeta1,bf2->zeta2);
    V += normConst(*bf2) * Overlap(*bf1,bf22) / normConst(bf22);
  }
  return V;
}

double
V_ee(const GenHylleraasBasisFunction& bfi,
     const GenHylleraasBasisFunction& bfj)
{
  if (bfi.spin != bfj.spin)
    return 0.0;
  const GenHylleraasBasisFunction* bf1;
  const GenHylleraasBasisFunction* bf2;
  if (bfi.k < bfj.k) {
    bf1 = &bfi; bf2 = &bfj;
  }
  else {
    bf1 = &bfj; bf2 = &bfi;
  }
  GenHylleraasBasisFunction bf22(bf2->spin,bf2->i,bf2->j,bf2->k-1,bf2->zeta1,bf2->zeta2);
  return normConst(*bf2) * Overlap(*bf1,bf22) / normConst(bf22);
}

/////////////

/** Computes integral of function r1^i r2^j r12^k exp(-alpha*r1 - beta*r2 - gamma*r12) over r1, r2, and r12. Jacobian is not included. See J. Chem. Phys. 121, 6323 (2004). */
template<>
  double I<GenSlaterHylleraasBasisFunction>(const GenSlaterHylleraasBasisFunction& bf)
{
  const int i = bf.i;
  const int j = bf.j;
  const int k = bf.k;
  const double alpha = bf.alpha;
  const double beta = bf.beta;
  const double gamma = bf.gamma;

  const double pfac = 2.0 * fac(i) * fac(j) * fac(k);
  const double ab = alpha+beta;
  const double ag = alpha+gamma;
  const double bg = beta+gamma;
  
  double sum = 0.0;
  for(int l=0; l<=i; ++l) {
    for(int m=0; m<=j; ++m) {
      const int jml = j - m + l;
      for(int n=0; n<=k; ++n) {
	const int iln = i - l + n;
	const int knm = k - n + m;
	
	const double num = binomial(jml,l) * binomial(iln,n) * binomial(knm,m);
	const double denom = pow(ab,jml+1) * pow(ag,iln+1) * pow(bg,knm+1);
	sum += num/denom;
      }
    }
  }
  
  const double result = sum * pfac;
  return result;
}

};  // namespace hyller
