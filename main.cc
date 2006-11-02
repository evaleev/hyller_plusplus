/************************************************************
*                                                           *
*   Hylleraas treatment of two-electron system (S-states)   *
*                                                           *
*************************************************************/

#include <iostream>
#include <stdexcept>

#include "defines.h"
#include "includes.h"
#include "polynom.h"
#include "misc.h"
#include "hylleraas.h"
#include "slaterhylleraas.h"
#include "projector.h"

#include <matrix.h>
#include <matrix.timpl.h>
#include <wfn.h>
#include <solve.h>
#include <basis.h>
#include <gsh_basissets.h>
#include <basisfn.h>
#include <hamiltonian.h>
#include <energy.h>
#include <optimizer.h>

using namespace hyller;

/* Note: Normalization constant does NOT include factor 1/sqrt(8*Pi) 
         arising from integration over angles */

// These must be global -- required by Psi3
extern "C" {
  FILE *infile,*outfile;
  char *psi_file_prefix;
}

//
// Local functions
//
namespace {
  void compute_density_at_nucleus(const HylleraasWfn& wfn);
  void plot_angular_profile(const HylleraasWfn& wfn, int npoints, double radius, bool psinorm);
  void plot_pair_distribution(const HylleraasWfn& wfn, int npoints, double r12max);
  HylleraasWfn& solvCi(HylleraasBasisSet& basis, double Z, int root_num, bool opt, double threshold, int maxiter);
  OrbitalWfn& solvHF(OrbitalBasisSet& basis, HylleraasBasisSet& hbs, double Z, double threshold, int maxiter, bool opt, int maxiter_zeta);
  void eigensolve(int nbf, int indDim, double** H, double** X, int root, double& E, std::vector<double>& coeff);
  void usage();

  void test_gsh();
  void test_gen_basisfn();
  void test_gen_basis();
  void test_hf_2body(const OrbitalWfn& hfwfn);
  void test_gen_energy();
}

int main(int argc, char **argv)
{
  
  /*------------------------------------------------------------------------
  Declarations
  -------------------------------------------------------------------------*/
  int mult = 1;    // do singlet by default
  double Z = 2.0;  // do helium by default
  double default_zeta1S = 2*1.8149; /* default value for zeta/2 */
  double default_zeta3S = 2*1.1351; /* default value for zeta/2 */
  double zeta;
  int root = 0;           /* number of the root we're shooting for */
  int opt = 0;
  const double threshold = 1.0e-12;  // how accurately to solve equations
  int maxiter = 100;  // maximum number of iterations
  int nlm_min = -1;       /* minimum value of the sum of indices n,l and m 
                             -1 in Kinoshita method or 0 in Hylleraas method */
  int nlm_max;   /* maximum n+l+m */
  int n_max;       /* maximum value of n */
  int l_max;
  int m_max;
  

  int errcode;     /* error code used in parsing section */
  char *string;
  
  
  /*-------------------------------------------------------------------------
    Parsing section
   --------------------------------------------------------------------------*/
  int num_extra_args = 0;
  char **extra_args;
  int nmax_given = 0,
    lmax_given = 0,
    mmax_given = 0,
    nlmmax_given = 0;
  if (argc < 3) {
    usage();
    exit(1);
  }
  for (int i=1; i<argc; i++) {
    /*--- read geometry from geom.dat file (in findif calculations) ---*/
    if (strcmp(argv[i], "--max") == 0) {
      nlm_max = atoi(argv[i+1]);  ++i;
      nlmmax_given = 1;
    }
    else if (strcmp(argv[i], "--nmax") == 0) {
      n_max = atoi(argv[i+1]);  ++i;
      nmax_given = 1;
    }
    else if (strcmp(argv[i], "--lmax") == 0) {
      l_max = atoi(argv[i+1]);  ++i;
      lmax_given = 1;
    }
    else if (strcmp(argv[i], "--mmax") == 0) {
      m_max = atoi(argv[i+1]);  ++i;
      mmax_given = 1;
    }
    else {
      extra_args[num_extra_args++] = argv[i];
    }
  }
  if (!nlmmax_given) { usage(); exit(1); }
  if (!nmax_given) { n_max = nlm_max; }
  if (!lmax_given) { l_max = nlm_max; }
  if (!mmax_given) { m_max = nlm_max; }
  
  errcode = psi_start(num_extra_args, extra_args, 0);
  if (errcode != PSI_RETURN_SUCCESS) {
    fprintf(stderr,"init_io -- psi_start failed");
    exit(1);
  }
  
  ip_cwk_add(":HYLLERAAS");
  tstart(outfile);
  fprintf(outfile,"-----------------------------------------------------------------------------\n\n");
  fprintf(outfile,"\t\t\tHYLLERAAS-TYPE CI EXPANSION\n\n");
  errcode = ip_data("MULT","%d",&mult,0);
  if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK)) {
    fprintf(stderr,"\tMultiplicity not specified or there's some other problem\n\n");
    exit(1);
  }
  errcode = ip_data("CHARGE","%lf",&Z,0);
  errcode = ip_string("CALC",&string,0);
  if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK) || !strcmp(string,"KINOSHITA"))
    nlm_min = -1;
  else {
    if (!strcmp(string,"HYLLERAAS"))
      nlm_min = 0;
    else
      fprintf(stderr,
	      "\tCalculation type %s is not recognized. Doing Kinoshita calculation\n",string);
  }
  errcode = ip_boolean("ZOPT",&opt,0); /* If there's some problem with parsing,
  there's no optimization by default */
  
  errcode = ip_data("ROOT","%d",&root,0);
  if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK)) {
    fprintf(stderr,"\tRoot number not found or there's some other problem\n\n");
    exit(1);
  }
  
  errcode = ip_data("ZETA","%lf",&zeta,0);
  if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK)) {
    zeta = (mult == 1) ? default_zeta1S : default_zeta3S;
  }
  
  errcode = ip_data("MAXITER","%d",&maxiter,0);
  
  
  /*---------------------------------------------------------------------------
    Compute Hylleraas wave function
   ---------------------------------------------------------------------------*/
  HylleraasBasisSet hyllbs(mult == 1, nlm_max, nlm_min, n_max, l_max, m_max, zeta);
  const int num_bf = hyllbs.num_bf();
  if (num_bf == 0)
    throw std::runtime_error("main -- no basis functions found which fit the constraints");
  if (num_bf < root)
    throw std::runtime_error("main -- root > num_bf");
  fprintf(outfile,"\tZ = %5.3lf\tn+l+m(max) = %d\tn+l+m(min) = %d\tnmax = %d\tlmax = %d\tmmax = %d\tN = %d\n",Z,nlm_max,nlm_min,n_max,l_max,m_max,hyllbs.num_bf());
  if (nlm_min == 0)
    fprintf(outfile,"\tPerforming Hylleraas-type calculation\n");
  else
    fprintf(outfile,"\tPerforming Kinoshita-type calculation\n");
  fflush(outfile);
  
#if !SKIP_HYLLERAAS
  const HylleraasWfn& wfn = solvCi(hyllbs,Z,root,opt,threshold,maxiter);

  /*----------------------------------------------------------------------------
    Compute properties
   ----------------------------------------------------------------------------*/

  compute_density_at_nucleus(wfn);

  int cusp_plot = 0;
  errcode = ip_boolean("CUSP_PLOT",&cusp_plot,0);
  if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK))
    cusp_plot = 0;  /* default - no plots */
  else {
    double radius;
    errcode = ip_data("RADIUS","%lf",&radius,0);
    if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK))
      radius = 1.0;
    if (radius < 0.0) {
      fprintf(stderr,"radius must be nonnegative");
      exit(1);
    }
    int npoints = 0;
    errcode = ip_data("NP_CUSP","%d",&npoints,0);
    if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK) || npoints <= 0)
      npoints = 100;
    int psinorm = 0;
    errcode = ip_boolean("PSINORM",&psinorm,0);
    if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK))
      psinorm = 0;

    plot_angular_profile(wfn, npoints, radius, psinorm);
  }

  int pd_plot = 0;
  errcode = ip_boolean("PD_PLOT",&pd_plot,0);
  if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK))
    pd_plot = 0;  /* default - no plots */
  else
  {
    double r12max = 0.0;
    errcode = ip_data("R12MAX","%lf",&r12max,0);
    if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK))
      r12max = 3.0;
    int npoints = 0;
    errcode = ip_data("NP_PD","%d",&npoints,0);
    if ((errcode == IPE_KEY_NOT_FOUND) || (errcode != IPE_OK) || npoints <= 0)
      npoints = 100;

    plot_pair_distribution(wfn, npoints, r12max);
  }
#endif
  
  /*---------------------------------------------------------------------------
    Compute Hartree-Fock wave function
   ---------------------------------------------------------------------------*/
#if !SKIP_HF
  OrbitalBasisSet orbs(nlm_max/2, 0, zeta/2);
  const int norbs = orbs.num_bf();
  const OrbitalWfn& hfwfn = solvHF(orbs,hyllbs,Z,threshold,100,opt,maxiter);
#endif

#if 0
  SDBasisSet sdbs(SpinAlphaBeta,orbs);
  SD_2_Hylleraas sd2hyll_1(sdbs,hyllbs);

  // Do triplet as well
  HylleraasBasisSet hyllbs3(false, nlm_max, nlm_min, n_max, l_max, m_max, zeta);
  SD_2_Hylleraas sd2hyll_3(sdbs,hyllbs3);

  // Now do CSFs
  {
    //OrbitalBasisSet orbs(nlm_max, 0, zeta/2);
    //SDBasisSet sdbs(SpinAlphaBeta,orbs);
    CSFBasisSet csfbs_1(SpinSinglet,sdbs);
    CSF_2_Hylleraas csf2hyll_1(csfbs_1,hyllbs);
    CSFBasisSet csfbs_3(SpinTriplet,sdbs);
    CSF_2_Hylleraas csf2hyll_3(csfbs_3,hyllbs3);
  }
#endif

  /*----------------------
    Test various features
   ----------------------*/
#if !SKIP_TESTS
  test_gsh();
  {
    const double alpha = zeta/2;
    fprintf(outfile,"\t-Solving CI problem in generalized Slater-Hylleraas basis:\n");
    typedef GenSlaterHylleraasBasisSet Basis;
    typedef Wavefunction<Basis> Wfn;
    Basis bs(nlm_max,nlm_min,n_max,m_max,alpha,alpha,0.0);
    const Wfn& wfn = solveCi(bs,Z,root,false,threshold,maxiter);
  }

  // test generic contracted functions
  test_gen_basisfn();

  /** Compute the energy using the wave function with following terms:
      e^{-\zeta (r_1 + r_2)}
      e^{-\zeta (r_1 + r_2)} e^{-\gamma r_{12}}
      e^{-\zeta (r_1 + r_2)} (r_1^2 + r_2^2 - r_{12}^2)
  */
  test_gen_basis();

#if !SKIP_HF
  /** Re-Compute the HF energy using a 2-body functionality
  */
  test_hf_2body(hfwfn);
#endif

#endif // end of tests
  // Test generic energy evaluator
  test_gen_energy();

  tstop(outfile);
  ip_done();
  fclose(outfile);
  fclose(infile);
}                          /* end of main() */

extern "C" char *gprgid()
{
  char *prgid = "hyller";
  return(prgid);
}

namespace {

void plot_angular_profile(const HylleraasWfn& wfn, int npoints, double radius, bool psinorm)
{
  FILE *fpcusp;
  
  /*---------------------------------------------------------------------------
    Plotting wave function vs. r12 on a circle
   ---------------------------------------------------------------------------*/
  fprintf(outfile,"\tThe wavefunction near the cusp is plotted\n");
  fprintf(outfile,"\tRadius = %lf\n",radius);
  double psi00 = 0.0;
  const HylleraasBasisSet& basis = wfn.basis();
  const std::vector<double>& coeff = wfn.coefs();
  const int nbf = basis.num_bf();
  const double zeta = basis.zeta();
  for(int i=0;i<nbf;i++) {
    const HylleraasBasisFunction& bf = basis.bf(i);
    if (((bf.l == 0) && (bf.m == 0)) ||
	((bf.l == 1) && (bf.m == -1))
	)
      psi00 += normConst(bf)*
	pow(zeta,bf.nlm()+3)*
	exp(-2.0*radius*zeta*0.5)*pow(radius*2.0,bf.n)*coeff[i];
  }
  fprintf(outfile,"\tPsi(r12=0) = %11.10e\n",psi00);
  fpcusp = fopen("cusp_plot.out","w");
  fprintf(fpcusp,"%d %d %d %d %d %d %lf %25.15lf %s\n",npoints,basis.nlm_max(),basis.nlm_min(),
	  basis.n_max(),basis.l_max(),basis.m_max(),
	  radius,psi00,(psinorm ? "true" : "false"));
  if (psinorm)
    fprintf(fpcusp,"0.000000\t%25.15lf\n",1.0);
  else
    fprintf(fpcusp,"0.000000\t%25.15lf\n",psi00);
  for(int i=1;i<npoints;i++) {
    double u = 2*i*radius/npoints;
    double psi = 0.0;
    for(int j=0;j<nbf;j++) {
      const HylleraasBasisFunction& bf = basis.bf(j);
      if (bf.m == 0)
	psi += normConst(bf)*
	  pow(zeta,bf.nlm()+3)*
	  exp(-2.0*radius*zeta*0.5)*pow(radius*2.0,bf.n)*
	  pow(u,bf.m)*coeff[j];
    }
    if (psinorm)
      psi /= psi00;
    fprintf(fpcusp,"%lf\t%25.15lflf\n",u,psi);
  }
  fclose(fpcusp);
}

void
plot_pair_distribution(const HylleraasWfn& wfn, int npoints, double r12max)
{
  FILE* fppd = fopen("pd_plot.out","w");
  fprintf(fppd,"%d\t%lf\n",npoints,r12max);
  fprintf(outfile,"\tThe pair distribution function is plotted\n\n");

  const HylleraasBasisSet& basis = wfn.basis();
  const std::vector<double>& coeff = wfn.coefs();
  const int nbf = basis.num_bf();
  const double zeta = basis.zeta();
  Polynom* PD = init_Polynom(2*basis.nlm_max()+4);

  for(int i=0;i<nbf;i++) {
    const HylleraasBasisFunction& bfi = basis.bf(i);

    const double nc1 = normConst(bfi)*
      pow(zeta,bfi.nlm()+3);

    for(int j=0;j<nbf;j++) {
      const HylleraasBasisFunction& bfj = basis.bf(j);

      const double nc2 = normConst(bfj)*
	pow(zeta,bfj.nlm()+3);

      double tmp = 0.125*nc1*nc2*(1+pow(-1,bfi.l+bfj.l))*
	coeff[i]*coeff[j]/(bfi.l+bfj.l+1);
      int kmax = bfi.n + bfj.n + 2;
      for(int k=0;k<=kmax;k++)
	(*PD).body[bfi.nlm()+bfj.nlm()+3-k]
	  += facdfac(kmax,kmax-k)*
	  tmp/pow(zeta,k+1);

      tmp = 0.125*nc1*nc2*(1+pow(-1,bfi.l+bfj.l))*
	coeff[i]*coeff[j]/(bfi.l+bfj.l+3);
      kmax = bfi.n + bfj.n;
      for(int k=0;k<=kmax;k++)
	(*PD).body[bfi.nlm()+bfj.nlm()+3-k]
	  -= facdfac(kmax,kmax-k)*
	  tmp/pow(zeta,k+1);

    }
  }

  fprintf(fppd,"0.000000\t%25.15lf\n",0.0);
  for(int i=1;i<npoints;i++) {
    const double x = i*r12max/npoints;
    fprintf(fppd,"%lf\t%25.15lf\n",x,exp(-zeta*x)*eval_Polynom(PD,x));
  }

  destr_Polynom(PD);
  fclose(fppd);
}

void
compute_density_at_nucleus(const HylleraasWfn& wfn)
{
  const HylleraasBasisSet& basis = wfn.basis();
  const std::vector<double>& coeff = wfn.coefs();
  const int nbf = basis.num_bf();
  /*---------------------------------------------------------------------------
  Computing electron density at nucleus in Hylleraas-type calculation
  -----------------------------------------------------------------------------*/
  if (basis.nlm_min() == 0) {
    double nuc_dens = 0.0;
    for(int i=0;i<nbf;i++) {
      const HylleraasBasisFunction& bfi = basis.bf(i);
      for(int j=0;j<nbf;j++) {
	const HylleraasBasisFunction& bfj = basis.bf(j);
	const double tmp = density_Nuc(bfi,bfj);
	nuc_dens += coeff[i]*coeff[j]*tmp;
      }
    }
    fprintf(outfile,"\tDensity at the nucleus = %3.12lf\n",nuc_dens);
  }
}

HylleraasWfn&
solvCi(HylleraasBasisSet& basis, double Z, int root_num, bool opt, double threshold, int maxiter)
{
  double **H,**SS,**X,**temp,*evals;
  double E = 1.0;
  double E_old = 0.0;
  double Ed = 0.0;
  double Ed2 = 0.0;
  double tmp,nuc_dens,zeta,u,x;
  int iter = 1;
  int indDim = 0;

  double zeta_new = basis.zeta();
  const int nbf = basis.num_bf();
  std::vector<double> coeff(basis.num_bf());
 
  while(iter == 1 || fabs(E-E_old) > threshold) {
    E_old = E;
    zeta = zeta_new;
    if (opt) {
      fprintf(outfile,"\t Test message : iteration #%d\n",iter);
      fprintf(outfile,"\t Test message : zeta = %3.12lf\n",zeta);
      fflush(outfile);
    }
    if (iter == 1) {
      SS = block_matrix(nbf,nbf);
      for(int i=0;i<nbf;i++) {
	const HylleraasBasisFunction& bfi = basis.bf(i);
        for(int j=0;j<=i;j++) {
	  const HylleraasBasisFunction& bfj = basis.bf(j);
	  SS[i][j] = SS[j][i] = S(bfi,bfj);
	}
      }
      fprintf(outfile,"\t Overlap matrix\n");
      print_mat(SS,nbf,nbf,outfile);
      indDim = nbf - sq_sqrt(SS,&X,nbf);
      fprintf(outfile,"\t%d basis functions are linearly independent\n",indDim);
      free_block(SS);
      fflush(outfile);
    }
    H = block_matrix(nbf,nbf);
    for(int i=0;i<nbf;i++) {
      const HylleraasBasisFunction& bfi = basis.bf(i);
      for(int j=0;j<=i;j++) {
	const HylleraasBasisFunction& bfj = basis.bf(j);

	double Hij = T(bfi,bfj);
	
	Hij += Z*V_en(bfi,bfj)/2.0;
	Hij += V_ee(bfi,bfj);
	
	H[j][i] = H[i][j] = Hij;
      }
    }
    fprintf(outfile,"\tHamiltonian matrix\n");
    print_mat(H,nbf,nbf,outfile);
    
    temp = block_matrix(nbf,indDim);
    
    mmult(H,0,X,0,temp,0,nbf,nbf,indDim,0);
    free_block(H);
    H = block_matrix(indDim,indDim);
    mmult(X,1,temp,0,H,0,indDim,nbf,indDim,0);
    free_block(temp);
    temp = block_matrix(indDim,indDim);
    
    evals = init_array(indDim);
    sq_rsp(indDim,indDim,H,evals,1,temp,1.0E-20);
    free_block(H);
    /* Test - prints out all roots */
    for(int i=0;i<=root_num;i++)
      fprintf(outfile,"\tState #%d  E = %3.12lf\n",i+1,evals[i]);
    E = evals[root_num];
    
    for(int i=0;i<indDim;i++) /* temporarily putting eigenvector to evals array */
      evals[i] = temp[i][root_num];
    free_block(temp);
    
    for(int i=0;i<nbf;i++) {
      double c = 0.0;
      for(int j=0;j<indDim;j++)
        c += X[i][j]*evals[j];
      coeff[i] = c;
    }
    
    if (!opt)
      free_block(X);
    
    tmp = 0;
    for(int i=0;i<nbf;i++) {
      const HylleraasBasisFunction& bfi = basis.bf(i);
      for(int j=0;j<nbf;j++) {
	const HylleraasBasisFunction& bfj = basis.bf(j);
	tmp =  T(bfi,bfj);
	Ed += tmp*coeff[i]*coeff[j]; 
      }
    }
    Ed2 = Ed*2/(zeta*zeta);
    Ed = (Ed+E)/zeta;
    //zeta_new = zeta-100000*Ed/Ed2;
    zeta_new = zeta-Ed/Ed2;
    if (opt) {
      fprintf(outfile,"\tdE/d(zeta) = %3.12lf\n",Ed);
      fprintf(outfile,"\t Test message : new zeta1S/2 = %3.12lf\n",zeta_new/2);
      fprintf(outfile,"\t Test message : delta(E) = %3.12lf\n",E-E_old);
    }
    iter++;
    
    if (!opt)
      break;     /* if optimization flag is not set - quit 'while' loop */
    else
      fflush(outfile);
  }
  
#if TEST_NORMALIZATION
  /* verify c^T . S . c = 1 */
  {
    double norm = 0.0;
    for(int i=0;i<nbf;i++) {
      const HylleraasBasisFunction& bfi = basis.bf(i);
      for(int j=0;j<nbf;j++) {
	const HylleraasBasisFunction& bfj = basis.bf(j);
	norm += coeff[i] * coeff[j] *
	  S(bfi,bfj);
      }
    }
    fprintf(outfile,"\tTEST: c^T * S * c = %10.9lf (should be 1.0)\n",norm);
  }
#endif

  HylleraasWfn* wfn = new HylleraasWfn(basis,coeff);
  return *wfn;
}

OrbitalWfn&
solvHF(OrbitalBasisSet& basis, HylleraasBasisSet& hbs, double Z, double threshold, int maxiter, bool opt, int maxiter_zeta)
{
  double **Hc, **SS,**X;
  double **Hee; // electron repulsion operator in basis of slater determinants
  double **F;   // Fock matrix
  double E = 1.0;
  double tmp,zeta,u,x;
  int ziter = 1;      // iteration counter for zeta optimization
  int iter = 1;       // SCF iteration counter
  int indDim = 0;

  const bool singlet = true;
  if (hbs.singlet() != singlet)
    throw std::runtime_error("Can only solve singlet HF problem");

  double zeta_new = basis.zeta();
  const int nbf = basis.num_bf();
  std::vector<double> coeff(basis.num_bf());
  F = block_matrix(nbf,nbf);
 
  // zeta-loop
  while(ziter == 1 || (opt && ziter > maxiter_zeta)) {
    zeta = zeta_new;
    if (opt) {
      fprintf(outfile,"\t Test message : zeta iteration #%d\n",ziter);
      fprintf(outfile,"\t Test message : zeta = %3.12lf\n",zeta);
      fflush(outfile);
    }

    SS = block_matrix(nbf,nbf);
    for(int i=0;i<nbf;i++) {
      const Orbital& bfi = basis.bf(i);
      for(int j=0;j<=i;j++) {
	const Orbital& bfj = basis.bf(j);
	SS[i][j] = SS[j][i] = S(bfi,bfj);
      }
    }
    fprintf(outfile,"\t Overlap matrix\n");
    print_mat(SS,nbf,nbf,outfile);
    indDim = nbf - sq_sqrt(SS,&X,nbf);
    fprintf(outfile,"\t%d basis functions are linearly independent\n",indDim);
    free_block(SS);
    fflush(outfile);

    //
    // One electron part of the Hamiltonian is computed in terms of Slater determinants
    //
    Hc = block_matrix(nbf,nbf);
    for(int i=0;i<nbf;i++) {
      const Orbital& bfi = basis.bf(i);
      for(int j=0;j<=i;j++) {
	const Orbital& bfj = basis.bf(j);

	double Hij = T(bfi,bfj);
	Hij -= Z * V_en(bfi,bfj);

	Hc[j][i] = Hc[i][j] = Hij;
      }
    }

    // get initial vector
    eigensolve(nbf,indDim,Hc,X,0,E,coeff);
    // compute density
    //double** D = density(coeff);

    //
    // Compute transformation matrices between bases
    //
    SDBasisSet sbs(SpinAlphaBeta,basis);
    const int ndet = sbs.num_bf();
    SD_2_Hylleraas sd_to_hyll(sbs, hbs);
    double** sXh = sd_to_hyll.X();
    CSFBasisSet cbs(SpinSinglet, sbs);
    const int ncsf = cbs.num_bf();
    const int nhyll = hbs.num_bf();
    double** hXs = block_matrix(nhyll,ndet);
    {
      SD_2_CSF sd_to_csf(sbs, cbs);
      double** cXs = sd_to_csf.X();
      CSF_2_Hylleraas csf_to_hyll(cbs,hbs);
      double** hXc = csf_to_hyll.X();
      mmult(hXc,0,cXs,0,hXs,0,nhyll,ncsf,ndet,0);
    }
    fprintf(outfile,"\tHylleraas functions in terms of determinants\n");
    print_mat(hXs,nhyll,ndet,outfile);
    double** hXh = block_matrix(nhyll,nhyll);
    mmult(hXs,0,sXh,0,hXh,0,nhyll,ndet,nhyll,0);
    fprintf(outfile,"\tShould be 1\n");
    print_mat(hXh,nhyll,nhyll,outfile);

    //
    // Two-electron part of H is computed in Hylleraas basis and transformed to CSF basis
    //
    {
      const int nhyll = hbs.num_bf();
      double** Hee_hyll = block_matrix(nhyll,nhyll);
      for(int i=0;i<nhyll;i++) {
	const HylleraasBasisFunction& bfi = hbs.bf(i);
	for(int j=0;j<=i;j++) {
	  const HylleraasBasisFunction& bfj = hbs.bf(j);
	  double Gij = V_ee(bfi,bfj);
	  Hee_hyll[j][i] = Hee_hyll[i][j] = Gij;
	}
      }

      // transform to SD basis
      Hee = block_matrix(ndet,ndet);
      double** temp = block_matrix(ndet,nhyll);
      mmult(sXh,0,Hee_hyll,0,temp,0,ndet,nhyll,nhyll,0);
      mmult(temp,0,sXh,1,Hee,0,ndet,nhyll,ndet,0);
      free_block(Hee_hyll);
      free_block(temp);
#if 0
      CSFBasisSet cbs(SpinSinglet, sbs);
      const int ncsf = cbs.num_bf();
      double** Hee_csf = block_matrix(ncsf,ncsf);
      CSF_2_hylleraas csf_to_hyll(cbs, hbs);
      double** X = csf_to_hyll.Xinv();
      double** temp = block_matrix(ncsf,nhyll);
      mmult(X,0,Hee_hyll,0,temp,0,ncsf,nhyll,nhyll,0);
      mmult(temp,0,X,1,Hee_csf,0,ncsf,nhyll,ncsf,0);
      free_block(Hee_hyll);
      free_block(temp);
#endif
    }

    //
    // SCF iterations
    // 
    double Escf = E;
    double Escf_old = 0.0;
    double dEscf = Escf - Escf_old;
    while(std::fabs(dEscf) > threshold && iter < maxiter) {
      Escf_old = Escf;

      // one electron contribution
      for(int i=0;i<nbf;i++) {
	for(int j=0;j<=i;j++) {
	  F[i][j] = F[j][i] = Hc[i][j];
	}
      }
      // two electron contribution
      int b12 = 0;
      for(int b1=0;b1<nbf;b1++) {
	for(int b2=0;b2<nbf;b2++,b12++) {
	  int k12 = 0;
	  for(int k1=0;k1<nbf;k1++) {
	    for(int k2=0;k2<nbf;k2++,k12++) {
	      F[b1][k1] += coeff[b2] * coeff[k2] * Hee[b12][k12];
	    }
	  }
	}
      }

      // get new vector
      eigensolve(nbf,indDim,F,X,0,E,coeff);

      // compute energy
      Escf = 2.0 * E;
      {
	int b12 = 0;
	for(int b1=0;b1<nbf;b1++) {
	  for(int b2=0;b2<nbf;b2++,b12++) {
	    int k12 = 0;
	    for(int k1=0;k1<nbf;k1++) {
	      for(int k2=0;k2<nbf;k2++,k12++) {
		Escf -= coeff[b1] * coeff[b2] * coeff[k1] * coeff[k2] * Hee[b12][k12];
	      }
	    }
	  }
	}
      }
      dEscf = Escf - Escf_old;
      fprintf(outfile,"\titer %d: E = %15.12lf  dE = %15.12lf  e = %15.12lf\n",iter,Escf,dEscf,E);

      ++iter;
    }


    //
    // Compute MPn wave functions
    //
    {
      // 1 = S
      double** I = block_matrix(nhyll,nhyll);
      {
	for(int i=0;i<nhyll;i++) {
	  const HylleraasBasisFunction& bfi = hbs.bf(i);
	  for(int j=0;j<=i;j++) {
	    const HylleraasBasisFunction& bfj = hbs.bf(j);
	    I[i][j] = I[j][i] = S(bfi,bfj);
	  }
	}
      }

      // Orthogonalizer (X = S^-1/2) and inverse (Xi = S^1/2 = SX)
      double** X;
      int rhyll;  // rank of Hylleraas basis
      {
	const int numDep = sq_sqrt(I,&X,nhyll);
	rhyll = nhyll - numDep;
      }
      double** Xi = block_matrix(nhyll,rhyll);
      mmult(I,0,X,0,Xi,0,nhyll,nhyll,rhyll,0);

      // 1h = XSX
      double** Ih = UtVU(X,I,nhyll,rhyll);

      // <sd|hyll>
      double** sSh = block_matrix(ndet,nhyll);
      mmult(sXh,0,I,0,sSh,0,ndet,nhyll,nhyll,0);

      // |0>
      double* V0;
      {
	double* V0s = init_array(ndet);
	int b12 = 0;
	for(int b1=0;b1<nbf;b1++) {
	  for(int b2=0;b2<nbf;b2++,b12++) {
	    V0s[b12] = coeff[b1] * coeff[b2];
	  }
	}
	V0 = Utv(sXh,V0s,ndet,nhyll);
	free(V0s);
      }
      double* V0h = Utv(Xi,V0,nhyll,rhyll);

      // |0><0|
      double** O;
      {
	double** Os = block_matrix(ndet,ndet);
	int b12 = 0;
	for(int b1=0;b1<nbf;b1++) {
	  for(int b2=0;b2<nbf;b2++,b12++) {
	    int k12 = 0;
	    for(int k1=0;k1<nbf;k1++) {
	      for(int k2=0;k2<nbf;k2++,k12++) {
		Os[b12][k12] = coeff[b1] * coeff[b2] * coeff[k1] * coeff[k2];
	      }
	    }
	  }
	}
	O = UtVU(sXh,Os,ndet,nhyll);
	free_block(Os);
      }
      fprintf(outfile,"\tO\n");
      print_mat(O,nhyll,nhyll,outfile);
      double** Oh = UtVU(Xi,O,nhyll,rhyll);
      
      // P = 1-O
      double** Ph = block_matrix(rhyll,rhyll);
      {
	for(int i=0;i<rhyll;i++) {
	  for(int j=0;j<=i;j++) {
	    const double p = Ih[i][j] - Oh[i][j];
	    Ph[i][j] = Ph[j][i] = p;
	  }
	}
      }
      
      // H0
      double** H0;
      {
	double** H0s = block_matrix(ndet,ndet);
	double** Ss = block_matrix(ndet,ndet);
	int b12 = 0;
	for(int b1=0;b1<nbf;b1++) {
	  for(int b2=0;b2<nbf;b2++,b12++) {
	    int k12 = 0;
	    for(int k1=0;k1<nbf;k1++) {
	      double S11 = S(basis.bf(b1),basis.bf(k1));
	      double F11 = F[b1][k1];
	      for(int k2=0;k2<nbf;k2++,k12++) {
		double S22 = S(basis.bf(b2),basis.bf(k2));
		double F22 = F[b2][k2];
		H0s[b12][k12] = F11*S22 + F22*S11;
		Ss[b12][k12] = S(sbs.bf(b12),sbs.bf(k12));
	      }
	    }
	  }
	}
	double** Sinv = sq_inverse(Ss,ndet);
	double** SHS = UVUt(Sinv,H0s,ndet,ndet);
	H0 = UtVU(sSh,SHS,ndet,nhyll);
	//H0 = UVUt(hXs,H0s,ndet,nhyll);
      }
      double** H0h = UtVU(X,H0,nhyll,rhyll);

      // H
      double** H = block_matrix(nhyll,nhyll);
      {
	for(int i=0;i<nhyll;i++) {
	  const HylleraasBasisFunction& bfi = hbs.bf(i);
	  for(int j=0;j<=i;j++) {
	    const HylleraasBasisFunction& bfj = hbs.bf(j);
	    double hij = T(bfi,bfj);
	    hij += Z*V_en(bfi,bfj)/2.0;
	    hij += V_ee(bfi,bfj);
	    
	    H[j][i] = H[i][j] = hij;
	  }
	}
      }
      double** Hh = UtVU(X,H,nhyll,rhyll);

      // H1 = H - H0
      double** H1 = block_matrix(nhyll,nhyll);
      {
	for(int i=0;i<nhyll;i++) {
	  for(int j=0;j<=i;j++) {
	    const double h1 = H[i][j] - H0[i][j];
	    H1[i][j] = H1[j][i] = h1;
	  }
	}
      }
      double** H1h = UtVU(X,H1,nhyll,rhyll);

      // H1|0>
      double* H10h = Utv(H1h,V0h,rhyll,rhyll);

      // H0-E0
      double** HE0h = block_matrix(rhyll,rhyll);
      {
	for(int i=0;i<rhyll;i++) {
	  for(int j=0;j<=i;j++) {
	    const double h = H0h[i][j] - 2.0*E*Ih[i][j];
	    HE0h[i][j] = HE0h[j][i] = h;
	  }
	}
      }
      fprintf(outfile,"\t(H0 - E0)h\n");
      print_mat(HE0h,rhyll,rhyll,outfile);

      // (H0-E0)^{-1}
      double**HE0h_i = sq_inverse(HE0h,rhyll);

      // P (H0-E0)^-1 P
      double** PHE0iP = UtVU(Ph,HE0h_i,rhyll,rhyll);

      // MP2 energy
      {
	double* Psi1 = Utv(PHE0iP,H10h,rhyll,rhyll);
	double emp2 = 0.0; dot_arr(H10h,Psi1,rhyll,&emp2); emp2 *= -1.0;
	fprintf(outfile,"\t E(MP2) = %15.10lf\n",emp2);
      }

      // Test: compute <0|H0-E0|0>
      double Etest0 = 0.0;
      for(int i=0;i<rhyll;i++) {
	for(int j=0;j<rhyll;j++) {
	  Etest0 += Oh[i][j] * HE0h[i][j];
	}
      }
      fprintf(outfile,"test 0: E = %15.10lf\n",Etest0);
      // Test: compute <0|H|0>
      double Etest1 = 0.0;
      for(int i=0;i<rhyll;i++) {
	for(int j=0;j<rhyll;j++) {
	  Etest1 += Oh[i][j] * Hh[i][j];
	}
      }
      fprintf(outfile,"test 1: E = %15.10lf\n",Etest1);
      // Test: compute <0|S|0>
      double Stest2 = 0.0;
      for(int i=0;i<rhyll;i++) {
	for(int j=0;j<rhyll;j++) {
	  Stest2 += Oh[i][j] * Ih[i][j];
	}
      }
      fprintf(outfile,"test 2: S = %15.10lf\n",Stest2);
      // Test: compute <0|H1|0>
      double Etest3 = 0.0;
      for(int i=0;i<rhyll;i++) {
	for(int j=0;j<rhyll;j++) {
	  Etest3 += Oh[i][j] * H1h[i][j];
	}
      }
      fprintf(outfile,"test 3: E = %15.10lf\n",Etest3);
      // Test: compute <0|H0|0>
      double Etest4 = 0.0;
      double* H00 = Utv(H0h,V0h,rhyll,rhyll);
      dot_arr(H00,V0h,rhyll,&Etest4);
      fprintf(outfile,"test 4: E = %15.10lf\n",Etest4);

    }

  
    if (!opt)
      free_block(X);
    
    ziter++;
    
    if (!opt)
      break;     /* if optimization flag is not set - quit 'while' loop */
    else
      fflush(outfile);
  }
  
#if TEST_NORMALIZATION
  /* verify c^T . S . c = 1 */
  {
    double norm = 0.0;
    for(int i=0;i<nbf;i++) {
      const Orbital& bfi = basis.bf(i);
      for(int j=0;j<nbf;j++) {
	const Orbital& bfj = basis.bf(j);
	norm += coeff[i] * coeff[j] *
	  S(bfi,bfj);
      }
    }
    fprintf(outfile,"\tTEST: c^T * S * c = %10.9lf (should be 1.0)\n",norm);
  }
#endif

  OrbitalWfn* wfn = new OrbitalWfn(basis,coeff);
  return *wfn;
}

void
eigensolve(int nbf, int indDim, double** H, double** X, int root, double& E, std::vector<double>& coeff)
{
  //fprintf(outfile,"\t Hamiltonian matrix\n");
  //print_mat(H,nbf,nbf,outfile);
    
  double** temp = block_matrix(nbf,indDim);
  mmult(H,0,X,0,temp,0,nbf,nbf,indDim,0);
  double** XtHX = block_matrix(indDim,indDim);
  mmult(X,1,temp,0,XtHX,0,indDim,nbf,indDim,0);
  free_block(temp);
  
  temp = block_matrix(indDim,indDim);
  double* evals = init_array(indDim);
  sq_rsp(indDim,indDim,XtHX,evals,1,temp,1.0E-20);
  free_block(XtHX);
  
  E = evals[root];
  //fprintf(outfile,"\tE = %3.12lf\n",E);
  
  for(int i=0;i<indDim;i++) /* temporarily putting eigenvector to evals array */
    evals[i] = temp[i][0];
  free_block(temp);
  
  for(int i=0;i<nbf;i++) {
    double c = 0.0;
    for(int j=0;j<indDim;j++)
      c += X[i][j]*evals[j];
    coeff[i] = c;
  }
  free(evals);
}

void
usage() {
  fprintf(stderr,"  USAGE: hyller --max <max> [--nmax <nmax> --lmax <lmax> --mmax <mmax>]\n");
  fprintf(stderr,"           nmax - upper limit of n (power of r1+r2)\n");
  fprintf(stderr,"           lmax - upper limit of l (power of r1-r2)\n");
  fprintf(stderr,"           mmax - upper limit of m (power of r12)\n");
  fprintf(stderr,"           max  - upper limit of n+l+m\n");
}


void test_gsh()
{
  // Test matrix elements over generalized Slater-Hylleraas functions
  {
    const int m = 1;
    const double zeta = 3.375;
    const double gamma = 0.0;
    std::cout << "Testing matrix elements of GenSlaterHylleraasBasisFunction (m=" << m <<")" << std::endl;
    HylleraasBasisFunction H000(0,0,0,zeta);
    HylleraasBasisFunction H00m(0,0,m,zeta);
    const double S_ref = S(H00m,H00m);
    const double uS_ref = S_ref / (normConst(H00m) * normConst(H00m));
    const double Ven_ref = V_en(H00m,H00m);
    const double Vee_ref = V_ee(H00m,H00m);
    const double T0m_ref = T(H000,H00m);
    const double Tm0_ref = T(H00m,H000);
    const double Tmm_ref = T(H00m,H00m);
    GenSlaterHylleraasBasisFunction GSH000(0,0,0,0.5*zeta,0.5*zeta,gamma);
    GenSlaterHylleraasBasisFunction GSH00m(0,0,m,0.5*zeta,0.5*zeta,gamma);
    const double us = S(GSH00m,GSH00m);
    const double s = us * NormConst(GSH00m) * NormConst(GSH00m);
    const double v_en = 2.0*V_ne(GSH00m,GSH00m) * NormConst(GSH00m) * NormConst(GSH00m);
    const double v_ee = V_ee(GSH00m,GSH00m) * NormConst(GSH00m) * NormConst(GSH00m);
    const double t0m = T(GSH000,GSH00m) * NormConst(GSH000) * NormConst(GSH00m);
    const double tm0 = T(GSH00m,GSH000) * NormConst(GSH00m) * NormConst(GSH000);
    const double tmm = T(GSH00m,GSH00m) * NormConst(GSH00m) * NormConst(GSH00m);
    std::cout << "Unnorm: S(ref) = " << uS_ref << "   S = " << us << std::endl;
    std::cout << "S(ref) = " << S_ref << "   S = " << s << std::endl;
    std::cout << "Ven(ref) = " << Ven_ref << "   Ven = " << v_en << std::endl;
    std::cout << "Vee(ref) = " << Vee_ref << "   Vee = " << v_ee << std::endl;
    std::cout << "mm:T(ref) = " << Tmm_ref << "   T = " << tmm << std::endl;
    std::cout << "0m:T(ref) = " << T0m_ref << "   T = " << t0m << std::endl;
    std::cout << "m0:T(ref) = " << Tm0_ref << "   T = " << tm0 << std::endl;
    std::cout << "E(ref) = " << (Ven_ref + Vee_ref + Tmm_ref) << "   E = " << (v_en + v_ee + tmm) << std::endl;
  }
  // Test symmetry of kinetic energy matrix elements over generalized Slater-Hylleraas functions
  {
    std::cout << "Checking that kinetic energy is symmetric with respect to particles 1 and 2" << std::endl;
    const int i = 2;
    const int j = 3;
    const int k = 1;
    const double alpha = 2.0;
    const double beta = 1.0;
    const double gamma = 1.0;

    GenSlaterHylleraasBasisFunction GSHijk(i,j,k,alpha,beta,gamma);
    GenSlaterHylleraasBasisFunction GSHjik(j,i,k,beta,alpha,gamma);
    const double T11 = T(GSHijk,GSHijk) * NormConst(GSHijk) * NormConst(GSHijk);
    const double T22 = T(GSHjik,GSHjik) * NormConst(GSHjik) * NormConst(GSHjik);
    const double T12 = T(GSHijk,GSHjik) * NormConst(GSHijk) * NormConst(GSHjik);
    const double T21 = T(GSHjik,GSHijk) * NormConst(GSHjik) * NormConst(GSHijk);
    std::cout << T11 << " " << T22 << std::endl;
    std::cout << T12 << " " << T21 << std::endl;
  }
  // Check matrix elements of T over generalized Slater-Hylleraas functions in terms of Hylleraas matrix elements
  {
    std::cout << "Testing kinetic energy over gen. Slater-Hylleraas" << std::endl;
    const int n = 0;
    const int l = 1;
    const int m = 7;
    const double zeta = 3.75;
    if (n + l != 1)
      throw std::runtime_error("This test only works when n+l == 1");
    const double alpha = 0.5*zeta;
    const double beta = 0.5*zeta;
    const double gamma = 0.0;

    GenSlaterHylleraasBasisFunction GSH10m(1,0,m,alpha,beta,gamma);
    GenSlaterHylleraasBasisFunction GSH01m(0,1,m,alpha,beta,gamma);
    const double sign = (l == 1 ? -1.0 : 1.0);
    const double t = (  T(GSH10m,GSH10m)*NormConst(GSH10m) * NormConst(GSH10m) + 
			T(GSH01m,GSH01m)*NormConst(GSH01m) * NormConst(GSH01m) + 
			sign*(T(GSH10m,GSH01m) + T(GSH01m,GSH10m))*NormConst(GSH10m)*NormConst(GSH01m)
		     );
    const double s = (  S(GSH10m,GSH10m)*NormConst(GSH10m) * NormConst(GSH10m) + 
			S(GSH01m,GSH01m)*NormConst(GSH01m) * NormConst(GSH01m) + 
			sign*(S(GSH10m,GSH01m) + S(GSH01m,GSH10m))*NormConst(GSH10m)*NormConst(GSH01m)
		     );

    HylleraasBasisFunction Hnlm(n,l,m,zeta);
    const double T_ref = T(Hnlm,Hnlm);
    std::cout << "n = " << n << " l = " << l << " m = " << m << "  T_ref = " << T_ref << " T = " << t/s << std::endl;
  }
  // Test matrix elements of delta functions
  {
    std::cout << "Testing delta function integrals over gen. Slater-Hylleraas" << std::endl;
    const int i = 0;
    const int j = 0;
    const int k = 0;
    const double alpha = 200.0;
    const double beta = 200.0;
    const double gamma = 0.0;

    GenSlaterHylleraasBasisFunction bf(i,j,k,alpha,beta,gamma);
    Orbital o2(j+k,beta);
    const double nc12 = NormConst(bf);
    const double nc2 = normConst(o2);
    const double s2 = S(o2,o2);
    const double delta_r1 = DeltaR1(bf,bf);
    const double s = S(bf,bf);
    const double t = T(bf,bf);
    std::cout << "nc(2) = " << nc2 << " nc(12) = " << nc12 << " s2 = " << s2 << " <delta(r1)> = " << delta_r1 << " <S> = " << s << " <T> = "<< t << std::endl;
    std::cout << "<T - 2 delta(r1)> = " << nc12 * nc12 * (t - 2.0*delta_r1) << std::endl;
  }
}

void test_gen_basisfn()
{
  typedef GenSlaterHylleraasBasisFunction GSH;
  typedef ContractedBasisFunction<GSH> CBF;

  double alpha = 1.0;
  double beta  = 2.0;
  double gamma = 3.0;

  std::cout << "Testing generic contracted basis functions" << std::endl;
  CBF c;
  c.add(GSH(0,1,0,alpha,beta,gamma),1.0);
  c.add(GSH(0,0,1,alpha,beta,gamma),-1.0);
  std::cout << "c:" << std::endl << c.to_string() << std::endl;
  CBF d(c*c);
  std::cout << "c*c:" << std::endl << d.to_string() << std::endl;
  CBF e(d*c);
  std::cout << "c*c*c:" << std::endl << e.to_string() << std::endl;
}

void test_gen_basis()
{
  typedef GenSlaterHylleraasBasisFunction GSH;
  BasisSet<GSH> bs;

  double alpha = 1.82;

#if 0
  // my 3-term wavefunction

  // Add e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,1,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.push_back(std::make_pair(GSH(2,0,0,alpha,alpha,0.0),+1.0));
    bf.push_back(std::make_pair(GSH(0,2,0,alpha,alpha,0.0),+1.0));
    bf.push_back(std::make_pair(GSH(0,0,2,alpha,alpha,0.0),-1.0));
    bs.add(bf);
  }
#endif

#if 0
  // my 3-term wavefunction uncontracted

  // Add e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,1,alpha,alpha,0.0));
  // Add r_1^2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(2,0,0,alpha,alpha,0.0));
  // Add r_2^2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,2,0,alpha,alpha,0.0));
  // Add r_{12}^2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,2,alpha,alpha,0.0));
#endif

#if 1
  // Hylleraas 3-term wave function

  // Add e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,0,alpha,alpha,0.0));
  // Add e^(- zeta r_1 - zeta r_2 - gamma r_{12})
  bs.add(GSH(0,0,1,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - 2 r_1 r_2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,0,alpha,alpha,0.0),-2.0);
    bs.add(bf);
  }
#endif

#if 0
  // my 4-term wavefunction

  // Add e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,1,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs.add(bf);
  }
  // Add (r_1^2 + r_2^2 - 2 r_1 r_2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,0,alpha,alpha,0.0),-2.0);
    bs.add(bf);
  }
#endif

#if 0
  // my 4-term wavefunction uncontracted

  // Add e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,1,alpha,alpha,0.0));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(1,1,0,alpha,alpha,0.0));
  // Add r_1^2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(2,0,0,alpha,alpha,0.0));
  // Add r_2^2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,2,0,alpha,alpha,0.0));
  // Add r_{12}^2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,2,alpha,alpha,0.0));
#endif

#if 0
  // Hartree-Fock-like wavefunction

  // Add e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_1 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(1,0,0,alpha,alpha,0.0));
  // Add r_2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(0,1,0,alpha,alpha,0.0));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs.add(GSH(1,1,0,alpha,alpha,0.0));
#endif


  double** U = bs.coefs();
  fprintf(outfile,"  -Contraction matrix for the test basis set\n");
  print_mat(U,bs.nprim(),bs.nbf(),outfile);

  solveCi<GSH>(bs,2.0,0);
}

void
test_hf_2body(const OrbitalWfn& hfwfn)
{
  const double gamma = 1.0;
  typedef GenSlaterHylleraasBasisFunction GSH;
  typedef ContractedBasisFunction<GSH> CGSH;

  const ContractedBasisFunction<Orbital> orb(contract(hfwfn));
  const CGSH phi0(operator^<GSH,Orbital,Orbital>(orb,orb));
  std::cout << "Hartree-Fock wave function (phi0):" << std::endl
	    << phi0.to_string() << std::endl;
  {
    BasisSet<GSH> bs;
    bs.add(phi0);
    // should give the Hartree-Fock energy back
    solveCi<GSH>(bs,2.0,0);
  }

  BasisSet<GSH> bs;
  // Now try wave function like phi0 * (some linear combination)
#if 0
  // my 3-term wavefunction

  // Add 1
  CGSH t0(gen_r1r2r12_oper(0,0,0));
  CGSH phi0t0(phi0*t0);
  std::cout << "t0:" << std::endl
	    << t0.to_string() << std::endl;
  std::cout << "phi0 * t0:" << std::endl
	    << phi0t0.to_string() << std::endl;
  bs.add(phi0t0);

  // Add r_{12}
  CGSH t1(GSH(0,0,1,0.0,0.0,0.0));
  CGSH phi0t1(phi0*t1);
  std::cout << "t1:" << std::endl
	    << t1.to_string() << std::endl;
  std::cout << "phi0 * t1:" << std::endl
	    << phi0t1.to_string() << std::endl;
  bs.add(phi0t1);

  // Add (r_1^2 + r_2^2 - r_{12}^2)
  {
    CGSH t2;
    t2.add(GSH(2,0,0,0.0,0.0,0.0),+1.0);
    t2.add(GSH(0,2,0,0.0,0.0,0.0),+1.0);
    t2.add(GSH(0,0,2,0.0,0.0,0.0),-1.0);
    CGSH phi0t2(phi0*t2);
    std::cout << "t2:" << std::endl
	      << t2.to_string() << std::endl;
    std::cout << "phi0 * t2:" << std::endl
	      << phi0t2.to_string() << std::endl;
    bs.add(phi0t2);
  }

#endif

#if 1
  // Hylleraas 3-term function

  // Add 1
  CGSH t0(gen_r1r2r12_oper(0,0,0));
  CGSH phi0t0(phi0*t0);
  std::cout << "t0:" << std::endl
	    << t0.to_string() << std::endl;
  std::cout << "phi0 * t0:" << std::endl
	    << phi0t0.to_string() << std::endl;
  bs.add(phi0t0);

  // Add r_{12}
  CGSH t1(GSH(0,0,1,0.0,0.0,0.0));
  CGSH phi0t1(phi0*t1);
  std::cout << "t1:" << std::endl
	    << t1.to_string() << std::endl;
  std::cout << "phi0 * t1:" << std::endl
	    << phi0t1.to_string() << std::endl;
  bs.add(phi0t1);

  // Add (r_1^2 + r_2^2 - 2 r_1 r_2)
  {
    CGSH t2;
    t2.add(GSH(2,0,0,0.0,0.0,0.0),+1.0);
    t2.add(GSH(0,2,0,0.0,0.0,0.0),+1.0);
    t2.add(GSH(1,1,0,0.0,0.0,0.0),-2.0);
    CGSH phi0t2(phi0*t2);
    std::cout << "t2:" << std::endl
	      << t2.to_string() << std::endl;
    std::cout << "phi0 * t2:" << std::endl
	      << phi0t2.to_string() << std::endl;
    bs.add(phi0t2);
  }

#endif


  solveCi<GSH>(bs,2.0,0);

}


void test_gen_energy()
{
  double alpha = 1.81607;
  double gamma = 0.00;
  double Z = 2.0;

  typedef GenSlaterHylleraasBasisFunction GSH;
  typedef SymmGSHBasisSet Basis;

#if 0
  // Hylleraas 3-term wave function
  Ptr<Basis> bs(new Basis(alpha,0.0,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,1,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - 2 r_1 r_2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,0,alpha,alpha,0.0),-2.0);
    bs->add(bf);
  }
#endif

#if 0
  // Hylleraas 3-term analog with en exponential instead of r12
  bool mutable_gamma = (gamma == 0.0)?false:true;
  Ptr<Basis> bs(new Basis(alpha,gamma,true,mutable_gamma));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add e^(- zeta r_1 - zeta r_2 - gamma r_{12})
  bs->add(GSH(0,0,0,alpha,alpha,gamma));
  // Add (r_1^2 + r_2^2 - 2 r_1 r_2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,0,alpha,alpha,0.0),-2.0);
    bs->add(bf);
  }
#endif

#if 0
  // my 3-term wavefunction
  bool mutable_gamma = (gamma == 0.0)?false:true;
  Ptr<Basis> bs(new Basis(alpha,gamma,true,mutable_gamma));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,1,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
#endif

#if 0
  // my trial wave function: (1 + r1r2) (1 + r12)
  Ptr<Basis> bs(new Basis(alpha,0.0,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,1,alpha,alpha,0.0));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,0,alpha,alpha,0.0));
  // Add r_1 r_2 r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,1,alpha,alpha,0.0));
#endif

#if 0
  // my trial wave function: (1 + r1r2) (1 + r1 \dot r2)
  Ptr<Basis> bs(new Basis(alpha,0.0,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,0,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,3,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
#endif

#if 0
  // my trial wave function: (1 + r1 \dot r2) (1 + r12)
  Ptr<Basis> bs(new Basis(alpha,0.0,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,1,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_{12} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,1,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,1,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,3,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
#endif

#if 0
  // my trial wave function: (1 + r1r2) (1 + r1 \dot r2) (1 + r12)
  Ptr<Basis> bs(new Basis(alpha,0.0,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,1,alpha,alpha,0.0));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,0,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_{12} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,1,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,1,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,3,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,3,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 r_{12} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,1,alpha,alpha,0.0));
  // Add r_1 r_2 r_{12} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,1,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,3,1,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,3,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
#endif

#if 0
  gamma = 0.05;
  // my trial wave function: (1 + r1r2) (1 + r1 \dot r2) (1 + e^{-\gamma r_{12}})
  Ptr<Basis> bs(new Basis(alpha,gamma,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add e^{-\gamma r_{12}} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,gamma));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,0,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add e^{-\gamma r_{12}} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,gamma),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,gamma),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,gamma),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,3,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 e^{-\gamma r_{12}} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,0,alpha,alpha,gamma));
  // Add r_1 r_2 e^{-\gamma r_{12}} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,0,alpha,alpha,gamma),+1.0);
    bf.add(GSH(1,3,0,alpha,alpha,gamma),+1.0);
    bf.add(GSH(1,1,2,alpha,alpha,gamma),-1.0);
    bs->add(bf);
  }
#endif

#if 1
  gamma = 0.10;
  // my trial wave function: (1 + r1r2) (1 + r1 \dot r2) (1 + r_{12} e^{-\gamma r_{12}})
  Ptr<Basis> bs(new Basis(alpha,gamma,true,false));

  // Add e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,0,alpha,alpha,0.0));
  // Add r_{12} e^{-\gamma r_{12}} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(0,0,1,alpha,alpha,gamma));
  // Add r_1 r_2 e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,0,alpha,alpha,0.0));
  // Add (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,2,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(0,0,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_{12} e^{-\gamma r_{12}} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(2,0,1,alpha,alpha,gamma),+1.0);
    bf.add(GSH(0,2,1,alpha,alpha,gamma),+1.0);
    bf.add(GSH(0,0,3,alpha,alpha,gamma),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,3,0,alpha,alpha,0.0),+1.0);
    bf.add(GSH(1,1,2,alpha,alpha,0.0),-1.0);
    bs->add(bf);
  }
  // Add r_1 r_2 r_{12} e^{-\gamma r_{12}} e^(- zeta r_1 - zeta r_2)
  bs->add(GSH(1,1,1,alpha,alpha,gamma));
  // Add r_1 r_2 r_{12} e^{-\gamma r_{12}} (r_1^2 + r_2^2 - r_{12}^2) e^(- zeta r_1 - zeta r_2)
  {
    typedef BasisSet<GSH>::ContrBF BF;
    BF bf;
    bf.add(GSH(3,1,1,alpha,alpha,gamma),+1.0);
    bf.add(GSH(1,3,1,alpha,alpha,gamma),+1.0);
    bf.add(GSH(1,1,3,alpha,alpha,gamma),-1.0);
    bs->add(bf);
  }
#endif

#if 0
  typedef Overlap<Basis> Overlap;
  Ptr<Overlap> s(new Overlap(bs));
  (*s)();
#endif

#if 1
  // we'll use this Hamiltonian ...
  typedef TwoBodyHamiltonian<Basis> Hamiltonian;
  Ptr<Hamiltonian> h(new Hamiltonian(bs,Z));
  //(*h)();
  // to get energy variationally ...
  typedef EigenEnergy<Hamiltonian> Energy;
  Ptr<Energy> energy(new Energy(0,h));
  (*energy)();
  // and optimize (basis set) parameters using Newton-Raphson method 
  typedef NewtonRaphsonOptimizer<Energy> Optimizer;
  Ptr<Optimizer> optimizer(new Optimizer(energy,1e-9,1e-3,1000));
  optimizer->optimize();
#endif
}

};
