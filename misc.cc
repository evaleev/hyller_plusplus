
#include <stdexcept>
#include "misc.h"
#include "includes.h"

extern "C" FILE* outfile;

namespace hyller {
#if CORRECT_RATHER_THAN_FAST
  Precomputed<true> cache(HYLLER_MAXINDEX);
#else
  Precomputed<false> cache(HYLLER_MAXINDEX);
#endif


/* This function computes Matrix^(-1/2),
   result stored in *Result (dim x dim-num_dep) (must not be allocated yet).
   dim is dimentionality, returns number of linearly dependent vectors */
  int sq_sqrt(double **Matrix, double ***Result, int dim)
  {
    int i,j;
    int num_dep = 0; /* number of linearly dependent vectors */
    double *eval;
    double **U;
    
    eval = init_array(dim);
    U = block_matrix(dim,dim);
    
    sq_rsp(dim,dim,Matrix,eval,1,U,1.0E-18);
    
    for(i=0;i<dim;i++)
      if (eval[i] < 1.0E-15)
	num_dep++;
    
    /* Allocating truncated X matrix */
    (*Result) = block_matrix(dim,dim - num_dep);
    
    for(i=0;i<dim;i++)
      for(j=num_dep;j < dim;j++)
	(*Result)[i][j-num_dep] = U[i][j]/sqrt(eval[j]);
    
    free_block(U);
    if (num_dep != 0) fprintf(outfile,"\tBasis set is truncated\n");
    return num_dep;
  }

  /* This function computes Matrix^(-1) */
  double** sq_inverse(double **Matrix, int dim)
  {
    double *eval;
    double **U;
    
    eval = init_array(dim);
    U = block_matrix(dim,dim);
    
    sq_rsp(dim,dim,Matrix,eval,1,U,1.0E-18);

    fprintf(outfile,"\t sq_inverse: evals\n");
    for(int i=0; i<dim; i++) {
      fprintf(outfile,"\t %10.6e\n",eval[i]);
    }

#if 1
    int rank = dim;
#else
    int rank = 0;
    for(int i=dim-1;i>=0;i--) {
      if (fabs(eval[i]) > 1.0E-15)
	++rank;
      else
	i = -1;
      }
#endif
    
    /* Allocating truncated X matrix */
    double** Inverse = block_matrix(dim,dim);
    
    for(int i=0;i<dim;i++)
      for(int j=0;j<dim;j++)
	for(int k=0;k<rank;k++)
	   Inverse[i][j] += U[i][dim-1-k]*U[j][dim-1-k]/eval[dim-1-k];
    
    free_block(U);
    free(eval);
    return Inverse;
  }

  double** UtVU(double** U, double**V, int i, int o)
  {
    double** utv = block_matrix(o,i);
    mmult(U,1,V,0,utv,0,o,i,i,0);
    double** utvu = block_matrix(o,o);
    mmult(utv,0,U,0,utvu,0,o,i,o,0);
    free_block(utv);
    return utvu;
  }

  double** UVUt(double** U, double**V, int i, int o)
  {
    double** uv = block_matrix(o,i);
    mmult(U,0,V,0,uv,0,o,i,i,0);
    double** uvut = block_matrix(o,o);
    mmult(uv,0,U,1,uvut,0,o,i,o,0);
    free_block(uv);
    return uvut;
  }

  double** UtVX(double** U, int i1, int o1, double** V, double** X, int i2, int o2)
  {
    double** utv = block_matrix(o1,i2);
    mmult(U,1,V,0,utv,0,o1,i1,i2,0);
    double** utvx = block_matrix(o1,o2);
    mmult(utv,0,X,0,utvx,0,o1,i2,o2,0);
    free_block(utv);
    return utvx;
  }

  double* Utv(double** U, double* v, int i, int o)
  {
    double* utv = init_array(o);
    for(int r=0; r<i; r++)
      for(int c=0; c<o; c++)
	utv[c] += U[r][c] * v[r];
    return utv;
  }

  double* Uv(double** U, double* v, int i, int o)
  {
    double* utv = init_array(o);
    for(int r=0; r<i; r++)
      for(int c=0; c<o; c++)
	utv[c] += U[c][r] * v[r];
    return utv;
  }

  double** add(int nrow, int ncol, double a, double** A, double b, double** B)
  {
    double** C = block_matrix(nrow,ncol);
    double* aptr = &(A[0][0]);
    double* bptr = &(B[0][0]);
    double* cptr = &(C[0][0]);
    for(unsigned int r=0; r<nrow; ++r) {
      for(unsigned int c=0; c<ncol; ++c, ++aptr, ++bptr, ++cptr) {
	*cptr += a*(*aptr) + b*(*bptr);
      }
    }
    return C;
  }

};
