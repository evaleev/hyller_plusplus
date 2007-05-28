
#ifndef _hyller_optimizer_h_
#define _hyller_optimizer_h_

#include <iostream>

#define DIAG_HESSIAN 0
#define USE_SYMM_MIX_DERIV 1

namespace hyller {

  /// Optimizes functions of type F
  template <typename F>
  class Optimizer {
  public:
    typedef F Function;
    typedef typename Function::Value Value;
    /// Convergence is measured as a norm of a residual, which usually has the same type as Value.
    typedef typename Traits<Value>::Norm Tolerance;
    typedef typename Function::PSet PSet;
    typedef typename Function::Parameter Parameter;

    Optimizer(const Ptr<Function>& f) : f_(f), converged_(false) {
    }
    virtual ~Optimizer() {
    }

    const Ptr<Function>& function() const { return f_; }

    const Tolerance& tolerance() const { return tolerance_; }
    void tolerance(const Tolerance& tol) { tolerance_ = tol; }
    unsigned int max_niter() const { return max_niter_; }
    void max_niter(unsigned int n) { max_niter_ = n; }
    bool converged() const { return converged_; }

    /// Details of the optimization are handled by derived classes
    virtual void optimize() =0;

  protected:
    // Call to indicate that convergence criteria are met
    void success() { converged_ = true; }
    // Call if comething changed and the optimization needs to be redone
    void obsolete() { converged_ = false; }

  private:
    Ptr<Function> f_;
    Tolerance tolerance_;
    unsigned int max_niter_;
    bool converged_;

  };

  template <typename F>
  class NewtonRaphsonOptimizer : public Optimizer<F> {
  public:
    typedef Optimizer<F> BaseOptimizer;
#include <baseoptimizer_typedefs.h>

    NewtonRaphsonOptimizer(const Ptr<Function>& f, const Tolerance& tol, const ParamType& disp, const unsigned int max_niter) : BaseOptimizer(f), disp_(disp) {
      BaseOptimizer::tolerance(tol);
      BaseOptimizer::max_niter(max_niter);
    }
    ~NewtonRaphsonOptimizer() {
    }

    /// Implements Optimizer::optimize()
    void optimize() {

      unsigned int iter = 0;

      while(iter < BaseOptimizer::max_niter() && !BaseOptimizer::converged()) {

	// This is very general and necessary but ...
#if 0
	// compute derivatives
	typedef typename NumericalDerivative12<Function> NDer12;
	Ptr<NDer12> der12(new NDer12(f()));
	der12->compute();

	// check convergence
	Tolerance der1_norm = norm(der12->der1());
	if (der1_norm <= tolerance())
	  success();

	// if not converged
	if (!converged()) {

	  // compute a (Newton-Raphson) step
	  std::vector<Parameter> step( (-1.0) * inverse(der12->der2()) * der12->der2());
	  
	  // massage the step
	  
	  // make a step

	  ++iter;
	}
#endif

	//
	// ... must cut corners a bit
	//

	Ptr<Function> func = BaseOptimizer::function();
	const unsigned int nparam = func->nparam();
	// figure out number of parameters to optimize 
	std::vector<unsigned int> map_param_opt_to_all;
	{
	  for(unsigned int p=0; p<nparam; ++p) {
	    if (func->param(p).mut_able()) {
	      map_param_opt_to_all.push_back(p);
	    }
	  }
	}
	const unsigned int nparam_to_opt = map_param_opt_to_all.size();

	// compute derivative numerically
	typedef std::vector<Value> Gradient;
	typedef std::vector<Gradient> Hessian;
	Gradient grad;
	Hessian hess;
	{
	  if (nparam_to_opt < 1 || nparam_to_opt > 2)
	    throw std::runtime_error("NewtonRaphsonOptimizer -- cannot yet optimize > 2 parameters");

	  if (nparam_to_opt == 1) {

	    Value e0 = func->operator()();
	    const unsigned int p = map_param_opt_to_all[0];
	    Parameter param = func->param(p);
	    const double disp = disp_;

	    param += disp;  func->param(p,param);
	    Value ep1 = func->operator()();
	    param += -disp;  func->param(p,param);

	    param += -disp;  func->param(p,param);
	    Value em1 = func->operator()();
	    param += disp;  func->param(p,param);

	    Value d1 = (ep1 - em1)/(2.0*disp);
	    grad.push_back(d1);

	    Value d2 = (ep1 + em1 - 2.0*e0)/(disp*disp);
	    Gradient hess_row(1,d2);
	    hess.push_back(hess_row);

	  }

	  if (nparam_to_opt == 2) {

	    const double disp = disp_;

	    Value e0 = func->operator()();

	    const unsigned int p0 = map_param_opt_to_all[0];
	    Parameter param0 = func->param(p0);

	    param0 += disp;  func->param(p0,param0);
	    Value e_0p1 = func->operator()();
	    param0 += -disp;  func->param(p0,param0);

	    param0 += -disp;  func->param(p0,param0);
	    Value e_0m1 = func->operator()();
	    param0 += disp;  func->param(p0,param0);

	    const unsigned int p1 = map_param_opt_to_all[1];
	    Parameter param1 = func->param(p1);

	    param1 += disp;  func->param(p1,param1);
	    Value e_1p1 = func->operator()();
	    param1 += -disp;  func->param(p1,param1);

	    param1 += -disp;  func->param(p1,param1);
	    Value e_1m1 = func->operator()();
	    param1 += disp;  func->param(p1,param1);

	    param0 += disp;  func->param(p0,param0);
	    param1 += disp;  func->param(p1,param1);
	    Value e_01p1 = func->operator()();
	    param0 += -disp;  func->param(p0,param0);
	    param1 += -disp;  func->param(p1,param1);

	    param0 += -disp;  func->param(p0,param0);
	    param1 += -disp;  func->param(p1,param1);
	    Value e_01m1 = func->operator()();
	    param0 += disp;  func->param(p0,param0);
	    param1 += disp;  func->param(p1,param1);

#if USE_SYMM_MIX_DERIV
	    param0 += disp;  func->param(p0,param0);
	    param1 += -disp;  func->param(p1,param1);
	    Value e_01pm = func->operator()();
	    param0 += -disp;  func->param(p0,param0);
	    param1 += disp;  func->param(p1,param1);

	    param0 += -disp;  func->param(p0,param0);
	    param1 += disp;  func->param(p1,param1);
	    Value e_01mp = func->operator()();
	    param0 += disp;  func->param(p0,param0);
	    param1 += -disp;  func->param(p1,param1);
#endif

	    Value d1_0 = (e_0p1 - e_0m1)/(2.0*disp);
	    grad.push_back(d1_0);
	    Value d1_1 = (e_1p1 - e_1m1)/(2.0*disp);
	    grad.push_back(d1_1);

	    Value d2_00 = (e_0p1 + e_0m1 - 2.0*e0)/(disp*disp);
	    Value d2_11 = (e_1p1 + e_1m1 - 2.0*e0)/(disp*disp);
#if !DIAG_HESSIAN
#if USE_SYMM_MIX_DERIV
	    Value d2_01 = (e_01p1 + e_01m1 - e_01pm - e_01mp)/(4.0*disp*disp);
	    Value d2_01_nonsymm = (e_01p1 + e_01m1 - 2.0*e0)/(2.0*disp*disp) - 0.5*(d2_00 + d2_11);
	    std::cout << "symm = " << d2_01 <<  "  nonsymm = " << d2_01_nonsymm << std::endl;
#else
	    Value d2_01 = (e_01p1 + e_01m1 - 2.0*e0)/(2.0*disp*disp) - 0.5*(d2_00 + d2_11);
#endif
	    Value d2_10 = d2_01;
#else
	    Value d2_01 = 0.0;
	    Value d2_10 = 0.0;
#endif

	    Gradient hess_row(2);
	    hess.push_back(hess_row);
	    hess.push_back(hess_row);
	    hess[0][0] = d2_00;
	    hess[1][1] = d2_11;
	    hess[0][1] = d2_01;
	    hess[1][0] = d2_10;

	  }

	}

	// converged?
	Value gradnorm = norm(grad);
	if (gradnorm <= BaseOptimizer::tolerance())
	  BaseOptimizer::success();

	if (BaseOptimizer::converged()) {

	  // print the final summary
	  fprintf(outfile,"  -Final Optimization Report\n");
	  fprintf(outfile,"  Param    Value    Gradient\n");
	  fprintf(outfile,"  -----  ---------  --------\n");
	  unsigned int p_opt = 0;
	  for(unsigned int p=0; p<nparam; ++p) {
	    Value value = func->param(p).value();
	    if (func->param(p).mut_able()) {
	      fprintf(outfile,"  %3d    %9.5lf  %8.4e\n",
		      p,value,grad[p_opt]);
	      ++p_opt;
	    }
	    else {
	      fprintf(outfile,"  %3d    %9.5lf  =========\n",
		      p,value);
	    }
	  }

	}
	else {

	  // compute a (Newton-Raphson) step
	  Value* gradvec = to_C_array(grad);
	  Value** hessmat = to_C_array(hess);
	  Value** hessmat_inv = sq_inverse(hessmat,nparam_to_opt);
	  Value* step = Utv(hessmat_inv,gradvec,nparam_to_opt,nparam_to_opt);
	  for(unsigned int i=0; i<nparam_to_opt; ++i)
	    step[i] *= -1.0;

	  // update parameters
	  for(unsigned int p=0; p<nparam_to_opt; ++p) {
	    const unsigned int pp = map_param_opt_to_all[p];
	    Parameter param = func->param(pp);
	    param += step[p];  func->param(pp,param);
	  }

	  // print a summary of a step
	  fprintf(outfile,"  -Summary of step %d\n", iter);
	  fprintf(outfile,"  Param  Old Value  Gradient  New Value\n");
	  fprintf(outfile,"  -----  ---------  --------  ---------\n");
	  unsigned int p_opt = 0;
	  for(unsigned int p=0; p<nparam; ++p) {
	    Value value = func->param(p).value();
	    if (func->param(p).mut_able()) {
	      fprintf(outfile,"  %3d    %9.5lf  %8.4e  %9.5lf\n",
		      p,value-step[p_opt],grad[p_opt],value);
	      ++p_opt;
	    }
	    else {
	      fprintf(outfile,"  %3d    %9.5lf  =========  =========\n",
		      p,value);
	    }
	  }

	  ++iter;
	}
      }
    }

  private:
    ParamType disp_;

  };

};

#endif
