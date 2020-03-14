#include <iostream>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiroots.h"


struct rosenbrock_params { double a, b; };


int rosenbrock_func(const gsl_vector *x, void *params, gsl_vector *f) {
	const double a = ((rosenbrock_params *) params)->a;
	const double b = ((rosenbrock_params *) params)->b;
	const double x0 = gsl_vector_get(x, 0);
	const double x1 = gsl_vector_get(x, 1);
	gsl_vector_set(f, 0, a*(1-x0));
	gsl_vector_set(f, 1, b*(x1-x0*x0));
	return GSL_SUCCESS;
}


void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3lu x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}


int main() {

	int return_code = EXIT_SUCCESS;
	gsl_multiroot_fsolver *s;
	gsl_vector *x;

	const size_t ndim = 2;
	rosenbrock_params params = { 1., 10. };
	gsl_multiroot_function f = {&rosenbrock_func, ndim, &params};

	double x_init[ndim] = {-10., -5.};
	x = gsl_vector_alloc(ndim);
	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x, 1, x_init[1]);

	s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, ndim);
	gsl_multiroot_fsolver_set(s, &f, x);


	size_t i=0;
	print_state(i, s);
	int status;
	const size_t max_iter = 200;
	for (; i<max_iter; ++i) {
		status = gsl_multiroot_fsolver_iterate(s);
		print_state(i+1, s);
		if (status != 0) { break; } // check if the solver is stuck
		status = gsl_multiroot_test_residual(s->f, 1e-7);
		if (status != GSL_CONTINUE) { break; }
	}
	std::cout << "[ LOG ] status = " << gsl_strerror(status) << std::endl;
	if (i>=max_iter) { 
		std::cerr << "[ERROR] Max iteration reached.\n";
		return_code = EXIT_FAILURE;
		goto free_and_return;
	}
	if (status != GSL_SUCCESS) {
  	if (status == GSL_ENOPROG) {
  		std::cerr << "[ERROR] No progress\n";
  		return_code = EXIT_FAILURE;
  		goto free_and_return;
  	}
  	if (status == GSL_EBADFUNC) {
  		std::cerr << "[ERROR] Bad function value\n";
  		return_code = EXIT_FAILURE;
  		goto free_and_return;
  	}
	}

	// Free allocated memory and return from main function
	// 
free_and_return:
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return return_code;
}

