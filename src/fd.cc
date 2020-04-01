#include "../include/fd.h"
#include "../include/lapack.hh"

#include <iostream>


/**
 * Evaluate an index of the first stencil (`is0`) for finite difference method,
 * of a equidistanced (i.e. equal-spacing) spatial grid.
 *
 * @param[in] position at which the corresponding first-stencil (`is0`)
 * @Nx_tot[in] total number of spatial grid points
 * @dx[in] spacing of the spatial grid
 * @x0[in] minimum value of the grid
 * @is0[out] the first-stencil (`is0`) that is found
 */

int eval_is0(const double x, const size_t Nx_tot, const double dx, 
		const double x0, size_t *is0) {

	// Check arguments
	// 
	const double xmin = x0, xmax = x0 + ( Nx_tot - 1 ) * dx;
	if ((xmin > x) || (x >= xmax)) { return EXIT_FAILURE; }


	// Determine the index of spatial point for finite difference stencils: `is0`
	// 
	const size_t Ns = FD_STENCIL_NUM;
	const size_t Nsr = Ns / 2; // number of stencils on right side of il
	const size_t Nsl = Nsr-((Ns-1) % 2); // number of stencils on left side of il
	const int _il_int = (int) ((x - xmin) / dx);
	if (_il_int < 0) { return EXIT_FAILURE; }
	const size_t il = (size_t) abs(_il_int);
	const size_t _is0 = (il-Nsl) \
											+ (il<Nsl)*(Nsl-il) \
											+ (il>Nx_tot-1-Nsr)*(Nx_tot-1-Nsr-il);
	*is0 = _is0;
	return EXIT_SUCCESS;
}



/**
 * Evaluate the value of wavefunction and its partial derivatives 
 * at given coordinate by using finite difference method
 */

int eval_f_and_derivs_with_is0(
		const double x, const std::complex<double> *wf_tot, 
		const double dx, const double x0,
		std::complex<double> *wf_derivs, size_t is0) 
{
	const size_t Ns = FD_STENCIL_NUM;
	if (dx < 0) { return EXIT_FAILURE; }

	//// Construct matrices for finite difference: `A1d`, `b1d`
	// 
	// Construct an array of displacements of each stencil: x_arr[is]-x
	double xs_minus_x[Ns];  // values of finite difference stencils minus x
	double *const pxsmax = xs_minus_x + Ns;
	for (double *pxs=xs_minus_x, val=x0+is0*dx-x; 
			pxs < pxsmax; ++pxs, val += dx) { *pxs = val; }

	double A1d[Ns*Ns];
	for (double *pA=A1d, *pAmax=A1d+Ns; pA<pAmax; ++pA) { *pA = 1.; }
	double itay = 1;  // for making n! under the derivative of f in taylor series
	for (double *pAprev=A1d, *pA=A1d+Ns, *pAmax=A1d+Ns*Ns; pA<pAmax; itay+=1.) {
		for (double *pxs=xs_minus_x; pxs<pxsmax; ++pxs, ++pAprev, ++pA) 
		{ *pA = (*pAprev) * (*pxs) / itay;	}
	}
	const size_t ncol_b = 2;
	double b1d[ncol_b*Ns];
	const std::complex<double> *pwf = wf_tot+is0;
	const std::complex<double> *pwfmax = pwf + Ns;
	std::complex<double> _wf;
	for (double *breal=b1d, *bimag=b1d+Ns; pwf < pwfmax; ++breal, ++bimag, ++pwf)
	{ 
		_wf = *pwf;
		*breal = _wf.real(), *bimag = _wf.imag();
	}

	// Solve finite difference equation Ax=b to get derivatives of wf
	// 
	int ipiv[Ns], info;
	dgesv_((int *)&Ns,(int *)&ncol_b,A1d,(int *)&Ns,ipiv,b1d,(int *)&Ns,&info);
	if (handle_gesv_info(info) != EXIT_SUCCESS) { return EXIT_FAILURE; }
	
	// Store the results to given array `wf_derivs`
	//
	std::complex<double> *pwf_derivs = wf_derivs, *pwf_derivs_max = wf_derivs+Ns;
	for (double *breal=b1d, *bimag=b1d+Ns; pwf_derivs < pwf_derivs_max; 
			++breal, ++bimag, ++pwf_derivs)
	{ *pwf_derivs = {*breal, *bimag}; }

	// Return
	return EXIT_SUCCESS;
}



int eval_f_and_derivs(
		const double x, const std::complex<double> *wf, 
		const size_t Nx_tot, const double dx, const double x0,
		std::complex<double> *wf_derivs, size_t *is0) 
{
	if (EXIT_SUCCESS != eval_is0(x, Nx_tot, dx, x0, is0)) {	return EXIT_FAILURE; }
	return eval_f_and_derivs_with_is0(x, wf, dx, x0, wf_derivs, *is0); 
}



