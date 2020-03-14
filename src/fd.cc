#include "../include/fd.h"
#include "../include/lapack.hh"

int eval_f_and_derivs(
		const double x, const std::complex<double> *wf, 
		const size_t Nx_tot, const double dx, const double x0,
		std::complex<double> *wf_derivs) 
{
	const double xmin = x0, xmax = x0 + ( Nx_tot - 1 ) * dx;
	if ((xmin > x) || (x >= xmax)) { return EXIT_FAILURE; }
	
	const size_t Ns = 4, il = (size_t) (x - xmin) / dx;
	const size_t is0 = (il-1) + (il<1)*(1-il) + (il>Nx_tot-3)*(Nx_tot-3-il);
	double xs_minus_x[Ns];  // finite difference stencils minus x
	double *const pxsmax = xs_minus_x + Ns;
	for (double *pxs=xs_minus_x, val=xmin+is0*dx-x; 
			pxs < pxsmax; ++pxs, val += dx) { *pxs = val; }

	double A1d[Ns*Ns];
	for (double *pA=A1d, *pAmax=A1d+Ns; pA<pAmax; ++pA) { *pA = 1.; }
	for (double *pAprev=A1d, *pA=A1d+Ns, *pAmax=A1d+Ns*Ns; pA<pAmax; ) {
		for (double *pxs=xs_minus_x; pxs<pxsmax; ++pxs, ++pAprev, ++pA) 
		{ *pA = (*pAprev) * (*pxs);	}
	}
	const size_t ncol_b = 2;
	double b1d[ncol_b*Ns];
	const std::complex<double> *pwf = wf+is0, *pwfmax = wf+is0+Ns;
	std::complex<double> _wf;
	for (double *breal=b1d, *bimag=b1d+Ns; pwf < pwfmax; ++breal, ++bimag, ++pwf)
	{ 
		_wf = *pwf;
		*breal = _wf.real(), *bimag = _wf.imag();
	}

	int ipiv[Ns], info;
	dgesv_((int *)&Ns,(int *)&ncol_b,A1d,(int *)&Ns,ipiv,b1d,(int *)&Ns,&info);
	if (handle_gesv_info(info) != EXIT_SUCCESS) { return EXIT_FAILURE; }
//	{ std::cerr << "[ERROR] Failed to solve linear system.\n"; }
	
	std::complex<double> *pwf_derivs = wf_derivs, *pwf_derivs_max = wf_derivs+Ns;
	for (double *breal=b1d, *bimag=b1d+Ns; pwf_derivs < pwf_derivs_max; 
			++breal, ++bimag, ++pwf_derivs)
	{ *pwf_derivs = {*breal, *bimag}; }

//	for (size_t i=0; i<Ns; ++i) { A1d[i] = 1.; }
//	for (size_t icol=1; icol<Ns; ++icol) {
//		for (size_t irow=0; irow<Ns; ++irow) {
//			A1d[icol*Ns+irow] = A1d[(icol-1)*Ns+irow] * xs_minus_x[irow];
//		}
//	}

	return EXIT_SUCCESS;
}
