#ifndef _FD_H_
#define _FD_H_

#include <complex>

#define FD_STENCIL_NUM 4

int eval_is0(const double x, const size_t Nx_tot, const double dx, 
		const double x0, size_t *is0);

int eval_f_and_derivs_with_is0(
		const double x, const std::complex<double> *wf, 
		const double dx, const double x0,
		std::complex<double> *wf_derivs, size_t is0);

int eval_f_and_derivs(
		const double x, const std::complex<double> *wf, 
		const size_t Nx_tot, const double dx, const double x0,
		std::complex<double> *wf_derivs, size_t *is0);

#endif // _FD_H_
