#ifndef _FD_H_
#define _FD_H_

#include <complex>

int eval_f_and_derivs(
		const double x, const std::complex<double> *wf, 
		const size_t Nx_tot, const double dx, const double x0,
		std::complex<double> *wf_derivs);

#endif // _FD_H_
