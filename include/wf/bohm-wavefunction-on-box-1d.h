#ifndef _BOHM_WAVEFUNCTION_ON_BOX_1D_H_
#define _BOHM_WAVEFUNCTION_ON_BOX_1D_H_

#include "wf/wavefunction-on-box-1d.h"


class Bohm_Wavefunction_on_Box_1D: public Wavefunction_on_Box_1D {
	double hbar, mass;
public:
	Bohm_Wavefunction_on_Box_1D(
			size_t Nx, double dx, double hbar=1., double mass=1.);
	static int wf_and_grad_q_wf(
			std::complex<double> *wf_q, std::complex<double> grad_q_wf[Ndim],
			double qvec[Ndim], std::complex<double> *wf_tot, 
			double dx_grid, double xmin, size_t is0, bool check_node=true); 
	static int sample_by_stdev(
			std::complex<double> *wf, size_t Nx, double dx, double xmin, 
			double *x_samples, size_t N_sample);
};



#endif // _BOHM_WAVEFUNCTION_ON_BOX_1D_H_
