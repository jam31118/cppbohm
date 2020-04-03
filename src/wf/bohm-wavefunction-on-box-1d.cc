#include "../../include/wf/bohm-wavefunction-on-box-1d.h"

#include <iostream>

#include "../../include/fd.h"


Bohm_Wavefunction_on_Box_1D::Bohm_Wavefunction_on_Box_1D(
		size_t Nx, double dx, double hbar, double mass): 
	Wavefunction_on_Box_1D(Nx, dx), hbar(hbar), mass(mass) {}


int Bohm_Wavefunction_on_Box_1D::wf_and_grad_q_wf(
		std::complex<double> *wf_q, std::complex<double> grad_q_wf[Ndim],
		double qvec[Ndim], std::complex<double> *wf_tot, 
		double dx_grid, double xmin, size_t is0, bool check_node) 
{

	const double xp = qvec[0];

	std::complex<double> wf_derivs[FD_STENCIL_NUM];
	int stat = eval_f_and_derivs_with_is0(
			xp, wf_tot, dx_grid, xmin, wf_derivs, is0);
	if (stat != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to evaluate wf value and derivatives\n";
		return EXIT_FAILURE;
	}

	std::complex<double> _wf_q = wf_derivs[0], _dx_wf_q = wf_derivs[1];

	if ((_wf_q == 0.) && check_node) { 
		std::cerr << "[ERROR] Wavefunction node encountered\n";
		std::cerr << "[ERROR] The particle position xp_next=" << xp << std::endl;
		return EXIT_FAILURE;
	}

	// Store results to output variables and return
	//
	*wf_q = _wf_q; grad_q_wf[0] = _dx_wf_q;
	return EXIT_SUCCESS;
}



int Bohm_Wavefunction_on_Box_1D::sample_by_stdev(
		std::complex<double> *wf, size_t Nx, double dx, double xmin, 
		double *x_samples, size_t N_sample) {
	double mean_x, stdev_x;
	int status = Wavefunction_on_Box_1D::mean_and_stdev_x(
			wf, Nx, dx, xmin, &mean_x, &stdev_x);
	if (status != EXIT_SUCCESS) { return status; }
	if (stdev_x <= 0) { return EXIT_FAILURE; }
	if (N_sample <= 0) { return EXIT_FAILURE; }
	if (N_sample == 1) { x_samples[0] = mean_x; }
	if (N_sample > 1) {
		double sample_spacing = (2.*stdev_x) / (N_sample - 1);
		double x_sample_min = mean_x - stdev_x;
		for (double xs=x_sample_min, *ps=x_samples, *psmax=x_samples+N_sample;
				ps < psmax; ++ps, xs += sample_spacing)
		{ *ps = xs; }
	}
	return EXIT_SUCCESS;
}


