#include "../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../include/fd.h"

#include <iostream>
#include <complex>


Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(): Propagator_on_Box_1D() {}

Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(
		size_t Nx, double dx, double *Vx, double hbar, double mass):
	Propagator_on_Box_1D(Nx, dx, Vx, hbar, mass) {	
}

Bohm_Propagator_on_Box_1D::~Bohm_Propagator_on_Box_1D() {}

int implicit_eq(const gsl_vector *dqvec, void *params, gsl_vector *eq) {

	struct implicit_eq_params *const pp = (struct implicit_eq_params *) params;	

	const double xp = pp->qvec[0];
	const double dxp = gsl_vector_get(dqvec, 0);
	const double xp_next = xp + dxp;

	std::complex<double> wf_derivs[FD_STENCIL_NUM];
	int stat = eval_f_and_derivs_with_is0(
			xp_next, pp->wf_tot, pp->dx_grid, pp->xmin, wf_derivs, pp->is0);
	if (stat != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to evaluate wf value and derivatives\n";
		throw "wf_derivs evaluation failure";
	}

	std::complex<double> wf_q = wf_derivs[0], dx_wf_q = wf_derivs[1];

	if (wf_q == 0.) { 
		std::cerr << "[ERROR] Wavefunction node encountered\n";
		std::cerr << "[ERROR] The particle position xp_next=" 
			<< xp_next << std::endl;
		throw "wavefunction node encountered"; 
	}

	double v = (pp->hbar / pp->mass) * std::imag( dx_wf_q / wf_q );
	double eq_x = - dxp + pp->dt * v;
	gsl_vector_set(eq, 0, eq_x);

	return GSL_SUCCESS;
}

