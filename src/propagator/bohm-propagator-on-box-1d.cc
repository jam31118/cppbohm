#include "../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../include/fd.h"

#include <iostream>
#include <complex>


Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(): Propagator_on_Box_1D() {}


Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(
		size_t Nx, double dx, double *Vx, double hbar, double mass):
	Propagator_on_Box_1D(Nx, dx, Vx, hbar, mass) {	

	s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, Ndim);
}


Bohm_Propagator_on_Box_1D::~Bohm_Propagator_on_Box_1D() {
	gsl_multiroot_fsolver_free(s);
}


// may come from Wavefunction
// Nx, dx, Ndim

int Bohm_Propagator_on_Box_1D::propagate(
		std::complex<double> *wf_tot, double dt, 
		int (*prop_wf)(std::complex<double> *wf, double dt, void *params),
		void *prop_wf_params,
		double *qarr, size_t Nq, double xmin) 
{

	const size_t Nx_tot = 1 + Nx + 1;
	const double xmax = xmin + (Nx_tot-1)*dx;

	struct implicit_eq_params eq_params = { 
		NULL, wf_tot, dx, xmin, dt, 1, 1, 0 
	}; // NULL and 0 should be replaced by appropriate one
	gsl_multiroot_function eq_f = {implicit_eq, Ndim, &eq_params};

	std::complex<double> *const wf = wf_tot + 1;
//	this->Propagator_on_Box_1D::propagate(wf, dt, 1);
	if (EXIT_SUCCESS != prop_wf(wf, dt, prop_wf_params)) {
		std::cerr << "[ERROR] Failed to propagate wavefunction\n";
	}

	gsl_vector *dqvec = gsl_vector_alloc(Ndim);
	for (size_t iq=0; iq<Nq; ++iq) {
		double xp = qarr[iq];
		gsl_vector_set(dqvec, 0, 0.);
		double qvec[Ndim]; qvec[0]	= { qarr[iq] };
		if ((xp<xmin) || (xp>=xmax)) { continue; }
		eq_params.qvec = qvec;
		if (EXIT_SUCCESS !=	eval_is0(xp+0., Nx_tot, dx, xmin, &eq_params.is0)) {
			std::cerr << "[ERROR] Failed to evaluate `is0`\n";
			return EXIT_FAILURE;
		}
		gsl_multiroot_fsolver_set(s, &eq_f, dqvec); 
		size_t i=0;
		const size_t max_iter = 200;
		int status;
		for (; i<max_iter; ++i) {
			status = gsl_multiroot_fsolver_iterate(s);
			if (status) { break; }
			status = gsl_multiroot_test_residual(s->f, 1e-7);
			if (status != GSL_CONTINUE) { break; }
		}
		if (i>=max_iter) {
			std::cerr << "[ERROR] Max iteration reached\n";
			return EXIT_FAILURE;
		}
		if (status != GSL_SUCCESS) {
			std::cerr << "[ERROR] The iteration failed with error: " 
				<< gsl_strerror(status) << std::endl;
			return EXIT_FAILURE;
		}
		qarr[iq] += gsl_vector_get(dqvec, 0);
	}
	gsl_vector_free(dqvec);

	return EXIT_SUCCESS;
}


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

