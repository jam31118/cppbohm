#include "../../include/propagator/bohm-propagator-on-box-1d-with-field.h"


Bohm_Propagator_on_Box_1D_with_field::Bohm_Propagator_on_Box_1D_with_field(
		size_t Nx, double dx, double hbar, double mass, double charge) 
	: Bohm_Propagator_on_Box_1D(Nx, dx, hbar, mass), charge(charge) {}


int _implicit_eq_with_field(
		const gsl_vector *dqvec, void *params, gsl_vector *eq) 
{

	struct _implicit_eq_with_field_params *const pp = 
		(struct _implicit_eq_with_field_params *) params;	

	const size_t Ndim = Bohm_Wavefunction_on_Box_1D::Ndim;

	// Evaluate qvec at next time step, i.e. `qvec_next = qvec + dqvec`
	// 
	double qvec_next[Ndim];
	qvec_next[0] = pp->qvec[0] + gsl_vector_get(dqvec, 0);

	// Evaluate wf and grad_q_wf
	// 
	std::complex<double> wf_q, grad_q_wf[Ndim];
	const bool check_wf_node = true;
	int stat = Bohm_Wavefunction_on_Box_1D::wf_and_grad_q_wf(&wf_q, grad_q_wf, 
			qvec_next, pp->wf_tot, pp->dx_grid, pp->xmin, pp->is0, check_wf_node);
	if (stat != EXIT_SUCCESS) { 
		throw "Failed to evaluate wf and grad_q_wf";
	}

	// Evaluate the implicit equation and return
	//
	const std::complex<double> dx_wf_q = grad_q_wf[0];
	const double v = (pp->hbar / pp->mass) * std::imag( dx_wf_q / wf_q ) \
									 - pp->charge * pp->vecpot;  
										// charge * A(t) in the presence of electromanetic field
	const double dxp = gsl_vector_get(dqvec, 0);
	const double eq_x = - dxp + pp->dt * v;
	gsl_vector_set(eq, 0, eq_x);
	return GSL_SUCCESS;
}


int Bohm_Propagator_on_Box_1D_with_field::propagate_with_field(
		std::complex<double> *wf_tot, double dt,
		int (*prop_wf_with_field)(
			std::complex<double> *wf, double dt, void *params),
		void *prop_wf_with_field_params, double *qarr, size_t Nq, double xmin, 
		double t, double (*p_A_func)(double t)) {

	//// Propagate wavefunction
	//
	std::complex<double> *const wf = wf_tot + 1;
	if (EXIT_SUCCESS != prop_wf_with_field(wf, dt, prop_wf_with_field_params)) {
		return BOHM_ERRNO_WF_PROP_FAILED;
	}
	
	//// Propagate particles
	//
	// [NOTE] the first NULL and 0 should be replaced by appropriate one
	const double At_next = (*p_A_func)(t+dt);
	struct _implicit_eq_with_field_params eq_params =
	{ NULL, wf_tot, p_wf->get_dx(), xmin, dt, hbar, mass, charge, At_next, 0 };

	gsl_multiroot_function eq_f = {
		&_implicit_eq_with_field, p_wf->Ndim, &eq_params};

	if (EXIT_SUCCESS != _propagate_core(qarr, Nq, &eq_f)) {
		return BOHM_ERRNO_PARTICLE_PROP_FAILED;
	}
	
	//// Return 
	// if propagations for both wavefunction and particles 
	// have ended without error
	return EXIT_SUCCESS;
}

