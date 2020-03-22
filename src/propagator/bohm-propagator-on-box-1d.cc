#include "../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../include/fd.h"

#include <iostream>
#include <complex>


struct _implicit_eq_params {
	double *qvec;
	std::complex<double> *wf_tot;
	double dx_grid;
	double xmin;
	double dt;
	double hbar, mass;
	size_t is0;
};


int _implicit_eq(const gsl_vector *dqvec, void *params, gsl_vector *eq) {

	struct _implicit_eq_params *const pp = (struct _implicit_eq_params *) params;	

	// Evaluate qvec at next time step, i.e. `qvec_next = qvec + dqvec`
	// 
	double qvec_next[Bohm_Wavefunction_on_Box_1D::Ndim];
	for (size_t idim=0; idim<Bohm_Wavefunction_on_Box_1D::Ndim; ++idim)
	{ qvec_next[idim] = pp->qvec[idim] + gsl_vector_get(dqvec, idim);	}
	
	// Evaluate wf and grad_q_wf
	// 
	std::complex<double> wf_q, grad_q_wf[Bohm_Wavefunction_on_Box_1D::Ndim];
	const bool check_wf_node = true;
	int stat = Bohm_Wavefunction_on_Box_1D::wf_and_grad_q_wf(&wf_q, grad_q_wf, 
			qvec_next, pp->wf_tot, pp->dx_grid, pp->xmin, pp->is0, check_wf_node);
	if (stat != EXIT_SUCCESS) { 
		std::cerr << "[ERROR] Failed to evaluate wf and grad_q_wf\n"; 
		throw "Failed to evaluate wf and grad_q_wf";
	}

	// Evaluate the implicit equation and return
	//
	const std::complex<double> dx_wf_q = grad_q_wf[0];
	const double v = (pp->hbar / pp->mass) * std::imag( dx_wf_q / wf_q );
	const double dxp = gsl_vector_get(dqvec, 0);
	const double eq_x = - dxp + pp->dt * v;
	gsl_vector_set(eq, 0, eq_x);
	return GSL_SUCCESS;
}


Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(
		size_t Nx, double dx, double hbar, double mass): hbar(hbar), mass(mass)
{	
	p_wf = new Bohm_Wavefunction_on_Box_1D(Nx, dx, hbar, mass);
	s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, p_wf->Ndim);
}


Bohm_Propagator_on_Box_1D::~Bohm_Propagator_on_Box_1D() {
	delete p_wf;
	gsl_multiroot_fsolver_free(s);
}


int Bohm_Propagator_on_Box_1D::propagate(
		std::complex<double> *wf_tot, double dt, 
		int (*prop_wf)(std::complex<double> *wf, double dt, void *params),
		void *prop_wf_params,
		double *qarr, size_t Nq, double xmin) 
{

	// [NOTE] NULL and 0 should be replaced by appropriate one
	struct _implicit_eq_params eq_params =
	{ NULL, wf_tot, p_wf->get_dx(), xmin, dt, 1, 1, 0 };
	
	gsl_multiroot_function eq_f = {_implicit_eq, p_wf->Ndim, &eq_params};

	std::complex<double> *const wf = wf_tot + 1;
	if (EXIT_SUCCESS != prop_wf(wf, dt, prop_wf_params)) {
		std::cerr << "[ERROR] Failed to propagate wavefunction\n";
		return EXIT_FAILURE;
	}
	
	if (EXIT_SUCCESS != _propagate_core(qarr, Nq, &eq_f)) {
		std::cerr << "[ERROR] Failed to propagate particles\n"; 
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}


int Bohm_Propagator_on_Box_1D::_propagate_core(
		double *qarr, size_t Nq, gsl_multiroot_function *p_eq_f) 
{
	struct _implicit_eq_params *p_eq_params = 
		(_implicit_eq_params *) p_eq_f->params;

	const double xmin = p_eq_params->xmin;
	const double xmax = p_wf->get_xmax(xmin);
	const size_t Nx_tot = p_wf->get_Nx_tot();
	const double dx = p_wf->get_dx();

	const double residual_tol = 1e-7;
	const double dxp_init = 0.;
	gsl_vector *const dqvec = gsl_vector_alloc(p_wf->Ndim);
	gsl_vector_set(dqvec, 0, dxp_init);
	for (size_t iq=0; iq<Nq; ++iq) {
		double xp = qarr[iq];
		if ((xp<xmin) || (xp>=xmax)) { 
			std::cout << "[ LOG ] particles out of bound: [" 
				<< xmin << "," << xmax << "]"; 
			continue; 
		}
		// Prepare `eq_params`: set member `qvec`
		p_eq_params->qvec = qarr + iq; // i.e. &qarr[iq];
		// Prepare `eq_params`: set member `is0`
		if (EXIT_SUCCESS !=	eval_is0(xp+dxp_init,Nx_tot,dx,xmin,&p_eq_params->is0))
		{
			std::cerr << "[ERROR] Failed to evaluate `is0`\n";
			return EXIT_FAILURE;
		}
		gsl_multiroot_fsolver_set(s, p_eq_f, dqvec); 
		size_t i=0; const size_t max_iter = 200;
		int status;
		for (; i<max_iter; ++i) {
			status = gsl_multiroot_fsolver_iterate(s);
			if (status) { break; }
			status = gsl_multiroot_test_residual(s->f, residual_tol);
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
		
		const double dxp = gsl_vector_get(s->x, 0);
		qarr[iq] += dxp;
	}
	gsl_vector_free(dqvec);

	return EXIT_SUCCESS;
}



