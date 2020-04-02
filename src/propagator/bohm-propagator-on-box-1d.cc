#include "../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../include/fd.h"

#include <iostream>
#include <complex>
#include <fstream>



/**
 * Evaluate the value of implicit equation for the backward Euler method
 *
 * @param[in] a vector of displacement `dq` by time propagation by `dt` 
 *     so that `q + dq` is the position at the next time point
 * @param[in] extra parameters for evaluating the implicit equation
 * @param[out] the vector-valued implicit equation
 */

int _implicit_eq(const gsl_vector *dqvec, void *params, gsl_vector *eq) {

	struct _implicit_eq_params *const pp = (struct _implicit_eq_params *) params;	

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
	//// Propagate wavefunction
	//
	std::complex<double> *const wf = wf_tot + 1;
	if (EXIT_SUCCESS != prop_wf(wf, dt, prop_wf_params)) {
		std::cerr << "[ERROR] Failed to propagate wavefunction\n";
		return EXIT_FAILURE;
	}
	
	//// Propagate particles
	//
	// [NOTE] the first NULL and the last 0 should be replaced by appropriate one
	struct _implicit_eq_params eq_params =
	{ NULL, wf_tot, p_wf->get_dx(), xmin, dt, hbar, mass, 0 };
	
	gsl_multiroot_function eq_f = {&_implicit_eq, p_wf->Ndim, &eq_params};

	if (EXIT_SUCCESS != _propagate_core(qarr, Nq, &eq_f)) {
		std::cerr << "[ERROR] Failed to propagate particles\n"; 
		return EXIT_FAILURE;
	}

	//// Return 
	// if propagations for both wavefunction and particles 
	// have ended without error
	return EXIT_SUCCESS;
}



#ifdef DEBUG
void print_state (size_t iq, double x, size_t iter, gsl_multiroot_fsolver *s) {
  printf ("iq = %2lu xp = % .8f iter = %3lu dxp = % .8f "
          "implicit_eq(dxp) = % .8e\n",
          iq, x, iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->f, 0));
}
#endif // DEBUG


/**
 * The core part of particle propagation
 */

int Bohm_Propagator_on_Box_1D::_propagate_core(
		double *qarr, size_t Nq, gsl_multiroot_function *p_eq_f) 
{
	struct _implicit_eq_params *p_eq_params = 
		(struct _implicit_eq_params *) p_eq_f->params;

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
		p_eq_params->qvec = qarr + iq; // i.e. &qarr[iq];  // [suspicious]
		// Prepare `eq_params`: set member `is0`
		if (EXIT_SUCCESS !=	eval_is0(xp+dxp_init,Nx_tot,dx,xmin,&p_eq_params->is0))
		{
			std::cerr << "[ERROR] Failed to evaluate `is0`\n";
			return EXIT_FAILURE;
		}
		gsl_multiroot_fsolver_set(s, p_eq_f, dqvec); 
		size_t i=0; size_t max_iter = 200;
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
			std::cerr << "[ERROR] .. at xp = " << xp 
				<< ", particle index (`iq`) = " << iq 
				<< ", is0 = " << p_eq_params->is0 
				<< ", x_tot_arr[is0] = " 
				<< xmin + p_eq_params->dx_grid * p_eq_params->is0 << std::endl;

#ifdef DEBUG

			double x_is0 = xmin + p_eq_params->dx_grid * p_eq_params->is0;
			double x_is0_max = x_is0 + p_eq_params->dx_grid * (FD_STENCIL_NUM-1);
			size_t N_dxp_debug = 201;
			double dx_dxp_debug = (x_is0_max - x_is0) / (N_dxp_debug - 1);
			double *dxp_arr_debug = new double[N_dxp_debug];
			for (size_t idxp=0; idxp<N_dxp_debug; ++idxp) { 
				dxp_arr_debug[idxp] = (x_is0 + dx_dxp_debug * idxp) - xp;
			} std::cout << std::endl;
			double *_impl_eq_dxp_arr_debug = new double[N_dxp_debug];

			const size_t Ndim = Bohm_Wavefunction_on_Box_1D::Ndim;
			gsl_vector *_dxp_temp = gsl_vector_alloc(Ndim);
			gsl_vector *_fq_temp = gsl_vector_alloc(Ndim);

			for (size_t idxp=0; idxp<N_dxp_debug; ++idxp) {
				gsl_vector_set(_dxp_temp, 0, dxp_arr_debug[idxp]);
				(*(p_eq_f->f))(_dxp_temp, p_eq_params, _fq_temp);
				_impl_eq_dxp_arr_debug[idxp] = gsl_vector_get(_fq_temp, 0);
			} std::cerr << "\n";


			std::ofstream dxp_arr_file("debug-dxp-arr.bin");
			dxp_arr_file.write((char *) dxp_arr_debug, N_dxp_debug * sizeof(double));
			dxp_arr_file.close();
			std::ofstream impl_eq_dxp_arr_file("debug-impl-eq-dxp-arr.bin");
			impl_eq_dxp_arr_file.write(
					(char *) _impl_eq_dxp_arr_debug, N_dxp_debug * sizeof(double));
			impl_eq_dxp_arr_file.close();


			std::complex<double> _wf_q_temp, _grad_q_wf_temp[Ndim];
			std::complex<double> *_wf_q_temp_arr = 
				new std::complex<double>[N_dxp_debug];

			double *_debug_velo_arr = new double[N_dxp_debug];

			double _qvec_temp[Bohm_Wavefunction_on_Box_1D::Ndim];
			for (size_t idxp=0; idxp<N_dxp_debug; ++idxp) {
				_qvec_temp[0] = xp + dxp_arr_debug[idxp];
				p_wf->wf_and_grad_q_wf(&_wf_q_temp, _grad_q_wf_temp, _qvec_temp, 
						p_eq_params->wf_tot, p_eq_params->dx_grid, p_eq_params->xmin, 
						p_eq_params->is0, true);
				_wf_q_temp_arr[idxp] = _wf_q_temp;
				_debug_velo_arr[idxp] = (p_eq_params->hbar / p_eq_params->mass) 
					* std::imag(_grad_q_wf_temp[0] / _wf_q_temp);
			}

			std::ofstream wf_q_temp_arr_file("debug-wf-dxp-arr.bin");
			wf_q_temp_arr_file.write(
					(char *) _wf_q_temp_arr, N_dxp_debug * sizeof(std::complex<double>));
			wf_q_temp_arr_file.close();
			delete [] _wf_q_temp_arr;

			std::ofstream velo_temp_arr_file("debug-velocity-dxp-arr.bin");
			velo_temp_arr_file.write(
					(char *) _debug_velo_arr, N_dxp_debug * sizeof(double));
			velo_temp_arr_file.close();
			delete [] _debug_velo_arr;


			(*(p_eq_f->f))(dqvec, p_eq_params, _fq_temp);
			std::cout << "implicit_eq(0,...): " << _fq_temp->data[0] 
				<< ", dqvec[0]: " << dqvec->data[0] << std::endl;
			(*(p_eq_f->f))(s->x, p_eq_params, _fq_temp);
			std::cout << "implicit_eq(s->x,...): " << _fq_temp->data[0] 
				<< ", s->x[0]: " << s->x->data[0] << std::endl;
			std::cout << "s->f: " << s->f->data[0] << std::endl;

			gsl_vector_free(_fq_temp);
			gsl_vector_free(_dxp_temp);
			delete [] dxp_arr_debug;
			delete [] _impl_eq_dxp_arr_debug;

#endif // DEBUG

			return EXIT_FAILURE;
		}
		
		const double dxp = gsl_vector_get(s->x, 0);
		qarr[iq] += dxp;
	}
	gsl_vector_free(dqvec);

	return EXIT_SUCCESS;
}



