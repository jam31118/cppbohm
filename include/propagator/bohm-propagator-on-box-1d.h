#ifndef _BOHM_PROPAGATOR_ON_BOX_1D_H_
#define _BOHM_PROPAGATOR_ON_BOX_1D_H_

// From GNU Scientific Library
// 
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiroots.h"

// From TDSE Library
//
#include "propagator/propagator-on-box-1d.h"
#include "wf/wavefunction-on-box-1d.h"

// From this library (BOHM)
//
#include "../../include/wf/bohm-wavefunction-on-box-1d.h"


class Bohm_Propagator_on_Box_1D {

	Bohm_Wavefunction_on_Box_1D *p_wf;
	double hbar, mass;
	gsl_multiroot_fsolver *s;

public:

	Bohm_Propagator_on_Box_1D(
			size_t Nx, double dx, double hbar=1., double mass=1.);

	~Bohm_Propagator_on_Box_1D();

	int propagate(std::complex<double> *wf_tot, double dt, 
			int (*prop_wf)(std::complex<double> *wf, double dt, void *params), 
			void *prop_wf_params,
			double *qarr, size_t Nq, double xmin); 

	int _propagate_core(
			double *qarr, size_t Nq, gsl_multiroot_function *p_eq_f);
};

#endif // _BOHM_PROPAGATOR_ON_BOX_1D_H_
