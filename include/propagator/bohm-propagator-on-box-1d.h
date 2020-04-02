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
#include "../wf/bohm-wavefunction-on-box-1d.h"


class Bohm_Propagator_on_Box_1D {

	gsl_multiroot_fsolver *s;

protected:

	double hbar, mass;
	Bohm_Wavefunction_on_Box_1D *p_wf;

public:

	Bohm_Propagator_on_Box_1D(
			size_t Nx, double dx, double hbar=1., double mass=1.);

	~Bohm_Propagator_on_Box_1D();

	int propagate(std::complex<double> *wf_tot, double dt, 
			int (*prop_wf)(std::complex<double> *wf, double dt, void *params), 
			void *prop_wf_params,
			double *qarr, size_t Nq, double xmin); 

protected:
	int _propagate_core(
			double *qarr, size_t Nq, gsl_multiroot_function *p_eq_f);
};



struct _implicit_eq_params {
	double *qvec;
	std::complex<double> *wf_tot;
	double dx_grid;
	double xmin;
	double dt;
	double hbar, mass;
	size_t is0;
	_implicit_eq_params(
			double *qvec, std::complex<double> *wf_tot, double dx_grid, 
			double xmin, double dt, double hbar, double mass, size_t is0)
		: qvec(qvec), wf_tot(wf_tot), dx_grid(dx_grid), 
		xmin(xmin), dt(dt), hbar(hbar), mass(mass), is0(is0) {}
};




#endif // _BOHM_PROPAGATOR_ON_BOX_1D_H_
