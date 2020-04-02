#ifndef _BOHM_PROPAGATOR_ON_BOX_1D_WITH_FIELD_H_
#define _BOHM_PROPAGATOR_ON_BOX_1D_WITH_FIELD_H_

#include "bohm-propagator-on-box-1d.h"
#include "../bohm-errno.h"


class Bohm_Propagator_on_Box_1D_with_field : public Bohm_Propagator_on_Box_1D {

	double charge;

public:

	Bohm_Propagator_on_Box_1D_with_field(
			size_t Nx, double dx, double hbar=1., double mass=1., double charge=-1.);

	int propagate_with_field(std::complex<double> *wf_tot, double dt,
			int (*prop_wf_with_field)(
				std::complex<double> *wf, double dt, void *params),
			void *prop_wf_with_field_params, double *qarr, size_t Nq, double xmin, 
			double t, double (*p_A_func)(double t), int wf_only=0);
};



struct _implicit_eq_with_field_params	: public _implicit_eq_params {
	double charge;
	double vecpot;
	_implicit_eq_with_field_params(
			double *qvec, std::complex<double> *wf_tot, double dx_grid, 
			double xmin, double dt, double hbar, double mass, double charge, 
			double vecpot, size_t is0)
		: _implicit_eq_params(qvec, wf_tot, dx_grid, xmin, dt, hbar, mass, is0), 
		charge(charge), vecpot(vecpot) {};
};

#endif // _BOHM_PROPAGATOR_ON_BOX_1D_WITH_FIELD_H_
