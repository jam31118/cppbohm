#ifndef _BOHM_PROPAGATOR_ON_BOX_1D_H_
#define _BOHM_PROPAGATOR_ON_BOX_1D_H_

// From GNU Scientific Library
// 
#include "gsl/gsl_vector.h"

// From TDSE Library
//
#include "propagator/propagator-on-box-1d.h"


class Bohm_Propagator_on_Box_1D: public Propagator_on_Box_1D {
public:
	Bohm_Propagator_on_Box_1D();
	Bohm_Propagator_on_Box_1D(
			size_t Nx, double dx, double *Vx, double hbar=1, double mass=1);
	~Bohm_Propagator_on_Box_1D();
};

struct implicit_eq_params {
	double *qvec;
	std::complex<double> *wf_tot;
	double dx_grid;
	double xmin;
	double dt;
	double hbar, mass;
	size_t is0;
};

int implicit_eq(const gsl_vector *dqvec, void *params, gsl_vector *eq);


#endif // _BOHM_PROPAGATOR_ON_BOX_1D_H_
