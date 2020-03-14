#ifndef _BOHM_PROPAGATOR_ON_BOX_1D_H_
#define _BOHM_PROPAGATOR_ON_BOX_1D_H_

#include "propagator/propagator-on-box-1d.h"

class Bohm_Propagator_on_Box_1D: public Propagator_on_Box_1D {
public:
	Bohm_Propagator_on_Box_1D();
	Bohm_Propagator_on_Box_1D(
			size_t Nx, double dx, double *Vx, double hbar=1, double mass=1);
	~Bohm_Propagator_on_Box_1D();
};

#endif // _BOHM_PROPAGATOR_ON_BOX_1D_H_
