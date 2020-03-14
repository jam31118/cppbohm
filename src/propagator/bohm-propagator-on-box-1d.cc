#include "../../include/propagator/bohm-propagator-on-box-1d.h"

Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(): Propagator_on_Box_1D() {}

Bohm_Propagator_on_Box_1D::Bohm_Propagator_on_Box_1D(
		size_t Nx, double dx, double *Vx, double hbar, double mass):
	Propagator_on_Box_1D(Nx, dx, Vx, hbar, mass) {	
}

Bohm_Propagator_on_Box_1D::~Bohm_Propagator_on_Box_1D() {}

