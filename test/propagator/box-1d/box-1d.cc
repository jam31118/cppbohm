#include <iostream>

#include "array.h"

#include "../../../include/propagator/bohm-propagator-on-box-1d.h"


int main() {
	const size_t Nx = 11;
	const double dx = 0.2;

	double *Vx = new double[Nx];
	set_to_zeros(Vx, Nx);

	Bohm_Propagator_on_Box_1D prop = Bohm_Propagator_on_Box_1D(Nx, dx, Vx);

	delete [] Vx;	
	return EXIT_SUCCESS;
}
