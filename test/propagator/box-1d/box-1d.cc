#include <iostream>

#include "array.h"

#include "../../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../../include/fd.h"


int main() {

	int return_code = EXIT_SUCCESS;

	double *Vx;
	Bohm_Propagator_on_Box_1D prop;
	std::complex<double> *wf;

	const size_t Nx = 11;
	const double dx = 0.2;
	const double dt = 0.05;

	Vx = new double[Nx];
	set_to_zeros(Vx, Nx);

	prop = Bohm_Propagator_on_Box_1D(Nx, dx, Vx);

	const size_t Nx_tot = 1 + Nx + 1;
	wf = new std::complex<double>[Nx_tot];
	set_to_randoms(wf, Nx_tot);

	if (prop.propagate_to_ground_state(wf+1, dt, 5000, 1e-10) != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to propagator to ground state\n";
		return_code = EXIT_FAILURE;
		goto free_and_return;
	}

free_and_return:
	prop.~Bohm_Propagator_on_Box_1D();
	delete [] Vx;	
	delete [] wf;
	return return_code;
}
