#include <iostream>
#include <complex>
#include <cmath>

#include "wf/wavefunction-on-box-1d.h"
#include "../../../include/wf/bohm-wavefunction-on-box-1d.h"


int main() {

	const size_t Nx = 21;
	const double xmin = 0., dx = 0.2;
	std::complex<double> *wf = new std::complex<double>[Nx];
	if (EXIT_SUCCESS != eval_ground_state_in_box_1d(wf, Nx, dx)) {
		std::cerr << "[ERROR] Failed to evaluate ground state of a box\n";
		return EXIT_FAILURE;
	}

	const size_t N_sample = 4;
	double *x_samples = new double[N_sample];
	int status = Bohm_Wavefunction_on_Box_1D::sample_by_stdev(
			wf, Nx, dx, xmin, x_samples, N_sample);
	if (status != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to get sample\n";
		return EXIT_FAILURE;
	}

	std::cout << "[ LOG ] Samples: ";
	for (size_t isamp = 0; isamp < N_sample; ++isamp) {
		std::cout << x_samples[isamp] << " ";
	} std::cout << std::endl;


	//// Get mean and stdard deviation
	//
	double mean_x, stdev_x;
	status = Wavefunction_on_Box_1D::mean_and_stdev_x(
			wf, Nx, dx, xmin, &mean_x, &stdev_x);
	if (status != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to evaluate mean and standard deviation\n";
	}
	std::cout << "[ LOG ] mean = " << mean_x << std::endl;
	std::cout << "[ LOG ] standard deviation = " << stdev_x << std::endl;
	std::cout << "[ LOG ] Sample range: [ mean - stdev, mean + stdev ] = [ " 
		<< mean_x - stdev_x << ", " << mean_x + stdev_x << " ]\n";

	delete [] wf;	
	delete [] x_samples;

	return EXIT_SUCCESS;
}
