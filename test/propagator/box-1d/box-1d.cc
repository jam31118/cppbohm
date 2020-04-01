// Headers from the standards
//
#include <iostream>
#include <fstream>

// Headers from GNU Scientific Library (GSL)
//
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multiroots.h"

// Headers from TDSE package
//
#include "array.h"
#include "wf/wavefunction-on-box-1d.h"
#include "propagator/propagator-on-box-1d.h"

// Headers from PARAM package
//
#include "param.h"

// Headers from this package (BOHM)
//
#include "../../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../../include/fd.h"


struct wf_prop_params {
	Propagator_on_Box_1D *p_wf_propagator;
};

int prop_wf(std::complex<double> *wf, double dt, void *params) {
	struct wf_prop_params *pparams = (wf_prop_params *) params;
	return pparams->p_wf_propagator->propagate(wf, dt, 1);
}





int main() {

	// Get parameters from input parameter file
	// 
	ParamFile param;
	try {	param = ParamFile("in.param"); }
	catch (std::exception& e) { 
		std::cerr << "[ERROR] During constructing param: " 
			<< e.what() << std::endl; 
		return EXIT_FAILURE;
	}
	catch (...) { 
		std::cerr << "[ERROR] Failed to construct param file with unknown reason"; 
		return EXIT_FAILURE;
	}
	const size_t Nx = param.get_int("Nx");
	const double dx = param.get_double("dx");
	const double xmin = param.get_double("xmin");
	const double dt = param.get_double("dt");
	const size_t Nq = param.get_int("Nq");  // number of particles
	const size_t Nt = param.get_int("Nt");

	const size_t Nx_tot = 1 + Nx + 1;


	// Prepare potential array
	//
	double *Vx = new double[Nx];
	set_to_zeros(Vx, Nx);


	// Construct propagator
	//
	Bohm_Propagator_on_Box_1D prop(Nx, dx);  // for particles
	Propagator_on_Box_1D wf_propagator(Nx, dx, Vx);  // for wavefuntion
	struct wf_prop_params pparams = { &wf_propagator };


	// Prepare initial state
	// 
	std::complex<double> *const wf_tot = new std::complex<double>[Nx_tot];
	std::complex<double> *const wf_tot_max = wf_tot + Nx_tot;
	std::complex<double> *const wf = wf_tot + 1;

	std::string wf_t0_fname;
	try { wf_t0_fname = param.get_string("wf-t0-file"); }
	catch (...) { 
		std::cout << "[ LOG ] No file for initial wavefunction found. "
			"Falling back to search for ground state\n";
	}

	if (!wf_t0_fname.empty()) {

		// Initialize wavefunction from input file
		// 
		try { 
			std::ifstream wf_t0_file(wf_t0_fname, std::ios::binary);
			wf_t0_file.read((char *) wf_tot, Nx_tot * sizeof(std::complex<double>));
			wf_t0_file.close();
		} catch (...) {
			std::cerr << "[ERROR] Failed to read initial wavefunction from: " 
				<< wf_t0_fname << std::endl;
			return EXIT_FAILURE;
		}

	} else {

		// Initialize wavefunction as that of lowest energy possible
		// 
		set_to_randoms(wf, Nx);
		wf_tot[0] = 0.; wf_tot[Nx_tot-1] = 0.;
		// Proagate into the state with lowest energy possible
		if (wf_propagator.propagate_to_ground_state(wf, dt, 20000, 1e-10) != EXIT_SUCCESS) {
			std::cerr << "[ERROR] Failed to propagator to ground state\n";
			return EXIT_FAILURE;
		}

	}

	// Set particle coordinates 
	// 
	double *const qarr = new double[Nq];
	double *const qarr_max = qarr + Nq;
	for (size_t i=0; i<Nq; ++i) 
	{ qarr[i] = (i+1) * (Nx_tot/2-1)*dx / (Nq-1+2) + (Nx_tot/4) * dx; }


	// Print particle coordinates
	//
	std::cout << "[ LOG ] Particle coordinates: ";
	for (size_t i=0; i<Nq; ++i) { std::cout << qarr[i] << " "; }
	std::cout << std::endl;


	// Prepare storage for propagation
	// 
	std::complex<double> *wf_t_1d = new std::complex<double>[Nt*Nx_tot];
	std::complex<double> **wf_t = new std::complex<double>*[Nt];
	set_2d_view_of_1d(wf_t, wf_t_1d, Nt, Nx_tot);

	double *qarr_t_1d = new double[Nt*Nq];
	double **qarr_t = new double*[Nt];
	set_2d_view_of_1d(qarr_t, qarr_t_1d, Nt, Nq);
	

	
	// Propagete
	//
	std::copy(wf_tot, wf_tot_max, wf_t[0]);
	std::copy(qarr, qarr_max, qarr_t[0]);
	int stat;
	for (size_t it=1; it<Nt; ++it) {
		stat = prop.propagate(wf_tot, dt, &prop_wf, &pparams, qarr, Nq, xmin);
		if (stat != EXIT_SUCCESS) {
			std::cerr << "[ERROR] Failed to propagate with particles\n";
			return EXIT_FAILURE;
		}
		std::copy(wf_tot, wf_tot_max, wf_t[it]);
		std::copy(qarr, qarr_max, qarr_t[it]);
	}
	

	// Store data to files
	// 	
	std::ofstream wf_file("wf.bin", std::ios::binary);
	wf_file.write((char *) wf_tot, Nx_tot * sizeof(std::complex<double>));
	wf_file.close();

	std::ofstream wf_t_file("wf_t.bin", std::ios::binary);
	wf_t_file.write((char *) wf_t_1d, Nt*Nx_tot*sizeof(std::complex<double>));
	wf_t_file.close();

	std::ofstream qarr_t_file("qarr_t.bin", std::ios::binary);
	qarr_t_file.write((char *) qarr_t_1d, Nt*Nq*sizeof(double));
	qarr_t_file.close();



	// Free memory and return
	//
	delete [] Vx;	
	delete [] wf_tot;

	delete [] qarr;

	delete [] wf_t;
	delete [] wf_t_1d;

	delete [] qarr_t;
	delete [] qarr_t_1d;

	return EXIT_SUCCESS;
}
