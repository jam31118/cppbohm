// Headers from the standards
//
#include <iostream>
#include <fstream>

// Headers from TDSE package
//
#include "array.h"

// Headers from PARAM package
//
#include "param.h"

// Headers from this package (BOHM)
//
#include "../../../include/propagator/bohm-propagator-on-box-1d.h"
#include "../../../include/fd.h"


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
		std::cerr << "[ERROR] Failed to construct param file with reason unknown"; 
		return EXIT_FAILURE;
	}
	const size_t Nx = param.get_int("Nx");
	const double dx = param.get_double("dx");
	const double xmin = param.get_double("xmin");
	const double dt = param.get_double("dt");
	const size_t Nq = param.get_int("Nq");  // number of particles


	// Prepare potential array
	//
	double *Vx = new double[Nx];
	set_to_zeros(Vx, Nx);


	// Construct propagator
	//
	Bohm_Propagator_on_Box_1D prop(Nx, dx, Vx);


	// Prepare initial state
	// 
	const size_t Nx_tot = 1 + Nx + 1;
	std::complex<double> *wf_tot = new std::complex<double>[Nx_tot];
	std::complex<double> *wf = wf_tot + 1;
	set_to_randoms(wf, Nx);
	wf_tot[0] = 0.; wf_tot[Nx_tot-1] = 0.;
	// Proagate into the state with lowest energy possible
	if (prop.propagate_to_ground_state(wf, dt, 5000, 1e-10) != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to propagator to ground state\n";
		return EXIT_FAILURE;
	}


	// Set particle coordinates 
	// 
	double *qarr = new double[Nq];
//	set_to_randoms(qarr, Nq, xmin, xmin+(Nx+1)*dx);
	for (size_t i=0; i<Nq; ++i) { qarr[i] = (i+1) * (Nx_tot-1)*dx / (Nq-1+2); }


	// Print particle coordinates
	//
	std::cout << "[ LOG ] Particle coordinates: ";
	for (size_t i=0; i<Nq; ++i) { std::cout << qarr[i] << " "; }
	std::cout << std::endl;


	std::complex<double> *wf_qarr = new std::complex<double>[Nq];
	std::complex<double> *dx_wf_qarr = new std::complex<double>[Nq];

	std::complex<double> wf_derivs[FD_STENCIL_NUM];
	int stat;
	size_t is0;
	for (size_t i=0; i<Nq; ++i) {
		stat = eval_f_and_derivs(qarr[i], wf_tot, Nx_tot, dx, xmin, wf_derivs, &is0);
		if (stat != EXIT_SUCCESS) {
			std::cout << "[ERROR] Failed to evaluate wf value and derivatives "
				"(i.e. wf, dx_wf, ...)\n";
			std::cout << "[ERROR] The particle position x=" << qarr[i] << std::endl;
			return EXIT_FAILURE;
		}
		wf_qarr[i] = wf_derivs[0];
		dx_wf_qarr[i] = wf_derivs[1];
	}


	// Construct spatial array
	//
	double *xarr = new double[Nx_tot];
	for (double *px=xarr, *pxmax=xarr+Nx_tot, val=xmin; px<pxmax; ++px, val+=dx)
	{ *px = val; }


	// Store data to files
	// 
	std::ofstream xarr_file("xarr.bin", std::ios::binary);
	xarr_file.write((char *) xarr, Nx_tot * sizeof(double));
	xarr_file.close();

	std::ofstream qarr_file("qarr.bin", std::ios::binary);
	qarr_file.write((char *) qarr, Nq * sizeof(double));
	qarr_file.close();

	std::ofstream wf_qarr_file("wf_qarr.bin", std::ios::binary);
	wf_qarr_file.write((char *) wf_qarr, Nq * sizeof(std::complex<double>));
	wf_qarr_file.close();	

	std::ofstream dx_wf_qarr_file("dx_wf_qarr.bin", std::ios::binary);
	dx_wf_qarr_file.write((char *) dx_wf_qarr, Nq * sizeof(std::complex<double>));
	dx_wf_qarr_file.close();
	
	std::ofstream wf_file("wf.bin", std::ios::binary);
	wf_file.write((char *) wf_tot, Nx_tot * sizeof(std::complex<double>));
	wf_file.close();


	// Free memory and return
	//
	delete [] Vx;	
	delete [] wf_tot;

	delete [] xarr;
	delete [] qarr;

	delete [] wf_qarr;
	delete [] dx_wf_qarr;

	return EXIT_SUCCESS;
}
