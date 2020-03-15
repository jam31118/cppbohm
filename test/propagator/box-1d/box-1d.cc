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

//	int return_code = EXIT_SUCCESS;

//	double *Vx;
//	Bohm_Propagator_on_Box_1D prop;
//	std::complex<double> *wf;

	ParamFile param;
	try {	param = ParamFile("in.param"); }
	catch (std::exception& e) { 
		std::cerr << "[ERROR] During constructing param: " 
			<< e.what() << std::endl; 
	}
	catch (...) { 
		std::cerr << "[ERROR] Failed to construct param file with reason unknown"; 
	}

	const size_t Nx = param.get_int("Nx");
	const double dx = param.get_double("dx");
	const double dt = param.get_double("dt");
	const double xmin = 0.0;
	
	const size_t Nq = param.get_int("Nq");  // number of particles


	double *Vx = new double[Nx];
	set_to_zeros(Vx, Nx);

	Bohm_Propagator_on_Box_1D prop(Nx, dx, Vx);

	const size_t Nx_tot = 1 + Nx + 1;
	std::complex<double> *wf = new std::complex<double>[Nx_tot];
	set_to_randoms(wf+1, Nx);
	wf[0] = 0.; wf[Nx_tot-1] = 0.;

	if (prop.propagate_to_ground_state(wf+1, dt, 5000, 1e-10) != EXIT_SUCCESS) {
		std::cerr << "[ERROR] Failed to propagator to ground state\n";
		return EXIT_FAILURE;
//		return_code = EXIT_FAILURE;
//		goto free_and_return;
	}

	// Test eval wf_derivs for given paritlce positions
	double *qarr = new double[Nq];

//	set_to_randoms(qarr, Nq, xmin, xmin+(Nx+1)*dx);
	for (size_t i=0; i<Nq; ++i) { qarr[i] = i * (Nx_tot-1)*dx / (Nq-1); }

	std::cout << "[ LOG ] Particle coordinates: ";
	for (size_t i=0; i<Nq; ++i) { std::cout << qarr[i] << " "; }
	std::cout << std::endl;


	

	std::complex<double> *wf_qarr = new std::complex<double>[Nq];
	std::complex<double> *dx_wf_qarr = new std::complex<double>[Nq];



//	std::cout << "qarr: " << qarr << std::endl;
//	std::cout << "qarr[0]: " << qarr[0] << std::endl;
//	std::cout << "qarr[1]: " << qarr[1] << std::endl;



	std::complex<double> wf_derivs[FD_STENCIL_NUM];
	int stat;
	size_t is0;
	for (size_t i=0; i<Nq; ++i) {
//		std::cout << "qarr[i]: " << qarr[i] << std::endl;
		stat = eval_f_and_derivs(qarr[i], wf, Nx_tot, dx, xmin, wf_derivs, &is0);
//		std::cout << "is0=" << is0 << "x=" << qarr[i] << std::endl;
//		std::cout << "qarr[i] (after): " << qarr[i] << std::endl;

		if (stat != EXIT_SUCCESS) {
			std::cout << "[ERROR] Failed to evaluate eval_f_and_derivs\n";
			return EXIT_FAILURE;
//			return_code = EXIT_FAILURE;	
//			goto free_and_return;
		}

		wf_qarr[i] = wf_derivs[0];
		dx_wf_qarr[i] = wf_derivs[1];

//		std::cout << "[ LOG ] wf_derivs of particle x=" << qarr[i] << ": ";
//		for (size_t is=0; is<FD_STENCIL_NUM; ++is) {
//			std::cout << wf_derivs[is] << " ";
//		}	std::cout << std::endl;
	}


//	std::cout << "[ LOG ] Particle coordinates: ";
//	for (size_t i=0; i<Nq; ++i) { std::cout << qarr[i] << " "; }


//	std::cout << "qarr: " << qarr << std::endl;
//	std::cout << "qarr[0]: " << qarr[0] << std::endl;
//	std::cout << "qarr[1]: " << qarr[1] << std::endl;


	double *xarr = new double[Nx_tot];
	for (double *px=xarr, *pxmax=xarr+Nx_tot, val=xmin; px<pxmax; ++px, val+=dx) { *px = val; }


	std::ofstream xarr_file("xarr.bin", std::ios::binary);
	xarr_file.write((char *) xarr, Nx_tot * sizeof(double));
	xarr_file.close();


//	std::cout << "Nq: " << Nq << std::endl;

	std::ofstream qarr_file("qarr.bin", std::ios::binary);
	qarr_file.write((char *) qarr, Nq * sizeof(double));
//	if (!qarr_file.good()) { 
//		std::cerr << "[ERROR] Failed to write qarr\n"; return EXIT_FAILURE; }
	qarr_file.close();
//	if (! qarr_file) { 
//		std::cerr << "[ERROR] Something strange with qrr_file\n";
//		return EXIT_FAILURE;
//	} 

//	std::cout << "qarr: " << qarr << std::endl;

	std::ofstream wf_qarr_file("wf_qarr.bin", std::ios::binary);
	wf_qarr_file.write((char *) wf_qarr, Nq * sizeof(std::complex<double>));
	wf_qarr_file.close();	


	std::ofstream dx_wf_qarr_file("dx_wf_qarr.bin", std::ios::binary);
	dx_wf_qarr_file.write((char *) dx_wf_qarr, Nq * sizeof(std::complex<double>));
	dx_wf_qarr_file.close();


	
	std::ofstream wf_file("wf.bin", std::ios::binary);
	wf_file.write((char *) wf, Nx_tot * sizeof(std::complex<double>));
	wf_file.close();



//free_and_return:
//	prop.~Bohm_Propagator_on_Box_1D();
	delete [] Vx;	
	delete [] wf;

	delete [] xarr;
	delete [] qarr;

	delete [] wf_qarr;
	delete [] dx_wf_qarr;

	return EXIT_SUCCESS;
}
