#include <iostream>
#include <cstdlib>

#include "../../include/lapack.hh"

int main() {

	const int n = 2, nrhs = 1, lda = 2, idb = 2;
  double A1d[4] = {2, 1, 0, 1};
  double B1d[2] = {2, 3};
  int ipiv[2], info;

	std::cout << "Solve the following linear system:\n";
  for (int irow = 0; irow < 2; irow++) {
		for (int icol = 0; icol < 2; icol++) {
			std::cout << A1d[icol*2+irow] << " ";
		} std::cout << "| " << B1d[irow] << std::endl;
	}
	
	dgesv_(&n, &nrhs, A1d, &lda, ipiv, B1d, &idb, &info);
	if (handle_gesv_info(info) != EXIT_SUCCESS) 
	{ fprintf(stderr, "[ERROR] Failed to solve linear system.\n"); }

	std::cout << "Result:\n";
	for (int irow = 0; irow < 2; irow++) 
	{ std::cout << B1d[irow] << std::endl; }
	
	return EXIT_SUCCESS;
}

