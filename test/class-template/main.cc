#include <iostream>


class Wavefunction {
private:
	size_t Nx;
	double dx;
public:
	Wavefunction(size_t Nx, double dx): Nx(Nx), dx(dx) {}
	size_t get_Nx() { return Nx; }
};

class BohmWavefunction : public Wavefunction {
public:
	BohmWavefunction(size_t Nx, double dx): Wavefunction(Nx, dx) {}
};


template<class WF = Wavefunction>
class Propagator {
protected:
	size_t Nx;
	double dx;
public:
	WF *wf;
	Propagator(size_t Nx, double dx): Nx(Nx), dx(dx) {
		wf = new WF(Nx, dx);
	}
	~Propagator() {
		delete wf;
	}
};

class BohmPropagator : public Propagator<BohmWavefunction> {
public:
	BohmPropagator(size_t Nx, double dx): Propagator<BohmWavefunction>(Nx, dx) {}
};


int main() {

	const size_t Nx = 11;
	const double dx = 0.2;

	Propagator<> prop = Propagator<>(Nx, dx);
	std::cout << "prop.wf->get_Nx() = " << prop.wf->get_Nx() << std::endl;

	Propagator<BohmWavefunction> prop_with_bohm \
		= Propagator<BohmWavefunction>(Nx, dx);
	std::cout << "prop_with_bohm.wf->get_Nx() = " 
		<< prop_with_bohm.wf->get_Nx() << std::endl;

	BohmPropagator bohm_prop = BohmPropagator(Nx, dx);
	std::cout << "bohm_prop.wf->get_Nx() = " 
		<< bohm_prop.wf->get_Nx() << std::endl;

	return EXIT_SUCCESS;
}

