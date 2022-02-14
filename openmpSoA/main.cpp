#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "start_fdtd+pml.h"

int main(int argc, char *argv[]) {

	//setlocale(LC_ALL, "english");

	//omp_set_num_threads(1);
    if (argc>1) {
		omp_set_num_threads(std::atoi(argv[1]));
		std::cout<<"num_tread = " << std::atoi(argv[1]) << std::endl;
	}

	double n = 4;

	double T = 8.0 / 3.0 * M_PI;
	int delta_x = 8, delta_y = 8, delta_z = 8;
	int Nx = 256, Ny = 128, Nz = 32;
	std::cout << "Time: " << T << std::endl;
	run_fdtd<double, double>(n, Nx, Ny, Nz, (double)T, (double)0.005,
		delta_x, delta_y, delta_z, 46.5, 46.5);

	return 0;
}